#!/bin/bash
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)

# Input structure
#
#    [input]
#    ├── sub-01
#    │   ├── freesurfer
#    │   │   └─── aparc.DKTatlas+aseg.mgz
#    │   ├── sub-01__fa.nii.gz
#    │   ├── sub-01__fodf.nii.gz
#    │   └── sub-01__t1_warped.nii.gz
#    │
#    ├── S2
#    .
#    .

usage() { echo "$(basename $0) [-s path/to/subjects] [-m path/to/mni] [-o output_dir] [-t nb_threads] [-p step_size] [-e theta_deg] -g true" 1>&2; exit 1; }

while getopts "s:m:o:t:g:p:e::" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        o) o=${OPTARG};;
        t) t=${OPTARG};;
        g) g=${OPTARG};;
        p) step_size=${OPTARG};;  # step size for tracking
        e) theta=${OPTARG};;      # theta (deg) for tracking
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${m}" ] || [ -z "${o}" ]; then
    usage
fi

subject_dir=${s}
mni_dir=${m}
out_dir=${o}
nb_thread=${t}
step_size=${step_size}
theta=${theta}


if [ -n "${step_size}" ] && [ -n "${theta}" ]; then
    step_list=(${step_size})
    theta_list=(${theta})
else
    step_list=(0.1 0.5 1.0)
    theta_list=(20 30 40)
fi

fa_threshold=0.15 #before 0.2
npv_first_order=2222  #20000
opposite_side=leftright

gpu=""
if [ ! -z "${g}" ]; then
    gpu="--use_gpu"
fi

echo "Folder subjects: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Output folder: " ${out_dir}
echo "Use GPU: " ${gpu}
echo "Number of threads" ${nb_thread}
echo "Tracking grid: steps=${step_list[*]}  thetas=${theta_list[*]}"

#for nsub in ${subject_dir}/*/; do
    nsub=${subject_dir}
    nsub=$(basename "$nsub")

    
    mkdir -p ${out_dir}/orig_space/{rois,tracking_first_order,transfo}
    mkdir -p ${out_dir}/mni_space/{rois,tracking_first_order}

    # Per-combo dirs will be created later, inside the loops.
    orig_rois_dir=${out_dir}/orig_space/rois
    mni_rois_dir=${out_dir}/mni_space/rois
    orig_tracking_dir=${out_dir}/orig_space/tracking_first_order
    mni_tracking_root=${out_dir}/mni_space/tracking_first_order

    # =========================
    # Grid over (step, theta)
    # =========================
    for step_size in "${step_list[@]}"; do
      for theta in "${theta_list[@]}"; do
        combo_tag=step_${step_size}_theta_${theta}
        echo "|=== Running combo: ${combo_tag} ===|"

       
        mkdir -p ${orig_tracking_dir}/${combo_tag}/orig
        mkdir -p ${mni_tracking_root}/${combo_tag}/{orig,filtered,segmented,final}
        mni_tracking_dir=${mni_tracking_root}/${combo_tag}

        echo "|------------- 4) [ORIG-SPACE] Generate local tractography with inclusions ROI (${combo_tag}) -------------|"
        for nside in left right; do
            # Single run per side per combo 
            scil_tracking_local ${subject_dir}/tractoflow/${nsub}__fodf.nii.gz \
                ${orig_rois_dir}/${nsub}_cp_${nside}_orig.nii.gz \
                ${orig_rois_dir}/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
                ${orig_tracking_dir}/${combo_tag}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk \
                --npv $npv_first_order \
                --step ${step_size} \
                --theta ${theta} \
                ${gpu} -v -f
        done
        echo "|------------- 4) Done -------------|"
        echo ""

        echo "|------------- 5) Register Tracking in MNI space (${combo_tag}) -------------|"
        for nside in left right; do
            scil_tractogram_apply_transform \
                ${orig_tracking_dir}/${combo_tag}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk \
                ${mni_dir}/MNI/mni_masked.nii.gz \
                ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
                ${mni_tracking_dir}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk \
                --in_deformation ${out_dir}/orig_space/transfo/2orig_1Warp.nii.gz \
                --remove_invalid \
                --reverse_operation -f
        done
        echo "|------------- 5) Done -------------|"
        echo ""

        echo "|------------- 6) [MNI-SPACE] Filter tractography (${combo_tag}) -------------|"
        for nside in left right; do
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk \
                ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered_${combo_tag}.trk \
                --drawn_roi ${mni_rois_dir}/${nsub}_any_exclusion_roi_mni.nii.gz 'any' 'exclude' \
                --drawn_roi ${mni_rois_dir}/${nsub}_${nside/${nside}/${opposite_side/${nside}/}}_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
                --drawn_roi ${mni_rois_dir}/${nsub}_${nside}_cerebellum_wm_mni.nii.gz 'either_end' 'exclude' \
                --drawn_roi ${mni_dir}/MNI/midsagittal_plane.nii.gz 'any' 'exclude' -f
        done
        echo "|------------- 6) Done -------------|"
        echo ""

        echo "|------------- 7) [MNI-SPACE] Segmentation and cleaning tractography (${combo_tag}) -------------|"
        for nside in left right; do
            echo "|------------- 7) Segmentation - SIDE: ${nside} -------------|"
            ## Mesencephalic Tract (Top)
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered_${combo_tag}.trk \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_mesencephalic_${combo_tag}.trk  \
                --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'include' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane_for_mesencephalic.nii.gz 'any' 'include' \
                -f

            ## Spinal Tract (bottom)
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered_${combo_tag}.trk \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_spinal_${combo_tag}.trk  \
                --drawn_roi ${mni_dir}/MNI/lowest_cut_brainstem.nii.gz 'any' 'include' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
                -f

            ## Two remaining nucleus/tract
            scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered/${nsub}_${nside}_from_cp_filtered_${combo_tag}.trk \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_remaining_cp_${combo_tag}.trk  \
                --drawn_roi ${mni_dir}/MNI/lower_cut_brainstem.nii.gz 'any' 'exclude' \
                --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'exclude' \
                --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
                --bdo ${mni_dir}/MNI/sphere_exclusion_for_remaining_cp.bdo 'any' 'exclude' \
                -f

            echo "|------------- 8) Cleaning - SIDE: ${nside} -------------|"
            scil_bundle_reject_outliers \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_mesencephalic_${combo_tag}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_mesencephalic_${combo_tag}.trk -f

            scil_bundle_reject_outliers \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_spinal_${combo_tag}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_spinal_${combo_tag}.trk -f



# Length-filter ONLY the left spinal bundle: 61 < length < 66 mm
if [ "${nside}" = "left" ]; then
  scil_tractogram_filter_by_length \
    "${mni_tracking_dir}/final/${nsub}_${nside}_spinal_${combo_tag}.trk" \
    "${mni_tracking_dir}/final/${nsub}_${nside}_spinal_len61_66_${combo_tag}.trk" \
    --minL 61 --maxL 66 --display_counts

  # Overwrite the original filename so the rest of the pipeline (incl. merge) stays unchanged
  mv -f \
    "${mni_tracking_dir}/final/${nsub}_${nside}_spinal_len61_66_${combo_tag}.trk" \
    "${mni_tracking_dir}/final/${nsub}_${nside}_spinal_${combo_tag}.trk"
fi

    

            scil_bundle_reject_outliers \
                ${mni_tracking_dir}/segmented/${nsub}_${nside}_remaining_cp_${combo_tag}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp_${combo_tag}.trk \
                --alpha 0.97

            scil_tractogram_filter_by_orientation \
                ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp_${combo_tag}.trk \
                ${mni_tracking_dir}/final/${nsub}_${nside}_remaining_cp_${combo_tag}.trk \
                --max_z 7 --use_abs -f
        done
        echo "|------------- 7/8) Done (${combo_tag}) -------------|"
        echo ""
      done
    done

    # =========================
    # Merge/concatenate outputs across all combos into one place
    # =========================
    echo "|============= Concatenating final outputs across combos for ${nsub} =============|"
    merged_dir=${mni_tracking_root}/final_merged
    mkdir -p "${merged_dir}"

    for nside in left right; do
        # Collect per-combo files (only keep those that actually exist)
        mes_files=()
        spi_files=()
        rem_files=()
        for step_size in "${step_list[@]}"; do
          for theta in "${theta_list[@]}"; do
            combo_tag=step_${step_size}_theta_${theta}
            f_mes="${mni_tracking_root}/${combo_tag}/final/${nsub}_${nside}_mesencephalic_${combo_tag}.trk"
            f_spi="${mni_tracking_root}/${combo_tag}/final/${nsub}_${nside}_spinal_${combo_tag}.trk"
            f_rem="${mni_tracking_root}/${combo_tag}/final/${nsub}_${nside}_remaining_cp_${combo_tag}.trk"
            [[ -f "$f_mes" ]] && mes_files+=("$f_mes")
            [[ -f "$f_spi" ]] && spi_files+=("$f_spi")
            [[ -f "$f_rem" ]] && rem_files+=("$f_rem")
          done
        done

        # Concatenate only if we have at least 1 input (or 2+ ideally)
        if (( ${#mes_files[@]} )); then
          scil_tractogram_math -f concatenate "${mes_files[@]}" \
            "${merged_dir}/${nsub}_${nside}_mesencephalic_merged.trk"
        else
          echo "WARN: No mesencephalic files found for ${nsub} ${nside}; skipping."
        fi

        if (( ${#spi_files[@]} )); then
          scil_tractogram_math -f concatenate "${spi_files[@]}" \
            "${merged_dir}/${nsub}_${nside}_spinal_merged.trk"
        else
          echo "WARN: No spinal files found for ${nsub} ${nside}; skipping."
        fi

        if (( ${#rem_files[@]} )); then
          scil_tractogram_math -f concatenate "${rem_files[@]}" \
            "${merged_dir}/${nsub}_${nside}_remaining_cp_merged.trk"
        else
          echo "WARN: No remaining_cp files found for ${nsub} ${nside}; skipping."
        fi
    done
    echo "|============= Concatenation done for ${nsub} =============|"
#done
