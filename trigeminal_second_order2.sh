#!/bin/bash
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)
#
# Second-order script (PATH-COMPATIBLE with your first-order output structure):
#   first-order final bundles are expected in:
#     <out_dir>/<sub>/mni_space/tracking_first_order/final_merged/final/*.trk
#
# Usage:
#   bash trigeminal_second_order.sh -s <subjects_parent_or_single_subject> -m <ROIs_clean> -o <out_dir> [-t threads] [-g true]

set -euo pipefail

usage() {
  echo "Usage: $(basename "$0") -s <subjects_parent_or_single_subject_dir> -m <ROIs_clean_dir> -o <out_dir> [-t threads] [-g true]" 1>&2
  exit 1
}

s=""
m=""
o=""
g=""
t=""

while getopts "s:m:o:g:t:" args; do
  case "${args}" in
    s) s=${OPTARG} ;;
    m) m=${OPTARG} ;;
    o) o=${OPTARG} ;;
    g) g=${OPTARG} ;;
    t) t=${OPTARG} ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${m}" ] || [ -z "${o}" ]; then
  usage
fi

subject_dir="${s}"
mni_dir="${m}"
out_dir="${o}"

gpu=""
if [ -n "${g}" ]; then
  gpu="--use_gpu"
fi

if [ -n "${t}" ]; then
  export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${t}"
fi

# -------------------------------
# Tracking params (original)
# -------------------------------
npv_from_spinal_track_long=1000
npv_from_spinal_track_short=100
npv_from_thalamus_track=500

# -------------------------------
# FIRST-ORDER final folder (YOUR structure)
# -------------------------------
first_order_final_rel="mni_space/tracking_first_order/final_merged/final"
group_final_dir="${out_dir}/${first_order_final_rel}"
mkdir -p "${group_final_dir}"

echo "Folder subjects: ${subject_dir}"
echo "Folder MNI:      ${mni_dir}"
echo "Output folder:   ${out_dir}"
echo "Threads:         ${t:-"(default)"}"
echo "GPU:             ${gpu:-"(cpu)"}"
echo "First-order final expected at: ${out_dir}/*/${first_order_final_rel}/"
echo ""

# -------------------------------
# Detect single-subject vs parent folder
# -------------------------------
subjects=()
subjects_parent=""

if [ -d "${subject_dir}/tractoflow" ] && [ -d "${subject_dir}/freesurfer" ]; then
  sub_name="$(basename "${subject_dir}")"
  subjects=( "${sub_name}" )
  subjects_parent="$(dirname "${subject_dir}")"
else
  subjects_parent="${subject_dir}"
  shopt -s nullglob
  for d in "${subjects_parent}"/*/ ; do
    subjects+=( "$(basename "$d")" )
  done
  shopt -u nullglob
fi

if [ ${#subjects[@]} -eq 0 ]; then
  echo "ERROR: No subjects found under -s ${subject_dir}"
  exit 2
fi

# -------------------------------
# 0) Build group density masks from first-order bundles (spinal + remaining_cp)
# -------------------------------
echo "|------------- 0) Building group density masks from first-order bundles -------------|"

for nside in left right; do
  shopt -s nullglob
  spinal_inputs=( "${out_dir}"/*/"${first_order_final_rel}"/*_"${nside}"_spinal.trk )
  rcp_inputs=( "${out_dir}"/*/"${first_order_final_rel}"/*_"${nside}"_remaining_cp.trk )
  shopt -u nullglob

  if [ ${#spinal_inputs[@]} -eq 0 ]; then
    echo "ERROR: No spinal tracts found for side=${nside} at:"
    echo "  ${out_dir}/*/${first_order_final_rel}/*_${nside}_spinal.trk"
    exit 3
  fi
  if [ ${#rcp_inputs[@]} -eq 0 ]; then
    echo "ERROR: No remaining_cp tracts found for side=${nside} at:"
    echo "  ${out_dir}/*/${first_order_final_rel}/*_${nside}_remaining_cp.trk"
    exit 3
  fi

  echo "Union (${nside}) spinal: ${#spinal_inputs[@]} files -> ${group_final_dir}/all_${nside}_spinal.trk"
  scil_tractogram_math union "${spinal_inputs[@]}" \
    "${group_final_dir}/all_${nside}_spinal.trk" -f

  echo "Density (${nside}) spinal seed"
  scil_tractogram_compute_density_map \
    "${group_final_dir}/all_${nside}_spinal.trk" \
    "${group_final_dir}/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
    --binary -f

  echo "Union (${nside}) remaining_cp: ${#rcp_inputs[@]} files -> ${group_final_dir}/all_${nside}_remaining_cp.trk"
  scil_tractogram_math union "${rcp_inputs[@]}" \
    "${group_final_dir}/all_${nside}_remaining_cp.trk" -f

  echo "Density (${nside}) remaining_cp"
  scil_tractogram_compute_density_map \
    "${group_final_dir}/all_${nside}_remaining_cp.trk" \
    "${group_final_dir}/all_${nside}_remaining_cp_density_mni.nii.gz" \
    --binary -f
done

echo "|------------- 0) Done -------------|"
echo ""

# -------------------------------
# Subject loop
# -------------------------------
for nsub in "${subjects[@]}"; do
  subj_path="${subjects_parent}/${nsub}"

  orig_rois_dir="${out_dir}/${nsub}/orig_space/rois"
  mni_rois_dir="${out_dir}/${nsub}/mni_space/rois"
  orig_tracking_dir="${out_dir}/${nsub}/orig_space/tracking_second_order"
  mni_tracking_dir_second_order="${out_dir}/${nsub}/mni_space/tracking_second_order"

  mkdir -p "${out_dir}/${nsub}/orig_space/"{rois,tracking_second_order,transfo}
  mkdir -p "${out_dir}/${nsub}/orig_space/tracking_second_order/orig"
  mkdir -p "${out_dir}/${nsub}/mni_space/"{rois,tracking_second_order}
  mkdir -p "${out_dir}/${nsub}/mni_space/tracking_second_order/"{orig,filtered,segmented,final,cut}

  echo "|------------- PROCESSING SECOND-ORDER FIBERS TRACTOGRAPHY FOR TRIGEMINAL SYSTEM -------------|"
  echo "|------------- FOR DATASET ${nsub} -------------|"
  echo ""

  # Required first-order transforms (produced by your first-order)
  if [ ! -f "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" ] || \
     [ ! -f "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" ]; then
    echo "ERROR: Missing first-order transforms for subject ${nsub} in:"
    echo "  ${out_dir}/${nsub}/orig_space/transfo/"
    exit 4
  fi

  # Inputs
  if [ ! -f "${subj_path}/tractoflow/${nsub}__t1_warped.nii.gz" ]; then
    echo "ERROR: Missing: ${subj_path}/tractoflow/${nsub}__t1_warped.nii.gz"
    exit 5
  fi
  if [ ! -f "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" ]; then
    echo "ERROR: Missing: ${subj_path}/tractoflow/${nsub}__fodf.nii.gz"
    exit 5
  fi

  # WM mask: do NOT hardcode FA threshold; reuse whatever first-order created
  shopt -s nullglob
  wm_candidates=( "${orig_rois_dir}/${nsub}_wm_mask_"*"_orig.nii.gz" )
  shopt -u nullglob
  if [ ${#wm_candidates[@]} -eq 0 ]; then
    echo "ERROR: No WM mask found (expected from first-order) in:"
    echo "  ${orig_rois_dir}/${nsub}_wm_mask_*_orig.nii.gz"
    exit 6
  fi
  wm_mask="${wm_candidates[0]}"
  echo "Using WM mask: ${wm_mask}"
  echo ""

  echo "|------------- 1) Generate local tractography with spinal bundle -------------|"

  echo "|------------- 1.2) Seeding Mask + Registration in orig space -------------|"
  for nside in left right; do
    # copy group seed to subject MNI rois folder (original behavior)
    cp "${group_final_dir}/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz"

    # register density into orig space
    antsApplyTransforms \
      -d 3 \
      -i "${group_final_dir}/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      -r "${subj_path}/tractoflow/${nsub}__t1_warped.nii.gz" \
      -t "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
      -t "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
      -o "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz"
  done

  echo "|------------- 1.3) Extract thalamus -------------|"
  if [ ! -f "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz" ] || \
     [ ! -f "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz" ]; then
    echo "ERROR: Missing aparc volumes produced by first-order in:"
    echo "  ${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz"
    echo "  ${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz"
    exit 7
  fi

  Right_Thalamus=(49)
  Left_Thalamus=(10)

  scil_labels_combine "${orig_rois_dir}/${nsub}_right_thalamus_orig.nii.gz" \
    --volume_ids "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz" ${Right_Thalamus[*]} --merge_groups -f
  scil_labels_combine "${mni_rois_dir}/${nsub}_right_thalamus_mni.nii.gz" \
    --volume_ids "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz" ${Right_Thalamus[*]} --merge_groups -f

  scil_labels_combine "${orig_rois_dir}/${nsub}_left_thalamus_orig.nii.gz" \
    --volume_ids "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz" ${Left_Thalamus[*]} --merge_groups -f
  scil_labels_combine "${mni_rois_dir}/${nsub}_left_thalamus_mni.nii.gz" \
    --volume_ids "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz" ${Left_Thalamus[*]} --merge_groups -f

  echo "|------------- 1.4) Tracking -------------|"
  for nside in left right; do
    echo "|------------- 1.4a) Tracking from Spinal bundle - npv ${npv_from_spinal_track_long} -------------|"
    scil_tracking_local \
      "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" \
      "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz" \
      "${wm_mask}" \
      "${orig_tracking_dir}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk" \
      --npv ${npv_from_spinal_track_long} -v -f ${gpu}

    echo "|------------- 1.4b) Tracking from Spinal bundle - npv ${npv_from_spinal_track_short} -------------|"
    scil_tracking_local \
      "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" \
      "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz" \
      "${wm_mask}" \
      "${orig_tracking_dir}/orig/${nsub}_${nside}_from_spinal_track_npv100.trk" \
      --npv ${npv_from_spinal_track_short} -v -f ${gpu}

    echo "|------------- 1.4c) Tracking from Thalamus - npv ${npv_from_thalamus_track} -------------|"
    scil_tracking_local \
      "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" \
      "${orig_rois_dir}/${nsub}_${nside}_thalamus_orig.nii.gz" \
      "${wm_mask}" \
      "${orig_tracking_dir}/orig/${nsub}_${nside}_from_thalamus_npv500.trk" \
      --npv ${npv_from_thalamus_track} -v -f ${gpu}
  done

  echo "|------------- 1.5) Register Tracking in MNI space -------------|"
  for ntracking in from_thalamus_npv500.trk from_spinal_track_npv100.trk from_spinal_track_npv1000.trk; do
    for nside in left right; do
      scil_tractogram_apply_transform \
        "${orig_tracking_dir}/orig/${nsub}_${nside}_${ntracking}" \
        "${mni_dir}/MNI/mni_masked.nii.gz" \
        "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
        "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_${ntracking}" \
        --in_deformation "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
        --remove_invalid --reverse_operation -f
    done
  done

  echo "|------------- 1) Done -------------|"
  echo ""

  echo "|------------- 2) Generate ROIs -------------|"
  echo "|------------- 2.1) Registration VPM -------------|"

  for nside in left right; do
    cp "${mni_dir}/MNI/Distal/${nside}/VPM.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz"

    for nroi in "${mni_dir}/MNI/from_${nside}"/*; do
      ROI_basename=$(basename "$nroi")
      cp "$nroi" "${mni_rois_dir}/${nsub}_second_order_${ROI_basename/nii/_mni.nii}"
    done
  done

  echo "|------------- 2.2) Generate masks -------------|"
  for nside in left right; do
    if [ "$nside" == "left" ]; then
      contra_nside="right"
    else
      contra_nside="left"
    fi

    # DTTT ipsilat dPSN cuts: Remaining_CP to VPM
    scil_volume_math union \
      "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" \
      "${group_final_dir}/all_${nside}_remaining_cp_density_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz" \
      --data_type uint8 -f
    scil_labels_from_mask \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz" -f

    # DTTT ipsilat CS cuts: Spinal seed to VPM
    scil_volume_math union \
      "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz" \
      --data_type uint8 -f
    scil_labels_from_mask \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz" -f

    # DTTT controlat CS cuts: contra thalamus + spinal seed
    scil_volume_math union \
      "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz" \
      --data_type uint8 -f
    scil_labels_from_mask \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz" -f

    # VTTT controlat OS/IS cuts: contra thalamus + spinal seed
    scil_volume_math union \
      "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz" \
      --data_type uint8 -f
    scil_labels_from_mask \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz" -f

    # VTTT controlat vPSN cuts: contra VPM + remaining_cp density
    scil_volume_math union \
      "${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz" \
      "${group_final_dir}/all_${nside}_remaining_cp_density_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz" \
      --data_type uint8 -f
    scil_labels_from_mask \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz" -f
  done

  echo "|------------- 2) Done -------------|"
  echo ""

  echo "|------------- 3) Generate Second-order bundles for trigeminal system -------------|"
  for nside in left right; do
    if [ "$nside" == "left" ]; then
      contra_nside="right"
    else
      contra_nside="left"
    fi

    echo "|------------- 3.1) From ${nside} - VTTT (only controlateral) - OS/IS and vPSN -------------|"

    # OS and IS
    scil_tractogram_filter_by_roi \
      "${mni_tracking_dir_second_order}/orig/${nsub}_${contra_nside}_from_thalamus_npv500.trk" \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk" \
      --drawn_roi "${mni_rois_dir}/${nsub}_left_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_rois_dir}/${nsub}_right_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_rois_dir}/${nsub}_${contra_nside}_thalamus_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_Pons_Controlat.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_CaudalMedulla_Controlat.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Pons_Ipsilat.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/cs_plaque.nii.gz" 'any' 'exclude' -f

    # vPSN
    scil_tractogram_filter_by_roi \
      "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk" \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk" \
      --drawn_roi "${mni_rois_dir}/${nsub}_left_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_rois_dir}/${nsub}_right_cerebellum_wm_mni.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'either_end' 'include' \
      --drawn_roi "${mni_rois_dir}/${nsub}_${contra_nside}_VPM_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_CaudalMedulla_Controlat.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_INC_VTT_Area.nii.gz" 'any' 'include' -f

    echo "|------------- 3.2) ${nside} - DTTT (controlateral) - CS -------------|"
    scil_tractogram_filter_by_roi \
      "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_thalamus_npv500.trk" \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${contra_nside}_DTTT_Controlat_CS.trk" \
      --drawn_roi "${mni_rois_dir}/${nsub}_${contra_nside}_spinal_density_second_order_seed_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_thalamus_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_INC_CaudalMedulla_Ipsilat.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_INC_Medulla_Controlat.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${contra_nside}/DTTT_Controlat_EXC_Midbrain_Ipsilat.nii.gz" 'any' 'exclude' -f

    echo "|------------- 3.3) ${nside} - DTTT (ipsilateral) - dPSN and CS -------------|"

    # CS (ipsilat)
    scil_tractogram_filter_by_roi \
      "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv100.trk" \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
      --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz" 'either_end' 'include' \
      --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${mni_dir}/MNI/from_${nside}/VTTT_Controlat_EXC_Ventral_Brainstem.nii.gz" 'any' 'exclude' \
      --drawn_roi "${mni_dir}/MNI/midsagittal_plane.nii.gz" 'any' 'exclude' -f

    # dPSN (ipsilat)
    scil_tractogram_filter_by_roi \
      "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_from_spinal_track_npv1000.trk" \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk" \
      --drawn_roi "${mni_rois_dir}/${nsub}_${nside}_VPM_mni.nii.gz" 'any' 'include' \
      --drawn_roi "${group_final_dir}/all_${nside}_remaining_cp_density_mni.nii.gz" 'either_end' 'include' -f
  done

  echo "|------------- 3) Done -------------|"
  echo ""

  echo "|------------- 4) Cutting the Second-order bundles -------------|"
  for nside in left right; do
    if [ "$nside" == "left" ]; then
      contra_nside="right"
    else
      contra_nside="left"
    fi

    echo "|------------- 4.1) ${nside} - VTTT (only controlateral) - vPSN and OS/IS -------------|"
    scil_tractogram_cut_streamlines \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk" \
      --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_OSandIS_Cuts_labels_mni.nii.gz" \
      "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_VTTT_Controlat_OSandIS.trk" -f

    scil_tractogram_cut_streamlines \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk" \
      --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_VTTT_Controlat_vPSN_Cuts_labels_mni.nii.gz" \
      "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_VTTT_Controlat_vPSN.trk" -f

    echo "|------------- 4.2) ${nside} - DTTT (controlateral) - CS -------------|"
    scil_tractogram_cut_streamlines \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Controlat_CS.trk" \
      --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Controlat_CS_Cuts_labels_mni.nii.gz" \
      "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Controlat_CS.trk" -f

    echo "|------------- 4.3) ${nside} - DTTT (ipsilateral) - dPSN and CS -------------|"
    scil_tractogram_cut_streamlines \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
      --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_CS_Cuts_labels_mni.nii.gz" \
      "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" -f

    scil_tractogram_cut_streamlines \
      "${mni_tracking_dir_second_order}/filtered/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk" \
      --labels "${mni_rois_dir}/${nsub}_${nside}_second_order_DTTT_Ipsilat_dPSN_Cuts_labels_mni.nii.gz" \
      "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_dPSN.trk" -f
  done

  echo "|------------- 4) Done -------------|"
  echo ""

  echo "|------------- 5) Cleaning the Second-order bundles -------------|"
  for nside in left right; do
    scil_bundle_reject_outliers \
      "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
      "${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_DTTT_Ipsilat_CS.trk" \
      --alpha 0.30 -f

    for nbundle in DTTT_Ipsilat_dPSN DTTT_Controlat_CS VTTT_Controlat_OSandIS VTTT_Controlat_vPSN; do
      scil_bundle_reject_outliers \
        "${mni_tracking_dir_second_order}/cut/${nsub}_from_${nside}_${nbundle}.trk" \
        "${mni_tracking_dir_second_order}/final/${nsub}_from_${nside}_${nbundle}.trk" \
        --alpha 0.50 -f
    done
  done

  echo "|------------- 5) Done -------------|"
  echo ""
  echo "|------------- SECOND-ORDER FIBERS TRACTOGRAPHY FOR ${nsub} IS COMPLETED -------------|"
  echo ""
done
