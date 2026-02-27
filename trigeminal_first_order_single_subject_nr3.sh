#!/bin/bash
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Samir Akeb (2022-2023)
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Arnaud Bore (2023-2024)
# MODIFIED VERSION by Zineb - Oct 2025
# Added multiple step/theta runs and kept previous outputs

usage() { echo "$(basename $0) [-s path/to/subjects] [-m path/to/mni] [-o output_dir] [-t nb_threads] [-g true]" 1>&2; exit 1; }

while getopts "s:m:o:t:g:" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        o) o=${OPTARG};;
        t) t=${OPTARG};;
        g) g=${OPTARG};;
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

fa_threshold=0.20
npv_first_order=20000
opposite_side=leftright

gpu=""
if [ ! -z "${g}" ]; then
    gpu="--use_gpu"
fi

# 🔹 ADD parameter lists
step_sizes=(0.1 0.5 1.0)
thetas=(20 30 40)

echo "Folder subjects: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Output folder: " ${out_dir}
echo "Use GPU: " ${gpu}
echo "Number of threads" ${nb_thread}

for nsub in ${subject_dir}/*/
do
    nsub=`basename "$nsub"`

    # 🔹 Do NOT remove previous results anymore
    mkdir -p ${out_dir}/${nsub}/orig_space/{rois,tracking_first_order,transfo}
    mkdir -p ${out_dir}/${nsub}/orig_space/tracking_first_order/orig
    mkdir -p ${out_dir}/${nsub}/mni_space/{rois,tracking_first_order}
    mkdir -p ${out_dir}/${nsub}/mni_space/tracking_first_order/{orig,filtered,segmented,final}

    echo ""
    echo "|------------- PROCESSING FIRST ORDER TGN TRACTOGRAPHY FOR ${nsub} -------------|"
    echo ""

    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${nb_thread}"

    # 🔹 LOOP over all parameter combinations
    for step_size in "${step_sizes[@]}"; do
        for theta in "${thetas[@]}"; do
            echo ""
            echo "=== Running with step_size=${step_size}, theta=${theta} ==="

            # 🟢 Create a unique subfolder for each combination
            combo_tag="step${step_size}_theta${theta}"
            combo_dir=${out_dir}/${nsub}/mni_space/tracking_first_order/${combo_tag}
            mkdir -p ${combo_dir}/{orig,filtered,segmented,final}

            echo "|------------- 1) Registrations (if not done) -------------|"
            if [ ! -f "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" ]; then
                antsRegistrationSyN.sh \
                  -d 3 \
                  -f "${subject_dir}/${nsub}/tractoflow/${nsub}__t1_warped.nii.gz" \
                  -m "${mni_dir}/MNI/mni_masked.nii.gz" \
                  -t s \
                  -o "${out_dir}/${nsub}/orig_space/transfo/2orig_"
            fi

            ## The rest of your steps stay identical
            ## We only modify where outputs go and add combo tags
            ## To keep this short, I only show the key modified section (tracking)

            echo "|------------- 4) [ORIG-SPACE] Generate local tractography -------------|"
            for nside in left right
            do
                out_trk=${combo_dir}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk

                scil_tracking_local ${subject_dir}/${nsub}/tractoflow/${nsub}__fodf.nii.gz \
                    ${out_dir}/${nsub}/orig_space/rois/${nsub}_cp_${nside}_orig.nii.gz \
                    ${out_dir}/${nsub}/orig_space/rois/${nsub}_wm_mask_${fa_threshold}_orig.nii.gz \
                    ${out_trk} \
                    --npv $npv_first_order \
                    --step ${step_size} \
                    --theta ${theta} \
                    ${gpu} -v -f
            done

            echo "|------------- 5) Register Tracking in MNI space -------------|"
            for nside in left right
            do
                scil_tractogram_apply_transform \
                    ${combo_dir}/orig/${nsub}_${nside}_from_cp_${combo_tag}.trk \
                    ${mni_dir}/MNI/mni_masked.nii.gz \
                    ${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
                    ${combo_dir}/orig/${nsub}_${nside}_from_cp_${combo_tag}_mni.trk \
                    --in_deformation ${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
                    --remove_invalid \
                    --reverse_operation -f
            done

            echo "|------------- 9) Concatenate all runs for ${nsub} -------------|"
            for nside in left right
            do
                scil_tractogram_math concatenate ${out_dir}/${nsub}/mni_space/tracking_first_order/${nsub}_${nside}_from_cp_all.trk \
                    ${combo_dir}/orig/${nsub}_${nside}_from_cp_${combo_tag}_mni.trk -f
            done

            echo "=== Done for step=${step_size}, theta=${theta} ==="
        done
    done
done

