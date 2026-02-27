#!/usr/bin/env bash
# TRIGEMINAL SYSTEM TRACTOGRAPHY - Second-order (ensemble-ready, first-order-style CLI)
#
# Features:
# - Compatible with your first-order output structure:
#     <out_dir>/<sub>/mni_space/tracking_first_order/final_merged/final/*.trk
# - Builds group density masks from first-order (spinal + remaining_cp)
# - Second-order tracking runs in ORIG space over (step, theta) grid (ensemble)
# - Per-role seed budgets configurable (TOTAL npv), auto-split across combos
# - Merges all combos back into original naming convention:
#     *_from_spinal_track_npv1000.trk
#     *_from_spinal_track_npv100.trk
#     *_from_thalamus_npv500.trk
#
# Requirements:
# - scilpy CLI tools available in PATH
# - ANTs (antsApplyTransforms) available
# - First-order outputs already created (rois + transfo + first-order final bundles)

set -euo pipefail

# -------------------------------
# Defaults (override via flags)
# -------------------------------
STEP_DEFAULTS=(0.1 0.5 1.0)
THETA_DEFAULTS=(20 30 40)

NPV_SPINAL_LONG=1000
NPV_SPINAL_SHORT=100
NPV_THALAMUS=500

THREADS=""
USE_GPU="false"

# -------------------------------
# Usage
# -------------------------------
usage() {
  cat <<EOF 1>&2
Usage:
  $(basename "$0") -s <subjects_parent_or_single_subject_dir> -m <ROIs_clean_dir> -o <out_dir>
                   [-t threads] [-g true|false]
                   [-p step_size] [-e theta_deg]
                   [--npv_spinal_long N] [--npv_spinal_short N] [--npv_thalamus N]

Tracking params:
  -p step_size    If set with -e, runs SINGLE combo (no ensemble)
  -e theta_deg    If set with -p, runs SINGLE combo (no ensemble)
  If -p/-e not provided, ensemble is used:
    step  = ${STEP_DEFAULTS[*]}
    theta = ${THETA_DEFAULTS[*]}

Second-order seed budgets (TOTAL per role; auto-split across combos):
  --npv_spinal_long N   default: ${NPV_SPINAL_LONG}
  --npv_spinal_short N  default: ${NPV_SPINAL_SHORT}
  --npv_thalamus N      default: ${NPV_THALAMUS}

Notes:
  - First-order final expected at: <out_dir>/*/mni_space/tracking_first_order/final_merged/final/
  - This script expects first-order transfo in: <out_dir>/<sub>/orig_space/transfo/2orig_0GenericAffine.mat and 2orig_1Warp.nii.gz
EOF
  exit 1
}

die() { echo "ERROR: $*" 1>&2; exit 1; }
need_file() { [[ -f "$1" ]] || die "Missing file: $1"; }
need_dir()  { [[ -d "$1" ]] || die "Missing directory: $1"; }

ceil_div() { echo $(( ( $1 + $2 - 1 ) / $2 )); }

# -------------------------------
# Parse args (long opts first)
# -------------------------------
s=""; m=""; o=""
step_size=""; theta=""

args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --npv_spinal_long)  NPV_SPINAL_LONG="${2:?}"; shift 2 ;;
    --npv_spinal_short) NPV_SPINAL_SHORT="${2:?}"; shift 2 ;;
    --npv_thalamus)     NPV_THALAMUS="${2:?}"; shift 2 ;;
    -h|--help) usage ;;
    *) args+=("$1"); shift ;;
  esac
done
set -- "${args[@]}"

while getopts "s:m:o:t:g:p:e:" opt; do
  case "${opt}" in
    s) s="${OPTARG}" ;;
    m) m="${OPTARG}" ;;
    o) o="${OPTARG}" ;;
    t) THREADS="${OPTARG}" ;;
    g) USE_GPU="${OPTARG}" ;;
    p) step_size="${OPTARG}" ;;
    e) theta="${OPTARG}" ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

[[ -n "${s}" && -n "${m}" && -n "${o}" ]] || usage

subject_dir="${s}"
mni_dir="${m}"
out_dir="${o}"

need_dir "${subject_dir}"
need_dir "${mni_dir}"
mkdir -p "${out_dir}"

# -------------------------------
# Threads / GPU
# -------------------------------
if [[ -n "${THREADS}" ]]; then
  export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${THREADS}"
fi

gpu_flag=""
if [[ "${USE_GPU}" == "true" ]]; then
  gpu_flag="--use_gpu"
fi

# -------------------------------
# Ensemble grid
# -------------------------------
if [[ -n "${step_size}" || -n "${theta}" ]]; then
  [[ -n "${step_size}" && -n "${theta}" ]] || die "If using -p or -e, you must provide BOTH -p <step> and -e <theta>."
  step_list=("${step_size}")
  theta_list=("${theta}")
else
  step_list=("${STEP_DEFAULTS[@]}")
  theta_list=("${THETA_DEFAULTS[@]}")
fi

n_combos=$(( ${#step_list[@]} * ${#theta_list[@]} ))
(( n_combos > 0 )) || die "Internal error: n_combos=0"

npv_spinal_long_per_run="$(ceil_div "${NPV_SPINAL_LONG}" "${n_combos}")"
npv_spinal_short_per_run="$(ceil_div "${NPV_SPINAL_SHORT}" "${n_combos}")"
npv_thalamus_per_run="$(ceil_div "${NPV_THALAMUS}" "${n_combos}")"

# -------------------------------
# First-order final folder (group masks stored here)
# -------------------------------
first_order_final_rel="mni_space/tracking_first_order/final_merged/final"
group_final_dir="${out_dir}/${first_order_final_rel}"
mkdir -p "${group_final_dir}"

echo "Folder subjects: ${subject_dir}"
echo "Folder MNI:      ${mni_dir}"
echo "Output folder:   ${out_dir}"
echo "Threads:         ${THREADS:-'(default)'}"
echo "GPU:             ${gpu_flag:-'(cpu)'}"
echo "First-order final expected at: ${out_dir}/*/${first_order_final_rel}/"
echo "Tracking grid:   steps=${step_list[*]}  thetas=${theta_list[*]} (n_combos=${n_combos})"
echo "NPV totals:      spinal_long=${NPV_SPINAL_LONG}, spinal_short=${NPV_SPINAL_SHORT}, thalamus=${NPV_THALAMUS}"
echo "NPV per-run:     spinal_long=${npv_spinal_long_per_run}, spinal_short=${npv_spinal_short_per_run}, thalamus=${npv_thalamus_per_run}"
echo ""

# -------------------------------
# Detect single-subject vs parent folder
# -------------------------------
subjects=()
subjects_parent=""

if [[ -d "${subject_dir}/tractoflow" && -d "${subject_dir}/freesurfer" ]]; then
  subjects=( "$(basename "${subject_dir}")" )
  subjects_parent="$(dirname "${subject_dir}")"
else
  subjects_parent="${subject_dir}"
  shopt -s nullglob
  for d in "${subjects_parent}"/*/ ; do
    subjects+=( "$(basename "$d")" )
  done
  shopt -u nullglob
fi

(( ${#subjects[@]} > 0 )) || die "No subjects found under -s ${subject_dir}"

# -------------------------------
# 0) Build group density masks from first-order bundles (spinal + remaining_cp)
# -------------------------------
echo "|------------- 0) Building group density masks from first-order bundles -------------|"

for nside in left right; do
  shopt -s nullglob
  spinal_inputs=( "${out_dir}"/*/"${first_order_final_rel}"/*_"${nside}"_spinal.trk )
  rcp_inputs=( "${out_dir}"/*/"${first_order_final_rel}"/*_"${nside}"_remaining_cp.trk )
  shopt -u nullglob

  (( ${#spinal_inputs[@]} > 0 )) || die "No spinal tracts found for side=${nside} at ${out_dir}/*/${first_order_final_rel}/*_${nside}_spinal.trk"
  (( ${#rcp_inputs[@]} > 0 ))    || die "No remaining_cp tracts found for side=${nside} at ${out_dir}/*/${first_order_final_rel}/*_${nside}_remaining_cp.trk"

  echo "Union (${nside}) spinal: ${#spinal_inputs[@]} files -> ${group_final_dir}/all_${nside}_spinal.trk"
  scil_tractogram_math union "${spinal_inputs[@]}" "${group_final_dir}/all_${nside}_spinal.trk" -f

  echo "Density (${nside}) spinal seed"
  scil_tractogram_compute_density_map \
    "${group_final_dir}/all_${nside}_spinal.trk" \
    "${group_final_dir}/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
    --binary -f

  echo "Union (${nside}) remaining_cp: ${#rcp_inputs[@]} files -> ${group_final_dir}/all_${nside}_remaining_cp.trk"
  scil_tractogram_math union "${rcp_inputs[@]}" "${group_final_dir}/all_${nside}_remaining_cp.trk" -f

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
  need_dir "${subj_path}"

  orig_rois_dir="${out_dir}/${nsub}/orig_space/rois"
  mni_rois_dir="${out_dir}/${nsub}/mni_space/rois"
  orig_tracking_dir="${out_dir}/${nsub}/orig_space/tracking_second_order"
  mni_tracking_dir_second_order="${out_dir}/${nsub}/mni_space/tracking_second_order"

  mkdir -p "${out_dir}/${nsub}/orig_space/"{rois,tracking_second_order,transfo}
  mkdir -p "${out_dir}/${nsub}/orig_space/tracking_second_order/"{orig,trials}
  mkdir -p "${out_dir}/${nsub}/mni_space/"{rois,tracking_second_order}
  mkdir -p "${out_dir}/${nsub}/mni_space/tracking_second_order/"{orig,filtered,final,cut}

  trials_dir="${orig_tracking_dir}/trials"
  merged_orig_dir="${orig_tracking_dir}/orig"

  echo "|------------- PROCESSING SECOND-ORDER FIBERS TRACTOGRAPHY -------------|"
  echo "|------------- SUBJECT ${nsub} -------------|"
  echo ""

  # Required first-order transforms
  need_file "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat"
  need_file "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz"

  # TractoFlow inputs
  need_file "${subj_path}/tractoflow/${nsub}__t1_warped.nii.gz"
  need_file "${subj_path}/tractoflow/${nsub}__fodf.nii.gz"

  # WM mask from first-order
  shopt -s nullglob
  wm_candidates=( "${orig_rois_dir}/${nsub}_wm_mask_"*"_orig.nii.gz" )
  shopt -u nullglob
  (( ${#wm_candidates[@]} > 0 )) || die "No WM mask found in ${orig_rois_dir}/${nsub}_wm_mask_*_orig.nii.gz"
  wm_mask="${wm_candidates[0]}"
  echo "Using WM mask: ${wm_mask}"
  echo ""

  # -------------------------------
  # 1) Seed mask in orig space + extract thalamus + tracking + merge + transform to MNI
  # -------------------------------
  echo "|------------- 1) Tracking (orig) + Transform to MNI -------------|"
  echo "|------------- 1.1) Seeding masks -> orig space -------------|"

  for nside in left right; do
    cp "${group_final_dir}/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      "${mni_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_mni.nii.gz"

    antsApplyTransforms \
      -d 3 \
      -i "${group_final_dir}/all_${nside}_spinal_density_second_order_seed_mni.nii.gz" \
      -r "${subj_path}/tractoflow/${nsub}__t1_warped.nii.gz" \
      -t "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
      -t "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
      -o "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz"
  done

  echo "|------------- 1.2) Extract thalamus (orig + mni) -------------|"
  need_file "${orig_rois_dir}/${nsub}_aparc.DKTatlas+aseg_orig.nii.gz"
  need_file "${mni_rois_dir}/${nsub}_aparc.DKTatlas+aseg_mni.nii.gz"

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

  echo "|------------- 1.3) Tracking ENSEMBLE (per combo) -------------|"
  for step_val in "${step_list[@]}"; do
    for theta_val in "${theta_list[@]}"; do
      combo_tag="step_${step_val}_theta_${theta_val}"
      combo_dir="${trials_dir}/${combo_tag}"
      mkdir -p "${combo_dir}"

      for nside in left right; do
        scil_tracking_local \
          "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" \
          "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz" \
          "${wm_mask}" \
          "${combo_dir}/${nsub}_${nside}_from_spinal_long_${combo_tag}.trk" \
          --npv "${npv_spinal_long_per_run}" --step "${step_val}" --theta "${theta_val}" -v -f ${gpu_flag}

        scil_tracking_local \
          "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" \
          "${orig_rois_dir}/${nsub}_${nside}_spinal_density_second_order_seed_orig.nii.gz" \
          "${wm_mask}" \
          "${combo_dir}/${nsub}_${nside}_from_spinal_short_${combo_tag}.trk" \
          --npv "${npv_spinal_short_per_run}" --step "${step_val}" --theta "${theta_val}" -v -f ${gpu_flag}

        scil_tracking_local \
          "${subj_path}/tractoflow/${nsub}__fodf.nii.gz" \
          "${orig_rois_dir}/${nsub}_${nside}_thalamus_orig.nii.gz" \
          "${wm_mask}" \
          "${combo_dir}/${nsub}_${nside}_from_thalamus_${combo_tag}.trk" \
          --npv "${npv_thalamus_per_run}" --step "${step_val}" --theta "${theta_val}" -v -f ${gpu_flag}
      done
    done
  done

  echo "|------------- 1.4) Merge combos -> original naming -------------|"
  mkdir -p "${merged_orig_dir}"

  for nside in left right; do
    files=()
    for step_val in "${step_list[@]}"; do
      for theta_val in "${theta_list[@]}"; do
        combo_tag="step_${step_val}_theta_${theta_val}"
        f="${trials_dir}/${combo_tag}/${nsub}_${nside}_from_spinal_long_${combo_tag}.trk"
        [[ -f "$f" ]] && files+=("$f")
      done
    done
    (( ${#files[@]} > 0 )) || die "No spinal_long trial files for ${nsub} ${nside}"
    scil_tractogram_math -f concatenate "${files[@]}" \
      "${merged_orig_dir}/${nsub}_${nside}_from_spinal_track_npv1000.trk"

    files=()
    for step_val in "${step_list[@]}"; do
      for theta_val in "${theta_list[@]}"; do
        combo_tag="step_${step_val}_theta_${theta_val}"
        f="${trials_dir}/${combo_tag}/${nsub}_${nside}_from_spinal_short_${combo_tag}.trk"
        [[ -f "$f" ]] && files+=("$f")
      done
    done
    (( ${#files[@]} > 0 )) || die "No spinal_short trial files for ${nsub} ${nside}"
    scil_tractogram_math -f concatenate "${files[@]}" \
      "${merged_orig_dir}/${nsub}_${nside}_from_spinal_track_npv100.trk"

    files=()
    for step_val in "${step_list[@]}"; do
      for theta_val in "${theta_list[@]}"; do
        combo_tag="step_${step_val}_theta_${theta_val}"
        f="${trials_dir}/${combo_tag}/${nsub}_${nside}_from_thalamus_${combo_tag}.trk"
        [[ -f "$f" ]] && files+=("$f")
      done
    done
    (( ${#files[@]} > 0 )) || die "No thalamus trial files for ${nsub} ${nside}"
    scil_tractogram_math -f concatenate "${files[@]}" \
      "${merged_orig_dir}/${nsub}_${nside}_from_thalamus_npv500.trk"
  done

  echo "|------------- 1.5) Apply transform -> MNI space -------------|"
  need_file "${mni_dir}/MNI/mni_masked.nii.gz"

  for ntracking in from_thalamus_npv500.trk from_spinal_track_npv100.trk from_spinal_track_npv1000.trk; do
    for nside in left right; do
      scil_tractogram_apply_transform \
        "${merged_orig_dir}/${nsub}_${nside}_${ntracking}" \
        "${mni_dir}/MNI/mni_masked.nii.gz" \
        "${out_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat" \
        "${mni_tracking_dir_second_order}/orig/${nsub}_${nside}_${ntracking}" \
        --in_deformation "${out_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz" \
        --remove_invalid --reverse_operation -f
    done
  done

  echo "|------------- 1) Done -------------|"
  echo ""

  # -------------------------------
  # 2) Generate ROIs  (original)
  # -------------------------------
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
    if [[ "$nside" == "left" ]]; then
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

  # -------------------------------
  # 3) Generate Second-order bundles (original)
  # -------------------------------
  echo "|------------- 3) Generate Second-order bundles for trigeminal system -------------|"
  for nside in left right; do
    if [[ "$nside" == "left" ]]; then
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
    scil_tractogram_filter_by HOLDERS_by_roi \
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

  # -------------------------------
  # 4) Cutting (original)
  # -------------------------------
  echo "|------------- 4) Cutting the Second-order bundles -------------|"
  for nside in left right; do
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

  # -------------------------------
  # 5) Cleaning (original)
  # -------------------------------
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
