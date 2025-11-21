#!/bin/bash

#
# This version of the script performs a Bundle-Specific Tractography (BST)
# based on the a priori of the trigeminal_nerve atlas. That is, the atlas
# is used as prior to seed aggressively and define masks to help tractography.
# The FA threshold is lowered to 0.01 to track everywhere. 

# Input structure
#
#    [input]
#    ├── sub-01
#    │   ├── freesurfer
#    │   │   └─── aparc.DKTatlas+aseg.mgz
#    │   ├── tractoflow    
#    │   │   └─── sub-01__fa.nii.gz
#    │   │   └─── sub-01__fodf.nii.gz
#    │   │   └─── sub-01__t1_warped.nii.gz

#
# Example of running the script
# trigeminal_apply_atlas.sh -s input/S1/ -m ~/Research/Source/trigeminal/ROIs_clean/ -a ~/Research/Source/trigeminal/atlas/ -o output_atlas/S1/ -t 8 -g true


usage() { echo "$(basename $0) [-s path/to/subject] [-m path/to/trigeminal/ROIs_clean] [-a path/to/trigeminal/atlas] [-o output_dir] [-t nb_threads] [-g] (if you have a gpu)" 1>&2; exit 1; }

while getopts "s:m:a:o:t:g:" args; do
    case "${args}" in
        s) s=${OPTARG};;
        m) m=${OPTARG};;
        a) a=${OPTARG};;
        o) o=${OPTARG};;
        t) t=${OPTARG};;
        g) g=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${s}" ] || [ -z "${m}" ] || [ -z "${o}" ] || [ -z "${a}" ]; then
    usage
fi

gpu=""
if [ -n "${g}" ]; then
    gpu="--use_gpu"
fi

subject_dir=${s}
atlas_dir=${a}
mni_dir=${m}
out_dir=${o}
nb_thread=${t}

echo "Folder subject: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Atlas MNI: " ${atlas_dir}
echo "Output folder: " ${out_dir}
echo "GPU: " ${gpu}

npv="400" #TODO: add an npv option. No need for npy=100. Sub-optimal.
          #TODO: once Nasrine's PR is in, we mimick her ensemble tractography with 9 local tracking and npv/9
opposite_side=leftright
fa_threshold=0.01

mkdir -p ${out_dir}/orig_space/{bundles_mask,transfo}
mkdir -p ${out_dir}/{orig_space,mni_space}/rois
mkdir -p ${out_dir}/mni_space/tractograms/{orig,filtered,length,segmented,final}
mkdir -p ${out_dir}/orig_space/tractograms/{orig,final}

orig_rois_dir=${out_dir}/orig_space/rois
mni_rois_dir=${out_dir}/mni_space/rois
orig_tracking_dir=${out_dir}/orig_space/tractograms/
mni_tracking_dir=${out_dir}/mni_space/tractograms/

echo "|------------- 1) Register mni_masked in orig space -------------|"
antsRegistrationSyN.sh -d 3 \
		       -f ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
		       -m ${mni_dir}/MNI/mni_masked.nii.gz \
		       -t s \
		       -o ${out_dir}/orig_space/transfo/2orig_ \
		       -y 1 \
		       -n ${nb_thread} >> ${out_dir}/orig_space/transfo/2orig_log.txt

echo "|------------- 1.1) Reshape aparc.DKTatlas+aseg.mgz orig space -------------|"
## [ORIG-SPACE] Reshape aparc.DKTatlas+aseg.mgz
scil_volume_reslice_to_reference ${subject_dir}/freesurfer/aparc.DKTatlas+aseg.mgz \
				 ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
				 ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz \
				 --interpolation nearest --keep_dtype -f
echo "|------------- 1) Done -------------|"
echo ""

echo "|------------- 2) Generate exclusions and inclusions ROI -------------|"
## [ORIG-SPACE] Exclusions ROI labels 
Left_Cerebral_Cortex=($(seq 1000 1035)) # warning for missing labels is normal
Right_Cerebral_Cortex=($(seq 2000 2035))
Left_Cerebral_WM=(2)
Right_Cerebral_WM=(41)
Cerebellum_Cortex=(8 47)
Right_Cerebellum_WM=(46)
Left_Cerebellum_WM=(7)
Any_Exclusion_ROI=(${Left_Cerebral_Cortex[*]} ${Right_Cerebral_Cortex[*]} ${Left_Cerebral_WM[*]} ${Right_Cerebral_WM[*]} ${Cerebellum_Cortex[*]}) #Generate bilateral exclusion array 

## Exclusions ROIs masks
## Generate bilateral exclusion ROI
echo "|------------- 2.1) any_exclusion_roi_orig -------------|"
scil_labels_combine ${orig_rois_dir}/any_exclusion_roi_orig.nii.gz \
		    --volume_ids ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz ${Any_Exclusion_ROI[*]} \
		    --merge_groups -f

echo "|------------- 2.2) cerebellum_wm_right_orig -------------|"
## WM Cerebellum Right
scil_labels_combine ${orig_rois_dir}/right_cerebellum_wm_orig.nii.gz \
		    --volume_ids ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz ${Right_Cerebellum_WM[*]} \
		    --merge_groups -f

## WM Cerebellum Left
echo "|------------- 2.3) cerebellum_wm_left_orig -------------|"
scil_labels_combine ${orig_rois_dir}/left_cerebellum_wm_orig.nii.gz \
		    --volume_ids ${orig_rois_dir}/aparc.DKTatlas+aseg_orig.nii.gz ${Left_Cerebellum_WM[*]} \
		    --merge_groups -f

# WM mask
scil_volume_math lower_threshold ${subject_dir}/tractoflow/*__fa.nii.gz \
		 ${fa_threshold} \
		 ${orig_rois_dir}/wm_mask_${fa_threshold}_orig.nii.gz \
		 --data_type uint8 -f

# Register these ROI in MNI space
echo "|------------- 2.4) Registration in MNI space -------------|"
for roi in right_cerebellum_wm left_cerebellum_wm any_exclusion_roi
do
    antsApplyTransforms -d 3 \
			-i ${orig_rois_dir}/${roi}_orig.nii.gz \
			-r ${mni_dir}/MNI/mni_masked.nii.gz \
			-t [${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat, 1] \
			-t ${out_dir}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
			-o ${mni_rois_dir}/${roi}_mni.nii.gz \
			-n NearestNeighbor;
    
    scil_volume_math convert \
		     ${mni_rois_dir}/${roi}_mni.nii.gz \
		     ${mni_rois_dir}/${roi}_mni.nii.gz --data_type int16 -f
done
echo "|------------- 2) Done -------------|"
echo ""


echo "|------------- 3) Register atlas components  -------------|"
for mask_type in bundles_mask
do
    echo "|------------- 3.1) Component : Register atlas components from folder ${mask_type} in orig space -------------|"
    for component in ${atlas_dir}/${mask_type}/*;
    do
        atlas_component=`basename "$component"`
        echo "|------------- 3.2) Component : ${atlas_component} -------------|"
        antsApplyTransforms -d 3 \
			    -i ${component} \
			    -r ${subject_dir}/tractoflow/*__t1_warped.nii.gz \
			    -t ${out_dir}/orig_space/transfo/2orig_1Warp.nii.gz \
			    -t ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
			    -o ${out_dir}/orig_space/${mask_type}/${atlas_component} \
			    -n NearestNeighbor;
	
        scil_volume_math convert \
			 ${out_dir}/orig_space/${mask_type}/${atlas_component} \
			 ${out_dir}/orig_space/${mask_type}/${atlas_component} \
			 --data_type int16 -f;
    done
done
echo "|------------- 3) Done -------------|"
echo ""

# TODO: real BST
#     i) generate dilated endpoints mask
#     ii) TODI maps from tempalte
#     iii) e-fodf
#     iv) perform local tracking seeding from dilated union masks with endpoints, using e-fodf and a low FA threshold
echo "|------------- 4) Tracking from atlas component  -------------|"
for component in left_mesencephalic.nii.gz left_spinal.nii.gz left_remaining_cp.nii.gz right_mesencephalic.nii.gz right_spinal.nii.gz right_remaining_cp.nii.gz #${atlas_dir}/bundles_mask/*;
do
    atlas_component=$component

    # TODO: replace by the ensemble tractography strategy of Nasrin and trigeminal_first_order.sh
    echo "|------------- 4.1) Tracking from atlas component ${atlas_component} with npv=${npv} -------------|"
    scil_tracking_local ${subject_dir}/tractoflow/*__fodf.nii.gz \
			${out_dir}/orig_space/bundles_mask/${atlas_component} \
			${orig_rois_dir}/wm_mask_${fa_threshold}_orig.nii.gz \
			${out_dir}/orig_space/tractograms/orig__${atlas_component/.nii.gz/.trk} \
			--npv ${npv} \
			${gpu} -f;
	
    echo "|------------- 4.2) Register tracking to MNI space -------------|"
    scil_tractogram_apply_transform \
        ${out_dir}/orig_space/tractograms/orig__${atlas_component/.nii.gz/.trk} \
        ${mni_dir}/MNI/mni_masked.nii.gz \
        ${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
        ${mni_tracking_dir}/orig__${atlas_component/.nii.gz/.trk} \
        --in_deformation ${out_dir}/orig_space/transfo/2orig_1Warp.nii.gz \
        --remove_invalid \
        --reverse_operation -f
    
    # moving tractogram to the orig folder
    mv ${out_dir}/orig_space/tractograms/orig_*trk ${out_dir}/orig_space/tractograms/orig/
done
echo "|------------- 4) Done -------------|"
echo ""

echo "|------------- 5) Filtering & Segmenting -------------|"
echo "|------------- 5.1) Major filtering for mesencephalic, remaining cp, spinal -------------|"
# Main filter for mesencephalic, remaining cp, spinal
for atlas_component in mesencephalic spinal remaining_cp
do
    for nside in left right
    do
        echo "|------------- 5.1b) Major filtering for ${atlas_component} nside: ${nside} npv: ${npv} -------------|"
        scil_tractogram_filter_by_roi ${mni_tracking_dir}/orig__${nside}_${atlas_component}.trk \
				      ${mni_tracking_dir}/filtered_${nside}_${atlas_component}.trk \
				      --drawn_roi ${mni_rois_dir}/any_exclusion_roi_mni.nii.gz 'any' 'exclude' \
				      --drawn_roi ${mni_rois_dir}/${nside/${nside}/${opposite_side/${nside}/}}_cerebellum_wm_mni.nii.gz 'any' 'exclude' \
				      --drawn_roi ${mni_rois_dir}/${nside}_cerebellum_wm_mni.nii.gz 'either_end' 'exclude' \
				      --drawn_roi ${mni_dir}/MNI/midsagittal_plane.nii.gz 'any' 'exclude' \
				      --drawn_roi ${mni_dir}/MNI/cp_${nside}_bin.nii.gz any include -f -v # It was not needed when creating atlas since everything was tracked from cp
    done
done

# Filter mesencephalic
echo "|------------- 5.2) Segmentation for mesencephalic -------------|"
for nside in left right
do
    scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered_${nside}_mesencephalic.trk \
				  ${mni_tracking_dir}/segmented_${nside}_mesencephalic.trk  \
				  --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'include' \
				  --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
				  --drawn_roi ${mni_dir}/MNI/coronal_plane_for_mesencephalic.nii.gz 'any' 'include' \
				  -f

    # TODO: length threshold
    # 37 mm < length < 82 mm

    
    scil_bundle_reject_outliers \
        ${mni_tracking_dir}/segmented_${nside}_mesencephalic.trk \
        ${mni_tracking_dir}/final_${nside}_mesencephalic.trk -f
done

# Filter remaining_cp
echo "|------------- 5.3) Segmentation remaining cp -------------|"
for nside in left right
do
    scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered_${nside}_remaining_cp.trk \
				  ${mni_tracking_dir}/segmented_${nside}_remaining_cp.trk  \
				  --drawn_roi ${mni_dir}/MNI/lower_cut_brainstem.nii.gz 'any' 'exclude' \
				  --drawn_roi ${mni_dir}/MNI/upper_cut_brainstem.nii.gz 'any' 'exclude' \
				  --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
				  --bdo ${mni_dir}/MNI/sphere_exclusion_for_remaining_cp.bdo 'any' 'exclude' \
				  -f
    
    # TODO: length threshold
    # 9 mm < length < 38 mm
    
    # TODO: we apply the length and ROI rules of Nasrine in the future.
    scil_bundle_reject_outliers \
        ${mni_tracking_dir}/segmented_${nside}_remaining_cp.trk \
        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk  
    
    scil_tractogram_filter_by_orientation \
        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk  \
        ${mni_tracking_dir}/final_${nside}_remaining_cp.trk  \
        --max_z 7 --use_abs -f
done

# Filter spinal
echo "|------------- 5.4) Segmentation spinal -------------|"
for nside in left right
do
    scil_tractogram_filter_by_roi ${mni_tracking_dir}/filtered_${nside}_spinal.trk \
				  ${mni_tracking_dir}/segmented_${nside}_spinal.trk \
				  --drawn_roi ${mni_dir}/MNI/lower_cut_brainstem.nii.gz 'any' 'include' \
				  --drawn_roi ${mni_dir}/MNI/coronal_plane.nii.gz 'any' 'include' \
				  -f

    # TODO: length threshold
    # 36 mm < length < 70 mm
    
    scil_bundle_reject_outliers \
        ${mni_tracking_dir}/segmented_${nside}_spinal.trk \
        ${mni_tracking_dir}/final_${nside}_spinal.trk \
        -f
done


echo "|------------- 6 Transform back final MNI trk to orig space -------------|"
for atlas_component in mesencephalic spinal remaining_cp
do
    for nside in left right
    do
        scil_tractogram_apply_transform ${mni_tracking_dir}/final_${nside}_${atlas_component}.trk \
					${subject_dir}/tractoflow/*__t1_warped.nii.gz \
					${out_dir}/orig_space/transfo/2orig_0GenericAffine.mat \
					${orig_tracking_dir}/final/final_${nside}_${atlas_component}.trk \
					--inverse \
					--in_deformation ${out_dir}/orig_space/transfo/2orig_1InverseWarp.nii.gz \
					--remove_invalid -f
    done
done


# TODO: need better cleaning up
#  1) Length thresholds first - easy
#  2) Outlier, better alpha? - easy 
#  3) Recobundle with atlas? - later

# Last organization move of files in the proper output directories in mni_space
mv  ${mni_tracking_dir}/orig_*  ${mni_tracking_dir}/orig/
mv  ${mni_tracking_dir}/filtered_*  ${mni_tracking_dir}/filtered/
mv  ${mni_tracking_dir}/segmented_*  ${mni_tracking_dir}/segmented/
mv  ${mni_tracking_dir}/final_*  ${mni_tracking_dir}/final/


# TODO: Last organization move of files in the proper output directories in orig_space 

echo "|------------- Done -------------|"
echo ""


