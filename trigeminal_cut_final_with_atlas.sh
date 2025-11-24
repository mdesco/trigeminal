
# example run:
# trigeminal_cut_final_with_atlassh -f ~/Research/data_temp/re-move/RE-MOVE_MRI_protocol/RE-MOVE_001_processing/output_atlas/S1/mni_space/tractograms/final/ -a ~/Research/Source/trigeminal/atlas/bundles_mask/

while getopts "f:a:" args; do
    case "${args}" in
        f) f=${OPTARG};;
	a) a=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${f}" ] || [ -z "${a}" ]; then
    usage
fi

trk_dir=${f}
atlas_dir=${a}

echo "Final Folder with TRKs: " ${trk_dir}
echo "Atlas Folder with TRKs: " ${atlas_dir}

for t in ${trk_dir}/*.trk
do
    echo $t
    
    if [[ $t == *"mesencephalic"* ]]; then
	for nside in left right
	do
	    scil_tractogram_cut_streamlines final_${nside}_mesencephalic.trk final_${nside}_mesencephalic_cut.trk \
					    --mask ${atlas_dir}/${nside}_mesencephalic_mask_dil3.nii.gz \
					    --trim_endpoints
	done
    fi

    if [[ $t == *"remaining_cp"* ]]; then
	for nside in left right
	do
	    scil_tractogram_cut_streamlines final_${nside}_remaining_cp.trk final_${nside}_remaining_cp_cut.trk \
					    --mask ${atlas_dir}/${nside}_remaining_cp_mask_dil3.nii.gz \
					    --trim_endpoints
	done

    fi
    
    if [[ $t == *"spinal"* ]]; then
	for nside in left right
	do
	    scil_tractogram_cut_streamlines final_${nside}_spinal.trk final_${nside}_spinal_cut.trk \
					    --mask ${atlas_dir}/${nside}_spinal_mask_dil3.nii.gz \
					    --trim_endpoints
	done
    fi
done
