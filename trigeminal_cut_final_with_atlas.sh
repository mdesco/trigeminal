# example run:
# trigeminal_cut_final_with_atlassh -f ~/Research/data_temp/re-move/RE-MOVE_MRI_protocol/RE-MOVE_001_processing/output_atlas/S1/mni_space/tractograms/final/ -a ~/Research/Source/trigeminal/atlas/bundles_mask/

usage() { echo "$(basename $0) [-f path/to/final/dir/trks] [-a path/to/atlas/trks]" 1>&2; exit 1; }

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
    if [[ $t == *"mesencephalic"* ]]; then
	if [[ $t == *"left"* ]]; then
	    echo $t
	    scil_tractogram_cut_streamlines $t ${t/.trk/_cut.trk} \
					    --mask ${atlas_dir}/left_mesencephalic_mask_dil3.nii.gz \
					    --trim_endpoints --processes 8 -f
	else
	    echo $t
	    scil_tractogram_cut_streamlines $t ${t/.trk/_cut.trk} \
					    --mask ${atlas_dir}/right_mesencephalic_mask_dil3.nii.gz \
					    --trim_endpoints --processes 8 -f
	fi
    fi
    if [[ $t == *"remaining_cp"* ]]; then
	if [[ $t == *"left"* ]]; then
	    echo $t
	    scil_tractogram_cut_streamlines $t ${t/.trk/_cut.trk} \
					    --mask ${atlas_dir}/left_remaining_cp_mask_dil3.nii.gz \
					    --trim_endpoints --processes 8 -f
	else
	    echo $t
	    scil_tractogram_cut_streamlines $t ${t/.trk/_cut.trk} \
					    --mask ${atlas_dir}/right_remaining_cp_mask_dil3.nii.gz \
					    --trim_endpoints --processes 8 -f
	fi
    fi
    if [[ $t == *"spinal"* ]]; then
	if [[ $t == *"left"* ]]; then
	    echo $t
	    scil_tractogram_cut_streamlines $t ${t/.trk/_cut.trk} \
					    --mask ${atlas_dir}/left_spinal_mask_dil3.nii.gz \
					    --trim_endpoints --processes 8 -f
	else
	    echo $t
	    scil_tractogram_cut_streamlines $t ${t/.trk/_cut.trk} \
					    --mask ${atlas_dir}/left_spinal_mask_dil3.nii.gz \
					    --trim_endpoints --processes 8 -f
	fi
    fi
done
