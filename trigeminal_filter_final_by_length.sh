
# trigeminal_filter_final_by_length.sh -f path/to/final/dir/trks

# example run:
# trigeminal_filter_final_by_length.sh -f ~/Research/data_temp/re-move/RE-MOVE_MRI_protocol/RE-MOVE_001_processing/output_atlas/S1/mni_space/tractograms/final/

usage() { echo "$(basename $0) [-f path/to/final/dir/trks]" 1>&2; exit 1; }

while getopts "f:" args; do
    case "${args}" in
        f) f=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${f}" ] ; then
    usage
fi

trk_dir=${f}

echo "Final Folder with TRKs: " ${trk_dir}

for t in ${trk_dir}/*.trk
do
    echo $t
    
    if [[ $t == *"mesencephalic"* ]]; then
	echo "Length filtering $t - 37 mm < length < 82 mm"

	scil_tractogram_filter_by_length $t ${t/.trk/_length_37-82mm.trk} \
					 --minL 37 \
					 --maxL 82 \
					 --no_empty \
					 --display_counts 	
    fi

    if [[ $t == *"remaining_cp"* ]]; then
	echo "Length filtering $t - 9 mm < length < 38 mm"

	scil_tractogram_filter_by_length $t ${t/.trk/_length_9-38mm.trk} \
					 --minL 9 \
					 --maxL 38 \
					 --no_empty \
					 --display_counts 	
    fi
    
    if [[ $t == *"spinal"* ]]; then
	echo "Length filtering $t - 36 mm < length < 70 mm"
	
	scil_tractogram_filter_by_length $t ${t/.trk/_length_36-70mm.trk} \
					 --minL 36 \
					 --maxL 70 \
					 --no_empty \
					 --display_counts 	
    fi

done
