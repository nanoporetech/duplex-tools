#!/bin/bash
set -eo pipefail

usage="$(basename "$0") [-h] -i <pod5_dir> -D <dorado_executable> -m <small_model> -M <large_model>

Basecall with dorado. Do pairing for both stages.
    -h  show this help text.
    -D  dorado executable
    -i  input pod5 directory
    -M  large model (typically sup, used for simplex and stereo calling)
    -m  small model (typically fast, used for fast calling)
    -n  dry run. Print out all commands which would have been executed. Do not run them.
    -d  devices. (e.g. cuda:1,2,3,4).
    "

dflag=false
pflag=false
mflag=false
Mflag=false
dryflag=false
deviceflag=false
OUTPUT_DIR=dorado_calls

while getopts ':hD:o:i:M:m:nd:' option; do
  case "$option" in
    h  ) echo "$usage" >&2; exit;;
    D  ) DORADO=$OPTARG; dflag=true;;
    o  ) OUTPUT_DIR=$OPTARG; oflag=true;;
    i  ) POD5=$OPTARG; pflag=true;;
    M  ) MODEL_LARGE=$OPTARG; Mflag=true;;
    m  ) MODEL_SMALL=$OPTARG; mflag=true;;
    n  ) dryflag=true;;
    d  ) DEVICE_STRING="--device $OPTARG"; deviceflag=true;;
    \? ) echo "Invalid option: -${OPTARG}." >&2; exit 1;;
    :  ) echo "Option -$OPTARG requires an argument." >&2; exit 1;;
  esac
done
shift $(($OPTIND - 1))

PAIRS_DIR=${OUTPUT_DIR}/pairs_stage1

SIMPLEX_UNMAPPED_DIR=${OUTPUT_DIR}/simplex_stage1
SIMPLEX_UNMAPPED=${SIMPLEX_UNMAPPED_DIR}/simplex_unmapped.sam

SIMPLEX_UNMAPPED_MOVES_DIR=${OUTPUT_DIR}/simplex_stage2
SIMPLEX_UNMAPPED_MOVES=${SIMPLEX_UNMAPPED_MOVES_DIR}/simplex_unmapped_moves.sam

DUPLEX_DIR=${OUTPUT_DIR}/duplex
DUPLEX_DIR2=${OUTPUT_DIR}/duplex_stage2

SPLITDUPLEX_DIR=${OUTPUT_DIR}/pairs_stage2

warn="\033[0;33m"
good="\033[0;32m"
off="\033[0m"

if ! "$dflag"; then
    echo "$usage" >&2;
    echo -e "${warn}Need to specify dorado executable (with -D)${off}";
    exit 1;
fi

if ! "$pflag"; then
    echo "$usage" >&2;
    echo -e "${warn}Need to specify pod5 input (with -i)${off}";
    exit 1;
fi

if ! "$mflag"; then
    echo "$usage" >&2;
    echo -e "${warn}Need to specify small model (with -m)${off}";
    echo -e "${warn}Example: -m dna_r10.4.1_e8.2_400bps_fast@v4.0.0${off}";
    exit 1;
fi

if ! "$Mflag"; then
    echo "$usage" >&2;
    echo -e "${warn}Need to specify large model (with -M)${off}";
    echo -e "${warn}Example: -M dna_r10.4.1_e8.2_400bps_sup@v4.0.0${off}";
    exit 1;
fi

if "$dryflag"; then
    echo "$usage" >&2;
    echo -e "${warn}Will dry-run. Remove -n to execute the commands.${off}";
fi


if ! command -v duplex_tools pair > /dev/null; then
    echo -e "${warn}Need duplex_tools on path. ${off}";
    echo -e "${warn}To install: pip install duplex-tools.${off}";
    exit 1;
fi

if ! $dryflag; then
    mkdir -p ${DUPLEX_DIR}
    mkdir -p ${SIMPLEX_UNMAPPED_MOVES_DIR}
    mkdir -p ${DUPLEX_DIR2}
    mkdir -p ${SIMPLEX_UNMAPPED_DIR}
    mkdir -p ${SPLITDUPLEX_DIR}
fi

# Stage 1 basecalling
stage1_simplex="$DORADO basecaller $MODEL_LARGE $POD5 $DEVICE_STRING\
 > $SIMPLEX_UNMAPPED # Stage 1 simplex"

# Stage 1 pairs
stage1_pairing="duplex_tools pair ${SIMPLEX_UNMAPPED}\
 --output_dir ${PAIRS_DIR} # Stage 1 pairs"

# Stage 1 stereo
stage1_stereo="$DORADO duplex $MODEL_LARGE $POD5 $DEVICE_STRING\
 --pairs ${PAIRS_DIR}/pair_ids_filtered.txt\
 > ${DUPLEX_DIR}/duplex.sam # Stage 1 stereo"

# Stage 2 simplex
stage2_simplex="$DORADO basecaller\
 $MODEL_SMALL ${POD5} $DEVICE_STRING --emit-moves\
 > ${SIMPLEX_UNMAPPED_MOVES} # Stage 2 simplex"

# Stage 2 pairing
stage2_pairing="duplex_tools split_pairs --force_overwrite\
 ${SIMPLEX_UNMAPPED_MOVES} ${POD5} ${SPLITDUPLEX_DIR} # Stage 2 pairs"

cat="cat ${SPLITDUPLEX_DIR}/*_pair_ids.txt > ${SPLITDUPLEX_DIR}/pair_ids_all.txt"

# Stage 2 stereo
stage2_stereo="$DORADO duplex ${MODEL_LARGE} ${SPLITDUPLEX_DIR} $DEVICE_STRING\
 --pairs ${SPLITDUPLEX_DIR}/pair_ids_all.txt\
 > ${DUPLEX_DIR2}/duplex_splitduplex.sam # Stage 2 stereo"

printf "$stage1_simplex\n\n";
printf "$stage1_pairing\n\n";
printf "$stage1_stereo\n\n";
printf "$stage2_simplex\n\n";
printf "$stage2_pairing\n\n";
printf "$cat\n\n";
printf "$stage2_stereo\n\n";

if ! $dryflag; then
    echo -e "${good}Starting Stage1 - Simplex${off}";
    eval $stage1_simplex;
    echo -e "${good}Starting Stage1 - Pairing${off}";
    eval $stage1_pairing;
    echo -e "${good}Starting Stage1 - Stereo${off}";
    eval $stage1_stereo;

    echo -e "${good}Starting Stage2 - Simplex${off}";
    eval $stage2_simplex;
    echo -e "${good}Starting Stage2 - Pairing${off}";
    eval $stage2_pairing;
    eval $cat;
    echo -e "${good}Starting Stage2 - Stereo${off}";
    eval $stage2_stereo;
fi
