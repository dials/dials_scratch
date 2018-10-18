#!/bin/bash

set -e

usage() {
  echo "Usage: index_stills.sh [--submit | --collate] [options] datablock.json <reflections> [phil options]..."
}
help() {
  usage
  cat <<HERE

Options:
  -h, --help  Display this help
  --submit    Submit to SGE instead of running locally
  --collate   Collate the results of submitting to remote cluster
  -j JOBS     Submit to JOBS separate jobs. Only valid with --submit
  -n PROCS    Run PROCS instances in parallel. If GNU parallel is not available
              then this will be ignored.
  --python X  Use X as the python interpreter (defaults: libtbx.python)
  --parallel Y  Use GNU parallel from a specific location
  --smp N     When submitting, request a specific SMP count. Default: PROCS
HERE
}

# What mode to run in. index | submit | collate
MODE=index
# The number of jobs to submit. Only works in submit mode
NJOBS=1
# The number of processses to run in parallel
NPROCS=1
# Which python interpreter to use when launching sub-jobs
PYEXE=$(command -v libtbx.python 2>/dev/null || true)
# Do we have GNU parallel available?
PARALLEL=$(command -v parallel 2>/dev/null || true)
# Default arguments to always pass to dials.index
# In this case, we want to always fix the parameters so that we can combine sensibly
DEFAULT_PHIL="goniometer.fix=all detector.fix=all beam.fix=all"
# The specific processor count to request from qsub
PE_SMP_COUNT=""

# Find the full path to *this* script - we may need to call ourself
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null && pwd )"
SOURCE=$DIR/$(basename $SOURCE)

# Do the bash parsing
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    --submit)
    MODE=submit
    shift # past argument
    ;;
    --collate)
    MODE=collate
    shift # past argument
    ;;
    --combine)
    MODE=collate
    shift # past argument
    ;;
    -h|--help)
    shift # past value
    help
    exit 0
    ;;
    -n)
    NPROCS="$2"
    shift # past argument
    shift # past value
    ;;
    -j)
    NJOBS="$2"
    shift # past argument
    shift # past value
    ;;
    --python)
    PYEXE="$2"
    shift # past argument
    shift # past value
    ;;
    --parallel)
    PARALLEL="$2"
    shift # past argument
    shift # past value
    ;;
    --smp)
    PE_SMP_COUNT="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# MUST have at least two positionals - both files
if [[ $# -lt 2 ]]; then
    usage
    exit 1
fi

trap ctrl_c INT

function ctrl_c() {
  echo "** User cancelled"
  exit 1
}


# Validate that we have parallel if more than one job
if [[ -z $PARALLEL && $NPROCS -gt 1 ]]; then
  echo "$(tput bold)Warning: NPROCS=$NPROCS but GNU parallel not found. Running sequentially.$(tput sgr0)"
  NPROCS=1
fi

# Check that we have a python interpreter
if [[ -z "$PYEXE" ]]; then
  echo "Error: Could not find libtbx.python"
  exit 1
fi

# Validate that we didn't set NJOBS for non-submission
if [[ $MODE != "submit" && $NJOBS -gt 1 ]]; then
  echo "Error: Not submitting but asked for NJOBS=$NJOBS. Did you mean -n instead of -j?"
  exit 1
fi

# Validate that we have qsub available if submit mode
if [[ $MODE == "submit" ]]; then
  if ! command -v qsub > /dev/null 2>&1; then
    echo "$(tput bold)Error: Could not find qsub$(tput sgr0)"
    exit 1
  fi
fi


# Print out the loaded configuration
echo "Mode:   $MODE"
if [[ $MODE == "submit" ]]; then
  echo "NJOBS:  $NJOBS"
fi
echo "NPROCS: $NPROCS"
echo "Python: $PYEXE"

# Try to read the datablock to extract info e.g. number of images
datablock=$1
reflections=$2
shift; shift
num_images=$($PYEXE -c "from dxtbx.datablock import DataBlockFactory; a = DataBlockFactory.from_serialized_format('$datablock', check_format=False); print(a[0].num_images())")
NWIDTH=$($PYEXE -c "import math; print(int(math.ceil(math.log10($num_images))))")

# If this failed, try to give some info
if [[ ! $? ]]; then
  echo "Could not read datablock: $num_images"
  usage
  exit 1
fi

# Check if we're running in SGE mode and use this to specify a range of jobs
if [[ -n "$SGE_TASK_ID" && "$SGE_TASK_ID" != "undefined" ]]; then
  # Running in SGE task mode
  echo "Running in SGE Batch mode"
  id_start=$SGE_TASK_ID
  id_end=$(($SGE_TASK_ID+$SGE_TASK_STEPSIZE-1))
  # Make sure we don't go too far over the end
  if [[ $id_end -gt $SGE_TASK_LAST ]]; then
    id_end=$SGE_TASK_LAST
  fi
  SGE_MODE=1
else
  # Use the information from the datablock
  id_start=1
  id_end=$num_images
  SGE_MODE=""
fi
echo "Image range [1,$id_end]"

# # Get the info out of the datablock
# datablock="$1"
# reflections="$2"
# shift
# shift

# echo "Datablock:   $datablock"
# echo "Reflections: $reflections"
# num_images=$(libtbx.python -c "from dxtbx.datablock import DataBlockFactory; a = DataBlockFactory.from_serialized_format('$datablock', check_format=False); print(a[0].num_images())")
# num_len=$(python -c "import math; print(int(math.ceil(math.log10($num_images))))")

# Collate the results into a single experiment file
collate() {
  echo "Collating results"
  # Find every entry that passed
  indexed=$(find _proc -name "indexed.pickle" | sort)
  count=$(find _proc -name "indexed.pickle" | wc -l)
  echo "  Found $count indexed images"
  args=""
  # Build the argument list to combine_experiments
  for entry in $indexed; do
    args="$args $(dirname $entry)/experiments.json $entry "
  done
  echo "  Running dials.combine_experiments"
  set -x
  $PYEXE -mdials.command_line.combine_experiments \
    reference_from_experiment.beam=0 reference_from_experiment.goniometer=0 reference_from_experiment.detector=0 \
    $args
  set +x
  # Count the number of separate experiments (e.g crystals)
  ccount=$($PYEXE -c "from dxtbx.model.experiment_list import ExperimentListFactory; print(len(ExperimentListFactory.from_serialized_format('combined_experiments.json', check_format=False)))")
  echo "Found $ccount crystals in $count images"
}

# echo "$num_images still images in datablock ($num_len)"
if [[ $MODE == "index" ]]; then
  if [[ -n "$PARALLEL" ]]; then
    # Generate the list of folders
    folders=$($PYEXE -c "print(' '.join('_proc/{0:0{2}d}/{1:0{2}d}'.format((x//100)*100, x, $NWIDTH) for x in range($id_start,$id_end+1)))")
    # Generate the list of image numbers
    numbers=$(seq $id_start $id_end)
    # Generate the offset list of image numbers
    numplus=$(seq $(($id_start+1)) $(($id_end+1)))

    # echo "TESTING: Output to file par.txt"
    echo $folders ::: $numbers ::: $numplus | \
      $PARALLEL --bar -j $NPROCS --link \
      "mkdir -p {1} && cd {1} && echo \$(pwd) && $PYEXE -mdials.command_line.index ../../../$datablock ../../../$reflections scan_range={2},{3} $DEFAULT_PHIL $*" \
      ::: $folders ::: $numbers ::: $numplus 2>&1
  else
    # Run sequentially if no parallel available
    for i in $(seq $id_start $id_end); do
      folder=$($PYEXE -c "print('_proc/{0:0{2}d}/{1:0{2}d}'.format(($i//100)*100, $i, $NWIDTH))")
      ( cd $folder
        echo $(pwd)
        $PYEXE -mdials.command_line.index ../../../$datablock ../../../$reflections scan_range=$i,$(($i+1)) $DEFAULT_PHIL $* 2>&1
      ) || true
    done
  fi
  # If running interactively e.g. not on the grid, then collate now
  if [[ -z "$SGE_MODE" ]]; then
    collate
  else
    echo "Not collating output in SGE mode"
  fi
elif [[ $MODE == "submit" ]]; then
  if [[ -z "$PE_SMP_COUNT" && $NPROCS -gt 1 ]]; then
    PE_SMP_COUNT=$NPROCS
  fi
  if [[ -n "$PE_SMP_COUNT" ]]; then
    PE_SMP_ARGS="-pe smp $PE_SMP_COUNT"
  fi
  # Work out the rough step so that we get N jobs
  step=$(python -c "import math; print(int(math.ceil(float($id_end-$id_start)/$NJOBS)))")
  # Recursively submit ourselves in index mode
  set -x
  qsub -t $id_start-$id_end:$step -b y -e /dev/null $PE_SMP_ARGS \
      $SOURCE -n $NPROCS \
      --python "$PYEXE" --parallel "$PARALLEL" \
      $datablock $reflections $*
  set +x
elif [[ $MODE == "collate" ]]; then
  collate
fi
