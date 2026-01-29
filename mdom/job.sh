#!/usr/bin/bash

# grab the arguments
usage() { echo "Usage: $0 [-h] [-s script] [-w wait] [-j jobn] [-v venv]" 1>&2; exit 1; }

while getopts "hs:w:j:v:" opt; do
    case ${opt} in
        s)
            script=${OPTARG}
            ;;
        w)
            wait_t=${OPTARG}
            ;;
        j)
            jobn=${OPTARG}
            ;;
        v)
            venv=${OPTARG}
            ;;
        h:*)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

# Check for missing arguments, exit if we don't have everything
if [[ -z $script ]] || [[ -z $wait_t ]] || [[ -z $jobn ]] ; then
    usage
fi

set -e # exit on error
printf "Start time: "; /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job is running in directory: "; /bin/pwd

# choose your python from cvmfs
eval $(/cvmfs/icecube.opensciencegrid.org/py3-v4.4.0/setup.sh)
unset PYTHONPATH

# you may need to change these paths depending on where your virtual environment is
#source $venv/bin/activate
echo "Script: ${script}"
echo "Args: $wait_t $jobn"

#adding this to set up the environment for bulkice
cd /scratch/wkarim/bulkice_doumeki/mdom/build/
source env.sh

# if using icetray for your job, your env-shell.sh should be passed on the same line. For example:
# /cvmfs/icecube.opensciencegrid.org/users/jthwaites/icetray_updated/build/env-shell.sh python $script --wait $2 --jobn $3
#python $script --wait $wait_t --jobn $jobn
./bulkice_doumeki mdom ibd 88 output $jobn

echo "Job complete!"
printf "Finish time: "; /bin/date