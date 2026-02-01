#!/bin/bash
set -e

echo "Start time:"; date
echo "Node:"; hostname
echo "PWD:"; pwd

# unpack build directory into sandbox
#tar -xzf build.tar.gz

cd /home/wkarim/bulkice_doumeki/mdom/build/
source env.sh
OUTPUT_DIR="/data/user/wkarim/output/"
# define job number explicitly if needed
JOBN=${PROCESS:-0}

#defining dir for shuffled data
SHUFFLE_INPUT_DIR='/home/wkarim/bulkice_doumeki/mdom/InputFile/'
SHUFFLE_OUTPUT_DIR='/home/wkarim/bulkice_doumeki/mdom/InputFile/Positron'
python shuffle_data.py \
    --particle Positron \
    --input-dir "$SHUFFLE_INPUT_DIR" \
    --output-dir "SHUFFLE_OUTPUT_DIR" \
    --seed "$JOBN"
    
./bulkice_doumeki mdom ibd 88 "$OUTPUT_DIR" "$JOBN"

echo "Job complete!"

