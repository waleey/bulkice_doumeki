#!/bin/bash

SUFFIX="flat"
SENSOR="mdom"
SENSOR_DIR="mDOM"
MODE="radioactivity"
DEPTH="88"
OUTPUT_DIR="output"
OUTPUT_NR="0"

# Name of output file in output directory
OUTPUT_HITS_FILE=$(echo "./${OUTPUT_DIR}/${SENSOR}_${MODE}_${OUTPUT_NR}.dat")
echo $OUTPUT_HITS_FILE

# Source environment
source /home/jakob/.doumeki_ini.sh

# Set working directory
WORKING_DIR="/home/jakob/software/doumeki/bulkice_doumeki/mdom/build"
cd "$WORKING_DIR" || { echo "Failed to enter directory"; exit 1; }

# Destination directory for Geant4 output files
DEST_DIR=$(echo "/home/jakob/software/doumeki/bulkice_doumeki/analysis/files/output_geant4/background/${SENSOR_DIR}/")

# Name of decay chain files
FULL_CHAIN_FILE="../InputFile/U238DecayChain_copy.txt"
SIM_CHAIN_FILE="../InputFile/U238DecayChain.txt"
EXECUTABLE=$(echo "./bulkice_doumeki ${SENSOR} ${MODE} ${DEPTH} ${OUTPUT_DIR} ${OUTPUT_NR}")

# Get the number of lines in the FULL_CHAIN_FILE file
LINES=$(wc -l < "$FULL_CHAIN_FILE")
HEADER=$(sed -n "1p" "$FULL_CHAIN_FILE")
# Loop through each line and process it
for ((i=2; i<=LINES; i++)); do
    # Read the i-th line
    LINE=$(sed -n "${i}p" "$FULL_CHAIN_FILE")
    ISOTOPE=$(echo "$LINE" | cut -f1)

    echo $i, $LINES
    echo $LINE
    echo $ISOTOPE

    # Overwrite SIM_CHAIN_FILE with the current line
    echo -e "$HEADER\n$LINE" > "$SIM_CHAIN_FILE"
    
    # Name of simulation output files
    VERIFICATION_FILE=$(echo "veri_${ISOTOPE}_${SUFFIX}.dat" )
    PHOTON_HITS_FILE=$(echo "hits_${ISOTOPE}_${SUFFIX}.dat" )
    LOG_FILE=$(echo "log_${ISOTOPE}_${SUFFIX}.dat")

    echo $VERIFICATION_FILE, $PHOTON_HITS_FILE, $LOG_FILE
    echo "${EXECUTABLE} 2>&1 | tee ${LOG_FILE}"

    # Run the executable
    #${EXECUTABLE} 2>&1 | tee ${LOG_FILE}

    # Move output files into destination directory
    #mv ./radiodecay_output.csv ${DEST_DIR}${VERIFICATION_FILE}
    #mv ${OUTPUT_HITS_FILE} ${DEST_DIR}${PHOTON_HITS_FILE}
    #mv ${LOG_FILE} ${DEST_DIR}${LOG_FILE}
    
done

echo "Processing complete."

