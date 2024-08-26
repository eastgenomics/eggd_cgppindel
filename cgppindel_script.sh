#!/bin/bash

# Step 1: Set the working directory (where your input files are stored)
WORK_DIR=$(pwd)


# Step 2: Find input files in the working directory
# Paths to input files and directories
#REFERENCE="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/GRCh38.no_alt_analysis_set_chr_mask21.fa"
#SIMREP="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/simpleRepeats_sorted.bed.gz"
#GENES="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/coding_unrestricted_GRCh38_myeloid_v1.0.bed"
#UNMATCHED="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/normalPanel.gff3.gz"
#FILTER="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/targetedRules.lst"
TUMOUR="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/tumour_bamfile_markdup.bam"
#NORMAL="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/TA2_S59_L008_tumor_markdup.bam"
#DOCKER_IMAGE="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/cgppindel_image.tar"

REFERENCE="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/reference.fa"
SIMREP="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/genes.bed.gz"
GENES="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/simrep.bed.gz"
UNMATCHED="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/unmatched.gff3.gz"
FILTER="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/filter.lst"
#TUMOUR="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/131517686-24211K0015-24NGSHO35-8128-M-96527893_markdup.bam"
NORMAL="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/normal.bam"
DOCKER_IMAGE="/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/cgppindel_image.tar"

# Other parameters
ASSEMBLY="GRCh38"
SEQTYPE="TG"
SPECIES="Human"
CGPPINDEL_ID="cb44a611a143"


echo "Inputs loaded."

# Step 1: Load the Docker image
sudo docker load --input "${DOCKER_IMAGE}"

# Step 2: Extract the CGPPINDEL Docker ID
docker images --format="{{.Repository}} {{.ID}}" | grep "^quay.io" | cut -d' ' -f2 > cgppindel_id
CGPPINDEL_ID=$(cat cgppindel_id)

# Step 3: Run the CGPPINDEL Docker container
sudo docker run -v "/home/raymondmiles/Desktop/Software_Development/cgppindel_bashexperiment/":/"$(pwd)" \
                -w "$(pwd)" $CGPPINDEL_ID \
                pindel.pl \
                -reference "$REFERENCE" \
                -simrep "$SIMREP" \
                -genes "$GENES" \
                -unmatched "$UNMATCHED" \
                -assembly "$ASSEMBLY" \
                -species "$SPECIES" \
                -seqtype "$SEQTYPE" \
                -filter "$FILTER" \
                -tumour "$TUMOUR" \
                -normal "$NORMAL" \
                -outdir out/

# Step 4: Notify the user
if [[ $? -eq 0 ]]; then
    echo "pindel.pl ran successfully."
    echo "Results are saved in the 'out/' directory."
else
    echo "Error: pindel.pl failed to run."
fi