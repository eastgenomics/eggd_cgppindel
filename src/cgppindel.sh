#!/bin/bash
# cgppindel 1.1.0

set -exo pipefail

# Install all python packages from this app
sudo -H python3 -m pip install --no-index --no-deps packages/*

main() {

    echo "Value of reference: '$reference'"
    echo "Value of simrep: '$simrep'"
    echo "Value of genes: '$genes'"
    echo "Value of unmatched: '$unmatched'"
    echo "Value of assembly: '$assembly'"
    echo "Value of seqtype: '$seqtype'"
    echo "Value of filter: '$filter'"
    echo "Value of tumour: '$tumour'"
    echo "Value of normal: '$normal'"
    
    # Create input/output directories
    mkdir input
    mkdir -p out/cgppindel_output
    mkdir -p out/output_vcf
    mkdir -p out/vcf_index
    mkdir -p out/output_log
    mkdir -p out/output_vcf_with_vaf
    mkdir temp_logs

    # Make the output directory writeable for the app to work.
    chmod 777 out/cgppindel_output

    dx-download-all-inputs

    #Add all inputs in the same folder to enable indeces to be found.
    find ~/in -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/input

    # Unzip fasta reference
    gzip -d ~/input/$reference_name


    # add dnanexus user to docker group & start docker daemon
    sudo usermod -a -G docker dnanexus
    newgrp docker
    sudo systemctl start docker

    # load local container & get id
    sudo docker load --input ~/input/$docker_image_name 
    cgppindel_id=$(docker images --format="{{.Repository}} {{.ID}}" | grep "^quay.io" | cut -d' ' -f2) 
   
    # Run docker image using id
    sudo docker run -v `pwd`:/data -w "/data/out/cgppindel_output" $cgppindel_id \
    pindel.pl \
    -reference /data/input/${reference_prefix}fa \
    -simrep /data/input/$simrep_name \
    -genes /data/input/$genes_name \
    -unmatched /data/input/$unmatched_name \
    -assembly $assembly \
    -species Human \
    -seqtype $seqtype \
    -filter /data/input/$filter_name \
    -tumour /data/input/$tumour_name \
    -normal /data/input/$normal_name \
    -outdir /data/out/cgppindel_output

    # Add Allele frequency (AF) and Read depth (DP) onto cgppindel output file
    vcf=$(find out/cgppindel_output -type f -name "*.vcf.gz")
    basename=$(basename "$(basename "$vcf" .gz)" .vcf)
    python3 /tsv_file_generator.py -v "$vcf" 
    bgzip annots.tsv
    tabix -s1 -b2 -e2 annots.tsv.gz
    bcftools annotate -a annots.tsv.gz -h /annots.hdr -c CHROM,POS,ID,REF,ALT,+FORMAT/AF,+FORMAT/DP \
    "$vcf" > "$basename.tmp.af.vcf"

    # Remove variants that are smaller than or equal to 2 bp in length
    bcftools view -i 'INFO/LEN > 2' "$basename.tmp.af.vcf" > "$basename.af.vcf"
    bgzip "$basename.af.vcf"


    # Move vcf and index in specified runfolder to enable downstream use
    mv out/cgppindel_output/*.flagged.vcf.gz out/output_vcf
    mv out/cgppindel_output/*.flagged.vcf.gz.tbi out/vcf_index

    # Move vcf file with VAF in speified folder
    mv $basename.af.vcf.gz out/output_vcf_with_vaf

    # Move all .out and .err to a temporary folder
    mv out/cgppindel_output/logs/*.err temp_logs
    mv out/cgppindel_output/logs/*.out temp_logs

    # zip all .out & .err into logs.tar.gz
    tar -zcf out/output_log/logs.tar.gz --remove-files temp_logs

    # Upload output files
    dx-upload-all-outputs --parallel

    echo "Upload Complete"
}
