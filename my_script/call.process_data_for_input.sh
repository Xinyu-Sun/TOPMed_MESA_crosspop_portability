# Example usage of process_data_for_input.py
# python process_data_for_input.py -e ../my_input/raw_input/_expression_MAGENTA/AA_TPM_bed.cleaned.txt \
#     -g ../my_input/raw_input/_geno_plink_common/AFR \
#     -c ../my_input/raw_input/_cov_MAGENTA/magenta_aa_no_rel_no_pca_outlier.rnaseqc.low_expression_filtered.outlier_removed.tmm.expression.cov_pca.resid.PEER.cov \
#     -t whole_blood -a AFR -chr 22 \
#     -o ../my_input/processed_input

# Example code:
if [ $# -eq 0 ]; then
    echo "No arguments provided. Please provide the ancestry group (AA or NHW or HISP) and the chromosome number (1..22)"
    exit 1
fi

if [ "$1" != "AA" ] && [ "$1" != "NHW" ] && [ "$1" != "HISP" ]; then
    echo "Invalid ancestry group. Please provide the ancestry group (AA, NHW, or HISP)"
    exit 1
fi


if [ $2 -lt 1 -o $2 -gt 22 ]; then
    echo "Invalid chromosome number. Please provide the chromosome number (1..22)"
    exit 1
fi


if [ $1 == "AA" ]; then
  # call the python script with the appropriate arguments + chromosome number
  python process_data_for_input.py -e ../my_input/raw_input/_expression_MAGENTA/AA_TPM_bed.cleaned.txt \
    -g ../my_input/raw_input/_geno_plink_common/AFR \
    -c ../my_input/raw_input/_cov_MAGENTA/magenta_aa_no_rel_no_pca_outlier.rnaseqc.low_expression_filtered.outlier_removed.tmm.expression.cov_pca.resid.PEER.cov \
    -t whole_blood -a AFR -chr $2 \
    -o ../my_input/processed_input
fi


if [ $1 == "NHW" ]; then
  # call the python script with the appropriate arguments + chromosome number
  python process_data_for_input.py -e ../my_input/raw_input/_expression_MAGENTA/NHW_TPM_bed.cleaned.txt \
    -g ../my_input/raw_input/_geno_plink_common/EUR \
    -c ../my_input/raw_input/_cov_MAGENTA/nhw_no_rel.rnaseqc.low_expression_filtered.outlier_removed.tmm.expression.cov_pca.resid.PEER.cov \
    -t whole_blood -a EUR -chr $2 \
    -o ../my_input/processed_input
fi

if [ $1 == "HISP" ]; then
  # call the python script with the appropriate arguments + chromosome number
  python process_data_for_input.py -e ../my_input/raw_input/_expression_MAGENTA/HISP_TPM_bed.cleaned.txt \
    -g ../my_input/raw_input/_geno_plink_common/HISP \
    -c ../my_input/raw_input/_cov_MAGENTA/magenta_hisp_imp_fil.rnaseqc.low_expression_filtered.outlier_removed.tmm.expression.base_cov_pca.resid.PEER.cov \
    -t whole_blood -a HISP -chr $2 \
    -o ../my_input/processed_input
fi

