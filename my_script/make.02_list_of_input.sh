# List of inputs: it is a space- or tab-separated file, containing 3 columns and no header. 
# For each row, provide the condition code, the path to the first step's output, 
# and the path to the dosage file (the same used for step 1).

# Example:
# GBR	sample_data/cis_ES/GBR_chr22_cis.txt	sample_data/dosages/GEUVADIS_GBR_chr22_dosage_filtered.txt.gz
# YRI	sample_data/cis_ES/YRI_chr22_cis.txt	sample_data/dosages/GEUVADIS_YRI_chr22_dosage_filtered.txt.gz

mkdir -p ../pipeline_output/02_list_of_input

# remove the existing file
rm -f ../pipeline_output/02_list_of_input/chr*_input_files.txt

for pop in AFR EUR HISP
do
  for chr in {1..22}
  do
    echo "${pop} ../pipeline_output/MatrixeQTL/${pop}_chr${chr}_cis.txt ../my_input/processed_input/whole_blood_${pop}_chr${chr}/snp_dosage.txt" >> ../pipeline_output/02_list_of_input/chr${chr}_input_files.txt
  done
done

