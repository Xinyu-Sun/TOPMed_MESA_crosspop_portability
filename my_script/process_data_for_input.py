import pandas as pd
import numpy as np
import os
from pandas_plink import read_plink
import argparse
import statsmodels.api as sm

__author__ = 'Xinyu Sun'

# Three jobs in the script:
# 1. Read SNP dosage from PLINK files, and output in txt format
# 2. Read gene expression data, and output in txt format
# 3. Make gene annotation file that mainly records gene start and end position, and gene IDs

def parse_args():
    parser = argparse.ArgumentParser(description="Generate .sbams files for each gene with consistent sample order and ancestry.")
    parser.add_argument("-e", "--expr", required=True, help="Path to expression file")
    parser.add_argument("-g", "--geno", required=True, help="Prefix to PLINK genotype files. Note: prefix can contain path.")
    parser.add_argument("-c", "--covar", required=True, help="Path to covariate file")
    parser.add_argument("-t", "--tissue", default="tissue", help="Tissue type")
    parser.add_argument("-a", "--ancestry", required=True, help="Ancestry group (e.g., NHW, AA, HISP)")
    parser.add_argument("-chr", "--chromosome", required=True, type=str, help="Chromosome number")
    parser.add_argument("-o", "--output", default="../output", help="Output directory to save .sbams files")
    parser.add_argument("-anno", "--annotation", default="/mnt/pan/Data14/xxs410_xinyu/mashr_multiethnic_related/process_dapg_input/input/_expression_MAGENTA/MAGENTA_GENE_POS_filtered.txt", help="Path to gene annotation file")
    return parser.parse_args()

def load_data(expr_path, covar_path, plink_prefix, chromosome, annotation_path):
    # Load data
    df_expr = pd.read_csv(expr_path, sep="\t")
    df_expr = df_expr[df_expr['chr'].astype(str) == chromosome]
    df_covar = pd.read_csv(covar_path, sep="\t", comment=None)
    # df_gene's columns: ensg_ID, Chromosome, sbp, ebp, symbol
    df_gene = pd.read_csv(annotation_path, sep="\t")
    # rename and reorder the df_gene's columns to ["chr", "gene_id", "gene_name", "start", "end"]
    df_gene = df_gene.rename(columns={"ensg_ID": "gene_id", "Chromosome": "chr","sbp": "start", "ebp": "end", "symbol": "gene_name"})
    df_gene = df_gene[["chr", "gene_id", "gene_name", "start", "end"]]
    # add gene_type column to df_gene
    df_gene["gene_type"] = "protein_coding"
    df_gene = df_gene[df_gene['chr'].astype(str) == chromosome]
    
    # Ensure the order of gene_id in df_expr is the same as in df_gene
    df_expr = df_expr[df_expr['gene_id'].isin(df_gene['gene_id'])]
    # Filter df_gene to ensure it only contains gene_ids present in df_expr
    df_gene = df_gene[df_gene['gene_id'].isin(df_expr['gene_id'])]
    df_expr = df_expr.set_index('gene_id').loc[df_gene['gene_id']].reset_index()
    
    
    bim, fam, G = read_plink(f"{plink_prefix}_chr{chromosome}_common", verbose=False)
    G = 2 - G  # Convert genotype matrix to dosage format

    # Align samples across genotype, expression, and covariate data
    geno_samples = fam['iid'].tolist()
    expr_samples = df_expr.columns[4:].tolist()  # Assuming sample IDs start from the 5th column
    covar_samples = df_covar.columns[1:].tolist() # Assuming sample IDs start from the 2nd column

    # Find common samples and sort them to maintain consistency
    common_samples = sorted(set(geno_samples).intersection(expr_samples, covar_samples))

    # Filter and reorder genotype data
    fam_filtered = fam[fam['iid'].isin(common_samples)].sort_values(by='iid')
    G_filtered = G[:, fam_filtered.index]  # Aligning genotype matrix to the common sorted samples


    # Filter and reorder expression and covariate data to match the order in fam_filtered
    df_expr = df_expr[['chr', 'start', 'end', 'gene_id'] + common_samples]
    df_covar = df_covar[['#id'] + common_samples]

    
    return df_expr, df_covar, df_gene, bim, fam_filtered, G_filtered

def residual_expression(df_expr, df_covar):
    # Extract sample IDs from expression data
    sample_ids = df_expr.columns[4:]

    # Ensure the covariate data matches the sample IDs in expression data
    df_covar = df_covar.set_index(df_covar.columns[0])
    covar_data = df_covar.loc[:, sample_ids].T
    covar_data.columns = df_covar.index

    # Perform regression for each gene and get residuals
    residuals = pd.DataFrame(index=sample_ids, columns=df_expr.index)
    for idx, row in df_expr.iterrows():
        y = row[sample_ids].astype(float)
        X = covar_data.astype(float)
        X = sm.add_constant(X)
        model = sm.OLS(y, X, missing='drop').fit()
        residuals.loc[:, idx] = model.resid

    # Transpose back to original format and add required columns
    residuals = residuals.T
    residuals.insert(0, 'chr', df_expr['chr'].values[0])
    residuals.insert(1, 'start', df_expr['start'].values)
    residuals.insert(2, 'end', df_expr['end'].values)
    residuals.insert(3, 'gene_id', df_expr['gene_id'].values)

    return residuals

def write_files(df_expr, df_covar, df_gene, bim, fam, G, tissue, ancestry, output_path, chromosome):
    output_directory = os.path.join(output_path, f"{tissue}_{ancestry}_chr{chromosome}")
    os.makedirs(output_directory, exist_ok=True)

    # SNP Dosage File
    snp_file_path = os.path.join(output_directory, "snp_dosage.txt")
    with open(snp_file_path, 'w') as snp_file:
        snp_header = "chr\tsnp_ID\tpos\tref_allele\talt_allele\t" + "\t".join(fam['iid']) + "\n"
        snp_file.write(snp_header)
        for idx in bim.index:
            snp_info = (
                str(chromosome) + "\t" +
                bim.loc[idx, 'snp'] + "\t" +
                str(bim.loc[idx, 'pos']) + "\t" +
                bim.loc[idx, 'a1'] + "\t" +
                bim.loc[idx, 'a0'] + "\t" +
                "\t".join(map(str, G[idx].compute().astype(int))) + "\n"
            )
            snp_file.write(snp_info)

    # Gene Expression File
    expr_file_path = os.path.join(output_directory, "gene_expression.txt")
    with open(expr_file_path, 'w') as expr_file:
        expr_header_base = "gene_id\t"
        expr_header_samples = "\t".join(df_expr.columns[4:])
        expr_header = expr_header_base + expr_header_samples + "\n"
        expr_file.write(expr_header)
        for idx, row in df_expr.iterrows():
            if row['gene_id'] in df_gene['gene_id'].values:
                expr_values = '\t'.join(map(str, row[4:]))
                expr_line = f"{row['gene_id']}\t{expr_values}\n"
                expr_file.write(expr_line)

    # Gene Annotation File
    gene_anno_file_path = os.path.join(output_directory, "gene_annotation.txt")
    with open(gene_anno_file_path, 'w') as anno_file:
        anno_header = "chr\tgene_id\tgene_name\tstart\tend\tgene_type\n"
        anno_file.write(anno_header)
        for idx, row in df_gene.iterrows():
            gene_info = f"{chromosome}\t{row['gene_id']}\t{row['gene_name']}\t{row['start']}\t{row['end']}\t{row['gene_type']}\n"
            anno_file.write(gene_info)

def main():
    args = parse_args()
    df_expr, df_covar, df_gene, bim, fam, G = load_data(args.expr, args.covar, args.geno, args.chromosome, args.annotation)
    df_expr_residuals = residual_expression(df_expr, df_covar)
    write_files(df_expr_residuals, df_covar, df_gene, bim, fam, G, args.tissue, args.ancestry, args.output, args.chromosome)

if __name__ == "__main__":
    main()
