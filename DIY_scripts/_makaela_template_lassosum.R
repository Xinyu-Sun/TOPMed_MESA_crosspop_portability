library(data.table)
library(lassosum)
library(fdrtool)
suppressMessages(library(tidyverse))

args = commandArgs(trailing = TRUE)
PATHRef = args[1]
PATHTest = args[2]
PATHSumStats = args[3]
PATHPheno = args[4]
PATHPhenoTrain = args[5]
PATHCov = args[6]
GENE = args[7]
CHR = args[8]
TISSUE = args[9]
FILTER = args[10]

# CORTEX

bfile <- PATHRef
test.bfile <- PATHTest
sum.stat <- PATHSumStats

# We will need the EUR.hg38 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.

ld.file <- "EUR.hg38"

# output prefix
prefix <- "EUR"

# Read in the summary statistics:
## Create 2 versions: one for the lassosum standard method & other for the OTTERs method
ss <- read.delim(sum.stat, sep = "\t")

fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

#cat eqtlgen_aligned_snps_hd.txt
#chr     pos_b38 pos_b37 ref_allele      alt_allele      rsid    zscore  gene    hgnc    gene_chr        gene_pos_b37    n_cohorts  n_samples        p_value fdr     bonferroni_p_value

# Otters suggested correlation function

# Read the summary statistics of standardized beta in single variant test
# the standardized beta in single variant test = correlation
cor <- ss$Z

# lassosum only allow -1 < cor < 1
if (sum(abs(cor) >= 1) > 0){
  
  shrink_factor = max(abs(cor)) / 0.9999
  
  cor = cor / shrink_factor
  
}

ss$cor <- cor

# ADDED: new lambda values 29-09-2021

ori <- exp(seq(log(0.001), log(0.1), length.out=20))
vec <- c(0.05, 0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9, 0.95)
new_vec <- append(ori,vec)
#new_vec


# Run Lassosum

tryCatch({

	out <- lassosum.pipeline(
    		cor = ss$cor,
    		chr = ss$CHROM,
    		pos = ss$POS,
    		A1 = ss$A1 ,
    		A2 = ss$A2,
    		s = c(0.05, 0.1,0.2, 0.25, 0.33, 0.5, 0.75, 0.9, 0.95, 1),
    		lambda = new_vec ,
    		ref.bfile = bfile,
    		test.bfile = test.bfile,
    		LDblocks = ld.file,
    		exclude.ambiguous = FALSE, 
    		trace = 0,
    		destandardize = F)
}, error = function(e) {
	# Error handling code
	message('An Error Occurred when running the lassosum.pipeline')
 	print(GENE)
	print(TISSUE)
	print(FILTER)
	print(e)
	reason1 <- ifelse(grepl("bmerge", e) == 'TRUE', 'bmerge_no_ss', 'NA')
	reason2 <- ifelse(grepl("Not converging", e) == 'TRUE', 'not_converging', 'NA')
	res <- data.frame(gene = GENE, filter = FILTER, tissue = TISSUE, lambda_lsum = NA, s_lsum = NA, r_lsum = reason2, r2_lsum = NA, n_snps_in_win_ss = reason1, n_snps_in_mod_lsum = NA)
	name_results <- paste0("/mnt/pan/Data14/metabrain_lasso/lassosum_", TISSUE, "/results_table/", FILTER,"_fil/",TISSUE, "_", GENE, "_", FILTER , "_res_tab_lassosum_1kghighcov_ref_gtex_test.txt")
	write_tsv(res, file = name_results, col_names = TRUE)

})

		

## VALIDATE results: select optimal hyperparameters

# GTEx Cortex Expression Matrix
gtex_exp <- read.delim(file = PATHPheno, sep = "\t", header = FALSE)
# GTEx Cortex Covariates
gtex_cov <- read.delim(file = PATHCov, sep = "\t", header = FALSE)


gtex_exp_t <- data.frame(t(gtex_exp))
colnames(gtex_exp_t) <- gtex_exp_t[4,]
gtex_exp_t <- gtex_exp_t[-c(1:4),]
gtex_exp_t$IID <- gtex_exp_t$gene_id
gtex_exp_t <- gtex_exp_t[,c(1,3,2)]
colnames(gtex_exp_t) <- c('FID', 'IID', 'Exp')
gtex_exp_t$Exp <- as.numeric(gtex_exp_t$Exp)


gtex_cov_t <- data.frame(t(gtex_cov))
colnames(gtex_cov_t) <- gtex_cov_t[1,]
gtex_cov_t <- gtex_cov_t[-c(1),]
gtex_cov_t$FID <- gtex_cov_t$ID
gtex_cov_t$IID <- gtex_cov_t$ID
gtex_cov_t_pc <- gtex_cov_t[,-c(1,(dim(gtex_cov_t)[2] -1):dim(gtex_cov_t)[2])]
for (i in 1:dim(gtex_cov_t_pc)[2]){
    gtex_cov_t_pc[,i] <- as.numeric(gtex_cov_t_pc[,i])
}
gtex_cov_t_id <- gtex_cov_t[,c((dim(gtex_cov_t)[2] -1):dim(gtex_cov_t)[2])]
gtex_cov_t <- cbind(gtex_cov_t_id,gtex_cov_t_pc)


##### USING  META-ANALYSIS ######################
# Store the R2 results
v <- validate(out, pheno = gtex_exp_t, covar = gtex_cov_t)


rs <- data.frame(gene = GENE, filter = FILTER, tissue = TISSUE)

# Save Lassosum - Original Results
best_s <- v$best.s
best_lambda <- v$best.lambda
vt <- v$validation.table
r <- vt[vt$lambda == best_lambda & vt$s == best_s, ]
res <- cbind(rs,r)
res$r2_lsum <- (res$value)^2
res  <- res %>% rename(lambda_lsum = lambda,
                       s_lsum = s,
                       r_lsum = value)

# Add original number of SNPs

res$n_snps_in_win_ss <- dim(ss)[1]

# Add number of SNPs in each model! 

sus <- out$sumstats
sus$best_beta <- v$best.beta
res$n_snps_in_mod_lsum <- sum(sus$best_beta != 0)



# Save sumstats 


sumstats <- out$sumstats %>% select(-c(order)) 
sumstats$TargetID <- GENE
sumstats$ES <- v$best.beta

sumstats2 <- sumstats %>% rename(CHROM = chr, POS = pos) %>% filter(ES != 0) %>% select(CHROM, POS, A1, A2, TargetID, ES)


## VALIDATE results: select optimal hyperparameters

# GTEx Cortex Expression Matrix
gtex_exp <- read.delim(file = PATHPhenoTrain, sep = "\t", header = FALSE)


gtex_exp_t <- data.frame(t(gtex_exp))
colnames(gtex_exp_t) <- gtex_exp_t[4,]
gtex_exp_t <- gtex_exp_t[-c(1:4),]
gtex_exp_t$IID <- gtex_exp_t$gene_id
gtex_exp_t <- gtex_exp_t[,c(1,3,2)]
colnames(gtex_exp_t) <- c('FID', 'IID', 'Exp')
gtex_exp_t$Exp <- as.numeric(gtex_exp_t$Exp)




##### USING  META-ANALYSIS ######################
# Store the R2 results
v_train <- validate(out, pheno = gtex_exp_t, covar = gtex_cov_t)


# Save Lassosum - Original Results
best_s <- v_train$best.s
best_lambda <- v_train$best.lambda
v_train_t <- v_train$validation.table
r <- v_train_t[v_train_t$lambda == best_lambda & v_train_t$s == best_s, ]
res <- cbind(res,r)
res$r2_lsum_train <- (res$value)^2
res  <- res %>% rename(lambda_lsum_train = lambda,
                       s_lsum_train = s,
                       r_lsum_train = value)

# Add number of SNPs in each model!

sus$best_beta_train <- v_train$best.beta
res$n_snps_in_mod_lsum_train <- sum(sus$best_beta_train != 0)



# Save sumstats


sumstats_train <- out$sumstats %>% select(-c(order))
sumstats_train$TargetID <- GENE
sumstats_train$ES <- v_train$best.beta

sumstats_train <- sumstats_train %>% rename(CHROM = chr, POS = pos) %>% filter(ES != 0) %>% select(CHROM, POS, A1, A2, TargetID, ES)


# Goal: provide ES directly to OTTERS

# File path
#file_path_lsum <- file.path("/mnt/pan/Data14/metabrain_lasso/otters", TISSUE, "model_weights", paste0("chr_",CHR), paste0(FILTER, "_val"), "lsum.txt")

# Check if the file exists
#if (file.exists(file_path_lsum)) {
  # If file exists, read the existing data
#  existing_data <- fread(file_path_lsum, sep = "\t")

  # Append the new data to the existing data
#  updated_data <- rbind(existing_data, sumstats2)

  # Write the updated data back to the file
#  fwrite(updated_data, file_path_lsum, sep = "\t", quote = FALSE)
#  print("Credible SNPs appended to existing file in eQTL_ss")
#} else {
  # If file doesn't exist, create a new file with the new data
#  fwrite(sumstats2, file_path_lsum, sep = "\t", quote = FALSE)
#  print("New file created with credible SNPs in eQTL_ss")
#}


# Save Files 

name_results <- paste0("/mnt/pan/Data14/metabrain_lasso/lassosum_", TISSUE, "/results_table/", FILTER,"_fil/",TISSUE, "_", GENE, "_", FILTER , "_res_tab_lassosum_1kghighcov_ref_gtex_test.txt")

write_tsv(res, file = name_results, col_names = TRUE)

name_beta <- paste0("/mnt/pan/Data14/metabrain_lasso/lassosum_", TISSUE, "/best_beta/", FILTER,"_fil/chr_", CHR, "/", TISSUE, "_", GENE, "_", FILTER , "_best_beta_lassosum_1kghighcov_ref_gtex_test.txt")

write_tsv(sumstats2, file = name_beta, col_names = TRUE)

name_beta_train <- paste0("/mnt/pan/Data14/metabrain_lasso/lassosum_", TISSUE, "/best_beta/", FILTER,"_fil/chr_", CHR, "/", TISSUE, "_", GENE, "_", FILTER , "_best_beta_lassosum_1kghighcov_ref_gtex_test_split_validation.txt")

write_tsv(sumstats_train, file = name_beta_train, col_names = TRUE)


