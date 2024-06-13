library(lassosum)
library(fdrtool)
library(data.table)
library(argparse)
suppressMessages(library(tidyverse))


# TODO: in the future add HISP support, by reading HISP's LD block; rerun the script with HISP

# Arguments using argparse
parser <- ArgumentParser()
parser <- ArgumentParser(description = 
    "This is a script that read in the MASHR results, and perform LassoSum. 
    Currently it works only for MAGENTA data using in-sample LD.
    The script will output the LassoSum results for each gene on a given chromosome.
    Note!! HISP is not working yet, since the LD block is not general,due to admixture.")
parser$add_argument("--chr", help="Chromosome number of the gene", type="integer")
parser$add_argument("--genomebuild", help="whether the genome build is hg19 or hg38", type="character", default="hg38")
parser$add_argument("--pop", help="Population of interests. By default, only AFR and EUR is considered right now. HISP will be **ignored**", type="character", default="AFR-EUR")

args <- parser$parse_args()


# Define a list of populations, based on the `args$pop`
pop_list <- if (args$pop == "AFR-EUR") {
    c("AFR", "EUR")
} else {
    c(args$pop)
}

# chr	gene_id	gene_name	start	end	gene_type
# 1	ENSG00000000457	SCYL3	169849631	169894267	protein_coding
# 1	ENSG00000000460	C1orf112	169662007	169854080	protein_coding
# 1	ENSG00000000938	FGR	27612064	27635185	protein_coding
# 1	ENSG00000001460	STPG1	24356999	24416934	protein_coding
# 1	ENSG00000001461	NIPAL3	24415802	24472976	protein_coding
# 1	ENSG00000004455	AK2	33007986	33080996	protein_coding
gene_annotation <- fread(paste0("../my_input/processed_input/whole_blood_EUR", "_chr", args$chr, "/gene_annotation.txt"))


# Define an error log dataframe
error_log <- data.frame(gene = character(), pop = character(), error = character(), stringsAsFactors = FALSE)

# Function to log errors
log_error <- function(e, gene, pop) {
  error_log <<- rbind(error_log, data.frame(gene = gene, pop = pop, error = e$message, stringsAsFactors = FALSE))
}

# Function to log warnings
log_warning <- function(warning, gene, pop) {
  error_log <<- rbind(error_log, data.frame(gene = gene, pop = pop, error = warning, stringsAsFactors = FALSE))
}

# read in MASHR results
for (i in 1:nrow(gene_annotation)) {
# #just for testing
# for (i in 1:1) {

    gene <- gene_annotation$gene_id[i]
    gene_prefix <- paste0("../pipeline_output/MASHR_outputs/", gene, "_MASHR_")

    # gene snps snp_ID AFR_beta EUR_beta HISP_beta
    # ENSG00000000457 chr1:168849631 chr1:168849631_C_A -0.00807105705901084 -0.00825318958066658 -0.0091238403868618
    # ENSG00000000457 chr1:168851510 chr1:168851510_G_A -0.00890800483599222 0.00353274285522946 0.00563178893684715
    # ENSG00000000457 chr1:168851541 chr1:168851541_T_G -0.00765637691084664 0.00473959752054153 0.00428656768298717
    mashr_beta_ori <- fread(paste0(gene_prefix, "beta.txt"))
    # gene snps snp_ID AFR_SD EUR_SD HISP_SD
    # ENSG00000000457 chr1:168849631 chr1:168849631_C_A 0.0228784848562272 0.0232665592172763 0.0237501151481989
    # ENSG00000000457 chr1:168851510 chr1:168851510_G_A 0.0392085296751047 0.0278376038596746 0.0242421225809465
    # ENSG00000000457 chr1:168851541 chr1:168851541_T_G 0.0379419226099434 0.028505115076842 0.0219652869567535
    mashr_SD <- fread(paste0(gene_prefix, "SD.txt"))

    error_occurred <- FALSE  # Reset error flag for each gene

    # loop through the populations, AFR, EUR **(since HISP is not working yet, due to admixture)**
    # for each iteration update the mashr_beta_ori with the beta of the ancestry of interest
    for (pop in pop_list) {
        # Define plink file location
        plink_prefix <- paste0("../my_input/raw_input/_geno_plink_common/", pop, "_chr", args$chr, "_common")

        # select the beta & SD of the ancestry of interest
        cols <- c("gene", "snp_ID", paste0(pop, "_beta"))
        mashr_beta <- mashr_beta_ori[, ..cols]
        # print(head(mashr_beta_ori))
        # print(head(mashr_beta))

        # check is `paste0(pop, "_beta")` column is all zeros
        # if all zeros, skip the population
        if (sum(mashr_beta[[paste0(pop, "_beta")]] == 0) == nrow(mashr_beta)) {
            print(paste("Warning: All zeros for", gene, pop))
            # save the warning to the error log
            log_warning("Mashr beta estimates are all zeros", gene, pop)
            next
        }

        # exract A1: the alternative (effect) allele, A2: the reference allele
        # extract SNP's position and chromosome
        mashr_beta <- as_tibble(mashr_beta) %>%
            separate(snp_ID, into = c("chr", "pos", "a2", "a1"), sep = "[:_]", remove = FALSE) %>%
            mutate(chr = as.numeric(str_replace(chr, "chr", ""))) %>%
            mutate(pos = as.numeric(pos))

        # Lassosum template
        # out <- lassosum.pipeline(cor=cor,
        #                      chr=as.numeric(ss$Chrom),
        #                      pos=as.numeric(ss$SNPPos),
        #                      A1=ss$A1,
        #                      A2=ss$A2, # A2 is not required but advised
        #                      s = c(0.2, 0.5, 0.9, 1),
        #                      lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
        #                      ref.bfile = bfile, # The reference panel dataset
        #                      test.bfile = bfile, # We don't have test data here
        #                      LDblocks = argsL$LDblocks,
        #                      exclude.ambiguous = F,
        #                      destandardize = F,
        #                      trace = 0)

        # Read the summary statistics of standardized beta in single variant test
        # the standardized beta in single variant test = correlation
        z <- mashr_beta[[paste0(pop, "_beta")]] / mashr_SD[[paste0(pop, "_SD")]]
        cor <- z
        # lassosum only allow -1 < cor < 1
        if (sum(abs(cor) >= 1) > 0){
            shrink_factor = max(abs(cor)) / 0.9999
            cor = cor / shrink_factor
        }
        # My lambda values
        ori <- exp(seq(log(0.001), log(0.1), length.out=20))
        vec <- c(0.05, 0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9, 0.95)
        my_lambda <- append(ori,vec)
        my_lambda <- sort(my_lambda)

        # Run LassoSum within tryCatch
        tryCatch({
            out <- lassosum.pipeline(cor = cor,
                                     chr = as.numeric(mashr_beta$chr),
                                     pos = as.numeric(mashr_beta$pos),
                                     A1 = mashr_beta$a1,
                                     A2 = mashr_beta$a2,
                                     s = c(0.05, 0.1, 0.2, 0.25, 0.33, 0.5, 0.75, 0.9, 0.95, 1),
                                     lambda = my_lambda,
                                     ref.bfile = plink_prefix,
                                     LDblocks = paste0(pop, ".", args$genomebuild),
                                     exclude.ambiguous = FALSE,
                                     destandardize = FALSE,
                                     trace = 0)
            
        }, error = function(e) {
            log_error(e, gene, pop)
            error_occurred <<- TRUE  # Set the error flag
        })


        if (!error_occurred) {
            # perform pseudovalidation
            v <- pseudovalidate(out)
            lassosum_out <- subset(out, s = v$best.s, lambda = v$best.lambda)
            beta <- unlist(lassosum_out$beta)

            # Add the `beta` back to the original `mashr_beta_ori` df
            mashr_beta_ori[[paste0(pop, "_beta")]] <- beta
        } else {
            mashr_beta_ori[[paste0(pop, "_beta")]] <- 0
        }
    } 

    # TODO make directory and save the `mashr_beta_ori` df to a file under `pipeline_output/LassoSum_outputs/`
    # Create directory and save the `mashr_beta_ori` df to a file under `pipeline_output/LassoSum_outputs/`
    dir.create("../pipeline_output/LassoSum_outputs", showWarnings = FALSE)

    if (args$pop == "AFR-EUR") {
        mashr_beta_ori[, HISP_beta := NULL]
    }

    write.table(mashr_beta_ori, file = paste0("../pipeline_output/LassoSum_outputs/", gene, "_LassoSum_beta.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    
}

# save the error log for each chromosome, if any
if (nrow(error_log) > 0) {
    # Create directory to save the error log
    dir.create("lassosum_err_log", showWarnings = FALSE)
    write.table(error_log, file = paste0("lassosum_err_log/lassosum_error_log_chr", args$chr, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}
