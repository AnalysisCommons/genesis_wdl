# GENESIS: Statistical methods for analyzing genetic data from samples
Maintainer: Analysis Commons  
Version: 0.2

## Description:
Runs single variant and aggregate test for genetic data. Implements Single-variant, Burden, SKAT and SMMAT tests for Continuous or Dichotomous outcomes. All tests may account for familial relatedness through kinship matrices. Underlying functions adapted from: Conomos MP and Thornton T (2016). GENESIS: GENetic EStimation and Inference in Structured samples (GENESIS): Statistical methods for analyzing genetic data from samples with population structure and/or relatedness. http://bioconductor.org/packages/devel/bioc/html/GENESIS.htmlhttps://github.com/smgogarten/GENESIS

### What do these apps do?
Runs phenotype-genotype association analyses on sequence data. The workflow genesis_nullmodel will run the null model. This is fitting your model with your outcome, adjustments and kinship matrix, but does not use the genotypes. The second workflow genesis_tests takes the null model object from the first step and runs your association analysis. The same null model can be used for single-variant or aggregate tests.

### What are typical use cases for this app?
This workflow can sun single-variant, SKAT and burden tests. It will account for relatedness using a kinship matrix and implements GMMAT for logistic regression.

### What data are required for this app to run?
This workflow requires genotype files in GDS format (\*.gds) and a phenotype file. Rare variant tests require a gene-aggregation file and an annotation file.

### What does this app output?
genesis_null model will output a copy of your null model as an R object. The results of the genesis_tests workflow will be a compressed CSV file of association results

### How does this app work?
This workflow uses the GENESIS package developed by Matt Conomos, Tim Thornton and Stephanie Gogarten (TOPmed DCC). This workflow has been optimized to allow genotype random access from the GDS genotype format files and to parallelize over a number of cores.

## Inputs

-**this_outcome_name**: Column name of the outcome variable in the phenotype file \[string\]
-**this_outcome_type**: Continuous or Dichotomous \[string\] \[optional: Continuous\]
-**this_covariates_string**: Covariates, help: Comma separated list that match column names in the phenotype file. Leave blank for no adjustments \[string\] \[optional: \]
-**this_pheno_file**: Phenotypes for each sample in delimited text file \[file\]
-**these_genotype_files**: Genotypes stored in genomic data structure format, must only have a single chromosome in each file \[array\[file\]\] \[\*.gds, \*.GDS\]
-**this_results_file**: Prefix for output file name, no spaces \[string\]
-**this_kinship_matrix**: Kinship matrix with sample ids as the row and column names. Matricies saved as Rda will load faster, but csv is accepted as well. Rda files should contain a single numeric matrix object. \[file\] \[\*.Rda, \*.csv\] \[optional: \]
-**this_pheno_id**: Column name that contains the sample IDs. These IDs should match the genotype file IDs and the kinship file IDs \[string\] \[optional: ID\]
-**this_conditional**: Variants to condition on, chr pos ref alt format for the SNP that will be added to the model.  Multiple snps in a comma delimited list can be added. (e.g. '22:16425814:C:T' or '22:16425814:C:T,22:17808063:TA:T,22:18096610:G:T') \[string\] \[optional: NA\]
-**this_het_varsIn**: Grouping variable for heterogenous variances \[string\] \[optional: NA\]
-**this_transform**: Rank-normalize residuals and scale, and re-fit null model. Options none or transform \[string\] \[optional: none\]
-**this_transform_rankNorm**: Transform within het_vars groups or all samples together (e.g. rankNorm.option in updateNullModOutcome()), only used in conjuntion with 'transform', options are by.group or all \[string\] \[optional: all\]
-**this_transform_rescale**: Rescale residules  (e.g. rescale.option in updateNullModOutcome()) options are none, model, residSD.  Only used in conjuntion with 'transform' \[string\] \[optional: none\]
-**this_nullmodel_memory**: Memory desired for computation in GB \[int\] \[optional: 30\]
-**this_nullmodel_disk**: Disk space desired for computation in GB \[int\] \[optional: 100\]
-**this_agg_file**: Aggregation file contains lists of variants that should be aggregated into groups. Can also be used to filter variants in sliding window or single variant tests, or topass a variant weight column. File should be a CSV file with the headers: group_id, chromosome, position, ref and alt.  All variants for with the same group_idwill be combined into a single aggregate tests \[file\] \[\*.csv, \*.rda, \*.rdata, \*.Rda, \*.Rdata\] \[optional\]
-**this_top_maf**: Maximim minor allele frequency ( generally used for aggregate tests ) \[float\] \[optional: 1\]
-**this_test_stat**: Valid tests statistic types are: Score and Score.SPA \[string\] \[optional: Score\]
-**this_test_type**: Valid tests are one of the collapsing tests Burden, SKAT, fastSKAT, SMMAT,  SKATO or Single \[string\]
-**this_min_mac**: Minimum minor allele count for threshold \[int\] \[optional: 5\]
-**this_weights**: Beta weights set to flat weights (e.g. set to 'c(1,25)' for Wu weights or 'c(0.5,0.5)' for Madsen-Browning weights) \[string\] \[optional: FALSE\]
-**this_weights_col**: Name of column in aggfile that contains the variant weights \[string\] \[optional: FALSE\]
-**this_user_cores**: Number of cores to use in parallel \[int\] \[optional: min(machine cores - 1, 30)\]
-**this_window**: Runs a sliding window test based on this window size, in bp \[int\] \[optional: 0\]
-**this_step**: For use with 'window', indicates sliding window step size, in bp \[int\] \[optional: 0\]
-**this_genome_build**: hg38 or hg19 \[string\] \[optional: hg38\]
-**this_pass_only**: Filter variants to those with a PASS flag: TRUE or FALSE \[string\] \[optional: TRUE\]
-**this_imputed**: Input data is imputed: TRUE or FALSE \[string\] \[optional: FALSE\]
-**this_neig**: The number eigenvalues to approximate by using random projections for calculating p-values with fastSKAT \[int\] \[optional: 200\]
-**this_ntrace**: The number of vectors to sample when using random projections to estimate the trace needed for p-value calculation with fastSKAT \[int\] \[optional: 500\]
-**this_interaction**: Character string specifying the name of the variables for which a genotype interaction term should be included. Single variant only \[string\] \[optional: NULL\]
-**this_return_variants**: Returns single snp results for each aggregate test \[string\] \[optional: FALSE\]
-**this_tests_memory**: Memory desired for computation in GB \[int\] \[optional: 30\]
-**this_tests_disk**: Disk space desired for computation in GB \[int\] \[optional: 100\]
-**this_pval_threshold**: Threshold over association p-value for generating top results table from the summarize task \[float\] \[optional: 0.0001\]
-**this_summarize_memory**: Memory desired for computation in GB \[int\] \[optional: 30\]
-**this_summarize_disk**: memory desired for computation in GB \[int\] \[optional: 30\]