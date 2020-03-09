task genesis_nullmodel {
	String outcome_name
	String? outcome_type
	String? covariates_string
	File pheno_file
	File genotype_file
	String results_file
	File? kinship_matrix
	String? pheno_id
	String? conditional
	String? het_varsIn
	String? transform 
	String? transform_rankNorm
	String? transform_rescale

	Int memory
	Int disk

	command {
		R --vanilla --args ${outcome_name} ${default="Continuous" outcome_type} ${default="NA" covariates_string} ${pheno_file} ${genotype_file} ${results_file} ${default="NO_KINSHIP_FILE" kinship_matrix} ${default="ID" pheno_id} ${default="NA" conditional} ${default="NA" het_varsIn} ${default="none" transform} ${default="all" transform_rankNorm} ${default="none" transform_rescale} < /genesis_wdl/genesis_nullmodel.R
	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1_4_1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
	}

	output {
		File results = select_first(glob("${results_file}*"))
	}
}

task genesis_tests {
	File? agg_file
	Float? top_maf
	String? test_stat
	String test_type

	Int? min_mac
	String? weights
	String? weights_col
	Float? user_cores
	Float? window
	Float? step

	File genotype_file
	File null_model
	String results_file
	String? genome_build
	String? pass_only
	String? imputed
	Int? neig
	Float? ntrace
	String? interaction
	String? return_variants


	Int memory
	Int disk

	command {
		R --vanilla --args ${default="NONE" agg_file} ${default="1" top_maf} ${default="Score" test_stat} ${test_type} ${default="5" min_mac} ${default="FALSE" weights} ${default="FALSE" weights_col} ${default="30" user_cores} ${default="0" window} ${default="0" step} ${genotype_file} ${null_model} ${results_file} ${default="hg38" genome_build} ${default="T" pass_only} ${default="F" imputed} ${default="200" neig} ${default="500" ntrace} ${default="NULL" interaction} ${default="F" return_variants} < /genesis_wdl/genesis_tests.R
	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v1_4_1"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File results = select_first(glob("*.gz"))
		File varresults = select_first(glob("*.Rdata"))
	}
}

task summarize {
	Float? pval_threshold
	String results_file
	Array[File] assoc_files

	Int? memory
	Int? disk

	command {
		R --vanilla --args ${default="0.0001" pval_threshold} ${results_file} ${sep="," assoc_files} < /genesis_wdl/summarize_GWAS.R
	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v0.1"
		disks: "local-disk " + select_first([disk,"100"]) + " HDD"
		memory: select_first([memory,"20"]) + " GB"
	}

	output {
		File all_results = "${results_file}.all_variants.assoc.csv.gz"
		File top_results = "${results_file}.top_variants.assoc.csv.gz"
		File plots = "${results_file}.association.plots.png"
	}

}

workflow genesis_gwas_wf {
	# null model inputs
	String this_outcome_name
	String? this_outcome_type
	String? this_covariates_string
	File this_pheno_file
	Array[File] these_genotype_files
	String this_results_file
	File? this_kinship_matrix
	String? this_pheno_id
	String? this_conditional
	String? this_het_varsIn
	String? this_transform 
	String? this_transform_rankNorm
	String? this_transform_rescale
	Int? this_nullmodel_memory
	Int? this_nullmodel_disk

	# association test inputs
	File? this_agg_file
	Float? this_top_maf
	String? this_test_stat
	String this_test_type
	Int? this_min_mac
	String? this_weights
	String? this_weights_col
	Float? this_user_cores
	Float? this_window
	Float? this_step
	String? this_genome_build
	String? this_pass_only
	String? this_imputed
	Int? this_neig
	Float? this_ntrace
	String? this_interaction
	String? this_return_variants
	Int? this_tests_memory
	Int? this_tests_disk

	# summarization inputs
	Float? this_pval_threshold
	Int? this_summarize_memory
	Int? this_summarize_disk
	

	# Workflow metadata
	meta {
		description: "This workflow creates a null model from phenotype data with the GENESIS biostatistical package. This null model can then be used for association testing. This workflow also runs single variant and aggregate test for genetic data. Implements Single-variant, Burden, SKAT, SKAT-O and SMMAT tests for Continuous or Dichotomous outcomes. All tests account for familiar relatedness through kinship matrixes. Underlying functions adapted from: Conomos MP and Thornton T (2016). GENESIS: GENetic EStimation and Inference in Structured samples (GENESIS): Statistical methods for analyzing genetic data from samples with population structure and/or relatedness. R package version 2.3.4. http://bioconductor.org/packages/devel/bioc/html/GENESIS.html https://github.com/smgogarten/GENESIS"
		tags: "Statistics"
	    author: "Jen Brody (base code) & Tim Majarian (modifications and adaptations, WDL, docker)"
	    email: "tmajaria@broadinstitute.org"
	}

	# Parameter metadata
	parameter_meta {
		this_outcome_name: "name: outcome_name, label: outcome, help: Column name of the outcome variable in the phenotype file, class: string, optional: false"
		this_outcome_type: "name: outcome_type, label: Continuous or Dichotomous, class: string, optional: true, default: Continuous"
		this_covariates_string: "name: covariate_list, label: Covariates, help: Comma separated list that match column names in the phenotype file. Leave blank for no adjustments, class: string, optional: true, default: ''"
		this_pheno_file: "name: phenofile, class: file, patterns: [*.csv], optional: false"
		these_genotype_files: "name: genotypefile, class: array[file], patterns: [*.gds, *.GDS], optional: false"
		this_results_file: "name: outputfilename, label: prefix for output file name, no spaces, class: string, optional: false"
		this_kinship_matrix: "name: kinshipmatrix, label: kinship matrix with sample ids as the row and column names.  Matricies saved as Rda will load faster, but csv is accepted as well. Rda files should contain a single numeric matrix object., class: file,patterns: [*.Rda, *.csv], optional: true"
		this_pheno_id: "name: pheno_id, help: Column name that contains the sample IDs.  These IDs should match the genotype file IDs and the kinship file IDs., class: string, default: ID"
		this_conditional: "name: conditional, help: chr pos ref alt format for the SNP that will be added to the model.  Multiple snps in a comma delimited list can be added. (e.g. '22:16425814:C:T' or '22:16425814:C:T,22:17808063:TA:T,22:18096610:G:T'), class: string, optional: true, default: NA"
		this_het_varsIn: "name: het_vars, help: grouping variable for heterogenous variances, class: string, optional: true, default: NA"
		this_transform: "name: transform, label: flag for transforming, help: rank-normalize residuals and scale, and re-fit null model, class: string, optional: true, default: none"
		this_transform_rankNorm: "name:transform_rankNorm, label: transform within het_vars groups or all samples together ( e.g. rankNorm.option in updateNullModOutcome()), only used in conjuntion with 'transform', options are by.group or all, class: string, optional: true, default: all"
		this_transform_rescale: "name: transform_rescale, label: rescale residules  ( e.g. rescale.option in updateNullModOutcome() ) options are none, model, residSD.  Only used in conjuntion with 'transform', class: string, optional: true, default: none"
		this_nullmodel_memory: "help: memory desired for computation in GB, class: int, optional: true, default: 30"
		this_nullmodel_disk: "help: disk space desired for computation in GB, class:int, optional: true, default: 100"
		this_agg_file: "name: varaggfile, help: File contains lists of variants that should be aggregated into groups.  Can also be used to filter variants in sliding window or single variant tests, or topass a variant weight column. File should be a CSV file with the headers: group_id, chromosome, position, ref and alt.  All variants for with the same group_idwill be combined into a single aggregate tests.  , class: file,patterns: [*.csv, *.rda, *.rdata, *.Rda, *.Rdata], optional: true"
		this_top_maf: "name: top_maf, help: Maximim minor allele frequency ( generally used for aggregate tests ), class: float, optional: true, default: 1"
		this_test_stat: "name: test_stat, help: Valid tests statistic types are: Score, Wald. Firth can be used with Burden test only. , class: string, optional: true, default: Score"
		this_test_type: "name: test_type, help: Valid tests are one of the collapsing tests SKAT, Burden, SMMAT or Single, class: string, optional: false"
		this_min_mac: "name: min_mac, help: Minimum minor allele count for threshold ( only used for single variant tests ), class: int, optional: true, default: 5"
		this_weights: "name: weights, help: beta weights set to flat weights (e.g. set to 'c(1,25)' for Wu weights or 'c(0.5,0.5)' for Madsen-Browning weights), class: string, optional: true, default: FALSE"
		this_weights_col: "name: weights_col, help: Name of column in aggfile that contains the variant weights, class: string, default: FALSE"
		this_user_cores: "name: user_cores, class: int, optional: true, default: 30"
		this_window: "name: window, help: Runs a sliding window test based on this window size, in bp, class: int, optional: true, default: 0"
		this_step: "name: step, help: For use with 'window', indicates sliding window step size, in bp, class: int, optional: true, default: 0"
		this_genome_build: "name: genome_build, help: hg38 or hg19, class: string, optional: true, default: hg38"
		this_pass_only: "name: pass_only, help: Filter variants to those with a PASS flag: TRUE or FALSE, class: string, optional: true, default: TRUE"
		this_imputed: "name: imputed, help: Input data is imputed: TRUE or FALSE, class: string, optional: true, default: FALSE"
		this_neig: "name: neig, help: The number eigenvalues to approximate by using random projections for calculating p-values with fastSKAT, class: int, optional: true, default: 200"
		this_ntrace: "name: ntrace, help: The number of vectors to sample when using random projections to estimate the trace needed for p-value calculation with fastSKAT, class: int, optional: true, default: 500"
		this_interaction: "name: interaction, help: character string specifying the name of the variables for which a genotype interaction term should be included.  Single variant only, class: string, optional: true, default: NULL"
		this_return_variants: "name: return_variants, help: Returns single snp results for each aggregate test, class: string, optional: true, default: FALSE"
		this_tests_memory: "help: memory desired for computation in GB, class: int, optional: true, default: 30"
		this_tests_disk: "help: disk space desired for computation in GB, class:int, optional: true, default: 100"
		this_pval_threshold: "help: threshold over association p-value for generating top results table from the summarize task, class: float, optional: true, default: 0.0001"
		this_summarize_memory: "help: memory desired for computation in GB, class: int, optional: true, default: 30"
		this_summarize_disk: "help: memory desired for computation in GB, class: int, optional: true, default: 30"
	}

	call genesis_nullmodel {
		input: 
			outcome_name = this_outcome_name,
			outcome_type = this_outcome_type,
			covariates_string = this_covariates_string,
			pheno_file = this_pheno_file,
			genotype_file = these_genotype_files[0],
			results_file = this_results_file,
			kinship_matrix = this_kinship_matrix,
			pheno_id = this_pheno_id,
			conditional = this_conditional,
			het_varsIn = this_het_varsIn,
			transform = this_transform,
			transform_rankNorm = this_transform_rankNorm,
			transform_rescale = this_transform_rescale,
			memory = this_nullmodel_memory,
			disk = this_nullmodel_disk
	}

    scatter(this_genotype_file in these_genotype_files) {
		call genesis_tests {
			input:
				agg_file = this_agg_file,
				top_maf = this_top_maf,
				test_stat = this_test_stat,
				test_type = this_test_type,
				min_mac = this_min_mac,
				weights = this_weights,
				weights_col = this_weights_col,
				user_cores = this_user_cores,
				window = this_window,
				step = this_step,
				genotype_file = this_genotype_file,
				null_model = this_null_model,
				results_file = this_results_file,
				genome_build = this_genome_build,
				pass_only = this_pass_only,
				imputed = this_imputed,
				neig = this_neig,
				ntrace = this_ntrace,
				interaction = this_interaction,
				return_variants = this_return_variants,
				memory = this_tests_memory,
				disk = this_tests_disk

		}
	}

	call summarize {
		input:
			pval_threshold = this_pval_threshold,
			results_file = this_results_file,
			assoc_files = genesis_tests.results,
			memory = this_summarize_memory,
			disk = this_summarize_disk
	}

	output {
		File null_model = genesis_nullmodel.results
        Array[File] raw_association_files = genesis_tests.results
        File all_summary_statistics = summarize.all_results
        File top_summary_statistics = summarize.top_results
        File summary_plots = summarize.plots
    }
}