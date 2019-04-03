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

	Int memory
	Int disk

	command {
		R --vanilla --args ${default="NONE" agg_file} ${default="1" top_maf} ${default="Score" test_stat} ${test_type} ${default="5" min_mac} ${default="FALSE" weights} ${default="FALSE" weights_col} ${default="30" user_cores} ${default="0" window} ${default="0" step} ${genotype_file} ${null_model} ${results_file} < /genesis_dnanexus/genesis_tests.R
	}

	runtime {
		docker: "tmajarian/genesis_dnanexus:v0.1"
		disks: "local-disk ${disk} SSD"
		memory: "${memory} GB"
	}

	output {
		File results = select_first(glob("*.gz"))
	}
}

workflow genesis_tests_wf {
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

	Array[File] these_genotype_files
	File this_null_model
	String this_results_file

	Int this_memory
	Int this_disk

	# Workflow metadata
	meta {
		summary: "Runs single variant and aggregate test for genetic data. Implements Single-variant, Burden, SKAT, SKAT-O and SMMAT tests for Continuous or Dichotomous outcomes. All tests account for familiar relatedness through kinship matrixes. Underlying functions adapted from: Conomos MP and Thornton T (2016). GENESIS: GENetic EStimation and Inference in Structured samples (GENESIS): Statistical methods for analyzing genetic data from samples with population structure and/or relatedness. R package version 2.3.4. http://bioconductor.org/packages/devel/bioc/html/GENESIS.html https://github.com/smgogarten/GENESIS"
		tags: "Statistics"
	    author: "Jen Brody (base code) & Tim Majarian (modifications and adaptations, WDL, docker)"
	    email: "tmajaria@broadinstitute.org"
	}

	# Parameter metadata
	parameter_meta {
		this_agg_file: "name: varaggfile, help: File contains lists of variants that should be aggregated into groups.  Can also be used to filter variants in sliding window or single variant tests, or topass a variant weight column. File should be a CSV file with the headers: group_id, chromosome, position, ref and alt.  All variants for with the same group_idwill be combined into a single aggregate tests.  , class: file,patterns: [*.csv, *.rda, *.rdata, *.Rda, *.Rdata], optional: true"
		this_top_maf: "name: top_maf, help: Maximim minor allele frequency ( generally used for aggregate tests ), class: float, optional: true, default: 1"
		this_test_stat: "name: test_stat, help: Valid tests statistic types are: Score, Wald. Firth can be used with Burden test only. , class: string, optional: true, default: Score"
		this_test_type: "name: test_type, help: Valid tests are one of the collapsing tests SKAT, Burden, SMMAT or Single, class: string, optional: false"
		this_min_mac: "name: min_mac, help: Minimum minor allele count for threshold ( only used for single variant tests ), class: int, optional: true, default: 5"
		this_weights: "name: weights, help: beta weights set to flat weights (e.g. set to 'c(1,25)' for Wu weights or 'c(0.5,0.5)' for Madsen-Browning weights), class: string, optional: true, default: FALSE"
		this_weights_col: "name: weights_col, help: Name of column in aggfile that contains the variant weights, class: string, default: FALSE"
		this_user_cores: "name: user_cores, class: int, optional: false, default: 30"
		this_window: "name: window, help: Runs a sliding window test based on this window size, in bp, class: int, optional: true, default: 0"
		this_step: "name: step, help: For use with 'window', indicates sliding window step size, in bp, class: int, optional: true, default: 0"
		these_genotype_files: "name: genotypefile, class: array[file], patterns: [*.gds, *.GDS], optional: false"
		this_null_model: "name: null_model, class: file, patterns: [*.Rda, *.Rdata], optional: false"
		this_results_file: "name: outputfilename, label: prefix for output file name, no spaces, class: string, optional: false"
		this_memory: "help: memory desired for computation in GB, class: int, optional: false"
		this_disk: "help: disk space desired for computation in GB, class:int, optional: false"
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
				memory = this_memory,
				disk = this_disk

		}
	}

	output {
        Array[File] result = genesis_tests.results
    }
}