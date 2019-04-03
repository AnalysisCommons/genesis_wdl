task genesis_nullmodel {
	String outcome_name
	String? outcome_type
	String? covariates_string
	File pheno_file
	File? genotype_file
	String results_file
	File? kinship_matrix
	String? pheno_id
	String? test_stat
	String? conditional
	String? het_varsIn

	Int memory
	Int disk

	command {
		R --vanilla --args ${outcome_name} ${default="Continuous" outcome_type} ${default="NA" covariates_string} ${pheno_file} ${default="NA" genotype_file} ${results_file} ${default="NO_KINSHIP_FILE" kinship_matrix} ${default="ID" pheno_id} ${default="Score" test_stat} ${default="NA" conditional} ${default="NA" het_varsIn} < /genesis_dnanexus/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/genesis_dnanexus:v0.1"
		disks: "local-disk ${disk} HDD"
		memory: "${memory} GB"
	}

	output {
		File results = "${results_file}"
	}
}

workflow genesis_nullmodel_wf {
	String this_outcome_name
	String? this_outcome_type
	String? this_covariates_string
	File this_pheno_file
	File? this_genotype_file
	String this_results_file
	File? this_kinship_matrix
	String? this_pheno_id
	String? this_test_stat
	String? this_conditional
	String? this_het_varsIn

	Int this_memory
	Int this_disk

	# Workflow metadata
	meta {
		summary: "Runs the null model for use in the genesis app"
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
		this_genotype_file: "name: genotypefile, label: Genotypes, help: Use for providing genotypes for conditional analyses, class: file, optional: true, default: NA"
		this_results_file: "name: outputfilename, label: prefix for output file name, no spaces, class: string, optional: false"
		this_kinship_matrix: "name: kinshipmatrix, label: kinship matrix with sample ids as the row and column names.  Matricies saved as Rda will load faster, but csv is accepted as well. Rda files should contain a single numeric matrix object., class: file,patterns: [*.Rda, *.csv], optional: true"
		this_pheno_id: "name: pheno_id, help: Column name that contains the sample IDs.  These IDs should match the genotype file IDs and the kinship file IDs., class: string, default: ID"
		this_test_stat: "name: test_stat, help: Valid tests statistic types are: Score, Wald. Firth can be used with Burden test only. , class: string, optional: true, default: Score"
		this_conditional: "name: conditional, help: chr pos ref alt format for the SNP that will be added to the model.  Multiple snps in a comma delimited list can be added. (e.g. '22:16425814:C:T' or '22:16425814:C:T,22:17808063:TA:T,22:18096610:G:T'), class: string, optional: true, default: NA"
		this_het_varsIn: "name: het_vars, help: grouping variable for heterogenous variances, class: string, optional: true, default: NA"
		this_memory: "help: memory desired for computation in GB, class: int, optional: false"
		this_disk: "help: disk space desired for computation in GB, class:int, optional: false"
	}

	call genesis_nullmodel {
		input: 
			outcome_name = this_outcome_name,
			outcome_type = this_outcome_type,
			covariates_string = this_covariates_string,
			pheno_file = this_pheno_file,
			genotype_file = this_genotype_file,
			results_file = this_results_file,
			kinship_matrix = this_kinship_matrix,
			pheno_id = this_pheno_id,
			test_stat = this_test_stat,
			conditional = this_conditional,
			het_varsIn = this_het_varsIn,
			memory = this_memory,
			disk = this_disk
	}

	output {
        File null_model = genesis_nullmodel.results
    }
}