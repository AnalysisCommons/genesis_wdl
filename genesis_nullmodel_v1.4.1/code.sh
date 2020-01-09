#!/bin/bash

main() {

    echo "Value of phenofile: '$phenofile'"
    echo "Value of snpinfofile: '$snpinfofile'"
    echo "Value of genotypefile: '$genotypefile'"
    echo "Value of outputfilename: '$outputfilename'"

    echo "Value of kinshipmatrix: '$kinshipmatrix'"
    echo "Value of pheno_id: '$pheno_id'"
    echo "Value of snpNames: '$snpNames'"



    covariate_list=$(echo ${covariate_list} | sed 's/ //g') 
    echo "Value of covariate_list: '$covariate_list'"
    echo "value of conditional SNP: '$conditional'"    
    echo "Value of het_vars: '$het_vars'"


    dx download "$phenofile" -o phenofile &
    dx download "$genotypefile" -o genotypefile &

    if [[ "$kinshipmatrix" != "" ]] ; then
		echo 'downloading kinshipmatrix (s)'
	

		# Download input files and check in directory
		echo "dx-download-all-inputs"
		dx-download-all-inputs
		ls -sh in/kinshipmatrix/*/*

		# Check filetype extension inside first folder of in/datafile
		if [[ "${#kinshipmatrix[@]}" -lt "10" ]] ; then
    	    z="0"
		else 
            z="00"
		fi
		file0=( in/kinshipmatrix/$z/* )
    
		wait
    
		kinfiles=( in/kinshipmatrix/*/* )
		inkinshipfile="kinship=${kinfiles[@]}"

    else
		inkinshipfile="kinship=NO_KINSHIP_FILE"
    fi

    
    
    # wait if debug 
    if [ ${debug} -ne 0 ]
    then
       echo "DEBUG is on sleeping for ${debug}h"
       sleep ${debug}h
    fi
    wait


    
    NCORE=`getconf _NPROCESSORS_ONLN`
    export MKL_NUM_THREADS=$NCORE


    
    echo "Checking phenofile" 
    if [ -e phenofile ] 
    then
       head -n1 phenofile
    else
       echo "The phenofile is not ready"
    fi


    
    echo "Running code"
    dx-docker run --env MKL_NUM_THREADS=$NCORE -v /home/dnanexus/:/home/dnanexus/ -w /home/dnanexus/ analysiscommon/genesis_wdl:v0.2  Rscript genesis_nullmodel.R  --args $outcome_name $outcome_type "$covariate_list" $pheno_id $conditional $het_vars $transform $transform_rankNorm $transform_rescale "$inkinshipfile"
    echo "Finished running code"
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    dx mv ${results} ${outputfilename}.Rda
}
