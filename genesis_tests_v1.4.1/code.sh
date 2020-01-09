#!/bin/bash

main() {
	echo "App inputs"
    echo "Value of genotypefile: '$genotypefile'"
    echo "Value of outputfilename: '$outputfilename'"

    echo "Value of snpNames: '$snpNames'"
    echo "Value of varaggfile: '$varaggfile'"
    echo "Value of top_maf: '$top_maf'"

    echo "Value of min_mac: '$min_mac'"

    echo "Value of test_type: '$test_type'"
    echo "Value of test_stat: '$test_stat'"
    echo "Value of weights: '$weights'"

    echo "Value of weights_col: '$weights_col'"

    echo "Value of user_cores: '$user_cores'" 
    echo "Value of window: '$window'" 
    echo "Value of step: '$step'" 
    echo "Value of genome_build: '$genome_build'" 
    echo "Value of pass_only: '$pass_only'" 
    echo "Value of imputed: '$imputed'" 
    echo "Value of neig: '$neig'" 
    echo "Value of ntrace: '$ntrace'" 

    echo "Value of interaction: '$interaction'" 
    echo "Value of return_variants: '$return_variants'" 

    dx download "$null_model" -o null_model &
    dx download "$genotypefile" -o genotypefile &
    
    # varaggfile
    if [[ "$varaggfile" != "" ]] ; then
		echo 'downloading varaggfile'
		dx download "$varaggfile"  &
		varaggfile_filename=$( dx describe --name "$varaggfile" )
		invaraggfile=$varaggfile_filename
    else
		invaraggfile="NONE"
    fi
    
    export MKL_NUM_THREADS=1
    export JAVA_HOME=/usr/lib/jvm/default-java


    sudo chmod o+rw /myresources
    # wait if debug 
    if [ ${debug} -ne 0 ]
    then
       echo "DEBUG is on sleeping for ${debug}h"
       sleep ${debug}h
    fi
    wait
 
    echo "Launching R  code -"
    dx-docker run --env KMP_INIT_AT_FORK=FALSE --env MKL_NUM_THREADS=1 -v /home/dnanexus/:/home/dnanexus/ -w /home/dnanexus/ analysiscommon/genesis_wdl:v0.2  Rscript genesis_btest_combine.R  $invaraggfile $top_maf  $test_stat $test_type $min_mac $weights $weights_col $user_cores  $window $step $genome_build $pass_only $imputed $neig $ntrace $interaction $return_variants
    echo "Finished running R code"
    
    echo "Uploading Main Results"
    
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    dx mv ${results} ${outputfilename}.csv.gz
    
    if [[ "$return_variants" == "TRUE" && "$test_type" != "Single" ]] ; then
    	echo "Uploading Agg Test Single Var Results"
    	echo 'Returning Single Var Results from Agg Test'
    	varresults=$(dx upload varresults --brief)
    	dx-jobutil-add-output varresults "$varresults" --class=file
    	dx mv ${varresults} ${outputfilename}.variantInfo.Rda
    fi
    echo 'done'
}
