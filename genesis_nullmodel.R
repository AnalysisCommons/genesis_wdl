#== Args 
args<-commandArgs(TRUE)
#===mandatory parameters
outcome.name <- args[1]
outcome.type <-  args[2] # Continuous or Dichotomous
covariate.string <- args[3]
pheno.file <- args[4]
genotype.file <- args[5]
results.file <- args[6]

#==optional parameters
kinship.matrix <- args[7]
pheno.id <- args[8]


# added these to JSON
conditional <- args[9] # 1:237733935:G:A
het_varsIn <-  args[10]
transform <-  args[11]
transform.rankNorm <-  args[12]
transform.rescale <-  args[13]



GetFamilyDistribution <- function(response.type) {
               if (response.type == "Continuous"){
                      family = "gaussian"
               } else if (response.type == "Dichotomous"){
                      family = "binomial"
               } else {
                      msg = paste("Don't know how to deal with response type", response.type)
                      stop(msg)
               }
               return(family)
           }

GetKinshipMatrix <- function(kinship.matrix){
  cat('Loading Kinship Matrix:',kinship.matrix,'\n')
  if(grepl('Rda',kinship.matrix,ignore.case=TRUE)){
    kmatr = get(load(kinship.matrix))
  }
  else{
      warning('Reading matrix file from text - will be dense format\n')
    kmatr = as.matrix(read.csv(kinship.matrix,as.is=T,check.names=F,row.names=1))
  }

  cat('Loaded Kinship NROW:',NROW(kmatr),' NCOL:',NCOL(kmatr),'\n')
  kmatr
}

cat('kinship.matrix',kinship.matrix,'\n')
cat('outcome.name',outcome.name,'\n')
cat('covariate.string',covariate.string,'\n')
cat('outcome.type',outcome.type,'\n')
cat('het_vars',het_varsIn,'\n')
cat('conditional',conditional,'\n')
cat('pheno.id',pheno.id,'\n')
cat('transform',transform,'\n')
cat('transform.rankNorm',transform.rankNorm,'\n')
cat('transform.rescale',transform.rescale,'\n')


.libPaths()

suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
suppressMessages(library(gap))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))

## Setup
source("/genesis_wdl/pipelineFunctions.R")
sessionInfo()

if (covariate.string == "NA") {
  covariates <- NULL
} else {
  covariates <- split.by.comma(covariate.string)    
}

cat('covariates:',paste(covariates, collapse = ", "),'\n')
cat('all terms:', paste(unique(c(outcome.name, covariates, het_varsIn)), collapse = ", "), '\n')
 
## phenotype 
phenotype.data <- read.csv(pheno.file, header=TRUE, as.is=TRUE)

if(!transform %in% c('none','transform')){
  stop(paste0('transform must be "none" or "transform", got value ',transform))
}

# add sex if provided, needed for X chrom MAF calcs
if( 'sex' %in% names(phenotype.data)){
   
    gcol='sex'
    # if only contains 'F', sex gets read as T/F
    if (class(phenotype.data$sex) == 'logical') phenotype.data$sex = ifelse(phenotype.data$sex == FALSE,'F',phenotype.data$sex) 
    if(! all(phenotype.data$sex %in% c('F','M'))){
        msg = paste("sex column must be coded M/F")
        warning(msg)
        #stop(msg)
        
        
    }
}else{
    if(! 'sex' %in% names(phenotype.data)) warning("A column labeled 'sex' coded M/F must be included when performing analyses on the sex chromosomes")
    gcol = NULL
} 
if(NCOL(phenotype.data) < 2) warning("Is the phenotype file a CSV?  Too few columns from read.csv()") 



cat('Input pheno N=',nrow(phenotype.data),'\n')
if(het_varsIn != 'NA'){
    cat('prep pheno with het vars')
    het_vars = het_varsIn
    pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates,hetvars=het_vars, id=pheno.id, gender=gcol)
}else{
    cat('prep pheno without het vars\n')
    het_vars = NA
    pheno <- reducePheno(phenotype.data, outcome.name, covariates=covariates, id=pheno.id, gender=gcol)
}
cat('Output pheno N=',nrow(pheno),'\n')

## Report dropped individuals
dropped.ids.selector <- !(phenotype.data[[pheno.id]] %in% row.names(pheno))
dropped.ids <- phenotype.data[[pheno.id]][dropped.ids.selector] 
if (NROW(dropped.ids) != 0 ) {
  cat("Dropped because of incomplete cases:", length(dropped.ids) )
}

# For GDS files
all.terms <- unique(c(outcome.name, covariates, het_vars,gcol))
cat('Output pheno after mergeing with Genos N=',nrow(pheno),'\n')
if(nrow(pheno) == 0){
    msg = paste("Phenotype ID column doesn't match IDs in GDS")
    stop(msg)
}

## Conditional analaysis
if(conditional != 'NA'){
    
    f <- seqOpen(genotype.file)
    sample.ids <- seqGetData(f, "sample.id")
    # pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]
    cond_snps = strsplit(conditional,',')[[1]] 
    cond_snps = gsub(' ','',cond_snps)
    cond_snps = gsub('"','',cond_snps)

   
    chr_array = seqGetData(f, "chromosome")
    pos_array = seqGetData(f, "position")
    allele_array = seqGetData(f, "allele")
    allele_array = gsub(',',':',allele_array)
    snp_array = paste(chr_array,pos_array,allele_array,sep=':')

    
    cat('Conditioning on ',conditional,'...\n')
    cidx = which(snp_array %in% cond_snps)
    NCOND = length(cidx)
    cat('Found ',snp_array[cidx], ' in GDS file\n')
    if(NCOND < length(cond_snps)){
        warning('NOT ALL CONDITIONAL SNPS FOUND IN GDS')
    }

    
    if(any(cidx)){
        seqSetFilter(f,variant.sel=cidx, sample.id = row.names(pheno),verbose=TRUE)
        ncol = NCOL(pheno)
        cdoses = as.data.frame(altDosage(f,  use.names=TRUE)) 
        condheaders = paste0('csnp',1:NCOND)
        names(cdoses)  = condheaders
        pheno = merge(pheno,cdoses,by='row.names')
        row.names(pheno) = pheno$Row.names
        pheno = pheno[,names(pheno) != 'Row.names']
    }else{
        stop('Can not find snp ',conditional,' with position ',cpos,' to condition on in data file')
    }
  
    dropConditionalCases = NROW(pheno)-NROW(pheno[complete.cases(pheno),])
    if(dropConditionalCases > 0){
        cat('Warning: Dropping ',dropConditionalCases,' samples due to missing conditional genotype calls\n')
    }

    pheno = pheno[complete.cases(pheno),]
  
    covariates[(length(covariates) + 1):(length(covariates) +NCOND)] <- condheaders
    seqClose(f)

}




## Load KINSHIP matrix
## Kinship doesn't contain all samples
kinship.files = list()
if(kinship.matrix != 'NO_KINSHIP_FILE'){
    i = 0
    for(kfile in kinship.matrix){
        i = i+1
        kmatr = GetKinshipMatrix(kfile)
        pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
        kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
        cat('Output pheno in Kinship N=',nrow(pheno),'\n')
        kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]
        if(nrow(pheno) == 0){
            msg = paste("Phenotype ID column doesn't match IDs in Kinship Matrix")
            stop(msg)
        }
        if (!(identical(row.names(kmatr),row.names(pheno)))){
            stop("Something is off problem with re-ordering")
        }
        kinship.files[[i]]  =  kmatr
    }
    
} else {
  kmatr <- NULL
}



###################
## NULL MODEL
##################

## pheno does not need to be the same order as kmatr, but must have same dims at kmatr
cat('start fit....\n')

if (kinship.matrix == 'NO_KINSHIP_FILE'){
    kinship.files = NULL
}
  
if(het_varsIn == 'NA'){
    het_vars = NULL
}

cat('Fitting model with GRM and  het vars ...\n')
nullmod <- fitNullModel(pheno,
                        covars = covariates,
                        outcome = outcome.name,
                        group.var = het_vars,
                        family = GetFamilyDistribution(outcome.type),
                        cov.mat = kinship.files)
if(transform != 'none'){
  cat('Updating Null Model ...\n')
  nullmod <- nullModelInvNorm(nullmod,
                        cov.mat = kinship.files,norm.option = transform.rankNorm, rescale = transform.rescale)

}
save(nullmod,pheno,file= paste0(results.file, '.RData'))
