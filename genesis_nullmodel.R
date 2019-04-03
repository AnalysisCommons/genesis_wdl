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
test.stat <-  args[9] # Score, Wald, Firth
conditional <- args[10] # 1:237733935:G:A
het_varsIn <-  args[11]



# GLOBAL VARIABLES
collapsing.tests <- c("SKAT",  "Burden")
test.stat.vals <- c("Score", "Wald", "Firth")

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
cat('test.stat',test.stat,'\n')
cat('outcome.type',outcome.type,'\n')
cat('het_vars',het_varsIn,'\n')
cat('conditional',conditional,'\n')


if (!(test.stat %in% test.stat.vals)){
     msg = paste("The requested test statistic:", test.stat, "is not available (Use Firth, Score, Wald!")
     stop(msg)
}

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
# source("/home/dnanexus/pipelineFunctions.R")
source("/genesis_dnanexus/pipelineFunctions.R")
sessionInfo()

if (covariate.string == "NA") {
  covariates <- NULL
} else {
  covariates <- split.by.comma(covariate.string)    
}



## phenotype 
phenotype.data <- read.csv(pheno.file, header=TRUE, as.is=TRUE)

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
##f <- seqOpen('genotypefile')
#AAsample.ids <- seqGetData(f, "sample.id")
all.terms <- unique(c(outcome.name, covariates, het_vars,gcol))
#AApheno <- pheno[row.names(pheno) %in% sample.ids,na.omit(all.terms),drop=F]
cat('Output pheno after mergeing with Genos N=',nrow(pheno),'\n')
if(nrow(pheno) == 0){
    msg = paste("Phenotype ID column doesn't match IDs in GDS")
    stop(msg)
}


#subset to phenotyped samples
#AAseqSetFilter(f,sample.id = row.names(pheno))

# order pheno to the GDS subject order
#AAssample.ids <- seqGetData(f, "sample.id")
#AAspheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]



## Conditional analaysis
if(conditional != 'NA'){
    
    f <- seqOpen(genotype.file)
    sample.ids <- seqGetData(f, "sample.id")
    pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]
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
        seqSetFilter(f,variant.sel=cidx, sample.id = row.names(pheno),verbose=FALSE)
        ncol = NCOL(pheno)
        pheno = cbind(pheno,as.data.frame(altDosage(f,  use.names=FALSE)))
        condheaders  = paste0('csnp',1:NCOND)
        colnames(pheno)[(ncol+1):(ncol+NCOND)] = condheaders
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
if(kinship.matrix != 'NO_KINSHIP_FILE'){
    kmatr = GetKinshipMatrix(kinship.matrix)
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
    ## Get sample ids to check order 
    #AAseqSetFilter(f,sample.id = row.names(pheno))
    #AAsample.ids <- seqGetData(f, "sample.id")
    #AAif (!(identical(sample.ids,row.names(pheno)) && identical(row.names(kmatr),row.names(pheno)))){
    #AA    stop("Something is off problem with re-ordering")
    #AA}
} else {
  kmatr <- NULL
}



###################
## NULL MODEL
##################
print(class(kmatr))

## pheno does not need to be the same order as kmatr, but must have same dims at kmatr
cat('start fit....\n')
if (kinship.matrix == 'NO_KINSHIP_FILE'){
  cat('Fitting unrelated model')
  nullmod <- fitNullModel(pheno,
                        covars = covariates,
                     outcome = outcome.name,
                     family = GetFamilyDistribution(outcome.type))
}else{
#  kmatr = as.matrix(kmatr)
  if(het_varsIn == 'NA'){
      cat('Fitting model with GRM ... \n')
      nullmod <- fitNullModel(pheno,
                        covars = covariates,
                     outcome = outcome.name,
                     family = GetFamilyDistribution(outcome.type),
                     cov.mat = kmatr)
  }else{

      cat('Fitting model with GRM and  het vars ...\n')
      nullmod <- fitNullModel(pheno,
                        covars = covariates,
                     outcome = outcome.name,
                     group.var = het_vars,
                     family = GetFamilyDistribution(outcome.type),
                     cov.mat = kmatr)
      }
}
save(nullmod,pheno,file=results.file)
