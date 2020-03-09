###== Args 
args<-commandArgs(TRUE)

agg.file <- args[1] 
top.maf <- as.numeric(args[2]) 
test.stat <-  args[3] # Score, Wald
test.type  <-  args[4] # Burden, Single, SKAT, SMMAT

min.mac <- as.integer(args[5])
weights <- args[6]
weights.col <- args[7]
user_cores <-  as.numeric(args[8])
window = as.numeric(args[9])
step= as.numeric(args[10])
genotype.file <- args[11]
null.model <- args[12]
results.file <- args[13]

genome_build = args[14]
pass_only = ifelse(tolower(args[15]) == "F", F, T)
imputed = ifelse(tolower(args[16]) == "F", F, T)
neig = as.numeric(args[17])
ntrace = as.numeric(args[18])
interaction = args[19]
# interaction = ifelse(args[19] == "NA", NULL, args[19])
return_variants = ifelse(tolower(args[20]) == "F", F, T)


if (toupper(interaction) == 'NULL') interaction = NULL
# GLOBAL VARIABLES
collapsing.tests <- c("Burden", "SKAT", "fastSKAT", "SMMAT",  "SKATO")
test.type.vals <- c("Single",collapsing.tests)
test.stat.vals <- c("Score",  "Score.SPA")
genome_build.vals <- c("hg19",  "hg38")


cat("R script inputs\n")
cat('varaggfile',agg.file,'\n')
cat('top.maf',top.maf,'\n')
cat('test.stat',test.stat,'\n')
cat('test.type',test.type,'\n')
cat('weights.col',weights.col,'\n')
cat('window',window,'\n')
cat('step',step,'\n')
cat('genome_build ',genome_build,'\n')
cat('pass_only ',pass_only,'\n')
cat('imputed ',imputed,'\n')
cat('interaction ',!is.null(interaction), interaction,'\n')

cat('neig ',neig,'\n')
cat('ntrace ', ntrace,'\n')

if (!(test.stat %in% test.stat.vals)){
     msg = paste("The requested test statistic:", test.stat, "is not available Use Single, Burden, SKAT, fastSKAT, SMMAT or SKATO")
     stop(msg)
}

if(!test.type %in% test.type.vals){
    stop("Test type must be one of ",paste(test.type.vals,sep=','))
}

if(!genome_build %in% genome_build.vals){
    stop("genome_build must be one of ",paste(genome_build.vals,sep=','))
}


if (!is.null(interaction) & test.type %in% collapsing.tests){
    stop("Interaction must be FALSE for aggregate tests")

}

weights = eval(parse(text=weights))
cat('Weights',weights,class(weights),'\t',length(weights),'\n')
if(test.type %in% collapsing.tests){
	if(weights && weights.col != 'FALSE'){
		stop('You cannot use both weights ( beta function ) and weights_col ( user defined weights )')
	}
}
SINGLE_T = F # single var test
SW_T = F  # sliding window test
AGG_T = F # aggregate test defined by agg file
FILTER = F # aff file used for filtering.  Can filter both SW and Single as well
USE_AGG = F # use agg file for aggregation, filtering or weighting

if(agg.file != 'NONE') USE_AGG = T

if(!test.type %in% collapsing.tests){
	print("Selected Single Variant Test")
	SINGLE_T = T
	if(USE_AGG) FILTER = T
}else if(window > 0){
	warning("Selected Sliding Window Test")
	SW_T = T
	if(USE_AGG) FILTER = T
}else{

	print("Selected Aggregate File-defined Test")
	AGG_T = T
	if(!USE_AGG){
		warn("Aggregate testing requires an aggregation file, unless running a sliding window")
	}
}

if(return_variants & SINGLE_T){
	warning("'return_variants' is only valid for aggregate tests. Setting to FALSE")
	return_variants = F
}

if(length(weights) ==  1 & !SINGLE_T){
	weights = c(1,1)
	warning("No weights defined for aggregate tests, setting to flat weights")
}

suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
suppressMessages(library(gap))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(parallel))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))
suppressMessages(library(doMC))
suppressMessages(library(TopmedPipeline))
suppressMessages(library(survey))
suppressMessages(library(CompQuadForm))
suppressMessages(library(GenomicRanges))

cat("\n\nR Session Info:\n")
sessionInfo()

num_cores <- detectCores(logical=TRUE)
f <- seqOpen(genotype.file)

# get all samples in GDS
full.sample.ids <- seqGetData(f, "sample.id")
cat('\nload null model')

# loads null model and annotated data frame
load(null.model)
# get samples included in null model
nullmod$sample.id = row.names(nullmod$model.matrix)
sample_ids <- row.names(nullmod$model.matrix)

if(!exists("annotated.frame")){
    pheno$sample.id = row.names(pheno)
    annotated.frame = AnnotatedDataFrame(pheno)
    modified.pheno = pheno[full.sample.ids,,drop=FALSE]
    sample.data.for.annotated <- data.frame(sample.id = full.sample.ids, modified.pheno,stringsAsFactors=F)
    annotated.frame <- AnnotatedDataFrame(sample.data.for.annotated)
}

## Setup
source("/genesis_wdl/pipelineFunctions.R")
  
svd <- SeqVarData(f, sampleData=annotated.frame)
vi <- variantInfo(svd, alleles = FALSE, expanded=FALSE)
chr = vi$chr[1]
if(length(unique(chr)) > 1) stop(" Not expecting multiple chromosomes\n")
chr = chr[1] 

	
#
# chunks for splitting over cores
#
gr <- get(data(list = paste("chromosomes", genome_build, sep = "_"), 
        				package = "TopmedPipeline", envir = environment()))

gr = gr[seqnames(gr) == chr]    


nchunk = (user_cores*6)				
gr$seg.length <- ceiling(max(vi$pos)/(nchunk))
        			
dat = do.call(c, lapply(seq_along(gr), function(i) {
       					x <- gr[i]
        				window.start <- seq(BiocGenerics::start(x), BiocGenerics::end(x), 
            				x$seg.length)
        				GRanges(seqnames = seqnames(x), IRanges(window.start, 
            				width = x$seg.length))
}))
nchunk  = length(dat)

### Load Agg/Filter File
    
if(USE_AGG){
	if(grepl('\\.rda$',tolower(agg.file)) | grepl('\\.rdata$',tolower(agg.file)) ){
		agg = get(load(agg.file))
	}else{
            system.time({ agg = fread(agg.file, stringsAsFactors=F, sep=',', header=T,data.table=F) })
        }
    
	names(agg)[names(agg) %in%  c('CHR','CHROM','Chr','chrom')] = 'chr'	
        agg$chr = sub('chr','',agg$chr)
	names(agg) = tolower(names(agg))
	
	
	# merge in the weights, then pull back out with their variant.ids
	if(weights.col != 'FALSE'){
            wfile = agg[,c('chr','pos','ref','alt',weights.col)]
            names(wfile) = c('chr','pos','ref','alt','weights')
            wfile = unique(wfile)
    	
            svd <- addVariantData(svd, wfile)
            wdf = svd@variantData@data[,c('variant.id','weights')]
	}

	if(AGG_T){
		chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
                cat('N AggUnits=',length(unique(agg$group_id)),'done\n')
    	
                ## Turn Agg file into granges list
                aggVarList <- aggregateGRangesList(agg)
        }
        if(FILTER){
            filter_gr <- GRanges(seqnames = agg$chr, ranges = IRanges(start = agg$pos, 
                                                                      width = 1))
        }
}

seqClose(f)


doOne = function(idx,in_nullmod) {
				print(paste0(idx,' iter'))
	##############
				 ## Apply variant filters to get a list of variants in each gene
	###############
				 ## extract genotypes
				 f <- seqOpen(genotype.file)
				 svd = SeqVarData(f, sampleData=annotated.frame)
			
				 seqSetFilter(svd,sample.id = sample_ids,verbose=F)	
				 
				 # Sliding Window
				 if(SW_T){
				  	# Add a window width to make sure we get all possibile windows
    				dat <- resize(dat, width(dat) + window, fix = "start")
            		seg <- dat[idx]
    				seqSetFilter(svd, variant.sel = seg, verbose = TRUE)
    				if(FILTER) seqSetFilter(svd, variant.sel = filter_gr, verbose = TRUE, intersect=T)
    				
				 # Aggregating File 
   				 }else if(AGG_T){
   				 	if( length(names(aggVarList)) < idx) return(data.frame())
				 	avl = aggVarList[names(aggVarList) %in% chunk(names(aggVarList),nchunk)[[idx]]]
    				seqSetFilter(svd, variant.sel = avl, verbose = TRUE)
				 	
				 # Single Variant
				 }else if(SINGLE_T){
				 	seg <- dat[idx]
    				seqSetFilter(svd, variant.sel = seg, verbose = TRUE)
    				if(FILTER) seqSetFilter(svd, variant.sel = filter_gr, verbose = TRUE, intersect=T)
    				
					via <- variantInfo(svd, alleles = TRUE, expanded=FALSE)
					
				 }else{
				 	print('Warning -- I was not expecting to be here')
				 }
				 
				 
				 if(length(seqGetData(svd, "variant.id"))>0 & top.maf < 1){
				 	filterByRare(svd, sample_ids, af.max=top.maf,verbose=T)
				 }
				 
				 if(length(seqGetData(svd, "variant.id"))>0){
					## filter to maf and remove monomorphic or larger, if min.mac is set
					filterByMAC(svd, sample_ids,  mac.min=max(min.mac,1), build=genome_build,verbose=TRUE)
				 }

				 if (as.logical(pass_only)) {
					 print('Filtering to PASS variants only')
			    	 filterByPass(svd,verbose=TRUE)
				 }
			 
			 
			 
				 if(length(seqGetData(svd, "variant.id"))>0){
   
		if(SW_T){
			iterator <- SeqVarWindowIterator(svd, windowSize=window, windowShift=step,verbose = T)
		}else if(AGG_T){
			iterator <- SeqVarListIterator(svd, avl,verbose=T)
		}else if(SINGLE_T){
			iterator <- SeqVarBlockIterator(svd,variantBlock=4000,verbose = T)
		}
  
  
		## Collapse test

 		if (SW_T || AGG_T) {

			## beta function weights
			res <- assocTestAggregate( iterator, 
				in_nullmod, 
				weight.beta = weights,
				test=test.type,   genome.build = genome_build, imputed=imputed, verbose=T)
			if(SW_T){
				generes <- cbind(data.frame(gene=paste(res$results$chr,res$results$start,res$results$end,sep='_') ,  res$results,stringsAsFactors=F))
				res$results$windowname = paste(res$results$chr,res$results$start,res$results$end,sep='_')
				names(res$variantInfo) = paste(res$results$chr,res$results$start,res$results$end,sep='_')
 			}else{
 				generes <- cbind(data.frame(gene=row.names(res$results) ,  res$results,stringsAsFactors=F))
			}
 
 
		}else {  # Single variant test

 			generes <- assocTestSingle(iterator,  test =  test.stat, in_nullmod, genome.build = genome_build,imputed=imputed,GxE=interaction,verbose=T)
			#via <- variantInfo(svd, alleles = TRUE, expanded=FALSE)
			generes = merge(generes,via[,c('variant.id','ref','alt')],by='variant.id')
			cat(paste('Ran assoc',NROW(generes),'\n'))
			generes$snpID = paste0(generes$chr,':',generes$pos,':',generes$ref,':',generes$alt)  

		} # end single
	} # end is >0 variants to test
 
	seqClose(f)
 	if(!exists("generes")){
		generes= data.frame();
 		if (return_variants) res = list('results'=data.frame(),'variantInfo'=list())
 	}
 	if (return_variants) {
  		res
 	}else{
  		generes
  	}
}

combineAggRes = function(ares, varres.file){
	print('Combining Agg single vars')
	genereslist = list()
	varreslist = list()
	for(i in 1:length(ares)){
		genereslist[[i]] = ares[[i]]$result
		varreslist[[i]] = ares[[i]]$variantInfo
	}

	varres = do.call(c,varreslist)
	generes = as.data.frame(do.call(rbind,genereslist))

	# out <- gzfile('results', "w")
	# write.csv(generes, out, row.names=F)
	# close(out)
	save(varres,file=varres.file)
	return(generes)
}

runMainAnalysis = function(user_cores,in_nullmod,keepVars=NULL, keepGenes=NULL){
	registerDoMC(cores=min(c(user_cores,num_cores)))
	print('Start Env')
	cat('\n\n#########################\nRunning Analysis with ',min(c(as.numeric(user_cores),num_cores)),' cores of ',num_cores,'\n#########################\n\n')
	cat(paste('Running  - ',length(dat),' Units\n'))

	mcoptions <- list(preschedule=FALSE, set.seed=TRUE) #in_nullmod, 


	if (return_variants) {
		print('return vars - ')
		sm_obj <- 
			foreach (idx=1:length(dat),
 			.combine=function(...){c(list(...))},
 			.multicombine = TRUE,
			.inorder=FALSE,  .verbose = FALSE,
 			.options.multicore=mcoptions) %dopar% doOne(idx,in_nullmod)
 	}else{
		sm_obj <- 
			foreach (idx=1:length(dat),
 			.combine=function(...){rbindlist(list(...),fill=T)},
 			.multicombine = TRUE,
 			.inorder=FALSE,  .verbose = FALSE,
 			.options.multicore=mcoptions) %dopar% doOne(idx,in_nullmod)
	}
	cat("\nFinished Loop\n")
	sm_obj
}

print('Starting main analysis')
results = runMainAnalysis(user_cores=user_cores,in_nullmod=nullmod)

results.file <- sub(".csv", paste0(".", chr, ".csv"), results.file)
varresults.file <- sub(".csv", paste0(".", chr, ".Rdata"), results.file)
if (!endsWith(results.file, ".gz")) results.file <- paste0(results.file, ".gz")
out <- gzfile(results.file, "w")

if(return_variants){
	print('Returning single variants from aggregate tests...')
	comb.results <- combineAggRes(results, varresults.file)
	write.csv(comb.results, out, row.names=F)
	close(out)
}else{
	results = as.data.frame(results)
	write.csv(results, out, row.names=F)
	close(out)
	varres.save <- "Variant results not requested."
	save(varres.save, file = varresults.file)
}
