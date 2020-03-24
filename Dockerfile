FROM uwgac/topmed-master:latest

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN apt-get update && apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN Rscript -e "packageurl<-'https://cran.r-project.org/src/contrib/Archive/bdsmatrix/bdsmatrix_1.3-3.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"

RUN Rscript -e "packageurl<-'https://cran.r-project.org/src/contrib/Archive/gap/gap_1.2.1.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"

RUN Rscript -e "install.packages(c('doMC', 'qqman', 'RColorBrewer', 'foreach')); BiocManager::install('ramwas')"

RUN git clone https://github.com/AnalysisCommons/genesis_wdl.git

RUN cd genesis_wdl && git checkout v1_4_1  && git pull origin v1_4_1

RUN rm -rf /usr/local/analysis_pipeline/testdata/