FROM uwgac/topmed-master:2.10.0

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN apt-get update && apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# RUN sudo Rscript -e "packageurl<-'https://cran.r-project.org/src/contrib/Archive/bdsmatrix/bdsmatrix_1.3-3.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"

# RUN sudo Rscript -e "packageurl<-'https://cran.r-project.org/src/contrib/Archive/gap/gap_1.2.1.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"

RUN sudo Rscript -e "install.packages(c('dplyr', 'dbplyr', 'doMC', 'qqman', 'RColorBrewer', 'foreach')); BiocManager::install(c('biomaRt','ramwas'))"

RUN sudo git clone https://github.com/AnalysisCommons/genesis_wdl.git

RUN sudo cd genesis_wdl && git checkout v1_4_1 && git pull origin v1_4_1

RUN sudo rm -rf /usr/local/analysis_pipeline/testdata/