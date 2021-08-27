FROM uwgac/topmed-master:2.10.0

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN sudo apt-get update && sudo apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# RUN sudo Rscript -e "packageurl<-'https://cran.r-project.org/src/contrib/Archive/bdsmatrix/bdsmatrix_1.3-3.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"

# RUN sudo Rscript -e "packageurl<-'https://cran.r-project.org/src/contrib/Archive/gap/gap_1.2.1.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"

RUN sudo Rscript -e "install.packages(c('dplyr', 'dbplyr', 'doMC', 'qqman', 'RColorBrewer', 'foreach', 'gap', 'bdsmatrix')); BiocManager::install(c('biomaRt','ramwas'))"

RUN cd /usr/local

Run echo git

RUN sudo git clone https://github.com/AnalysisCommons/genesis_wdl.git

RUN cd genesis_wdl && sudo git checkout v1_5 && sudo git pull origin v1_5

RUN sudo rm -rf /usr/local/analysis_pipeline/testdata/

