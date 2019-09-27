FROM uwgac/topmed-master:latest

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN apt-get update && apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('gap', 'bdsmatrix', 'doMC', 'qqman', 'RColorBrewer')); BiocManager::install('ramwas')"

RUN git clone https://github.com/AnalysisCommons/genesis_wdl.git