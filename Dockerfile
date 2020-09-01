FROM uwgac/topmed-master:latest

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN sudo apt-get update && sudo apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('gap', 'bdsmatrix', 'doMC', 'qqman', 'RColorBrewer')); BiocManager::install('ramwas')"

RUN wget https://github.com/UW-GAC/GENESIS/archive/v2.15.6.tar.gz
RUN R CMD INSTALL v2.15.6.tar.gz

RUN git clone https://github.com/AnalysisCommons/genesis_wdl.git