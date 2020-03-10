FROM uwgac/topmed-master:latest

MAINTAINER Tim Majarian <tmajaria@broadinstitute.org>

RUN apt-get update && apt-get -y install git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('gap', 'bdsmatrix', 'doMC', 'qqman', 'RColorBrewer', 'cli', 'gh', 'usethis', 'devtools')); devtools::install_github('andreyshabalin/ramwas')"

RUN git clone https://github.com/AnalysisCommons/genesis_wdl.git

RUN cd genesis_wdl && git checkout v1_4_1  && git pull origin v1_4_1