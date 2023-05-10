FROM rocker/verse:4.0.3

RUN apt-get update

RUN Rscript -e "install.packages(c('optparse', 'jsonlite', 'tidyverse', 'RJSONIO', 'igraph', 'BiocManager'), dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install('edgeR')"

ADD predEngine.R /bin/
RUN chmod +x /bin/predEngine.R

CMD ["/bin/bash"]
