FROM rocker/shiny-verse:4.4.2
LABEL maintainer="Muhammad T. Suhail"

RUN apt-get update && apt-get install -y --no-install-recommends \
     sudo \
	 libcurl4-gnutls-dev \
  	 libcairo2-dev \
     libxt-dev \
     libssl-dev \
     libssh2-1-dev \
     && rm -rf /var/lib/apt/lists/*

	RUN apt-get -y update && apt-get install -y \
    default-jdk \
    r-cran-rjava \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/

# Install necessary R packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rhandsontable', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('rJava','corehunter'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('ggplot2', 'shinyFiles', 'DT','esquisse','scatterD3','shinycssloaders','shinythemes','plotly','shinyWidgets','bcrypt','RSQLite','DBI','shinyjs','RColorBrewer','shinyalert'), repos='http://cran.rstudio.com/')"

# Install BiocManager and then install SNPRelate
RUN R -e "BiocManager::install('SNPRelate')"

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site


RUN mkdir /root/lamington 
RUN mkdir /root/GDS
RUN mkdir /root/VCFs
RUN mkdir /root/DB

COPY R /root/lamington

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/root/lamington',host = '0.0.0.0', port=3838)"]
