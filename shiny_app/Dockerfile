# easyshiny

# Base image
FROM rocker/shiny:4.4.1

# General updates
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev libssl-dev libharfbuzz-dev libfribidi-dev libhdf5-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install the required packages
RUN Rscript -e 'reqPkg = c("bslib","data.table", "DT", "ggdendro", "ggplot2", "ggplotify", "ggrepel", "glue", "grid", "hdf5r", "magrittr", "Matrix", "patchwork", "RColorBrewer", "readr", "remotes", "reticulate", "Seurat", "shiny", "shinycssloaders", "shinyhelper", "showtext"); newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]; if(length(newPkg)){install.packages(newPkg)}'
RUN Rscript -e 'install.packages("DT")'

# Copy the app files (scripts, data, etc.)
RUN rm -rf /srv/shiny-server/*
COPY . /srv/shiny-server/

# Ensure that the expected user is present in the container
RUN if id shiny &>/dev/null && [ "$(id -u shiny)" -ne 999 ]; then \
        userdel -r shiny; \
        id -u 999 &>/dev/null && userdel -r "$(id -un 999)"; \
    fi; \
    useradd -u 999 -m -s /bin/bash shiny; \
    chown -R shiny:shiny /srv/shiny-server/ /var/lib/shiny-server/ /var/log/shiny-server/

# Other settings
USER shiny
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
