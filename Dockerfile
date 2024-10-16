FROM unigebsp/ngs

# Install Shiny server
RUN /rocker_scripts/install_shiny_server.sh

## Install R packages
RUN install2.r --error --skipinstalled -n 4 markdown kableExtra shinyjs openxlsx officer bsicons future
EXPOSE 3838