FROM unigebsp/ngs

## Install R packages
RUN install2.r --error --skipinstalled -n 4 markdown kableExtra shinyjs openxlsx officer

