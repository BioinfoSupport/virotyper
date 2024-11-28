FROM unigebsp/ngs

# Install Shiny server
RUN /rocker_scripts/install_shiny_server.sh

## Install R packages
RUN install2.r --error --skipinstalled -n 4 markdown kableExtra shinyjs openxlsx officer bsicons future

COPY --chmod=755 bin /srv/shiny-server/bin
COPY --chmod=755 notebooks /srv/shiny-server/notebooks
COPY --chmod=755 data/db/ /srv/shiny-server/data/db/
COPY --chmod=755 Makefile /srv/shiny-server/Makefile
COPY --chmod=755 app.R /srv/shiny-server/app.R

EXPOSE 3838
