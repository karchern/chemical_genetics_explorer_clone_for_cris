FROM rocker/shiny-verse
COPY ./packages.R packages.R
RUN Rscript packages.r
EXPOSE 3838
ADD chem_gen_explorer .
CMD ["R", "-e", "shiny::runApp('app.r', host = '0.0.0.0', port = 3838)"]