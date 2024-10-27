FROM rocker/shiny-verse
# RUN R -q -e 'install.packages(c("BiocManager", "GGally", "circlize", "GetoptLong", "shinydashboard"), lib="/usr/local/lib/R/site-library/")'
# RUN R -q -e 'BiocManager::install(c("ComplexHeatmap", "InteractiveComplexHeatmap"), lib="/usr/local/lib/R/site-library/")'
# RUN R -q -e 'install.packages(c("here"), lib="/usr/local/lib/R/site-library/")'
# RUN R -q -e 'install.packages(c("DT"), lib="/usr/local/lib/R/site-library/")'
RUN Rscript packages.r
EXPOSE 3838
ADD chem_gen_explorer .
CMD ["R", "-e", "shiny::runApp('app.r', host = '0.0.0.0', port = 3838)"]