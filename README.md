# Chemical Genetics Explorer

A shiny app allowing convenient exploration of chemical genetics data

## Setup

This shiny app can be run locally, either in your own `R` environment or via `Docker`. 

### Running via Docker

All you need for this to work is the `Docker` enginge on your system. See [here](https://docs.docker.com/engine/install/) on how to install it.

Once you have it installed, start the app by running this command from the command line

`bash run_app.sh`

### Running via local R installation

If you manage to install all the apps on your local machine (see `packages.r` for required packages and how to install them), you can start the shiny server by running

`Rscript app.r`

from the command line.

## Usage

To use the app, you'll need two files:
- File containing, for each pair of gene and mutant, effect size as well as q-value/FDR values (the last part is only necessary if you want to use empirical bayes-corrected effect size values, default is to use the (scaled) LFC values). The crucial column names are Name (==gene), scaledLFC (==effect size), contrast (==condition)
- File containing functional annotations. The file that contains functional annotations of your choice (for allowed column names see the app upon startup). The gene name in the file above (column Name there) should here be contained in a column called `locus tag`.

Access the shiny app by ging to `http://0.0.0.0:3838` in your browser.

TODO

## Contributions

If you find a bug or have any suggestions, please open up a ticket :)

