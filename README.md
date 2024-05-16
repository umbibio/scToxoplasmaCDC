# scToxoplasmaCDC

## This page provides instructions for running statistical analysis on scRNA-Seq and scATAC-Seq data of replicating tachyzoites. The analysis is as follows:

* Single-cell RNA-seq data processing
* Single-Cell ATAC-seq data processing
* Pseudo-time analysis
* Fitting gene curves(smoothing splines)
* Marker analysis
* CUT&RN data processing
* Reproducing Figures in the manuscript: https://www.biorxiv.org/content/10.1101/2023.10.06.561197v1

## Running the code:

To install the code please clone the library by copy-pasting the following in your terminal:

git clone https://github.com/umbibio/scToxoplasmaCDC.git

## Requirements
The code runs on a standard desktop/laptop computer with enough RMA (>16Gb). The code is tested on Mac and uses standard R packages. 

## Run time
Run time for most parts of the code is very short. Fitting Splines will take up to 1 hour. Multiple cores can be specified for faster run.

## Data availability
The processed data is available [here](https://umbibio.math.umb.edu/toxosc/assets/public-data/preprocessed/rds_ME49_59.zip) and can be used to run the analysis pipeline smoothly and generate the figures in the paper. We suggest you save the rds folder within the same directory where the code is cloned. Alternatively, change the relative path in the code.

Additional downloadable files related to this project are available [here](https://umbibio.math.umb.edu/toxosc/data).
