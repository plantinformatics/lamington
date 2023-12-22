# LAMINGTON
## _PCA & Core Diversity Set Explorer_

The primary objective of this tool is to facilitate the creation of core diversity sets from a genetically characterized population. Leveraging Core Hunter, a widely accepted tool for defining core subsets, this application serves as a wrapper. It introduces a user-friendly graphical interface (GUI) and extends its capabilities to include features like genotype and sample filtering, PCA plot visualization, and outlier detection.

image

## Features

- **Developed in R**: A widely used statistical programming language within the research and breeding community.
- **Flexibility of Deployment**: Designed for both local and server deployment, accommodating the analysis of large datasets.
- **Accessibility**: The source code is open source and available on Github under a GPL license. Docker builds are provided for easy deployment on any operating system.
- **Visualisation Functions**: Users can interactively plot statistics such as minor allele frequency (MAF) and call rate, defining suitable cutoffs.
- **PCA Plotting Functions**: Users can generate PCA plots for the population and color data points based on sample information from an associated spreadsheet.
- **Core Hunter Functions**: Exposes the main functions and parameters of Core Hunter through a graphical interface.
- **Detection of Outliers**: Supports a cyclical workflow for outlier detection. Users can identify outliers, remove them, and rerun the analysis without outliers.
- **Samples to Include/Exclude in Core Set**: Users have the ability to specify samples that should always be included or excluded in the core set, facilitating customization based on existing core sets or sample availability.

## Installation

To run **Lamington** on your local machine, follow these initial setup steps by downloading R or RStudio. Execute the following commands once to configure the environment:
**Install the dependencies**
```sh
install.packages(c('shiny', 'ggplot2', 'shinyFiles', 'dplyr','rjava', 'corehunter','SNPRelate',
                'DT','esquisse','scatterD3','shinycssloaders','shinythemes','rJava','corehunter'))
                
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
```

**You may now run the shiny app with just one command in R:**
```sh
shiny::runGitHub("lamington", "plantinformatics")
```
**Or**
```sh
# First clone the repository with git. If you have cloned it into
# ~/lamington, first change your working directory to ~/lamington, then use runApp('R') to start the app.
setwd("~/lamington") # change to match where you downloaded this repo to
runApp('R') # runs the app 
```

## Docker

**Lamington** is very easy to install and deploy in a Docker container.

By default, the Docker will expose port 8080, so change this within the
Dockerfile if necessary. When ready, simply use the Dockerfile to
build the image.

```sh
cd lamington
docker build -t  lamington  --progress=plain  .
```

This will create the **Lamington** image and pull in the necessary dependencies.
This process might take 10 to 15 minutes.

Once done, run the Docker image, you may change the port mapping from the DockerFile 
and map the port to whatever you wish on your host. In this example, we map port 3838 
of the host to port 3838 of the Docker (or whatever port was exposed in the Dockerfile).
The data folder needs to be also mapped, that contains the saved GDS files.
```sh
docker run -v /lamington/Data:/root/Data  -d -p 127.0.0.1:3838:3838 lamington
```

Verify the deployment by navigating to your server address in
your preferred browser.

```sh
127.0.0.1:3838
```

## Workflow

##### 1. Data Ingestion:
Genotypes (VCF Format) - Accepts genotype data in VCF format. Locally, users can 
load files from the local filesystem. On a server, a designated folder lists VCF
files for user selection.
Sample Passport/Metadata - Users can input metadata/passport data for VCF samples 
in CSV/TSV and/or Excel formats.
##### 2. Filtering Visualisation:
Using the slider users can visualise the change in Missing Rate and 
MAF on the histogram on the right.
##### 3. Genotype Data Filtering:

Prior to analysis, users can filter genotype data based on 
Minor Allele Frequency (MAF) and call rate.
SNP numbers are reducible through Linkage Disequilibrium (LD) pruning.
Users define a list of SNPs for downstream analyses.

##### 4. PCA Calculation:
After defining the set of SNPs, the app calculates Principal Component Analysis (PCA).

##### 5. Calculation of Core Sets:

Users define core sets by utilizing Core Hunter.
Main Core Hunter options are accessible to users, allowing core set 
definition in various sizes simultaneously.
Core sets are marked in the PCA data frame for visualization.

##### 6. Adding Metadata/Passport Data:

Metadata/passport data from the original VCF file, if available,
is added to the PCA data frame for future visualization.

##### 7. PCA Plot Visualization and Outlier Definition:
Core sets from step 5 are visualized within the PCA plot from step 4.
The calculated PCA is visualized with an interactive plot allowing zooming
and sample selection.Sample names are displayed on cursor hover.
After outlier removal, users define core sets by utilizing Core Hunter.
Users can add samples to an exclusion list, regenerating the PCA by revisiting step 3.
Core sets from step 5 are visualized within the PCA plot from step 4.
Core sets are exportable as a csv file containing the list of samples,PCA and 
Population data and the core set.

##### 8. Final Plot
User can use the addin to create final plot and customise it according to their needs. 
Selecting Tab2 from the list of data frames.
## Required Dependencies:
- shiny 
- corehunter
- SNPRelate
- ggplot2
- shinyFiles
- dplyr
- DT
- esquisse 
- scatterD3 
- shinycssloaders 
- shinythemes 

## License
This code is licensed under the GPLv3. Please see the file LICENSE.txt for information.

