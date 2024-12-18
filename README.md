# LAMINGTON

## *PCA & Core Diversity Set Explorer*

The primary objective of this tool is to facilitate the creation of core
diversity sets from a genetically characterized population. Leveraging
Core Hunter, a widely accepted tool for defining core subsets, this
application serves as a wrapper. It introduces a user-friendly graphical
interface (GUI) and extends its capabilities to include features like
genotype and sample filtering, PCA plot visualization, and outlier
detection.

![image](https://github.com/user-attachments/assets/aefd2946-51ce-4040-8733-03fb027b298f)

The initial development of Lamington was carried out by La Trobe
University student Muhammad Tahaa Suhail as part of a Work-based
Learning placement with Agriculture Victoria in the context of the
Australian Grains Genebank Strategic Partnership, a \$30m 5-year joint
investment between the Victorian State Government and Grains Research
and Development Corporation (GRDC) that aims to unlock the genetic
potential of plant genetic resources for the benefit of Australian grain
growers.

## Features

-   **Developed in R**: A widely used statistical programming language
    within the research and breeding community.
-   **Flexibility of Deployment**: Designed for both local and server
    deployment, accommodating the analysis of large datasets.
-   **Accessibility**: The source code is open source and available on
    Github under a GPL license. Docker builds are provided for easy
    deployment on any operating system.
-   **Visualisation Functions**: Users can interactively plot statistics
    such as minor allele frequency (MAF) and call rate, defining
    suitable cutoffs.
-   **PCA Plotting Functions**: Users can generate PCA plots for the
    population and color data points based on sample information from an
    associated spreadsheet.
-   **Core Hunter Functions**: Exposes the main functions and parameters
    of Core Hunter through a graphical interface.
-   **Detection of Outliers**: Supports a cyclical workflow for outlier
    detection. Users can identify outliers, remove them, and rerun the
    analysis without outliers.
-   **Samples to Include/Exclude in Core Set**: Users have the ability
    to specify samples that should always be included or excluded in the
    core set, facilitating customization based on existing core sets or
    sample availability.

## Installation

To run **Lamington** on your local machine, follow these initial setup
steps by downloading R or RStudio. Execute the following commands once
to configure the environment:

**Install the dependencies**

``` sh
install.packages(c('shiny', 'ggplot2', 'shinyFiles', 'dplyr','rjava', 'corehunter','SNPRelate',
                'DT','esquisse','scatterD3','shinycssloaders','shinythemes','rJava','corehunter'))
                
if (!requireNamespace("BiocManager",quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
```

**You may now run the shiny app with just one command in R:**

``` sh
#This command runs Lamington from Github without explicitly cloning lamington
shiny::runGitHub("lamington", "plantinformatics",ref="main",subdir = "R")
```

**Or**

``` sh
# First clone the repository with git.
git clone git@github.com:plantinformatics/lamington.git
setwd("~/lamington") # change to match where you downloaded this repo to
shiny::runApp('R') # runs the app 
```

## Docker

**Lamington** is very easy to install and deploy in a Docker container.

A Docker container for Lamington can be built using the
[*Dockerfile*](Dockerfile). This file contains all the instructions
required to assemble the container image, including the installation of
necessary R packages, system dependencies, and the Lamington software
itself.

The Docker image is configured to expose port 3838 by default. To change
this, edit the `EXPOSE` instruction in the [*Dockerfile*](Dockerfile).
Once you have configured the [*Dockerfile*](Dockerfile) as needed, build
the image using the docker build command.

``` sh
#Clone lamington if not already done so
git clone git@github.com:plantinformatics/lamington.git
cd lamington
docker build -t  lamington  --progress=plain  .
```

Building the **Lamington** Docker image will download all required
dependencies. This process may take 10 to 15 minutes.

After the build completes, run the Docker image. You can map the
container's exposed port (defined in the [*Dockerfile*](Dockerfile)) to
a port on your host machine. For example, to map host port 3838 to the
container's port 3838, use the appropriate docker run command.

To enable data access, mount the VCF directory (for direct uploads) and
the GDS directory (containing saved GDS files) to corresponding
locations within the container.

``` sh
docker run -d -v 'path2localdirectory/VCFs:/root/VCFs' -v 'path2localdirectory/GDS:/root/GDS' -p 127.0.0.1:3838:3838 lamington
```

Verify the deployment by navigating to your server address in your
preferred browser.

``` sh
127.0.0.1:3838
```

## Workflow

##### 1. Adding Metadata/Passport Data:

Metadata/passport data from the original VCF file, if available, is
added to the PCA data frame for future visualization.

##### 2. Data Ingestion:

Genotypes (VCF Format) - Accepts genotype data in VCF format. Locally,
users can load files from the local filesystem. On a server, a
designated folder lists VCF files for user selection. Sample
Passport/Metadata - Users can input metadata/passport data for VCF
samples in CSV/TSV and/or Excel formats. \

##### 3. Filtering
Visualisation: Using the slider users can visualise the change in
Missing Rate and MAF on the histogram on the right.

##### 4. Genotype Data Filtering:

Prior to analysis, users can filter genotype data based on Minor Allele
Frequency (MAF) and call rate. SNP numbers are reducible through Linkage
Disequilibrium (LD) pruning. Users define a list of SNPs for downstream
analyses.

##### 5. PCA Calculation:

After defining the set of SNPs, the app calculates Principal Component
Analysis (PCA).

##### 6. Calculation of Core Sets:

Users define core sets by utilizing Core Hunter. Main Core Hunter
options are accessible to users, allowing core set definition in various
sizes simultaneously. Core sets are marked in the PCA data frame for
visualization.


##### 7. PCA Plot Visualization and Outlier Definition:

Core sets from step 5 are visualized within the PCA plot from step 4.
The calculated PCA is visualized with an interactive plot allowing
zooming and sample selection.Sample names are displayed on cursor hover.
After outlier removal, users define core sets by utilizing Core Hunter.
Users can add samples to an exclusion list, regenerating the PCA by
revisiting step 3. Core sets from step 5 are visualized within the PCA
plot from step 4. Core sets are exportable as a csv file containing the
list of samples,PCA and Population data and the core set.

##### 8. Final Plot

User can use the addin to create final plot and customise it according
to their needs. Selecting Tab2 from the list of data frames. \##
Required Dependencies: - shiny - corehunter - SNPRelate - ggplot2 -
shinyFiles - dplyr - DT - esquisse - scatterD3 - shinycssloaders -
shinythemes

## License

This code is licensed under the GPLv3. Please see the file LICENSE.txt
for information.
