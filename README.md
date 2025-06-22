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
    GitHub is under a GPL license. Docker builds are provided for easy
    deployment on any operating system.
-   **Visualisation Functions**: Users can interactively plot statistics
    such as minor allele frequency (MAF) and call rate, defining
    suitable cutoffs.
-   **PCA Plotting Functions**: Users can generate PCA plots for the
    population and colour data points based on sample information from an
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
steps by first downloading [*R*](https://cran.r-project.org/), [*Java*](https://www.java.com/en/download/) and [*RStudio*](https://posit.co/download/rstudio-desktop/)(optional). 

Execute the following commands once installed to configure the R environment:

**Install the dependencies**

``` sh
install.packages(c('shiny', 'ggplot2', 'shinyFiles', 'dplyr', 'corehunter', 'rhandsontable', 'DBI', 'shinyjs', 'bcrypt', 'RSQLite', 'shinyalert', 'plotly', 'DT','esquisse','scatterD3','shinycssloaders','shinythemes','rJava','corehunter', 'BiocManager'))
if (!requireNamespace("BiocManager",quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
```

**You may now run the shiny app with just one command in R:**

``` sh
#This command runs Lamington from Github without explicitly cloning lamington
shiny::runGitHub("lamington", "plantinformatics",ref="main",subdir = "R")
```

**Alternatively**

You can clone **lamington** directly to your computer as follows;

``` sh
# First clone the repository with git.
git clone git@github.com:plantinformatics/lamington.git

#Then within an R terminal 
setwd("pathtodownload/lamington") # change to match where you downloaded this repo to run lamington
shiny::runApp('R') # runs the app 
```

## Using Docker

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

After the build is complete, run the Docker image. You can map the
container's exposed port (defined in the [*Dockerfile*](Dockerfile)) to
a port on your host machine. For example, you can use the appropriate docker run command to map host port 3838 to the container's port 3838.

To enable data access, mount the VCF directory (for direct uploads), the GDS directory (containing saved GDS files) and DB directory (containing an SQLITE database for user and file managemet) to the corresponding
locations within the container.

``` sh
docker run -d -v 'path2localdirectory/DB:/root/DB' -v 'path2localdirectory/VCFs:/root/VCFs' -v 'path2localdirectory/GDS:/root/GDS' -p 127.0.0.1:3838:3838 lamington
```

Verify the deployment by navigating to your server address in your
preferred browser.

``` sh
127.0.0.1:3838
```

## Workflow

##### 0. Creating User/Logging in
Lamington uses a simple user profile management to help manage data uploads. If using Lamington for the first time, you will be required to create a user account as shown on the screenshot below
![image](https://github.com/user-attachments/assets/73819e17-74cf-48bb-b616-9b2245a87ecf)

Once registered, you will be prompted to log in as in the screenshot below
![image](https://github.com/user-attachments/assets/8ce6be6c-d4a7-4c09-a45a-594843c79983)


##### 1. Importing Metadata/Passport Data:

Metadata/passport data can be imported into **Lamington** on the *'Add POP Data'* tab.

![image](https://github.com/user-attachments/assets/c775c6d2-5f1b-4946-84b5-acaf63e850b9)

As an example, the [*passport data*](Data/241126_AGG_Chickpea_Lamington_Input.fixID.txt) extracted using [*Genolink*](https://github.com/plantinformatics/genolink) for the [AGG Chickpea - Release 241203](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SQFKJW) is imported into **Lamington** as shown in the figure below

![image](https://github.com/user-attachments/assets/a8b042dd-ae2a-4f19-a9ac-c5ef138cd503)


##### 2. Genotype Data Ingestion:

In the *'Convert VCF File tab'*, provide your genotype data in VCF ([*Variant Call Format*](https://samtools.github.io/hts-specs/VCFv4.2.pdf)). 

You can download the following [*VCF*](https://dataverse.harvard.edu/file.xhtml?fileId=10744648&version=1.0) file from the [AGG Chickpea - Release 241203](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SQFKJW) as an example to testing **Lamington**

You then have two options:
* Upload directly: Select the VCF file from your local computer.
* Choose from the server: Browse and select the VCF file from a pre-populated list on the Server.
![image](https://github.com/user-attachments/assets/a995780d-4b6b-4fe7-8bbd-7c55329e1c7b)


* After entering a GDS file name (without the .gds extension which will be appended automatically), the 'Convert and Display' button is displayed which you can click to convert the input VCF file to GDS.
* Once the VCF file has been uploaded with a meaningful name provided for the GDS file, the data can then be loaded on the *'Select GDS file' tab*

![image](https://github.com/user-attachments/assets/7cd72a0b-caf9-4289-915b-f4615e77c942)



##### 3. Filtering

Visualisation: Using the slider users can visualise the change in Missing Rate and MAF on the histogram on the right.

![image](https://github.com/user-attachments/assets/83f68ebb-e3a8-4937-836e-df02537d2df2)


In addition, if the metadata/passport data has been uploaded, users can subset and compare different sets within the metadata/passport data.

![image](https://github.com/user-attachments/assets/6f1c57f6-0f33-4b9c-86b7-9a4d89961f44)



##### 4. Genotype Data Filtering:

Before analysis, you have the option to use the genotype data as is or filter under the *'Genotype Matrix'* tab based on the following; 
* Minor Allele Frequency (MAF)
* Call rate (CR).
* Linkage Disequilibrium (LD) pruning.
* Select specific samples by providing a list of sample IDs or selecting from the metadata/passport data.
![image](https://github.com/user-attachments/assets/26997117-608e-434a-83e9-ef688efee92b)




##### 5. PCA Calculation:

After defining the set of SNPs, the Principal Component Analysis (PCA) can be performed under the *PCA* tab.
![image](https://github.com/user-attachments/assets/5e419ff9-fbc3-4eeb-b2ff-447535c06dfd)


To gain deeper insights from the PCA results, you can include metadata/passport information. 

This allows you to:
* Explore population-specific patterns.
* Visualise relationships between population groups.
![image](https://github.com/user-attachments/assets/afcdf041-85e8-474d-8507-3872cb4f5d47)


##### 6. Calculation of Core Sets:
Lamington utilises the [**CoreHunter package**](https://cran.r-project.org/web/packages/corehunter/) to compute core sets, smaller representative subsets of your data. Lamington provides access to the main CoreHunter options, enabling you to define multiple core sets with varying sizes.  These core sets are then integrated into the PCA data frame for visualisation and analysis.

![image](https://github.com/user-attachments/assets/ca1aa767-5661-4119-841b-22ed76064b9e)



##### 7. PCA Plot Visualization and Outlier Definition:
* Core sets from step 6 are visualized using a PCA plot based on the PCA components from step 5. The calculated PCA is visualised with an interactive plot allowing for zooming and sample selection.
* Sample names are displayed on cursor hover.
* You can select and remove outliers and rerun Steps 5-7.
* You can add samples to an exclusion list, and rerun steps 5-7.
* Core sets are exportable as a CSV file containing the list of samples, PCA and population data and the core set.
![image](https://github.com/user-attachments/assets/1809390f-8c60-4411-8459-627535149ddd)


##### 8. Final Plot

User can use the addin to create final plot and customise it according
to their needs. Selecting Tab2 from the list of data frames. \##
Required Dependencies: - shiny - core hunter - SNPRelate - ggplot2 -
shinyFiles - dplyr - DT - esquisse - scatterD3 - shinycssloaders -
shinythemes

## License

This code is licensed under the GPLv3. Please see the file LICENSE.txt
for information.

