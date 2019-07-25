# MEM v3

Marker Enrichment Modeling (MEM) is a tool designed to calculate enrichment scores.  MEM generates human and machine readable labels that quantify the features enriched in a sample.  The classic use of MEM is to identify multiple populations of cells and to compare each population to all of the other remaining cells from the original sample.  MEM enrichment scores range from +10 (meaning greatly enriched) through 0 (meaning not enriched) to -10 (meaning greatly lacking).  MEM scores are built form two fundamental statistics, the median and interquartile range, and the output of MEM can be represented as a heatmap of the scores where the rows are each population and the columns are measured features.  This information can also be represented in a compact label where the most enriched features are listed first.

## Getting Started

MEM was designed for biology users who may not have much experience in R, and this code base has and install script and example scripts showing different uses of MEM with three example datasets.  These instructions will focus on helping biologist users.  

To get started, make sure you have the latest R and RStudio installed.  Then download this repository and open the 00_install_tools.rmd R markdown file in RStudio.  You can then press the green triangle play button on the right to step through each chunk of code.  The first chunk of code will simply print text explaining how to get started installing files.  The second file will list the EULA, which generally indicates that MEM can be used free in a not-profit / academic setting, but must be licensed from Vanderbilt University for any commercial use (http://mem.vueinnovations.com/licensing).  Continuing on, later lines of code will install the various packages used in the MEM examples.

Experts may wish to skip all this and simply install MEM and test out the vignettes within it.  This is fine, but half the fun is in how MEM can compare the results of different analysis strategies for the same dataset, so we recommend clicking through the examples (or knitting them).

### Prerequisites

Everything needed to install and run the example code is listed in the 00_install_tools.rmd file.  This will include the Bioconductor package and flowCore, which enable reading of Flow Cytometry Standard (FCS) files of single cell data.  Additional packages used in the example scripts include FlowSOM, t-SNE, and UMAP, in addition to MEM.

### Installing

Use 00_install_tools.rmd

## Running the cGVHD dataset

The first example and test is within the 01_cGVHD_example.rmd R markdown file.  

This dataset is 339 chronic graft-vs-host disease patients who have been scored according to 8 features representing involvement of 8 organ domains.  There are some clear clusters of patients in this dataset, and some less clear clusters.  The goal of this exercise is to generate patient clusters through FlowSOM clustering following t-SNE or UMAP analysis, calculate MEM labels from each analysis, and then compare the MEM labels from the two analyses using the root mean square deviation (RMSD) in the MEM labels.

1) After completing installation with 00_install_tools.rmd, load 01_cGVHD_example.rmd in RStudio and press play on the first chunk to load the libraries.  

2) Then play the next chunk, which will generate a t-SNE plot of the 339 patients in this dataset.  

3) The next chunk clusters the patients on the density across the t-SNE embedding axes using FlowSOM.  (Yes, we're aware this is a non-standard use of t-SNE and FlowSOM, but it works surprisingly well.)

4) The next chunk calculates MEM labels as a heatmap and text output to the console.  For example, the MEM label for FlowSOM cluster 1 was: 
UP Joint+4 Sclerosis+2 Fascia+1 . DN Mouth-1 Liver-1

This label means that Joint involvement was positively enriched (+4) and mouth involvement was negatively enriched (-1).  This means that patients in this group had more joint involvement and less mouth involvement than patients in other groups.  Cluster 6 gets the label:
UP Liver+10 BSA+2 GI+1 Eye+1 . DN None

This means that liver involvement was greatly enriched in this group.  The None following DN means that nothing was negatively enriched.

5) The next chunk is similar to #2 and will generate a UMAP plot of the data.

6) The next chunk is similar to #3 and will FlowSOM cluster on the density across the UMAP axes.

7) The next chucnk calculates the MEM labels for the FlowSOM clusters on the UMAP axes.  Two of them are:
5 : UP Joint+5 Sclerosis+3 Fascia+1 . DN Mouth-1 Liver-1
6 : UP Liver+10 . DN None

These two are somewhat similar to the ones we saw above.  

8) The next chunk shows the exact FlowSOM clustering on the t-SNE embedding axes shown in the original paper (Gandelman et al., Haematologica 2018).

9) The next chucnk calculates the MEM labels for the populations from #8, the published clusters from FlowSOM on t-SNE.  Two of them are:
2 : UP Joint+5 Sclerosis+3 Fascia+1 . DN Mouth-1 Liver-1
7 : UP Liver+10 . DN None

Note these are similar to the clusters 5 and 6 from the UMAP and 1 and 6 from the other t-SNE.  

10) The last chunk calculates a similarity score for the clusters from the three different analyses.  Ideally, if all analyses were finding the same number of clusters and they all had the same MEM phenotype, this heatmap would have several blocks of 3 red squares (indicating high similarity).

## Authors

* **Kirsten Diggins** - *Version 1.0* 
* **Sierra Barone** - *Versions 2.0 to 3.0* 
* **Jonathan Irish** - *Versions 1.0 to 3.0* 

## License

This project is licensed according to Vanderbilt University policy - see the [LICENSE.MD](LICENSE.MD) file for details
