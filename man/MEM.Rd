\name{MEM}
\alias{MEM}
\title{
Marker Enrichment Modeling
}
\description{
The MEM function takes pre-clustered, single-cell data as input and calculates relative enrichment scores for each marker on each population.
}
\usage{
MEM_values <- MEM(exp_data, transform=FALSE, cofactor=1, choose.markers=FALSE, markers="all",
choose.ref=FALSE, zero.ref=FALSE, rename.markers=FALSE, new.marker.names="none",
file.is.clust=FALSE, add.fileID=FALSE, IQR.thresh=NULL, output.prescaled.MEM=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp_data}{
list of file names or a \code{matrix} or \code{data.frame} object where the last column contains a numeric cluster ID for each cell (row). If \code{exp_data} is a list of files, either each file must be one cluster or each file must contain a \code{cluster} channel (column) that specifies a numeric cluster ID for each cell (row). See \code{details} for more information.
}
  \item{transform}{
\code{TRUE} or \code{FALSE}; whether or not to apply \code{\link{asinh}} transformation to the data. Default is \code{FALSE}.
}
  \item{cofactor}{
numeric; if \code{transform} is \code{TRUE}, what cofactor should be applied. Default is 1. Arcsinh transformed value = arcsinh(raw value/cofactor)
}
  \item{choose.markers}{
\code{TRUE} or \code{FALSE}; whether or not the user wants to choose the markers (columns) for analysis. If data contains markers that will not be used in the analysis (e.g. SSC or FSC channels in flow data), should be set to \code{TRUE}. If \code{FALSE}, either all of the markers in the experiment data will be used in MEM or the user can pass a character string of the markers to be used in the analysis using the function call (markers) below.
}
  \item{markers}{
\code{"all"} or ex.\code{"1:2,7,11:12,25"}; if the user wants to choose markers to be used in the MEM analysis without having the console ask for a user input, enter a character string similar to the one shown in the example. The markers chosen should be separated with colons or commas, without spaces spaces between. If \code{"all"}, all of the markers will be used in MEM.
}
  \item{choose.ref}{
\code{TRUE} or \code{FALSE}; Default reference for each population is all other populations in the dataset. For example, in a dataset containing 7 clusters, reference for population 1 would include clusters 2-7. If set to \code{TRUE}, user will be prompted to enter which cluster(s) should be used as reference instead of the default bulk non-population reference.
}
  \item{zero.ref}{
\code{TRUE} or \code{FALSE}; If set to \code{TRUE}, a zero, or synthetic negative, reference will be used for all populations. MAGref therefore is 0 and IQRref is the median IQR across all markers chosen.
}
  \item{rename.markers}{
\code{TRUE} or \code{FALSE}; if \code{TRUE}, user will be prompted to enter new column names. Default \code{FALSE}. If \code{FALSE}, either the column names will not be changed or the user can pass a character string of the new column names using the function call (new.marker.names) below.
}
  \item{new.marker.names}{
\code{"none"} or ex.\code{"CD4,CD19,HLA-DR,CD8,CD14,CD16"}; if user wants new column names for channels without having the console ask for a user input, enter a character string like the one shown in the example. Each new column name should be separated by a comma, without spaces between names. If \code{"none"}, the column names will not be changed.
}
  \item{file.is.clust}{
\code{TRUE} or \code{FALSE}; if multiple files are entered as input and each file contains cells from only one cluster, should be set to \code{TRUE}. This prompts function to merge data into one \code{matrix} for analysis and to add a file ID for each file that will stand in as the cluster ID. A text file indicating which file corresponds to which cluster number will be written to the \code{output files} folder (by default will be created as a subdirectory in your working directory).
}
  \item{add.fileID}{
\code{TRUE} or \code{FALSE}; if multiple files are entered but \code{file.is.clust} is \code{FALSE}, this indicates that there are multiple files but each contains cells from multiple clusters and that there is already a cluster channel included as the last column in each file. If \code{add.fileID} is \code{TRUE}, a file ID will be appended to the cluster ID so user can identify the file as well as cluster from which each population came.
}
  \item{IQR.thresh}{
Default \code{NULL}. Optionally can be set to a numeric value. See \code{Details} for more information.
}
  \item{output.prescaled.MEM}{
Default \code{FALSE}. If TRUE, creates folder in working directory called "output files" containing a TXT file with pre-scaled MEM values. The MEM matrix output by \code{\link{build.heatmaps}} contains post-scaled (-10 to +10 scale) MEM values.
}
}
\details{
For each population and its reference, MEM first calculates \code{\link{median}} marker levels and marker interquartile ranges (\code{\link{IQR}}), and then calculates MEM scores according to the equation

\code{MEM = |Median_Pop - Median_Ref| + IQR_Ref/IQR_Pop -1 ; if Median_Pop - Median_ref < 0, -MEM}

A dataset is provided as an example to be used with \code{\link{MEM}} and \code{\link{build.heatmaps}}. Please see dataset \code{\link{PBMC}} for more details.

Input data can be file type .txt, .fcs, or .csv. A matrix or data.frame object where the last column contains cluster identy per cell is also accepted. In all cases, the expected data structure is cells (datapoints) in rows and measured markers (i.e. features, parameters) in columns of the input data.

\code{IQR threshold}: The MEM equation takes the ratio of population and reference IQRs and adds this value to the difference in medians. Low IQR values below 1, like those resulting from background noise level measurements, can therefore artificially inflate the overall MEM score. In order to correct this, a threshold of 0.5 is automatically applied.
However, the function can calculate an IQR threshold using the input data. If IQR_thresh is set to "auto", the threshold will be calculated as the IQR associated with the 2nd quartile median value across all populations and corresponding reference populations. This should be used if the user anticipates that 0.5 will not be an adequate threshold for the particular dataset.
}
\value{
\item{MAGpop }{Matrix; Median expression level of markers on each population}
\item{MAGref }{Matrix; Median expression on each population's corresponding reference population}
\item{IQRpop }{Matrix; IQR of markers on each population}
\item{IQRref }{Matrix; IQR on each population's corresponding reference population}

}
\references{
Diggins et al., Nature Methods, 2017
}
\author{
Kirsten Diggins, Sierra Barone, and Jonathan Irish
}
\note{
The object generated from \code{\link{MEM}} is meant to be passed to \code{\link{build.heatmaps}} which will generate MEM labels and heatmaps.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{build.heatmaps}}
}
\examples{
## For multiple file input, set working directory to folder containing files, then
## infiles <- dir()

## For single file or object input (e.g. PBMC), input data directly into MEM function

## User inputs
MEM_values = MEM(PBMC, transform=TRUE, cofactor=15, choose.markers=FALSE, 
markers="all", choose.ref=FALSE, zero.ref = FALSE, rename.markers=FALSE, 
new.marker.names="none", IQR.thresh=NULL)}

