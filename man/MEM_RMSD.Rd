\name{MEM_RMSD}
\alias{MEM_RMSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MEM RMSD similarity between populations
}
\description{
MEM_RMSD calculates a normalized average RMSD score pairwise between populations given their MEM scores as input. This is meant to serve as a metric of similarity between populations.

The function calculates the sum of squares for all shared markers between two populations, then takes the square root of the average.

For "a" through n markers, the sum of squares is calculated as:
sum of squares = (a2-a1)^2 + (b2-b1)^2 ...(n2-n1)^2

Root-mean-square distance (RMSD) is calculated as:
RMSD = sqrt(sum of squares/number of markers)

The RMSD values are then converted to percentages with the maximum RMSD in the matrix set as 100 percent, so that the final RMSD score is the percent of the maximum RMSD.

Percent_max_RMSD = 100-RMSD/max_RMSD*100

The function then outputs a clustered heatmap of Percent_max_RMSD values and the matrix of numerical values used to build the heatmap.

}
\usage{
MEM_RMSD(data, format=NULL, newWindow.heatmaps=FALSE, output.matrix=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    The data input to MEM_RMSD can be either 1) a matrix of values, where populations are in rows and their MEM scores are in columns, 2) the list of matrices output by \code{\link{MEM}}, or 3) a file path pointing to a folder containing tab-delimited text files, one file for each population, where each file lists  marker names in the first column and the corresponding MEM scores in the second column.
}
 \item{format}{
 Default is NULL. When \code{format} is equal to "pop files", the function expects a file path as input where the designated folder contains one file for each population's set of MEM scores.

 }
  \item{newWindow.heatmaps}{
  If \code{TRUE}, the clustered heatmap will be output to a new R window. Default is \code{FALSE}.

 }
  \item{output.matrix}{
  If \code{TRUE}, the matrix of values calculated by the function will be output to a tab-delimited text file located in a folder called "output files" that is generated in the working directory.

 }

}
\details{
If you are calculating \code{MEM_RMSD} on population files, populations do not have to include all of the same markers. The function will determine which markers each pair of populations has in common and will use those common markers to calculate RMSD. If the populations have no markers in common, the function will terminate with an error. Note that population names must match exactly between files in order for them to be considered the same.
}
\value{
    \item{RMSD_vals }{Matrix of the calculated pairwise percent max RMSD scores}
    \item{RMSD heatmap}{Hierarchically clustered heatmap of RMSD_vals}
%% ...
}
\references{
Diggins et al., Nature Methods, 2017
}
\author{
Kirsten Diggins
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
