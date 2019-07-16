\name{build.heatmaps}
\alias{build.heatmaps}
\title{
    Build heatmaps
}
\description{
    Takes matrix list generated from \code{\link{MEM}} as input and outputs MEM labels, heatmap of population median values, and heatmap of MEM scores.
}
\usage{
    build.heatmaps(exp_data, cluster.MEM = "none", cluster.medians = "none", 
    cluster.IQRs = "none", newWindow.heatmaps=FALSE, display.thresh = 0, 
    output.files = FALSE, labels = FALSE, only.MEMheatmap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{exp_data}{
        List of matrices generated as output from \code{\link{MEM} aka MEM_values}
    }
    \item{cluster.MEM}{
        \code{"none"},\code{"row"},\code{"col"},or \code{"both"}; which dimension(s) of the MEM heatmap to hierarchically cluster. See \code{Details} for more information.
    }
    \item{cluster.medians}{
        \code{"none"},\code{"row"},\code{"col"},or \code{"both"}; which dimension(s) of the median heatmap to cluster. If \code{"none"}, the median row and column order will match order of the MEM heatmap.
    }
     \item{cluster.IQRs}{
        \code{"none"},\code{"row"},\code{"col"},or \code{"both"}; which dimension(s) of the IQR heatmap to cluster. If \code{"none"}, the IQR row and column order will match order of the MEM heatmap.
    }
    \item{display.thresh}{
        Numeric; 0-10. Markers with enrichment scores that are equal to or greater than \code{display.thresh} will be displayed as row names on \code{MEM heatmap}.
    }
    \item{newWindow.heatmaps}{
        \code{TRUE} or \code{FALSE}; whether or not to pop out new windows to display the MEM, median, and IQR heatmaps
    }
    \item{output.files}{
        \code{TRUE} or \code{FALSE}; If \code{TRUE}, output of \code{build.heatmaps} will be written to folder \code{output folder} that is created as a subdirectory of the working directory. Written files include txt file containing MEM labels (row names), and txt files containing MEM scores, medians, and IQR values.
    }
    \item{labels}{
         \code{TRUE} or \code{FALSE}; whether or not to include MEM labels in the MEM heatmap
    }
    \item{only.MEMheatmap}{
         \code{TRUE} or \code{FALSE}; whether to only show MEM heatmap or all relevant heatmaps (MEM, median, IQR)
    }
}
\details{
    Heatmaps are clustered using the default complete linkage hierarchical clustering in the \code{\link{hclust}} function. See \code{\link{heatmap.2}} and \code{\link{hclust}} for more information.
}
\value{
    Heatmaps of median, IQR, and MEM values on each population; optionally written to file. \code{\link{X11}} is used to open plots in a new window. When windows open, user must click within the window to permit R to print the heatmap to that window (the window will remain blank until user clicks on it).
}
\references{
Diggins et al., Nature Methods, 2017
}
\author{
Kirsten Diggins, Sierra Barone, Jonathan Irish
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
        \code{\link{MEM}}, \code{\link{heatmap.2}}, \code{\link{hclust}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

build.heatmaps(MEM_values, cluster.MEM = "both", cluster.medians = "none", 
    cluster.IQRs = "none", newWindow.heatmaps=FALSE, display.thresh = 1, 
    output.files = TRUE, labels = FALSE, only.MEMheatmap = FALSE)}
