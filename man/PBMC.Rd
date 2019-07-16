\name{PBMC}
\alias{PBMC}
\docType{data}
\title{
Normal Human Peripheral Blood Mononuclear Cells (PBMCs)
}
\description{
This dataset contains 49651 events from a 25-marker panel CyTOF analysis of normal human PBMCs. Expression values are the raw (pre-transformation) median intensity (MI) values.

Data has been pre-gated to include only DNA-intercalator (Iridium) positive and CD45-high events. Expert biaxial gating was used to separate these events into 7 major blood cell populations: CD4+ T cells, CD8+ T cells, IgM+ B cells, IgM- B cells, dendritic cells (DCs), natural killer (NK) cells, and monocytes. Per-cell population identity is specified in the \code{cluster} channel (variable).

This dataset is meant to be used as an example with the \code{\link{MEM}} package.

See refs for experimental protocol and further details.
}
\usage{data("PBMC")}
\format{
  A data frame with 49651 observations on the following 26 variables.

  The following 25 surface markers were measured by CyTOF.

  \describe{
    \item{\code{CD19}}{B cell receptor (BCR)}
    \item{\code{CD117}}{c-Kit; RTK expressed by stem and progenitor cells}
    \item{\code{CD11b}}{ITGAM, macrophage-1 antigen; complement receptor 3}
    \item{\code{CD4}}{T cell receptor (TCR) co-receptor; binds antigens presented by MHC II}
    \item{\code{CD8}}{T cell receptor (TCR) co-receptor; binds antigens presented by MHC I}
    \item{\code{CD20}}{B cell surface protein}
    \item{\code{CD34}}{surface protein expressed by hemotopoeitic stem cells and lost over course of differentiation}
    \item{\code{CD61}}{surface marker expressed on platelets}
    \item{\code{CD123}}{Interleukin-3 receptor; expressed by progenitor cells}
    \item{\code{CD45RA}}{CD45 isoform expressed by Naive T lymphocytes}
    \item{\code{CD45}}{protein tyrosine phosphatase; expressed by all mature leukocytes}
    \item{\code{CD10}}{membrane metallo-endopeptidase expressed by common lymphoid progenitors}
    \item{\code{CD33}}{Siglec-3; expressed by myeloid cells}
    \item{\code{CD11c}}{complement receptor; highly expressed on dendritic cells and myeloid cells}
    \item{\code{CD14}}{pattern recognition receptor expressed by innate lymphoid cells}
    \item{\code{CD69}}{involved in signaling and proliferation of activated t-lymphocytes and natural killer cells}
    \item{\code{CD15}}{plays role in phagocytosis and chemotaxis; expressed in multiple blood cell malignancies}
    \item{\code{CD16}}{low affinity Fc receptor for IgG; expressed by natural killer cells, neutrophils, and myeloid cells}
    \item{\code{CD44}}{cell adhesion molecule and hyaluronic acid receptor}
    \item{\code{CD38}}{highly expressed on germinal center B cells and plasma cells}
    \item{\code{CD25}}{IL-2 receptor; expressed by activated T cells}
    \item{\code{CD3}}{T cell receptor (TCR)}
    \item{\code{IgM}}{heavy chain isoform of BCR}
    \item{\code{HLADR}}{MHC class II receptor}
    \item{\code{CD56}}{NCAM; expressed by natural killer cells}
    \item{\code{cluster}}{
    1: CD4+ T cells
    2: CD8+ T cells
    3: Dendritic cells (DCs)
    4: IgM- B cells
    5: IgM+ B cells
    6: Monocytes
    7: Natural killer (NK) cells
    }
  }
}
\details{
The dataset should be arcsinh transformed with cofactor of 15. See \code{\link{MEM}} for more details.

}
\source{
Leelatian et al., Methods Mol Biol, 2015.
}
\references{
Diggins et al., Methods, 2016.
Diggins et al., Nature Methods, 2017
}
\examples{
data(PBMC)
## maybe str(PBMC) ; plot(PBMC) ...
}
\keyword{datasets}
