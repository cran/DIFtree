% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.DIFtree.R
\name{plot.DIFtree}
\alias{plot.DIFtree}
\title{Plotting of Item focussed Trees}
\usage{
\method{plot}{DIFtree}(x, item, component = "intercept", cex.lines = 2,
  cex.branches = 1, cex.coefs = 1, cex.main = 1, title = NULL, ...)
}
\arguments{
\item{x}{Object of class \code{\link[DIFtree]{DIFtree}}}

\item{item}{Number of the item, for which the tree shall be plotted}

\item{component}{Component of the model for which the tree shall be plotted; 
can be \code{"intercept"} or \code{"slope"}. For \code{"Rasch"} and \code{"PCM"} only one tree of item parameters
is available for each DIF item and therefore \code{component} will be ignored.}

\item{cex.lines}{Width of branches of the tree}

\item{cex.branches}{Size of the labels of branches of the tree}

\item{cex.coefs}{Size of coefficients in the terminal nodes of the tree}

\item{cex.main}{Size of the title of the tree}

\item{title}{Optional title, which is added to the tree;
if \code{title=NULL} the title is the number of the plotted item.}

\item{...}{Further arguments passed to or from other methods}
}
\description{
Visualization of trees for items with DIF identified by item focussed recursive partitioning 
in dichotomous or polytomous items.
}
\examples{
data(data_sim_Rasch)
 
Y <- data_sim_Rasch[,1]
X <- data_sim_Rasch[,-1]
 
\dontrun{
 
mod <- DIFtree(Y=Y,X=X,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
 
plot(mod,item=1)
}
}
\references{
Berger, Moritz and Tutz, Gerhard (2016): Detection of Uniform and Non-Uniform Differential Item Functioning 
by Item Focussed Trees, Journal of Educational and Behavioral Statistics 41(6), 559-592.

Bollmann, Stella, Berger, Moritz & Tutz, Gerhard (2018): Item-Focussed Trees for the Detection 
of Differential Item Functioning in Partial Credit Models, Educational and Psychological Measurement 78(5), 781-804.

Tutz, Gerhard and Berger, Moritz (2016): Item focussed Trees for the Identification of Items
in Differential Item Functioning, Psychometrika 81(3), 727-750.
}
\seealso{
\code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}, \code{\link[DIFtree]{summary.DIFtree}}
}
\author{
Moritz Berger <moritz.berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
}
