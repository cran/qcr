% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdqcs.depth.r
\name{fdqcs.depth}
\alias{fdqcs.depth}
\alias{fdqcs.depth.default}
\alias{fdqcs.depth.fdqcd}
\title{Function to plot depth functional data (DFD) - chart}
\usage{
fdqcs.depth(x, ...)

\method{fdqcs.depth}{default}(
  x,
  data.name = NULL,
  func.depth = depth.mode,
  nb = 200,
  type = c("trim", "pond"),
  ns = 0.01,
  plot = TRUE,
  trim = 0.025,
  smo = 0.05,
  draw.control = NULL,
  ...
)

\method{fdqcs.depth}{fdqcd}(
  x,
  func.depth = depth.mode,
  nb = 200,
  type = c("trim", "pond"),
  ns = 0.01,
  plot = TRUE,
  trim = 0.025,
  smo = 0.05,
  draw.control = NULL,
  ...
)
}
\arguments{
\item{x}{an R object (used to select the method). See details.}

\item{...}{arguments passed to or from methods.}

\item{data.name}{a string that specifies the title displayed on the plots. 
If not provided it is taken from the name of the object's data.}

\item{func.depth}{Type of depth measure, by default depth.mode.}

\item{nb}{The number of bootstrap samples.}

\item{type}{the method used to trim the data (trim or pond).}

\item{ns}{Quantile to determine the cutoff from the Bootstrap procedure}

\item{plot}{a logical value indicating that it should be plotted.}

\item{trim}{The porcentage of the trimming.}

\item{smo}{The smoothing parameter for the bootstrap samples.}

\item{draw.control}{ist that it specifies the col, lty and lwd for objects: fdataobj, statistic, IN and OUT.}
}
\description{
This function is used to compute statistics required by the DFD chart.
}
\examples{
\dontrun{
library(qcr)
m <- 30
tt<-seq(0,1,len=m)
mu<-30 * tt * (1 - tt)^(3/2)
n0 <- 100
set.seed(12345)
mdata<-matrix(NA,ncol=m,nrow=n0)
sigma <- exp(-3*as.matrix(dist(tt))/0.9)
for (i in 1:n0) mdata[i,]<- mu+0.5*mvrnorm(mu = mu,Sigma = sigma )
fdchart <- fdqcd(mdata)
plot.fdqcd(fdchart,type="l",col="gray")
set.seed(1234)
fddep <- fdqcs.depth(fdchart,plot = T)
plot(fddep,title.fdata = "Fdata",title.depth = "Depth")
summary(fddep)
}
}
\references{
Flores, M.; Naya, S.; Fernández-Casal,R.; Zaragoza, S.; Raña, P.; Tarrío-Saavedra, J. 
Constructing a Control Chart Using Functional Data. Mathematics 2020, 8, 58.
}
