% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcs.hat.cpm.r
\name{qcs.hat.cpm}
\alias{qcs.hat.cpm}
\title{Process capability index (estimate cpm)}
\usage{
qcs.hat.cpm(
  object,
  limits = c(lsl = -3, usl = 3),
  target = NULL,
  mu = 0,
  std.dev = 1,
  nsigmas = 3,
  k0 = 1,
  alpha = 0.05,
  n = 50,
  contour = TRUE,
  ylim = NULL,
  ...
)
}
\arguments{
\item{object}{qcs object of type \code{"qcs.xbar"} or \code{"qcs.one"}.}

\item{limits}{A vector specifying the lower and upper specification limits.}

\item{target}{A value specifying the target of the process. 
If is \code{NULL}, the target is set at the middle value bewteen specification limits.}

\item{mu}{A value specifying the mean of data.}

\item{std.dev}{A value specifying the within-group standard deviation.}

\item{nsigmas}{A numeric value specifying the number of sigmas to use.}

\item{k0}{A numeric value. If the capacity index exceeds the \code{k} value, then the process is capable.}

\item{alpha}{The significance level (0.05 for default)}

\item{n}{Size of the sample.}

\item{contour}{Logical value indicating whether contour graph should be plotted.}

\item{ylim}{The y limits of the plot.}

\item{...}{Arguments to be passed to or from methods.}
}
\description{
Estimate  \code{"cpm"} using the method described by Kerstin Vannman(2001).
}
\examples{
library(qcr)
data(pistonrings) 
xbar <- qcs.xbar(pistonrings[1:125,],plot = TRUE)
mu <-xbar$center
std.dev <-xbar$std.dev
LSL=73.99; USL=74.01
qcs.hat.cpm(limits = c(LSL,USL),
           mu = mu,std.dev = std.dev,ylim=c(0,1))
qcs.hat.cpm(object = xbar, limits = c(LSL,USL),ylim=c(0,1))
}
\references{
Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd
ed, New York, John Wiley & Sons. \cr
Vannman, K. (2001). \emph{A Graphical Method to Control Process Capability}. Frontiers in Statistical Quality Control, 
No 6, Editors: H-J Lenz and P-TH Wilrich. Physica-Verlag, Heidelberg, 290-311.\cr
Hubele and Vannman (2004). \emph{The E???ect of Pooled and Un-pooled Variance Estimators on Cpm When Using Subsamples}.
Journal Quality Technology, 36, 207-222.\cr
}
