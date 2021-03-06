% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.npqcs.R
\name{plot.npqcs}
\alias{plot.npqcs}
\alias{plot.npqcs.r}
\alias{plot.npqcs.Q}
\alias{plot.npqcs.S}
\title{Plot method for 'npqcs' objects}
\usage{
\method{plot}{npqcs}(x, title, subtitle, xlab, ylab, ylim, lim = TRUE, ...)

\method{plot}{npqcs.r}(
  x,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ...
)

\method{plot}{npqcs.Q}(
  x,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ...
)

\method{plot}{npqcs.S}(
  x,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ...
)
}
\arguments{
\item{x}{Object npqcs (Multivarite Quality Control Statical)}

\item{title}{an overall title for the plot}

\item{subtitle}{a sub title for the plot}

\item{xlab}{a title for the x axis}

\item{ylab}{a title for the y axis}

\item{ylim}{the y limits of the plot}

\item{lim}{a logical value indicating that limits should be constant.}

\item{...}{arguments to be passed to or from methods.}
}
\description{
Generic function for plotting Multivarite charts of object of class 'npqcs' to perform statistical 
quality control.
}
