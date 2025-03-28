% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splicejam-shiny-server.R
\name{sashimiAppServer}
\alias{sashimiAppServer}
\title{Sashimi Shiny app server}
\usage{
sashimiAppServer(input, output, session)
}
\arguments{
\item{input}{provided by shiny}

\item{output}{provided by shiny}

\item{session}{provided by shiny}
}
\description{
Sashimi Shiny app server
}
\details{
This function contains the server logic for the R-shiny
Splicejam Sashimi viewer.

The R-shiny app is started by \code{launchSashimiApp()}, which
calls \code{shiny::shinyApp()}, using arguments \code{server}, \code{ui},
\code{onStart}, and \code{options}. This function fulfills the
argument \code{server}.
}
\concept{splicejam R-shiny functions
\subsection{@import jamba}{
}

\subsection{@import dplyr}{
}

\subsection{@import ggplot2}{
}

\subsection{@import plotly}{
}

\subsection{@importFrom plotly subplot}{
}

\subsection{@import GenomicRanges}{
}}
