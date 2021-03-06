% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialcontrast.R
\name{plot_map}
\alias{plot_map}
\title{Plot a map comparing a focal area to every other area}
\usage{
plot_map(latlon, summaries, focal_rate = "net.div", focal_cluster = 1)
}
\arguments{
\item{latlon}{Data.frame with fields for taxon, lon, lat and perhaps other fields. Expects a cluster field giving the cluster id.}

\item{summaries}{Output list from summarize_cluster_results()}

\item{focal_rate}{Character string the names for the rate column to plot}

\item{focal_cluster}{Integer for the cluster to plot as the focal one}
}
\value{
Nothing, though it uses print to plot the ggplot2 object
}
\description{
Plot a map comparing a focal area to every other area
}
