% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{square_spots_vertices}
\alias{square_spots_vertices}
\title{Make vertices for each square spot}
\usage{
square_spots_vertices(
  sce,
  x = "col",
  y = "row",
  fill = "spatial.cluster",
  scale.factor = 1
)
}
\value{
Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
  vertices outlining the spot's border
}
\description{
Make vertices for each square spot
}
\keyword{internal}
