
#'  Make vertices outlining spots/subspots for geom_polygon()
#'
#' @param spot_positions Dataframe with row and col of each spots or subspots
#' @param vertex vertices of each spots or subspots

#' @return Data frame with vertices outlining spots/subspot

make_spot_vertices <- function(spot_positions, vertex) {

  spot_vertices <- merge(spot_positions, vertex)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  as.data.frame(spot_vertices)
}


#' Make vertices for each triangle subspot of a hex
#'
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' @keywords  internal

tri_spots_vertices <- function(sce, x = "sub.col", y = "sub.row", fill = "spatial.cluster"){

  cdata <- colData(sce)

  if(is.character(fill) && length(fill) == 1){
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }else if(is.vector(fill)){
    spot_positions <- cdata[, c(x, y)]
    spot_positions$fill <- fill
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }

  spot_positions$spot <- rownames(spot_positions)

  ## r = inradius, distance from center to edge midpoint
  r <- 1 / 2
  ## R = circumradius, distance from center to vertex
  R <- (2 / sqrt(3)) * r

  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1

  spot_positions$y.pos <- spot_positions$y.pos * sqrt(3) / 2
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2

  vertex <- do.call(rbind, list(
    data.frame(x.offset=c(0, 0, r), y.offset=c(0, -R, -R/2), subspot.idx=3),
    data.frame(x.offset=c(0, r, r), y.offset=c(0, -R/2, R/2), subspot.idx=5),
    data.frame(x.offset=c(0, r, 0), y.offset=c(0, R/2, R), subspot.idx=1),
    data.frame(x.offset=c(0, 0, -r), y.offset=c(0, R, R/2), subspot.idx=2),
    data.frame(x.offset=c(0, -r, -r), y.offset=c(0, R/2, -R/2), subspot.idx=6),
    data.frame(x.offset=c(0, -r, 0), y.offset=c(0, -R/2, -R), subspot.idx=4)
  ))

  spot_vertices <- make_spot_vertices(spot_positions, vertex)

  spot_vertices
}


#'  Make vertices for each hexagonal spot
#'
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' @keywords  internal
#'
hex_spots_vertices <- function(sce, x = "col", y = "row", fill = "spatial.cluster"){

  cdata <- colData(sce)

  if(is.character(fill) && length(fill) == 1){
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }else if(is.vector(fill)){
    spot_positions <- cdata[, c(x, y)]
    spot_positions$fill <- fill
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }

  spot_positions$spot <- rownames(spot_positions)

  ## r = inradius, distance from center to edge midpoint
  r <- 1 / 2
  ## R = circumradius, distance from center to vertex
  R <- (2 / sqrt(3)) * r

  ## Start at (1-indexed origin)
  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1

  spot_positions$y.pos <- spot_positions$y.pos * sqrt(3) / 2
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2

  ## Vertices of each hex
  vertex <- data.frame(x.offset = c(r, 0, -r, -r, 0, r), y.offset = c(R/2, R, R/2, -R/2, -R, -R/2))

  spot_vertices <- make_spot_vertices(spot_positions, vertex)
  spot_vertices

}


#'  Make vertices for each square spot
#'
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' @keywords  internal
square_spots_vertices  <- function(sce, x = "col", y = "row", fill = "spatial.cluster", scale.factor = 1){

  cdata <- colData(sce)


  if(is.character(fill) && length(fill) == 1){
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }else if(is.vector(fill) || is.factor(fill)){
    spot_positions <- cdata[, c(x, y)]
    spot_positions$fill <- fill
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
  }

  spot_positions$spot <- rownames(spot_positions)

  vertex <- data.frame(data.frame(x.offset=c(0, 1, 1, 0), y.offset=c(0, 0, 1, 1)))
  vertex <- vertex * scale.factor

  spot_vertices <- make_spot_vertices(spot_positions, vertex)
  return(spot_vertices)
}

#' Make vertices outlining spots/subspots for geom_polygon()
#'
#' @param sce SingleCellExperiment with row/col in colData
#' @param fill Name of a column in \code{colData(sce)} or a vector of values to
#'   use as fill for each spot
#' @param platform "Visium" or "ST", used to determine spot layout
#' @param is.enhanced If true, \code{sce} contains enhanced subspot data instead
#'   of spot-level expression. Used to determine spot layout.
#'
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#'
#' @keywords internal
make_vertices <- function(sce, fill, platform, is.enhanced) {

  if (platform == "Visium") {
    if (!is.enhanced) {
      vertices <- hex_spots_vertices(sce, x = "col", y = "row", fill = fill)
    }
    else{
      vertices <- tri_spots_vertices(sce, x = "sub.col", y = "sub.row", fill = fill)
    }
  }

  if (platform == "ST") {
    if (!is.enhanced) {
      vertices <- square_spots_vertices(sce, x = "col", y = "row", fill = fill, scale.factor = 1)
    }
    else {
      vertices <- square_spots_vertices(sce, x = "sub.col", y = "sub.row", fill = fill, scale.factor = 1 / 3)
    }
  }

  return(vertices)
}

#'  Plot spatial cluster assignments.
#'
#' @param sce SingleCellExperiment. If \code{label} is specified and is a string,
#'   it must exist as a column in \code{colData(sce)}.
#' @param label Labels used to color each spot. May be the name of a column in
#'   \code{colData(sce)}, or a vector of discrete values.
#' @param palette  Optional vector of hex codes to use for discrete spot values.
#' @param platform Spatial transcriptomic platform.
#' @param is.enhanced Used to determine spot layout.
#'
#' @return Returns a ggplot object.
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_manual coord_equal labs theme_void
ClusterPlot <- function(sce, label = "spatial.cluster", palette = NULL, platform=NULL, is.enhanced=NULL){

  if(is.null(platform))
    platform <- metadata(sce)$BayesSpace.data$platform

  if(is.null(is.enhanced))
    is.enhanced<- metadata(sce)$BayesSpace.data$is.enhanced

  vertices <- make_vertices(sce, fill = label, platform, is.enhanced)

  cplot <- ggplot(data = vertices, aes_(x = ~x.vertex, y = ~y.vertex, group = ~spot, fill = ~factor(fill))) + geom_polygon() + labs(fill = "cluster") + coord_equal() + theme_void()

  if (!is.null(palette))
    cplot <- cplot + scale_fill_manual(values=palette)

  cplot

}

#' Plot spatial gene expression
#'
#' @param sce SingleCellExperiment. If \code{label} is specified and is a string,
#'   it must exist as a column in \code{colData(sce)}.
#' @param feature Feature vector used to color each spot. May be the name of a
#'   gene/row in an assay of \code{sce}, or a vector of continuous values.
#' @param platform Spatial transcriptomic platform.
#' @param is.enhanced Used to determine spot layout.
#' @param diverging If true, use a diverging color gradient in
#'   \code{featurePlot()} (e.g. when plotting a fold change) instead of a
#'   sequential gradient (e.g. when plotting expression).
#' @param low,mid,high Optional hex codes for low, mid, and high values of the
#'   color gradient used for continuous spot values.
#'
#' @return Returns a ggplot object.
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_gradient scale_fill_gradient2 coord_equal labs theme_void
#' @importFrom scales muted

FeaturePlot <- function(sce, feature, platform  =NULL, is.enhanced=NULL, diverging=FALSE,
                        low=NULL, high=NULL, mid=NULL){

  if(is.null(platform))
    platform <- metadata(sce)$BayesSpace.data$platform

  if(is.null(is.enhanced))
    is.enhanced<- metadata(sce)$BayesSpace.data$is.enhanced

  if (is.character(feature)) {
    fill <- assay(sce, "logcounts")[feature, ]
    fill.name <- feature
  } else {
    fill <- feature
    fill.name <- "Expression"
  }

  vertices <- make_vertices(sce, fill = fill, platform, is.enhanced)

  fplot <- ggplot(data=vertices,
                  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) + geom_polygon() + labs(fill=fill.name) + coord_equal() + theme_void()

  if (diverging) {
    low <- ifelse(is.null(low), "#F0F0F0", low)
    high <- ifelse(is.null(high), muted("red"), high)
    fplot <- fplot + scale_fill_gradient(low=low, high=high)
  } else {
    low <- ifelse(is.null(low), muted("blue"), low)
    mid <- ifelse(is.null(mid), "#F0F0F0", mid)
    high <- ifelse(is.null(high), muted("red"), high)
    fplot <- fplot + scale_fill_gradient2(low=low, mid=mid, high=high)
  }

  fplot
}
