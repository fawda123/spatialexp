contig_plo <- function(ctch_in, flow_in, nlcd_in, subbasin = '0510', lu = 'PctUrbOp2011ws'){
  
  library(rgdal)
  library(scales)
  library(spdep)
  library(foreign)
  
  # import catchment, rename FEATUREID
  ctch <- readOGR(dirname(ctch_in), gsub('\\.shp', '', basename(ctch_in)), verbose = FALSE)
  names(ctch)[names(ctch) %in% 'FEATUREID'] <- 'COMID'
  
  # import flow lines, join with catchment
  flow <- read.dbf(flow_in)
  ctch$REACHCODE <- flow$REACHCODE[match(ctch$COMID, flow$COMID)]
  
  # subset catchments by basin in flow
  ctch$BASIN <- substr(ctch$REACHCODE, 1,4)
  ctch <- ctch[ctch$BASIN == subbasin & ! is.na(ctch$BASIN),]
  
  # import nlcd, select lu
  if(grepl('\\.RData$', nlcd_in)){
    load(nlcd_in)
    assign('nlcd', get(ls()[grep('^NLCD', ls())]))
  } else {
    nlcd <- read.csv(nlcd_in)
  }
  nlcds_sub <- nlcd[, c('COMID', lu)]
  
  # merge catchment with land use
  # create area summaries
  ctch <- merge(ctch, nlcd)

  # create neighbors and spatial weights
  neighb <- poly2nb(ctch, row.names = ctch$COMID, queen = TRUE)
  neighb <- nb2listw(neighb, style = "W")

  moran.plot.dot(ctch, neighb, lu)
  
}

moran.plot.dot <- function (x, listw, lu, zero.policy = NULL, spChk = NULL, labels = NULL, 
    xlab = NULL, ylab = NULL, quiet = NULL, ...){ 

    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (is.null(quiet)) 
        quiet <- !get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(quiet))
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))
    areasqkm <- data.frame(x)[, 'AreaSqKM']
    x <- data.frame(x)[, lu]
    xname <- lu
    if (!is.numeric(x)) 
        stop(paste(xname, "is not a numeric vector"))
    if (any(is.na(x))) 
        stop("NA in X")
    n <- length(listw$neighbours)
    if (n != length(x)) 
        stop("objects of different length")
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw)) 
        stop("Check of data and weights ID integrity failed")
    labs <- TRUE
    if (is.logical(labels) && !labels) 
        labs <- FALSE
    if (is.null(labels) || length(labels) != n) 
        labels <- as.character(attr(listw, "region.id"))
    wx <- lag.listw(listw, x, zero.policy = zero.policy)
    if (is.null(xlab)) 
        xlab <- xname
    if (is.null(ylab)) 
        ylab <- paste("spatially lagged", xname)
    plot(x, wx, pch=19, col= alpha("black",0.4),cex=2.1+(log10(areasqkm/max(areasqkm))), xlab = xlab, ylab = ylab, ...)
    if (zero.policy) {
        n0 <- wx == 0
        if (any(n0)) {
            symbols(x[n0], wx[n0], inches = FALSE, circles = rep(diff(range(x))/50, 
                length(which(n0))), bg = "grey", add = TRUE)
        }
    }
    xwx.lm <- lm(wx ~ x)
    abline(xwx.lm,col="red",lwd=2)
    abline(h = mean(wx), lty = 2)
    abline(v = mean(x), lty = 2)
    infl.xwx <- influence.measures(xwx.lm)
    is.inf <- which(apply(infl.xwx$is.inf, 1, any))
    points(x[is.inf], wx[is.inf], pch = 1, cex=1,col="blue")
    if (labs) 
        text(x[is.inf], wx[is.inf], labels = labels[is.inf], 
            pos = 2, cex = 0.7)
    rownames(infl.xwx$infmat) <- labels
    if (!quiet) 
        summary(infl.xwx)
    invisible(infl.xwx)
}

environment(moran.plot.dot) <- environment(moran.plot)