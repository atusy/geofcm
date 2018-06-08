#' GEOFCM based on Yoshida et al. (2018)
#' 
#' GEOFCM is an extention of spatial fuzzy c-means (sFCM), which is an extention of FCM,
#' 
#' @param x matrix or data.frame containing observed data other than coordinates
#' @param coord matrix or data.frame containing coordinate data
#' @param C number of clusters
#' @param m Controls fuzziness basedd on fuzzy c-means (default = 1.5)
#' @param q Controls fuzziness based on distance (default = 1.5)
#' @param r Distance in km to find neighbors of each data point. Default value is NULL and use median distance among data points.
#' @param preprocess Either 'whiten', 'scale', or NA. NA avoids preprocessing.
#' @param iter.max an integer indicating maximum iteration for clustering  (default = 100L).
#' @param th a numeric indicating threshold to determine convergence (default = 1e-5).
#' 
#' @importFrom geosphere distm
#' @importFrom ForeCA whiten
#' @importFrom stats median
#' @importFrom stats runif
#' 
#' @export
geofcm <- function(
  x,
  coord, 
  C, 
  m = 1.5, 
  q = 1.5,
  r = NULL, 
  preprocess = c('whiten', 'scale', NA),
  iter.max = 100L, 
  th = 1e-5 
) {
  x <- as.matrix(x)
  preprocess <- match.arg(preprocess)
  x <- if(preprocess == 'whiten') {
    ForeCA::whiten(x)$U
  } else if(preprocess == 'scale') {
    scale(x)
  } else {
    x
  }
  # return(x)
  N <- nrow(x) # Number of observations
  d <- h <- matrix( # initializing d and h (can be any)
    runif(N * C, 0, 1),
    N,
    C
  ) # initialization is used to prepare initial u
  u <- (d / rowSums(d)) ^ m #initializing u randomly
  D <- geosphere::distm(coord) / 1000# Distance matrix
  D_med <- median(D[upper.tri(D)]) # Median distance among locations
  if(is.null(r)) r <- D_med
  NB <- D < r # Neighborhood matrix
  rm(D) # remove D
  mu <- mu_old <- matrix(0, C, ncol(x)) # initializing center of clusters
  for(j in 1:iter.max) { # iteration for clustering
    for(k in 1:C) { # find new centers
      mu[k, ] <- colSums(u[, k] * x) / sum(u[, k])
    }
    if(all(abs(1 - mu_old / mu) < th)) break() # finish calculation if relative difference mu and mu_old are small enough
    mu_old <- mu # assign mu to mu_old
    for(k in 1:C) { # distance among cluster centers and data points
      d[, k] <- colSums((t(x) - mu[k, ]) ^ 2)
    }
    for(k in 1:C) { # membership degree in FCM
      u[, k] <- 1 / rowSums((d / d[, k]) ^ (1 / (1 - m)))
    }
    for(k in 1:C) { # 
      for(i in 1:N) {
        h[i, k] = sum(u[NB[, i], k])
      }
    }
    h <- (h / colSums(NB))^ q
    deno <- rowSums(h * u) #denominator for next iteration
    for(k in 1:C) { # modify u based on neighborhood 
      u[, k] <- (u[, k] * h[, k]) / deno
    }
  }
  
  structure( # return GEOFCM class object
    list(
      hard_clusters = apply(u, 1, which.max),
      centers = mu,
      membership = u,
      convergence = j < iter.max,
      median_distance = D_med
    ),
    class = c('GEOFCM', 'list')
  )
}
#' print for GEOFCM
#' 
#' @param x GEOFCM class object
#' @param n number of rows to view matrix of membership degrees
# #' @param plot TRUE or FALSE to plot the result
#' @param ... other arguments passed to plot
# #' @importFrom graphics plot
#' @importFrom utils head
# #' @importFrom colorspace plot
#' 
#' @export
print.GEOFCM <- function(
  x, n = 5, 
  # plot = FALSE, 
  ...
) {
  cat('hard clusters\n')
  print(x$hard_clusters)
  cat('center of clusters\n')
  print(x$centers)
  cat('\n')

  cat(paste(
    'median distance among data points are', 
    x$median_distance,
    'km\n\n'
  ))
  
  cat(paste(
    'first',
    n,
    'rows of membership matrix\n'
  ))
  print(head(x$membership, n))
  
  # if(plot) plot(x)
  
  invisible(x)
}

#' #' plot for GEOFCM
#' #' 
#' #' @param x GEOFCM class object
#' #' @param ... other arguments passed to print
#' #' @importFrom stats prcomp
#' #' 
#' #' @export
#' plot.GEOFCM <- function(x, ...) {
#'   plot(
#'     prcomp(x$x, center = TRUE, scale = TRUE)$x[, 1:2],
#'     col = x$hard_clusters,
#'     ...
#'   )
#'   invisible(x)
#' }
