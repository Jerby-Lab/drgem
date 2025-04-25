#' Negated version of `%in%`
#'
#' This function returns a logical vector indicating if there are no matches or not.
#'
#' @param x vector or NULL: the values to be matched.
#' @param table vector or NULL: the values to be matched against.
#' @return A logical vector indicating if a match was found (negated).
#' @export
#'
`%!in%` <- Negate(`%in%`)

#' Function for computing "transcripts per million"
#'
#' Further normalizes and transforms data to
#' common expression data output: "transcripts per million"
#'
#' @param r data list
#' @return r, your data list, with an r$tpm field added.
#' @export
#'
get.tpm<-function(X, comp.reads){
  X <-as.matrix(X) # making sure the count data is a matrix
  tpm<-1000000*sweep(X,2,comp.reads,FUN = '/') # calculating transcripts per million
  tpm<-log2((tpm/10)+1) # log transform
  colnames(tpm)<- colnames(X)
  row.names(tpm) <- row.names(tpm)
  return(tpm)
}

#' Calculate Mean Squared Error
#'
#' Calculates mean squared error for the the
#'
#' @param X input matrix
#' @param X_hat reconstructed matrix
#' @return vector of mean squared error loss
#' @export
#'
calc_mse <- function(X, X_hat){
  err = apply((X - X_hat)**2, 1, mean)
  return(err)
}

#' Capping the values of an object
#'
#' This function caps the values of a numeric object to specific quantiles.
#'
#' @param X Numeric vector or matrix.
#' @param q Numeric value between 0 and 1 indicating the quantile threshold for capping.
#' @return The input object with values capped to the specified quantiles.
#' @export
#'
cap_object <- function (X, q) {
  ceil_q <- 1 - q
  ceil <- quantile(X, ceil_q)
  floor <- quantile(X, q)
  X[X > ceil] <- ceil
  X[X < floor] <- floor
  return(X)
}

#' Prepocessing Function for Overall Expression
#'
#' prepares data to compute overall expression
#' signature
#'
#' @param r data list
#' @param n.cat (default = 50)
#' @return r, your data list
#' @export
#'
prep4OE <- function(r, n.cat = 50) {
  r$zscores <- center.matrix(r$tpm, dim = 1, sd.flag = T)
  X <- 10 * ((2^r$tpm) - 1)
  r$genes.dist <- log2(rowMeans(X, na.rm = T) + 1)
  r$genes.dist.q <- discretize.prvt(r$genes.dist, n.cat = n.cat)
  b <- rowSums(is.na(r$zscores)) == 0
  if (any(!b)) {
    r <- set.list(r, b)
  }
  r$binZ <- average.mat.rows(r$zscores, r$genes.dist.q, f = colMeans)
  return(r)
}

#' Overall Expression Wrapper
#'
#' wrapper around code that computes overall expression
#' from the given gene signature
#'
#' @param r data list
#' @param sig signature list
#' @return overall expression matrix
#' @export
get.OE <- function (r, sig) {
  scores <- get.OE1(r, sig)
  names(sig) <- gsub(" ", ".", names(sig))
  two.sided <- unique(gsub(".up", "", gsub(".down", "", names(sig))))
  b <- is.element(paste0(two.sided, ".up"), names(sig)) & is.element(paste0(two.sided,
                                                                            ".down"), names(sig))
  if (any(b)) {
    two.sided <- two.sided[b]
    scores2 <- as.matrix(scores[, paste0(two.sided, ".up")] -
                           scores[, paste0(two.sided, ".down")])
    colnames(scores2) <- two.sided
    scores <- cbind(scores2, scores)
  }
  if (!is.null(r$cells)) {
    rownames(scores) <- r$cells
  }
  else {
    if (!is.null(r$samples)) {
      rownames(scores) <- r$samples
    }
  }
  return(scores)
}

#' Calculate Reconstruction Loss
#'
#' Computes the reconstruction loss of a dimensionality reduction technique,
#' such as PCA, using mean squared
#' error between the original and reconstructed data.
#'
#' @param X_scaled Scaled input matrix.
#' @param U Matrix of principal components or loadings.
#' @param n_dims Number of dimensions to use for reconstruction.
#' @return Vector of reconstruction loss values (mean squared error).
#' @export
#'
calc_recon_loss <- function(X_scaled, U, n_dims) {
  feats = intersect(row.names(X_scaled), row.names(U))
  X_hat = t(X_scaled)[,feats] %*% U[feats,1:n_dims] %*% t(U[feats,1:n_dims])
  recon_loss = calc_mse(t(X_scaled)[,feats], X_hat)
  return(recon_loss)
}


#' Matrix centering
#'
#' Helper function to center your matrix
#'
#' @param m matrix object
#' @param dim which dimension, rows = 1, columns = 2
#' @param sd.flag boolean for standardizing data to sd = 1 (default = F)
#' @return a matrix of zscores
#' @export
#'
center.matrix <- function (m, dim = 1, sd.flag = F) {
  if (dim == 1) {
    zscores <- sweep(m, 1, rowMeans(m, na.rm = T), FUN = "-")
  }
  else {
    zscores <- sweep(m, 2, colMeans(m, na.rm = T), FUN = "-")
  }
  if (sd.flag) {
    zscores <- sweep(zscores, dim, apply(m, dim, function(x) 
      (sd(x, na.rm = T))), FUN = "/")
  }
  return(zscores)
}

#' Get Matrix
#'
#' helper function for casting data into a matrix
#'
#' @param m.rows matrix rows
#' @param m.cols matrix columns
#' @param data your data
#' @return your data in matrix format
#'
get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

#' Discretize data into Quantiles
#'
#' discretizes a vector of values
#'
#' @param v vector of values
#' @param n.cat number of bins, usually G/20 if G is the number of genes in your data
#' @param q1 quantiles
#' @return discretized data
#' @export
#'
discretize.prvt <- function(v, n.cat, q1) {
  q1 <- quantile(v, seq(from = (1/n.cat), to = 1, by = (1/n.cat)), 
                 na.rm = T)
  u <- matrix(data = 1, nrow = length(v))
  for (i in 2:n.cat) {
    u[(v >= q1[i - 1]) & (v <= q1[i])] <- i
  }
  u <- paste0("Q", u)
  return(u)
}

#' Get the average of each matrix
#'
#' Function that gets the average of each column
#' of a matrix
#'
#' @param m matrix
#' @param ids identifiers
#' @param f averaging function (colMeans)
#' @export
#' 
average.mat.rows <- function (m, ids, f = colMeans) {
  ids.u <- sort(unique(ids))
  m1 <- get.mat(ids.u, colnames(m))
  for (x in ids.u) {
    b <- is.element(ids, x)
    if (sum(b) == 1) {
      m1[x, ] <- m[b, ]
    }
    else {
      m1[x, ] <- f(m[b, ])
    }
  }
  return(m1)
}

#' Cap Data
#'
#' Cap the values of data (matrix, data.frame, or list)
#'
#' @param X a list-type object: matrix, data.frame, or simple list
#' @param q quantile to cap at (e.g. 0.05 is capping at the bottom 5th and top 95th quantiles)
#' @export
#' 
cap_object <- function (X, q) {
  ceil_q <- 1 - q
  ceil <- quantile(X, ceil_q)
  floor <- quantile(X, q)
  X[X > ceil] <- ceil
  X[X < floor] <- floor
  return(X)
}
