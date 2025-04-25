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

#' Overall Expression Subroutine
#'
#' code that computes overall expression
#' signature for each cell-type
#'
#' @param r data list
#' @param sig signature list
#' @return overall expression scores
#' @export
#'
get.OE1 <- function (r, sig) {
  if (is.list(sig)) {
    scores <- t(plyr::laply(sig, function(g) get.OE1(r, g)))
    rownames(scores) <- r$cells
    colnames(scores) <- names(sig)
    return(scores)
  }
  g <- sig
  b <- is.element(r$genes, g)
  assertthat::is.string(rownames(r$binZ)[1])
  n1 <- plyr::laply(rownames(r$binZ), function(x) sum(b[r$genes.dist.q ==
                                                          x]))
  rand.scores <- t(r$binZ) %*% n1
  if (sum(b) == 1) {
    raw.scores <- r$zscores[b, ]
  }
  else {
    raw.scores <- colSums(r$zscores[b, ])
  }
  scores <- (raw.scores - rand.scores)/sum(b)
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


#' t-test for a matrix
#'
#' performing t-tests on a per-row basis
#'
#' @param m matrix
#' @param b boolean
#' @param two.sided whether to perform two-sided or not
#' @param rankf ??
#' @param fold.changeF ??
#' @returns pvalues in a table
#' @export
t.test.mat<-function(m,b,two.sided=F,rankf = F,fold.changeF = T){
  if(length(b)!=ncol(m)){
    print("Error. Inconsistent no. of samples.")
    return()
  }
  if(sum(b)<2||sum(!b)<2){
    return(get.mat(rownames(m),c('more','less',"zscores")))
  }
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) t.test(x[b],x[!b])$p.value))
  }else{
    p<-t(apply(m,1,function(x) c(t.test(x[b],x[!b],alternative = 'greater')$p.value,
                                 t.test(x[b],x[!b],alternative = 'less')$p.value)))
    colnames(p)<-c('more','less')
    p<-cbind(p,get.p.zscores(p))
    colnames(p)[3]<-"zscores"
  }
  if(rankf){
    p<-cbind(p,rank(p[,1]),rank(p[,2]))
    colnames(p)[4:5]<-c("rank.more","rank.less")
  }
  if(fold.changeF){
    p<-cbind.data.frame(p,pos.mean = rowMeans(m[,b]),neg.mean = rowMeans(m[,!b]))
    p$logFC<-log2(p$pos.mean/p$neg.mean)
  }

  return(p)
}

#' Convert p-values to z-scores
#'
#' samples _size_ number of cells from each
#' labelled population
#'
#' @param p a table of p values where p[,1] is the greater than hypothesis, and the p[,2] is the less than hypothesis
#' @return _z_-scores
#' @export
get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}
