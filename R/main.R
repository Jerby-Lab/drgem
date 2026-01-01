#' Distributionally robust initialization
#'
#' Use PCA as diagnostic tool for a distributionally robust initialization
#'
#' @param X Seurat object to be processed.
#' @param meta Vector of class labels to subsample.
#' @param preprocess Logical, whether to reprocess the Seurat object during each iteration.
#' @param scale.data Number of iterations for subsampling.
#' @param n_dims Number of subsamples per class.
#' @param beta Number of principal components to use for PCA.
#' @param viz Logical, whether to visualize data (default = T)
#' @param var_features Number of features to use for variable feature selection.
#' @param align Logical, whether to align clusters during each iteration (default = T)
#' @param params List of additional parameters for clustering and UMAP.
#' @return seurat object with weak labels
#' @export
drgem_phaseI <- function(X,
                         meta,
                         preprocess = T,
                         scale.data = NULL,
                         n_dims = 20,
                         beta = 0.3,
                         viz = T,
                         var_features = NULL,
                         align = T,
                         field1 = "labels",
                         field2 = "seurat_clusters",
                         field_aligned = "weak_clusters",
                         params1 = list("nfeatures" = 200),
                         params2 = list("res" = 0.12, "min.dist" = 0.1,
                                        "k.param" = 20, "return.model" = F)
) {
  # diagnostic PCA
  obj <- initialize_pca(X,
                        meta = meta,
                        preprocess = preprocess,
                        scale.data = scale.data,
                        n_dims = n_dims,
                        params = params1,
                        var_features = var_features)

  # fit model on high risk cells to obtain E(init)
  high_recon_cells = row.names(
    dplyr::filter(obj@meta.data,
                  recon_loss > quantile(obj@meta.data$recon_loss, 1-beta)))
  hrlobj <- initialize_pca(X = X[,high_recon_cells],
                           meta = meta[high_recon_cells,],
                           preprocess = F,
                           scale.data = obj@assays$RNA$scale.data[,high_recon_cells],
                           n_dims = n_dims,
                           var_features = var_features,
                           params = params1)

  # transfer full dataset to E(init) to obtain E(1)
  iobj <- transfer_embedding(hrlobj, obj)
  # obtain c(1)
  iobj <- cluster_and_viz(iobj,
                          n_dims=n_dims,
                          reduction = "hrlpca",
                          viz = viz,
                          align = align,
                          field1 = field1,
                          field2 = field2,
                          field_aligned = field_aligned,
                          params = params2)
  return(iobj)
}

#' Calculate Sample Confidence
#'
#' Iteratively performs sub-sampling and PCA, then assigns confidence scores to clusters.
#'
#' @param iobj Seurat object to be processed.
#' @param query Vector of class labels to subsample.
#' @param reprocess Logical, whether to reprocess the Seurat object during each iteration.
#' @param n_iter Number of iterations for subsampling.
#' @param n_sub Number of subsamples per class.
#' @param n_dims Number of principal components to use for PCA.
#' @param n_features Number of features to use for variable feature selection.
#' @param field Metadata field to use for alignment.
#' @param align Logical, whether to align clusters during each iteration.
#' @param params List of additional parameters for clustering and UMAP.
#' @param autoconvergence Boolean to decide whether to run autoconvergence (default = F),
#' @param r_prime Increment at which to check convergence (default = 100)
#' @param delta_target Minimum change in cell type assignments (default = 0.05)
#' @return A list of two objects: (1) A dataframe with assigned clusters and
#' confidence scores for each cell. and (2) A dataframe with all the labels per
#' iteration.
#' @export
drgem_phaseII <- function(iobj,
                          query = c(1, 2, 3, 4, 5),
                          assign = c(1, 2, 3, 4, 5),
                          # confounding = c(),
                          sigs = NULL,
                          markers_field = "markers.scores",
                          reprocess = F,
                          n_iter = 800,
                          n_sub = 600,
                          n_dims = 30,
                          n_features = 155,
                          field = "predicted_clusters",
                          align = T,
                          params = list("res" = 1.0,
                                        "min.dist" = 0.1,
                                        "k.param"=10,
                                        "return.model"= F),
                          var_features = row.names(iobj),
                          autoconvergence = F,
                          r_prime = 100,
                          delta_target = 0.05
)
{
  # iterative sub-sampling based on equal density of weak labels
  if (autoconvergence){
    print("Running Autoconvergence...")
    n_batch = round(c(n_iter/r_prime))
    d = rep(1, length(query)) # initialize delta
    batch = 1 # start batch
    df = data.frame(row.names = colnames(iobj)) # initialize iteration data.frame
    while (sum(d >= delta_target) != 0 & batch <= n_batch) {
      print(batch)
      # run batch
      iters = 1:r_prime + r_prime*(batch - 1)
      df_tmp = parallel::mclapply(iters, mc.cores = 10, function(i){
        set.seed(i)
        # run dimension reduction on balanced sub-sample
        subsample <- subsample_by_class(query, iobj, n_sub, field = field, seed =i)
        temp_obj <- initialize_pca(X = iobj[["RNA"]]$counts[,subsample],
                                   meta = iobj@meta.data[subsample,],
                                   preprocess = F,
                                   scale.data = iobj[["RNA"]]$scale.data[,subsample],
                                   n_dims = n_dims,
                                   params = list("nfeatures" = n_features),
                                   var_features = var_features)

        # run clustering on balanced sub-sample
        temp_obj = cluster_and_viz(temp_obj,
                                   n_dims = n_dims,
                                   reduction = "pca",
                                   viz = F,
                                   align = align,
                                   field1 = field,
                                   field2 = "seurat_clusters",
                                   field_aligned = "predicted_clusters",
                                   params = params)

        final = rep(NA, dim(iobj)[2])
        names(final) = colnames(iobj)
        final[colnames(temp_obj)] <- as.character(temp_obj@meta.data[["predicted_clusters"]])
        return(final)
      })
      df_tmp1 = do.call("cbind.data.frame",
                        df_tmp[which(unlist(lapply(df_tmp, length)) == dim(iobj)[2])])
      df = cbind.data.frame(df, df_tmp1)

      # assess batch
      assignments <- t(apply(df, 1, label_consensus, query=query))
      assignments = data.frame(assignments)
      assignments$confidence <- as.numeric(assignments$confidence)
      assignments$n <- as.numeric(assignments$n)
      row.names(assignments) <- colnames(iobj)

      # if each cell has been sampled at least 5 times
      if (min(assignments$n) > 5) {
        # get last batch results
        last_iter = r_prime + r_prime*(batch - 2)
        last_ass <- t(apply(df[,1:last_iter], 1, label_consensus, query=query))
        last_ass = data.frame(last_ass)
        last_ass$confidence <- as.numeric(last_ass$confidence)
        last_ass$n <- as.numeric(last_ass$n)
        row.names(last_ass) <- colnames(iobj)

        # get d
        d = unlist(lapply(query, function(x){
          tmp = dplyr::filter(assignments, assignment == x)
          tmp1 = last_ass[row.names(tmp),]
          sum(tmp$assignment != tmp1$assignment)/length(tmp$assignment)
        }))
      }
      batch = batch + 1
    }
    # convergence statement
    if (batch > n_batch) {
      print(paste0("DR-GEM Phase II ran the specified ", n_iter, " iterations"))
    } else {
      print(paste0("DR-GEM Phase II converged in ", batch*r_prime, "iterations"))
    }
    colnames(df) <- NULL
  } else {
    df <- parallel::mclapply(1:n_iter, mc.cores=10, function(i){
      set.seed(i)

      # run dimension reduction on balanced sub-sample
      subsample <- subsample_by_class(query, iobj, n_sub, field = field, seed =i)
      temp_obj <- initialize_pca(X = iobj[["RNA"]]$counts[,subsample],
                                 meta = iobj@meta.data[subsample,],
                                 preprocess = F,
                                 scale.data = iobj[["RNA"]]$scale.data[,subsample],
                                 n_dims = n_dims,
                                 params = list("nfeatures" = n_features),
                                 var_features = var_features)

      # run clustering on balanced sub-sample
      temp_obj = cluster_and_viz(temp_obj,
                                 n_dims = n_dims,
                                 reduction = "pca",
                                 viz = F,
                                 align = align,
                                 field1 = field,
                                 field2 = "seurat_clusters",
                                 field_aligned = "predicted_clusters",
                                 params = params)

      final = rep(NA, dim(iobj)[2])
      names(final) = colnames(iobj)
      final[colnames(temp_obj)] <- as.character(temp_obj@meta.data[["predicted_clusters"]])
      return(final)
    })
    df = do.call("cbind.data.frame", df[which(unlist(lapply(df, length)) == dim(iobj)[2])])
    colnames(df) <- NULL
  }

  # do score aggregation
  assignments <- t(apply(df, 1, label_consensus, query=query))
  assignments = data.frame(assignments)
  assignments$confidence <- as.numeric(assignments$confidence)
  assignments$n <- as.numeric(assignments$n)
  row.names(assignments) <- colnames(iobj)

  # get counts
  counts <- data.frame(t(apply(df, 1, count_assignments, query=query)))
  row.names(counts) <- colnames(iobj)

  #merge
  assignments = cbind(assignments, counts)

  #return
  return(list(assigments = assignments,
              iterations = df))
}

#' Perform Final Embedding
#'
#' Performs final embedding on high-confidence subsampled cells.
#'
#' @param iobj Seurat object to be processed.
#' @param assignments Dataframe with assignment labels and confidence scores.
#' @param thres Confidence threshold for selecting samples.
#' @param thres_by_class Class-specific confidence thresholds.
#' @param n_dims Number of principal components to use for the final embedding.
#' @param subsample_density List indicating the number of samples to subsample per class.
#' @return The Seurat object with the final embedding and cluster visualization.
#' @export
drgem_phaseIII <- function(iobj,
                           assignments,
                           thres = 0.95,
                           thres_by_class = NULL,
                           n_dims = 20,
                           n_features = 200,
                           subsample_density = c("1" = 800,
                                                 "2" = 500,
                                                 "3" = 450,
                                                 "4" = 400,
                                                 "5" = 50),
                           params=list("res" = 0.18, "min.dist" = 0.1, "k.param" = 10),
                           return_ref=F){
  # threshold on confidence
  if (is.null(thres) & is.null(thres_by_class)) {
    cells = row.names(assignments)
  } else if (!is.null(thres) & !is.null(thres_by_class)) {
    print("ambiguous thresholding input...")
    return()
  } else if (!is.null(thres) & is.null(thres_by_class)) {
    cells = row.names(dplyr::filter(assignments, confidence > thres))
  } else {
    cells = do.call("c", lapply(names(thres_by_class), function(x){
      xcells = row.names(dplyr::filter(assignments,
                                       assignment == x,
                                       confidence >= thres_by_class[x]))
    }))
  }

  # subsample
  high_conf_cells <- subsample_per_class(assignments,
                                         cells,
                                         subsample_density)

  # fit on high confidence subsample
  refobj <- initialize_pca(X = iobj[["RNA"]]$counts[,high_conf_cells],
                           meta = iobj@meta.data[high_conf_cells,],
                           preprocess = T,
                           scale.data = NULL,
                           n_dims = n_dims,
                           params = list("nfeatures" = n_features))

  X_scaled = refobj[["RNA"]]$scale.data
  U = Seurat::Loadings(refobj, reduction = "pca")
  feats = intersect(row.names(X_scaled), row.names(U))
  recon_loss <- calc_recon_loss(X_scaled, U, n_dims)
  refobj@meta.data$recon_loss <- recon_loss
  refobj@meta.data$final_predictions <- assignments[colnames(refobj),]$assignment
  refobj@meta.data$confidence_scores <- assignments[colnames(refobj),]$confidence_scores

  refobj <- Seurat::RunUMAP(refobj, dims =1:n_dims,return.model = T)

  anchors <- Seurat::FindTransferAnchors(
    reference = refobj,
    features = row.names(iobj),
    query = iobj,
    reference.reduction = "pca",
    dims = 1:n_dims
  )

  query <- Seurat::MapQuery(
    anchorset = anchors,
    query = iobj,
    reference = refobj,
    reference.reduction = "pca",
    reduction.model = "umap"
  )

  query@meta.data$final_predictions <- assignments$assignment
  query@meta.data$confidence_scores <- assignments$confidence_scores
  out <- list(query, refobj)
  return(out)
}

#' scRNA Marker Based Cell Type Assignments
#'
#' Using gene expression markers, unsupervised clusters,
#' and gene expression, assign the likely celltype to those
#' clusters
#'
#' @param r data list including the gene expression data
#' @param cell.clusters a _n_-length list to where _n_ is the number of cells
#' @param cell.sig cell type signatures (list of lists of markers associated with each cell-type, default exists)
#' @param immune.markers immune markers (default exists)
#' @param non.immune.cell.types which cell types are not immune cell types?
#' @param EM.flag run EM? (default = F)
#' @param OE.type which version of overall-expression to calculate (default = "V1")
#' @param test.type which kind of test to test marker significance between clusters
#' @param minZ minimum significance cutoff to assign a cell-type (default = 10)
#' @return a list with objects named "c" and "m" corresponding to your subsampled counts and metadata
#' @export
drgem_phaseIV <- function(obj,
                          cell.clusters,
                          cell.sig,
                          EM.flag = F,
                          OE.type = "V1",
                          test.type = "ttest",
                          minZ = 10,
                          subsample = T){

  if(missing(cell.sig)){
    print("cell signatures missing")
    return()
  }

  # get marker scores
  r <- c()
  r$cd <- obj@assays$RNA$counts
  r$cells <- colnames(obj)
  r$genes <- row.names(obj)
  r$comp.reads <- apply(r$cd, 2, sum)
  r$tpm <- get.tpm(X = r$cd, comp.reads = r$comp.reads)
  r$clusters <- cell.clusters
  r <- prep4OE(r, n.cat = 20)
  r$markers.scores <- get.OE(r, sig = cell.sig)

  # run pairwise t-test analysis
  r<-scRNA_cluster.annotation.ttest(r,cell.clusters = cell.clusters,minZ = minZ)

  # save results
  obj@meta.data <- cbind(obj@meta.data, apply(r$markers.scores, 2, cap_object, q = 0.05))
  obj@meta.data$final_annotations <- r$cell.types
  out <- list(obj, r)
  return(out)
}

#' Add _r_'' iterations to DR-GEM Phase II
#'
#' Functionality for adding _r_'' iterations
#'
#' @param df assignments across iterations
#' @param r_pp number of additional iterations to run
#' @param iobj Seurat object to be processed.
#' @param query Vector of class labels to subsample.
#' @param reprocess Logical, whether to reprocess the Seurat object during each iteration.
#' @param n_iter Number of iterations for subsampling.
#' @param n_sub Number of subsamples per class.
#' @param n_dims Number of principal components to use for PCA.
#' @param n_features Number of features to use for variable feature selection.
#' @param field Metadata field to use for alignment.
#' @param align Logical, whether to align clusters during each iteration.
#' @param params List of additional parameters for clustering and UMAP.
#' @return A list of two objects: (1) A dataframe with assigned clusters and
#' confidence scores for each cell. and (2) A list with the proportion of
#' change of each cell cluster after the additional _r_'' iterations.
#' @export
add_iterations_drgemII <- function(df, r_pp, iobj,
                             query = c(1, 2, 3, 4, 5),
                             assign = c(1, 2, 3, 4, 5),
                             # confounding = c(),
                             sigs = NULL,
                             markers_field = "markers.scores",
                             reprocess = F,
                             n_iter = 800,
                             n_sub = 600,
                             n_dims = 30,
                             n_features = 155,
                             field = "predicted_clusters",
                             align = T,
                             params = list("res" = 1.0,
                                           "min.dist" = 0.1,
                                           "k.param"=10,
                                           "return.model"= F),
                             var_features = row.names(iobj)) {
  prev = dim(df)[2]
  iters = 1:r_pp + prev
  df_add <- parallel::mclapply(iters, mc.cores=10, function(i){
    set.seed(i)

    # run dimension reduction on balanced sub-sample
    subsample <- subsample_by_class(query, iobj, n_sub, field = field, seed =i)
    temp_obj <- initialize_pca(X = iobj[["RNA"]]$counts[,subsample],
                               meta = iobj@meta.data[subsample,],
                               preprocess = F,
                               scale.data = iobj[["RNA"]]$scale.data[,subsample],
                               n_dims = n_dims,
                               params = list("nfeatures" = n_features),
                               var_features = var_features)

    # run clustering on balanced sub-sample
    temp_obj = cluster_and_viz(temp_obj,
                               n_dims = n_dims,
                               reduction = "pca",
                               viz = F,
                               align = align,
                               field1 = field,
                               field2 = "seurat_clusters",
                               field_aligned = "predicted_clusters",
                               params = params)

    final = rep(NA, dim(iobj)[2])
    names(final) = colnames(iobj)
    final[colnames(temp_obj)] <- as.character(temp_obj@meta.data[["predicted_clusters"]])
    return(final)
  })
  df_add = do.call("cbind.data.frame",
                   df_add[which(unlist(lapply(df_add, length)) == dim(iobj)[2])])
  colnames(df_add) <- NULL
  df = cbind(df, df_add)

  assignments <- t(apply(df, 1, label_consensus, query=query))
  assignments = data.frame(assignments)
  assignments$confidence <- as.numeric(assignments$confidence)
  assignments$n <- as.numeric(assignments$n)
  row.names(assignments) <- colnames(iobj)

  last_ass <- t(apply(df[,1:prev], 1, label_consensus, query=query))
  last_ass = data.frame(last_ass)
  last_ass$confidence <- as.numeric(last_ass$confidence)
  last_ass$n <- as.numeric(last_ass$n)
  row.names(last_ass) <- colnames(iobj)

  # get d
  d = unlist(lapply(query, function(x){
    tmp = dplyr::filter(assignments, assignment == x)
    tmp1 = last_ass[row.names(tmp),]
    sum(tmp$assignment != tmp1$assignment)/length(tmp$assignment)
  }))

  return(list(assignments = assignments,
              deltas = d))
}

#' Compute _t_-tests
#'
#' function for computing two sample _t_-test
#' on single cell clustered data
#'
#' @param r data list including the gene expression data
#' @param cell.clusters a _n_-length list of cluster assignments, where _n_ = number of cells
#' @param minZ cutoff for significance in the -log10(p) space (default = 10)
#' @return data list, but now including `cell.types`, and `assignments.ttests` slots
#'
scRNA_cluster.annotation.ttest <- function(r,
                                           cell.clusters,
                                           minZ = 10,
                                           subsample = F) {
  if (subsample) {
    b<-sample.per.label(r$clusters,500)
  }

  b <- rep(TRUE, length(r$clusters))

  z<-plyr::laply(unique(r$clusters),function(x) t.test.mat(t(r$markers.scores[b,]),r$clusters[b]==x)[,3])
  rownames(z)<-unique(r$clusters)
  colnames(z) <- colnames(r$markers.scores)
  z<-cbind.data.frame(z,cell.type1 = colnames(z)[apply(z,1,function(x) which(x==max(x))[1])],
                      n = rowSums(z>minZ))
  b2<-z$n>1
  z$diff<-apply(z[,1:(ncol(z)-2)],1,function(x) sort(x,decreasing = T)[1])-apply(z[,1:(ncol(z)-2)],1,function(x) sort(x,decreasing = T)[2])
  z$unique<-z$n==1|z$diff>10
  z$cell.type<-z$cell.type1
  z$cell.type[!z$unique]<-"UD"
  r$cell.types<-z[r$clusters,"cell.type"]

  r$assignments.ttests<-z
  return(r)
}


#' Initialize PCA and Seurat Object
#'
#' This function initializes a Seurat object from input data, optionally preprocesses it,
#' performs PCA, and calculates the reconstruction loss.
#'
#' @param X Matrix of gene expression counts.
#' @param meta Metadata for the cells. If NULL, only expression data will be used.
#' @param preprocess Logical, if TRUE, data will be normalized, variable features selected, and scaled.
#' @param scale.data Scaled data to be used directly if preprocessing is FALSE.
#' @param n_dims Number of principal components to compute.
#' @param params List of additional parameters, such as the number of features for variable feature selection (`nfeatures`).
#' @param var_features List of variable features, if predetermined
#' @return A Seurat object with PCA results and reconstruction loss added to the metadata.
#' @export
initialize_pca <- function(X,
                           meta,
                           preprocess = T,
                           scale.data = NULL,
                           n_dims = 30,
                           params = list("nfeatures" = 200),
                           var_features = NULL){
  # Create Seurat Object
  if (is.null(meta)) {
    obj <- Seurat::CreateSeuratObject(X)
  } else {
    obj <- Seurat::CreateSeuratObject(X, meta.data = meta)
  }

  # Preprocess Data
  if (preprocess) {
    obj <- Seurat::NormalizeData(obj, verbose = F)
    if (params[["nfeatures"]] < dim(X)[1]) {
      if (is.null(var_features)) {
        tmp = Seurat::FindVariableFeatures(obj, nfeatures = params[["nfeatures"]])
        var_features = VariableFeatures(tmp)
        print(paste0("n var features: ", length(var_features)))
      }
    }
    obj <- Seurat::ScaleData(obj, verbose = F)
  } else {
    obj@assays$RNA$scale.data = scale.data
  }

  # Run PCA
  if (params[["nfeatures"]] < dim(X)[1]) {
    obj <- Seurat::RunPCA(obj, features = var_features, verbose = F)
  } else {
    obj <- Seurat::RunPCA(obj, features = row.names(obj), verbose = F)
  }

  # Calculate Reconstruction Error
  X_scaled = obj@assays$RNA$scale.data
  U = Seurat::Loadings(obj, reduction = "pca")
  recon_loss <- calc_recon_loss(X_scaled[row.names(U),], U, n_dims)
  obj@meta.data$recon_loss = recon_loss

  return(obj)
}

#' Cluster and Visualize Data
#'
#' This function clusters cells using PCA or other dimensional reduction methods and visualizes the results using UMAP.
#' Optionally aligns clusters between fields.
#'
#' @param obj Seurat object to be processed.
#' @param n_dims Number of dimensions to use for clustering and UMAP.
#' @param reduction Dimensionality reduction method to use (default is "pca").
#' @param clus Logical, if TRUE, clusters will be found using `FindNeighbors` and `FindClusters`.
#' @param viz Logical, if TRUE, UMAP will be run and visualized.
#' @param align Logical, if TRUE, clusters will be aligned between two fields.
#' @param field1 First metadata field for alignment (default is "labels").
#' @param field2 Second metadata field for alignment (default is "seurat_clusters").
#' @param field_aligned Field where aligned clusters are stored (default is "predicted_clusters").
#' @param params List of additional parameters, such as resolution for clustering and minimum distance for UMAP.
#' @return The Seurat object with clustering and alignment information added.
#' @export
cluster_and_viz <- function(obj,
                            n_dims = 10,
                            reduction = "pca",
                            clus = T,
                            viz = T,
                            align = T,
                            field1 = "labels",
                            field2 = "seurat_clusters",
                            field_aligned = "weak_clusters",
                            params = list("res" = 0.12, "min.dist" = 0.1,
                                          "k.param" = 20, "return.model" = F)
){
  if(clus){
    if (!("RNA_snn" %in% names(obj@graphs))) {
      obj <- Seurat::FindNeighbors(obj,
                                   reduction = reduction,
                                   dims = 1:(n_dims),
                                   verbose = F,
                                   k.param = params[["k.param"]])
    }
    obj <- Seurat::FindClusters(obj, resolution = params[["res"]])
  }

  if (viz){obj <- Seurat::RunUMAP(obj,
                                  dims = 1:(n_dims),
                                  reduction = reduction,
                                  min.dist = params[["min.dist"]],
                                  return.model = params[["return.model"]])}

  if (align){
    mapping = map_clusters1(obj@meta.data[[field1]], obj@meta.data[[field2]])
    obj@meta.data[[field_aligned]] = factor((mapping[obj@meta.data$seurat_clusters]),
                                            levels(obj@meta.data[[field1]]))
  }

  return(obj)
}

#' Transfer Embeddings Between Seurat Objects
#'
#' Transfers embeddings from a source Seurat object to a target Seurat object.
#'
#' @param sobj Source Seurat object from which to transfer PCA embeddings.
#' @param tobj Target Seurat object to which the embeddings are transferred.
#' @param n_dims Number of dimensions to use in the reconstruction.
#' @param semb Embedding reduction in the source object (default is "pca").
#' @param temb Embedding reduction in the target object (default is "hrlpca").
#' @param tkey Key for storing transferred embeddings in the target object.
#' @param trecon Field name for storing reconstruction loss in the target object.
#' @return The updated target Seurat object with transferred embeddings and reconstruction loss.
#' @export
transfer_embedding <- function(sobj,
                               tobj,
                               n_dims = 10,
                               semb = "pca",
                               temb = "hrlpca",
                               tkey = "hRLPC_",
                               trecon = "recon_loss_hrlpca"){
  sU = Seurat::Loadings(sobj, reduction = semb)
  tX_scaled = tobj[["RNA"]]$scale.data
  feats = intersect(row.names(sU), row.names(tX_scaled))
  tobj[[temb]] <- Seurat::CreateDimReducObject(embeddings = t(tX_scaled)[,feats] %*% sU[feats,],
                                               key = tkey,
                                               assay = 'RNA')
  recon_loss <- calc_recon_loss(tX_scaled, sU, n_dims)
  tobj@meta.data[[trecon]] = recon_loss
  return(tobj)
}

calculate_ct_confidence <- function(iobj,
                                    query = c(1, 2, 3, 4, 5),
                                    assign = c(1, 2, 3, 4, 5),
                                    sigs = NULL,
                                    markers_field = "markers.scores",
                                    reprocess = F,
                                    n_iter = 10,
                                    n_sub = 500,
                                    n_dims = 10,
                                    n_features = 200,
                                    field = "ct_init",
                                    align = F,
                                    params = list("res" = 0.12, "min.dist" = 0.1, "k.param" = 20)
){
  # iterative sub-sampling based on equal density and reconstruction loss
  df <- do.call("cbind", parallel::mclapply(1:n_iter, mc.cores = 10, function(i){
    set.seed(i)
    print(i)

    # Run dimension reduction on sub-sample
    subsample <- subsample_by_class(query, iobj, n_sub, field = field, seed =i)
    temp_obj <- subset(iobj, cellids %in% subsample)
    temp_obj <- FindVariableFeatures(temp_obj, nfeatures = n_features)
    temp_obj <- ScaleData(temp_obj, features = VariableFeatures(temp_obj), verbose = F)
    temp_obj <- RunPCA(temp_obj, features = VariableFeatures(temp_obj), verbose = F)

    # extract recon_loss
    X_scaled = temp_obj[["RNA"]]@scale.data
    U = Seurat::Loadings(temp_obj, reduction = "pca")
    recon_loss <- calc_recon_loss(X_scaled[row.names(U),], U, n_dims)
    temp_obj@meta.data$recon_loss = recon_loss

    temp_obj = cluster_and_viz(temp_obj,
                               n_dims = n_dims,
                               reduction = "pca",
                               viz = F,
                               align = align,
                               params = params)

    ##################################
    # Annotation based on signatures #
    ##################################

    temp_obj$seurat_clusters <- paste0("C", temp_obj$seurat_clusters)
    s <- c()
    s$cd <- as.matrix(temp_obj[["RNA"]]@counts)
    s$genes <- row.names(temp_obj)
    s$cells <- colnames(temp_obj)
    s$comp.reads <- apply(as.matrix(temp_obj[["RNA"]]@counts), 2, sum)
    s$tpm <- as.matrix(temp_obj[["RNA"]]@counts)
    s$clusters <- temp_obj$seurat_clusters

    s <- scRNA_markers_subtype_assign(s,
                                      cell.clusters = s$clusters,
                                      cell.sig = sigs[assign], resolve_ties = T,
                                      minZ = 10)
    print(table(temp_obj$targets, s$cell.types))
    final = rep(NA, dim(iobj)[2])
    names(final) = colnames(iobj)
    final[colnames(temp_obj)] <- s$cell.types
    return(final)
  }))
  # saveRDS(df, "~/df_tmp.rds")
  # do score aggregation

  assignments <- t(apply(df, 1, label_consensus_UD_NTC, query=assign))
  # assignments <- t(apply(df, 1, label_consensus, query=assign))
  assignments = data.frame(assignments)
  assignments$confidence <- as.numeric(assignments$confidence)
  assignments$n <- as.numeric(assignments$n)
  # row.names(assignments) <- colnames(iobj)

  # get counts
  counts <- data.frame(t(apply(df, 1, count_assignments, query=assign)))
  # row.names(counts) <- colnames(iobj)

  #merge
  assignments = cbind(assignments, counts)
  return(assignments)
}

#' Map Clusters
#'
#' Custom function that maps clusters between two cluster assignments based on
#' a contingency table; includes some edge cases.
#'
#' @param a vector of cluster assignments from the first method
#' @param b vector of cluster assignments from the second method
#' @return named vector indicating the mapped cluster assignments
#' @export
map_clusters <- function(a, b){
  contingency_matrix = table(a, b)
  cluster_mapping <- apply(contingency_matrix, 2, function(column) {
    rownames(contingency_matrix)[which.max(column)]})

  # If there are duplicates:
  dups = cluster_mapping[which(duplicated(cluster_mapping))]
  for (dup in dups){
    dups_tmp = cluster_mapping[cluster_mapping == dup]
    maximum = max(contingency_matrix[dup, names(dups_tmp)])
    losers = contingency_matrix[dup, names(dups_tmp)][
      contingency_matrix[dup, names(dups_tmp)] != maximum]
    pad = as.numeric(max(cluster_mapping))
    for (loser in names(losers)) {
      pad = pad + 1
      cluster_mapping[loser] = pad
    }
  }
  return(cluster_mapping)
}

#' Map Clusters (Kuhn-Munkres)
#'
#' Function that maps clusters between two cluster assignments using Kuhn-Munkres
#'
#' @param a vector of cluster assignments from the first method
#' @param b vector of cluster assignments from the second method
#' @return named vector indicating the mapped cluster assignments
#' @export
map_clusters1 <- function(a, b){
  contingency_matrix = table(b, a)
  solve = lpSolve::lp.assign(contingency_matrix, direction="max")
  map = solve$solution
  colnames(map) = colnames(contingency_matrix)
  row.names(map) = row.names(contingency_matrix)
  key = apply(map, 1, function(x){
    names(x)[x == "1"]
  })
  return(key)
}

#' Visualize Simple Plot
#'
#' Function that visualizes clustering results using multiple plots including PCA, UMAP, and reconstruction loss heatmap.
#'
#' @param obj Seurat object
#' @param true_labels string indicating the column name for reconstruction loss (default = "recon_loss")
#' @param predicted_labels string indicating the column name for cluster assignments (default = "seurat_clusters")
#' @param custom.colors list of custom color mappings (default = NULL)
#' @return None, prints plots to console
#' @export
visualize_simpleplot <- function(obj,
                                 true_labels = "labels",
                                 weak_labels = "seurat_clusters",
                                 custom.colors = NULL,
                                 type1 = "discrete",
                                 type2 = "discrete",
                                 reduction = "umap",
                                 midpoint = NULL){
  if (type1 == "cont") {
    a = Seurat::FeaturePlot(obj,
                                features = true_labels, reduction = "umap")
    b = b+ggplot2::scale_color_gradient2(low = '#A01C00',
                                         mid = '#f6edbd',
                                         high = '#009392',
                                         midpoint = midpoint)
  } else {
    a = Seurat::DimPlot(obj, group.by = true_labels, reduction = reduction)
    if (!is.null(custom.colors)) {
      a = a + ggplot2::scale_color_manual(values = custom.colors)
    }
  }

  if (type2 == "cont") {
    b = Seurat::FeaturePlot(obj,
    features = weak_labels, reduction = "umap")
    b = b+ggplot2::scale_color_gradient2(low = '#A01C00',
                                     mid = '#f6edbd',
                                     high = '#009392',
                                     midpoint = midpoint)
  } else {
    b = Seurat::DimPlot(obj, group.by = weak_labels, reduction = reduction)
    if (!is.null(custom.colors)) {
      b = b + ggplot2::scale_color_manual(values = custom.colors)
    }
  }
  print(a+b)
}

#' Visualize Simple Plot
#'
#' Function that visualizes clustering results using multiple plots including PCA, UMAP, and reconstruction loss heatmap.
#'
#' @param iobj Seurat object
#' @param assignments data frame of consensus matrix from phase II!
#' @param field1 string indicating the column name for the true labels
#' @param cluster_field string indicating the column name final predicted label
#' @return None, saves the multiplot visualization to a file
#' @export
plot_cm <- function(iobj,
                    assignments,
                    field1 = "labels",
                    field2 = "assignment") {
  count = table(iobj@meta.data[[field1]],
                assignments[[field2]])
  freq = table(iobj@meta.data[[field1]])
  prop <- sweep(count, 1, freq, "/")

  d = ggplot2::ggplot(data.frame(prop),
                      ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = '#009392',
                                  mid = '#f6edbd',
                                  high = '#A01C00',
                                  midpoint = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Freq, 2)), color = "white", size = 4) +
    ggplot2::ylab("Predicted") +
    ggplot2::xlab("True") +
    ggplot2::coord_flip()

  print(d)
}

#' Visualize Simple Plot
#'
#' Function that visualizes clustering results using multiple plots including PCA, UMAP, and reconstruction loss heatmap.
#'
#' @param fobj Seurat object
#' @param features features to plot
#' @param ncol No. of columns of plots (default = 3)
#' @return None, prints to console
#' @export
plot_sigs <- function(fobj,
                      features,
                      ncol = 3) {

  a = Seurat::FeaturePlot(fobj, features = features, ncol = ncol)
  a = a & ggplot2::scale_color_gradient2(low = '#009392',
                                         mid = '#f6edbd',
                                         high = '#A01C00',
                                         midpoint = 0)
  print(a)
}


#' Visualize Multiplot
#'
#' Function that visualizes clustering results using multiple plots including PCA, UMAP, and reconstruction loss heatmap.
#'
#' @param obj Seurat object
#' @param name string for naming the output file
#' @param recon_field string indicating the column name for reconstruction loss (default = "recon_loss")
#' @param cluster_field string indicating the column name for cluster assignments (default = "seurat_clusters")
#' @return None, saves the multiplot visualization to a file
#' @export
visualize_multiplot <- function(obj,
                                name = "",
                                recon_field = "recon_loss",
                                cluster_field = "seurat_clusters",
                                pub = T,
                                custom.colors = NULL,
                                height = 8,
                                width = 12){
  obj@meta.data$labels <- factor(obj@meta.data$labels)
  a = Seurat::DimPlot(obj, group.by = "labels", reduction = "pca")
  b = Seurat::DimPlot(obj, group.by = "labels", reduction = "umap")
  c = Seurat::DimPlot(obj, group.by = cluster_field, reduction = "umap")
  if (!is.null(custom.colors)) {
    a = a + ggplot2::scale_color_manual(values = custom.colors)
    b = b + ggplot2::scale_color_manual(values = custom.colors)
    c = c + ggplot2::scale_color_manual(values = custom.colors)
  }
  count = table(obj@meta.data$labels, obj@meta.data[[cluster_field]])
  freq = table(obj@meta.data$labels)
  prop <- sweep(count, 1, freq, "/")

  d = ggplot2::ggplot(data.frame(prop), aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile() +
    ggplot2::scale_fill_gradient2(low = '#009392',
                                  mid = '#f6edbd',
                                  high = '#A01C00',
                                  midpoint = 0.5) +
    geom_text(aes(label = round(Freq, 2)), color = "white", size = 4) +
    ylab("Predicted") +
    xlab("True") +
    coord_flip()

  obj@meta.data$recon_loss_cap <- cap_object(obj@meta.data[[recon_field]], 0.02)
  e = Seurat::FeaturePlot(obj, features = "recon_loss_cap", reduction = "umap")
  e = e+ggplot2::scale_color_gradient2(low = '#009392',
                                       mid = '#f6edbd',
                                       high = '#A01C00',
                                       midpoint =
                                         median(obj@meta.data$recon_loss_cap))
  f = ggplot2::ggplot(obj@meta.data, aes_string(x = "labels",
                                                y = recon_field,
                                                fill = "labels")) +
    ggplot2::geom_boxplot()

  if (!is.null(custom.colors)) {
    a = a + scale_color_manual(values = custom.colors)
    b = b + scale_color_manual(values = custom.colors)
    c = c + scale_color_manual(values = custom.colors)
    f = f + scale_fill_manual(values = custom.colors)
  }
  p = (a + b + c)/(d + e + f)
  png(paste0(name, ".png"), units = "in", height = height, width = width, res = 500)
  print(p)
  dev.off()

  if (pub) {
    pdf(paste0(name, ".pdf"), height = height, width = width)
    print(p)
    dev.off()
  }
}

#' Sub-sample Cells by Class Label
#'
#' This function sub-samples cells from a Seurat object by a specified class label, ensuring a consistent number of samples per class.
#'
#' @param query A vector of class labels to be sampled.
#' @param obj A Seurat object from which to sub-sample cells.
#' @param n_sub The number of cells to sub-sample per class.
#' @param seed A random seed to ensure reproducibility (default is 888).
#' @return A vector of cell names that have been sub-sampled from the Seurat object.
#' @export
subsample_by_class <- function(query, obj, n_sub, field, seed = 888){
  set.seed(seed)
  subsample <- do.call("c", lapply(query, function(x){
    if (sum(obj@meta.data[[field]] == x) < n_sub) {
      n <- sum(obj@meta.data[[field]] == x)
    } else {
      n <- n_sub
    }

    sample(colnames(obj)[obj@meta.data[[field]] == x], n)
  }))
  return(subsample )
}

#' Sub-sample Cells per Class Based on Density
#'
#' This function sub-samples cells based on a specified density for each class, considering assignments and cell counts.
#'
#' @param assignments A dataframe containing cell assignment labels.
#' @param cells A vector of cell names to be considered for sub-sampling.
#' @param subsample_density A named list indicating the number of cells to sub-sample for each class.
#' @param seed A random seed to ensure reproducibility (default is 888).
#' @return A vector of sub-sampled cell names.
#' @export
subsample_per_class <- function(assignments,
                                cells,
                                subsample_density,
                                seed = 888){
  set.seed(seed)
  tmp = assignments[cells,]
  subsample <- do.call("c", lapply(names(subsample_density), function(x){
    if (sum(tmp$assignment == x) < subsample_density[x]) {
      n <- sum(tmp$assignment == x)
    } else {
      n <- subsample_density[x]
    }
    sample(row.names(tmp)[tmp$assignment== x], n)
  }))

  return(subsample)
}

#' Compute Consensus Label and Confidence
#'
#' This function computes the consensus label and confidence score for a set of cluster assignments.
#'
#' @param x A vector of cluster assignments, with possible values from the `query` or "UD" for unassigned.
#' @return A named vector containing the consensus assignment, confidence score, and the number of votes.
#' @export
label_consensus <- function(x, query){
  x <- factor(x, levels = c(query, "UD"))
  calc <- table(x)/sum(table(x))
  ct <- calc[which(calc == max(calc))][1]
  conf <- as.numeric(ct)[1]
  if (!is.na(ct) & names(ct) == "UD") {
    calc1 <- calc[which(names(calc) != "UD")]
    ct <- calc1[which(calc1 == max(calc1))][1]
    conf <- as.numeric(calc[names(ct)])[1]
  }
  return(c(assignment = names(ct), confidence = conf, n = sum(table(x))))
}

#' Compute Consensus Label and Confidence
#'
#' This function computes the consensus label and confidence score for a set of cluster assignments.
#'
#' @param x A vector of cluster assignments, with possible values from the `query` or "UD" for unassigned.
#' @return A named vector containing the consensus assignment, confidence score, and the number of votes.
#' @export
label_consensus_UD_NTC <- function(x, query){
  x[x == "UD"] <- "NTC"
  x <- factor(x, levels = c(query, "UD"))
  calc <- table(x)/sum(table(x))
  ct <- calc[which(calc == max(calc))][1]
  conf <- as.numeric(ct)[1]
  if (!is.na(ct) & (names(ct) == "UD" | names(ct) == "NTC")) {
    calc1 <- calc
    ct <- sum(calc1[c("UD", "NTC")])
    names(ct) <- "NTC"
    conf <- as.numeric(ct)[1]
  }
  return(c(assignment = names(ct), confidence = conf, n = sum(table(x))))
}

count_assignments <- function(x, query){
  x <- factor(x, levels = c(query, "UD"))
  return(table(x))
}
