library(drgem, quietly = T)
cd = readRDS("../data/tubo-ovarian_st_sub_cd.rds")
meta.data = readRDS("../data/tubo-ovarian_st_sub_meta.rds")

iobj <- drgem_phaseI(cd,
                     meta.data,
                     preprocess = T,
                     scale.data = NULL,
                     n_dims = 20,
                     params1 = list("nfeatures" = 960),
                     params2 = list("res" = 0.12, "min.dist" = 0.1,
                                    "k.param" = 20, "return.model" = F))
visualize_simpleplot(iobj,
                     true_labels = "cell.types", 
                     weak_labels = "weak_clusters",
                     custom.colors = list("#56B4E9", "#E69F00",
                                          "#D55E00", "#009E73","#A869DB",
                                          "#000000"))
# saveRDS(iobj, "../../faircon/tubo-ovarian_st_sub_iobj.rds")
iobj = readRDS("../../faircon/tubo-ovarian_st_sub_iobj.rds")
# iobj
query = c(1, 2, 3, 4, 5, 6)
reprocess = F
n_iter = 100
n_sub = 1000
n_dims = 20
field = "weak_clusters"
n_features = 960
params=list("res" = 0.12,
            "min.dist" = 0.1,
            "k.param" =10,
            "return.model" = F)
align = T
autoconvergence = T
r_prime = 10
delta_target = 0.05

out = drgem_phaseII(iobj, 
              query = c(1, 2, 3, 4, 5, 6),
              reprocess = F,
              n_iter = 100,
              n_sub = 1000,
              n_dims = 20,
              field = "weak_clusters",
              n_features = 960,
              params=list("res" = 0.12,
                          "min.dist" = 0.1,
                          "k.param" =10,
                          "return.model" = F),
              align = T,
              autoconvergence = T,
              r_prime = 10,
              delta_target = 0.05)
