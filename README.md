# DR-GEM: Distributionally Robust and latent Group-awarE consensus Machine learning for single cell genomics

DR-GEM is a  self-supervised meta-algorithm that combines principles in distributionally robust optimization with balanced consensus machine learning to overcome the challenges of latent class imbalance and non-gaussian technical noise in analysis pipelines featuring dimensionality reduction and clustering. 

# **Quick start**

To install DR-GEM you can eiher run: 
`devtools::install_github(repo = "https://github.com/Jerby-Lab/drgem")` or you can download all the files in this R package and use: 
`devtools::install("drgem")`. Installation time is estimated to be < 1 minute.

We have provided a tutorial that will be available via the following link:
[DR-GEM tutorial on tubo-ovarian cancer spatial transcriptomics dataset](https://htmlpreview.github.io/?https://github.com/Jerby-Lab/drgem/blob/main/vignettes/tutorial_tubo-ovarian.html) or you can download the file in `vignettes/tutorial_tubo-ovarian.html` and view in any web browswer. The toy dataset for the tutorial and all synthetic datasets tested in this manuscript can be found on Zenodo [here](https://zenodo.org/records/15285190). 

Please note that this repository is associated with a study currently under peer review and may undergo update in the near future

# **Citation**
Please cite this work: Yeh _et al_. [**Robust self-supervised machine learning for single cell embeddings and annotations.**](https://www.biorxiv.org/content/10.1101/2025.06.05.658097v1) _bioRxiv_ (2025)

# **Software Requirements**

* R (tested on version 4.3.2 -- "Eye Holes")
* R library dependencies: lpSolve (5.6.23), plyr (1.8.9), Seurat (5.1.0), SeuratObject (5.0.2), sp (2.1-4), dplyr (1.1.4)
  
# **License** 

BSD 3-Clause License provided ([here](https://github.com/Jerby-Lab/drgem/blob/main/LICENSE)).

# **Contact**

For any inquiries on this repository please feel free to post these under "issues" or reach out to Christine Yiwen Yeh ([cyyeh@stanford.edu](cyyeh@stanford.edu)) and Livnat Jerby ([ljerby@stanford.edu](ljerby@stanford.edu)).
