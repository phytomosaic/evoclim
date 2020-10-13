# evoclim

Climatic niche evolution of lichens.

## Contributors

Matt Nelsen  
Rob Smith


## Motivation

Working repository for phylogenetic reconstruction of climate tolerances of lichenized fungi.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/evoclim')
```


## Approach:
1.  Construct phylogeny (https://doi.org/10.1073/pnas.2001913117)
1.  Obtain geo distribution info (GBIF)
1.  Obtain climate info (MERRAClim)
1.  For each taxon, and each climate predictor, estimate:
    * Niche minima (5th percentile)
    * Niche centroid (50th percentile)
    * Niche maxima (95th percentile)
    * Niche breadth (IQR = 75 minus 25 percentile)
1.  Ancestral state reconstruction of climate niches mapped onto the tree 
1.  Estimate rates of climatic niche evolution  


## Questions

1. Which traits more phylogenetically conserved: thermal, or moisture, or interactive predictors?
    * how: compare phylogenetic signal (K) among each predictor

1. Which clades are narrow-niche “climate specialists” that may be extinction-prone?
    * how: map "niche breadth (mean abs deviation)" onto phylogeny, identify clades with narrowest breadth, test for phylo signal (K)

1. Which clades are “climate losers” (low 95th percentiles) expected to decline with warming/drying (and “winners” for vice versa)
    * how: map "niche maxima (95th percentiles)" onto phylogeny, identify clades where low maxima are approaching modern conditions and therefore "vulnerable", test for phylo signal (K)

1. Which clades are faster-evolving, with faster rates of climatic niche evolution, possibly better able to adapt to imminent climate changes?
    * how: BAMM if single trait; l1ou for multivariate.  w/l1ou, could also examine whether shifts into novel adaptive regimes concentrated early in the phylogeny (early-burst).

1. Green-algal faster/slower climate evolution than cyanolichens?
    * how: OUwie. or BAMM (see Folk et al. 2018 PNAS)

1. To what degree is rate of diversification correlated with rate of climatic niche evolution?
    * how: BAMM (see Folk et al. 2018 PNAS)


## Questions/comments

phytomosaic@gmail.com
