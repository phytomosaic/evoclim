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

Then answer:

* Which clades are narrow-niche “climate specialists” that may be extinction-prone?
* Which clades are “climate losers” (low 95th percentiles) expected to decline with warming/drying (and “winners” for vice versa)
* Which clades are faster-evolving, with faster rates of niche evolution, possibly better able to adapt to climate changes?
* Green-algal faster/slower climate evolution than cyanolichens?
* Which more sensitive: thermal, or moisture, or interactive predictors?


## Load data

```r
require(evoclim)
?veg
data(veg)
d   <- veg
xy  <- d$xy    # spatial
spe <- d$spe   # species
env <- d$env   # environment
tra <- d$tra   # traits
phy <- d$phy   # phylogeny
```
