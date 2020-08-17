# evoclim

Climatic niche evolution of lichens.

## Contributors

Matt Nelsen  
Rob Smith


## Motivation

Working repository for phylogenetic reconstruction of ancestral climate tolerances of lichenized fungi.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/evoclim')
```

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

## Visualize

```r
# spatial
plot(xy, pch=19, col='#00000050')

#species
plot_heatmap(spe, logbase=10)

# environment
plot_heatmap(sapply(env, function(x)100+scale(x)))

# traits
plot_heatmap(sapply(tra, function(x)100+scale(x)))

# phylogeny
plot(phy, cex=0.6, no.margin=TRUE)

```