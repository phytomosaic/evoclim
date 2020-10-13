######################################################################
#
#  Phylogenetic signal for all climate traits
#
#    Rob Smith, phytomosaic@gmail.com, 14 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################


rm(list=ls())
gc()
require(ecole)      # for plotting and convenience functions
require(viridis)    # for plotting color scales
require(phytools)   # for phylogenetic comparative tasks
load(file='./data/tr.rda', verbose=T)
print(object.size(tr), units='MB')
head(tr) # traits, where rows = species and cols = niche positions/breadth

### read in the phylogeny (Nelsen et al 2020 PNAS)
p <- read.tree('./data_raw/Lecanoromycetes_ML_Tree_Timescaled')
# plot(p, cex=0.2, type='fan')
ptax <- p$tip.label  # taxa names in phylogeny
otax <- rownames(tr) # taxa names in GBIF occurrences
length(ptax)         # full phylogeny has 3373 taxa
length(otax)         # occurrence data has 2577 taxa
all(otax %in% ptax)  # expect TRUE
all(rownames(tr) %in% p$tip.label)

### proceed with pruning tree and occurrences to exact matches
p    <- keep.tip(p, tip = otax) # keep only taxa in GBIF occurrences
ptax <- p$tip.label  # (reduced) taxa names in phylogeny
otax <- rownames(tr) # (reduced) taxa names in GBIF occurrences
length(ptax)         # 2577 exact matches
length(otax)         # 2577 exact matches
rm(ptax,otax)

### phylogenetic signal of climatic niches (climate niche conservatism)
#     best metric is lambda, according to doi:10.1111/j.2041-210X.2012.00196.x
lam <- apply(tr, 2, FUN=function(j) { # ! ! ! TIMEWARN ! ! ! ~5 min per trait
  names(j) <- rownames(tr)
  x <- phytools::phylosig(p, j, method='lambda', test=T, nsim=999)
  round(c(lambda=x$lambda, pval_l=x$P), 4)})
lam

### visualize ALL climate niche values across phylogeny
png('./fig/fig_01_phy_heatmap.png', wid=9.5, hei=16.5, units='in',
    bg='transparent', res=500)
phylo.heatmap(p, tr, 0.2, viridis::inferno(99), standardize=T)
dev.off()

### pick a SINGLE climate niche trait
x  <- setNames(tr[,'bio1_q50'], rownames(tr)) # bio1_q50 = 'thermal centroid'
cc <- ecole::colvec(x, begin=0.1, end=0.95, alpha=1)
(rng <- range(x))

### color tips by climate niche value
png('./fig/fig_02_phy_heattips.png', wid=9.5, hei=16.5, units='in',
    bg='transparent', res=500)
plot(p, cex=0.2, label.offset=5, no.margin=T, tip.col=cc) 
# tiplabels(pch = 21, bg = cc, cex = 0.8, adj = 4)
dev.off()

### continuous trait reconstruction ! ! ! TIMEWARN ! ! ! ~1-2 min per biovar
trm <- contMap(p, x, lims=c(-20,30), plot=FALSE) # _across branches_
trm$cols[1:length(trm$cols)] <- viridis::inferno(length(trm$cols))
# fa <- fastAnc(p, x) # _at each node_

### plot ancestral trait reconstruction
png('./fig/fig_03_phy_traitreconstruction.png', wid=6.5, hei=16.5, units='in',
    bg='transparent', res=1000)
plot(trm, fsize=0.1, lwd=1, outline=F, leg.txt='Thermal centroid (\u0b00C)')
# nodelabels(pch=19, col = colvec(fa), cex = 0.9)
# tiplabels(pch=19, col=colvec(tr,alpha=1), cex = 0.1)
dev.off()

### plot ancestral trait reconstruction (fan)
png('./fig/fig_03_phy_traitreconstruction_fan.png', wid=8.5, hei=8.5, units='in',
    bg='transparent', res=1000)
plot(trm, fsize=0.1, lwd=1, outline=F, leg.txt='Thermal centroid (\u0b00C)', type='fan')
# nodelabels(pch=19, col = colvec(fa), cex = 0.9)
# tiplabels(pch=19, col=colvec(tr,alpha=1), cex = 0.1)
dev.off()

# ### traitgram: trait (y) vs time (x)
# phenogram(p, tr, fsize=0.6, spread.costs=c(1,0.5))
# phenogram(p, tr, fsize=0.6, spread.costs=c(1,0), colors=trm$cols)

### TODO rates of trait evolution: 
?ratebytree
ratebytree(trees, x, type='continuous', ...)
# `trees`	=
#   an object of class "multiPhylo". If x consists of a list of different traits
#   to be compared, then trees could also be a simple set of duplicates of the
#   same tree, e.g., rep(tree,length(x)).

### human-readable names, for publication
nm <- c('Ann T','T diurnal rng','Isothermality','T seasnlty','Max T warm mo',
        'Min T cold mo','T ann rng','T wet qtr','T dry qtr','T warm qtr',
        'T cold qtr','Ann P','P wet mo','P dry mo','P seasnlty','P wet qtr',
        'P dry qtr','P warm qtr','P cold qtr')

# ### variable names, just for reference:
#   bio1  = Mean annual temperature
#   bio2  = Mean diurnal range (mean of max temp - min temp)
#   bio3  = Isothermality (bio2/bio7) (* 100)
#   bio4  = Temperature seasonality (standard deviation *100)
#   bio5  = Max temperature of warmest month
#   bio6  = Min temperature of coldest month
#   bio7  = Temperature annual range (bio5-bio6)
#   bio8  = Mean temperature of the wettest quarter
#   bio9  = Mean temperature of driest quarter
#   bio10 = Mean temperature of warmest quarter
#   bio11 = Mean temperature of coldest quarter
#   bio12 = Total (annual) precipitation
#   bio13 = Precipitation of wettest month
#   bio14 = Precipitation of driest month
#   bio15 = Precipitation seasonality (coefficient of variation)
#   bio16 = Precipitation of wettest quarter
#   bio17 = Precipitation of driest quarter
#   bio18 = Precipitation of warmest quarter