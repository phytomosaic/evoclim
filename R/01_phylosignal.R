######################################################################
#
#  Phylogenetic signal for all climate traits
#
#    Rob Smith, phytomosaic@gmail.com, 14 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################

rm(list=ls())
require(ecole)       # for plotting and convenience functions
require(viridis)     # for plotting color scales
require(phytools)    # for phylogenetic comparative tasks
require(phylosignal) # ?phylosignal::phyloSignal()
require(phylobase)   # to convert `phylo` --> `phylo4d` object

### load traits data
load(file='./data/tr.rda', verbose=T) # traits matrix
print(object.size(tr), units='MB')    # size = 1.7 MB

### load the phylogeny (Nelsen et al 2020 PNAS)
p <- read.tree('./data/Lecanoromycetes_ML_Tree_Timescaled')
# plot(p, cex=0.2, type='fan')
ptax <- p$tip.label  # taxa names in phylogeny
otax <- rownames(tr) # taxa names in GBIF occurrences
length(ptax)         # full phylogeny has 3373 taxa
length(otax)         # occurrence data has 2577 taxa
all(otax %in% ptax)                # expect TRUE
all(rownames(tr) %in% p$tip.label) # expect TRUE

### proceed with pruning tree and occurrences to exact matches
p    <- keep.tip(p, tip = otax) # keep only taxa in GBIF occurrences
ptax <- p$tip.label  # (reduced) taxa names in phylogeny
otax <- rownames(tr) # (reduced) taxa names in GBIF occurrences
length(ptax)         # 2577 exact matches
length(otax)         # 2577 exact matches
rm(ptax,otax)

### phylogenetic signal of climatic niches (Blombergs K)
pp <- phylobase::phylo4d(p, tip.data=tr) # temporarily format to phylo4d
k  <- phylosignal::phyloSignal(pp, methods='K', reps=999) # 436.5 sec = 7.3 min
save(k, file='./res/k.rda')
rm(pp)

### quick visualize ALL climate niche values across phylogeny
# png('./fig/fig_01_phy_heatmap.png', wid=9.5, hei=16.5, units='in',
#     bg='transparent', res=500)
# phylo.heatmap(p, tr, 0.2, viridis::inferno(99), standardize=T)
# dev.off()

### human-readable names, for figures
nm <- c('Ann T','T diurnal rng','Isothermality','T seasnlty','Max T warm mo',
        'Min T cold mo','T ann rng','T wet qtr','T dry qtr','T warm qtr',
        'T cold qtr','Ann P','P wet mo','P dry mo','P seasnlty','P wet qtr',
        'P dry qtr','P warm qtr','P cold qtr')
nm <- paste0(rep(nm, ea=4), gsub('.*_', ' ', dimnames(tr)[[2]]))

### ancestral state reconstructions (save outputs and plots)
lapply(1:NCOL(tr), function(j) { # ! ! ! TIMEWARN ! ! ! ~1-2 min per trait
  cat(paste0(j, ' of ', NCOL(tr), '... '))
  jnm <- j
  x   <- tr[,j]
  lim <- c(floor(min(x, na.rm=T)),ceiling(max(x, na.rm=T)))
  # continuous trait reconstruction across branches 
  trm <- phytools::contMap(p, x, lims=lim, plot=F)
  trm$cols[1:length(trm$cols)] <- viridis::inferno(length(trm$cols))
  # save output
  fnm <- paste0('./res/asr/asr_',sprintf("%02d", jnm),'__',
                dimnames(tr)[[2]][jnm],'.rda')
  save(trm, file=fnm)
  # plot output (fan)
  fnm <- paste0('./fig/asr/fig_',sprintf("%02d", jnm),'_asr__',
                dimnames(tr)[[2]][jnm],'.png')
  png(fnm, wid=8.5, hei=8.5, units='in', bg='transparent', res=500)
  plot(trm, lwd=1, outline=F, leg.txt=nm[jnm], type='fan', fsize=c(0.1,1))
  dev.off()
})

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

####    END    ####