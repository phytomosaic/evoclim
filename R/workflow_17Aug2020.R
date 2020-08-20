######################################################################
#
#  Workflow of analysis
#
#    Rob Smith, phytomosaic@gmail.com, 20 Aug 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
# devtools::install_github('phytomosaic/ecole')
require(ecole)      # for plotting and convenience functions
require(viridis)    # for plotting color scales
require(phytools)   # for phylogenetic comparative tasks
require(sp)         # for spatial tasks
require(raster)     # for spatial tasks
require(data.table) # for memory-safe handling of large files
require(CoordinateCleaner) # for obvious

### read GBIF individual species files, binding them together
`f` <- function(x, keepcols=NULL, ...) {
  cat(round(file.info(x)$size/1024^2,1),', ') # file size, in MB
  return(data.table::fread(x, header = T, dec = '.', fill = T, data.table = T, 
                           select = keepcols, integer64='character', ...))
}
pth <- './data_raw/gbif/data_raw_gbif/'
j   <- c('decimalLatitude','decimalLongitude') # 'species'
fnm <- list.files(pth, pattern='.csv', full.names=T)
nm  <- gsub('_gbif_bioclim_records.csv', '', fnm)
nm  <- gsub('./data_raw/gbif/data_raw_gbif/', '', nm)
fnm <- setNames(fnm, nm)  ;  rm(nm)
xy  <- rbindlist(lapply(fnm, f, keepcols=j), 
                 fill=T, idcol='species') # ! ! ! TIMEWARN ! ! ! ~3-4 min
setnames(xy, names(xy), tolower(names(xy))) # cleanup column names
setnames(xy, c('decimallatitude','decimallongitude'), c('lat','lon')) # simplify
dim(xy) # 878681 observations
head(xy)

### filter occurrences (duplicates, artificial coords, etc)
sum(duplicated(xy))         # 245479 duplicates flagged for removal
xy  <- xy[!duplicated(xy),] # remove duplicates
flg <- CoordinateCleaner::clean_coordinates(
  x = xy, lon='lon', lat='lat',
  tests = c('capitals','centroids','equal','gbif','institutions','zeros'),
  capitals_rad=250, centroids_rad=250, centroids_detail='country', inst_rad=100)
sum(flg$.summary == 0)      # 1621 were flagged for removal
dim(flg)[[1]] - sapply(flg[,grep('\\.', names(flg))], function(i) sum(i)) 
xy <- xy[flg$.summary,]     # exclude 1621 records flagged by ANY test
dim(xy)                     # 633202 observations
rm(f,flg,fnm,j,pth)         # cleanup environment

### remove species with fewer than 5 occurrences
frq <- table(xy$species)    # species frequencies
set_par(2) ; hist(frq, breaks=55) ; hist(log10(1+frq), breaks=55)
length(unique(xy$species))        # 2906 taxa total
sum(frq < 5)                      # 329 taxa are too infrequent to analyze
length(j <- names(frq[frq >= 5])) # 2577 taxa to keep
xy <- xy[xy$species %in% j]
dim(xy)                     # 632474 observations

# ### obtain MERRAclim climate data
# yr <- c('80s','90s','00s')
# bc <- 1:19
# pth <- './data_raw/merraclim/'
# for(x in 1:length(bc)) {
#   s80s <- raster(paste0(pth,'2_5m_mean_',yr[1],'_bio',bc[x],'.tif'))
#   s90s <- raster(paste0(pth,'2_5m_mean_',yr[2],'_bio',bc[x],'.tif'))
#   s00s <- raster(paste0(pth,'2_5m_mean_',yr[3],'_bio',bc[x],'.tif'))
#   stk <- stack(s80s,s90s,s00s)
#   stkmean <- mean(stk)
#   writeRaster(stkmean, paste0(pth,'2_5_mean_mean80sto00s_bio',bc[x],'.tiff'),
#               type='GTiff', overwrite=T)
# }

### read in 1980-2010 30-y climate normals (from MERRAclim), ~5-km resolution
pth     <- './data_raw/merraclim/'
fnm     <- list.files(pth, pattern ='\\.tif$', full.names=T)
r       <- stack(fnm)
trg_prj <- projection(r)                     # raster projection
# png('./fig/fig_00_climatevars.png', wid=8.5, hei=8.5, units='in',
#     bg='transparent', res=500)
# set_par(1)
# plot(r, col=viridis::inferno(99, begin=0.1, end=0.95), axes=F, box=F)
# data(wrld_simpl,package='maptools'); plot(wrld_simpl,add=T); rm(wrld_simpl)
# dev.off()

### convert coordinates to spatial object
`make_xy` <- function(xy, crs, ...) {          # reproject points
  coordinates(xy) <- ~lon+lat
  proj4string(xy) <- sp::CRS('+init=epsg:4326')    # WGS84 dec deg
  return(sp::spTransform(xy, crs))
}
xy <- make_xy(xy, crs=trg_prj)

### extract env values at query coordinates
vals <- raster::extract(r,xy)      # ! ! ! TIMEWARN ! ! ! ~10 min
d    <- cbind.data.frame(xy, vals) # values and coordinates in new object
rm(xy, vals, r, make_xy, pth, fnm, trg_prj) # cleanup
d$optional <- NULL
sum(is.na(d$X2_5_mean_mean80sto00s_bio1))      # 43 NA values
d <- d[!is.na(d$X2_5_mean_mean80sto00s_bio1),] # remove 43 NA values
names(d) <- gsub('X2_5_mean_mean80sto00s_','', names(d)) # clean names
d[,paste0('bio',c(1,2,4:11))] <- 
  d[,paste0('bio',c(1,2,4:11))] * 0.10 # rescale temperature units to degrees C
dim(d)  # 632431 observations, beautiful clean complete data
save(d, file='./data/d.rda')

### clean restart here ???????????????????????????????
rm(list=ls())
gc()
require(ecole)      # for plotting and convenience functions
require(viridis)    # for plotting color scales
require(phytools)   # for phylogenetic comparative tasks
require(data.table) # for memory-safe handling of large files
load(file='./data/d.rda', verbose=T)
print(object.size(d), units='MB')
head(d)

# TODO memory overrun!
# ### bin occurrences to grid (better than thinning, Smith et al 2020 J Biogeogr)
# ###   quantiles of grid-cell means per species = niche traits
# `bingrid_multi` <- function (x, field='', nr=100, nc=150, 
#                              xmn=-180, xmx=180, ymn=-90, ymx=90, ...) {
#   t_start <- Sys.time()
#   cat(paste0('Began at: ', t_start, '\n'))
#   `bingrid_q` <- function (...) {
#     hascol <- all(c('lat','lon','species',field) %in% names(x))
#     if (!hascol) stop('need columns `lat`,`lon`,`species` and `field`')
#     # rasterize: grid-cell means BY SPECIES
#     b <- by(data = x[,!colnames(x) %in% 'species'],
#             INDICES = x[,'species'],
#             FUN = function(a) {
#               raster::rasterize(
#                 a[,c('lon','lat')],
#                 raster::raster(nrows=nr,ncols=nc,xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx),
#                 field = as.matrix(a[,field]),
#                 fun   = mean,
#                 na.rm = TRUE)}, simplify = FALSE)
#     b <- sapply(b, raster::getValues) # grid-cell-means (rows) by species (cols)
#     # for each species, get 'traits' = niche positions and breadth
#     tr <- apply(b, 2,
#                 function(i) {
#                   i   <- na.omit(i)
#                   out <- c(quantile(i, probs=c(0.05,0.50,0.95)), stats::mad(i))
#                   names(out) <- c('q05','q50','q95','mad')
#                   out
#                 })
#     return(tr)
#   }
#   # apply to each field
#   lst <- lapply(field, FUN=function(i) bingrid_q(x, field=i))
#   res <- do.call(rbind, lst)
#   rownames(res) <- paste0(rep(field,each=4),'_',rownames(res))
#   res <- t(res)
#   # timing
#   t_end  <- Sys.time()
#   t_diff <- round(difftime(t_end, t_start, units="mins"),2)
#   cat(paste0('Completed at: ', t_end),
#       '\n   after time elapsed of:', t_diff, 'minutes')
#   return(res) # rows = species, cols = climate traits
# }
# b <- bingrid_multi(d, field=paste0('bio',1:19))
# # save(b, file='./data/b.rda')
# load(file='./data/b.rda', verbose=T)
# print(object.size(b), units='MB')

### for each species and each biovar, get 'traits' = niche positions and breadth
tr <- lapply(d[,paste0('bio',1:19)], function(var) {
  do.call('rbind', 
          tapply(var, d$species,
                 function(i) {
                   i   <- na.omit(i)
                   out <- c(stats::quantile(i, probs=c(0.05,0.50,0.95), na.rm=T),
                            stats::mad(i, na.rm=T))
                   names(out) <- c('q05','q50','q95','mad')
                   out
                 }))})
tr <- do.call('cbind', tr)
colnames(tr) <- paste0(rep(paste0('bio',1:19),each=4),'_',colnames(tr))
head(tr) # traits, where rows = species and cols = niche positions/breadth
rm(d)

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
?phytools::phylosig
lam <- apply(tr, 2, FUN=function(j) {
  names(j) <- rownames(tr)
  x <- phytools::phylosig(p, j, method='lambda', test=T, nsim=999)
  round(c(lambda=x$lambda, pval_l=x$P), 4)})
lam

### visualize ALL climate niche values across phylogeny
png('./fig/fig_01_phy_heatmap.png', wid=9.5, hei=16.5, units='in',
    bg='transparent', res=500)
phylo.heatmap(p, tr, 0.2, viridis::inferno(99), standardize=T)
dev.off()

### color tips by climate tolerance value
png('./fig/fig_02_phy_heattips.png', wid=9.5, hei=16.5, units='in',
    bg='transparent', res=500)
x  <- setNames(tr[,'bio12_q05'], rownames(tr))        # lowest extreme of MAP
cc <- ecole::colvec(x, begin=0.1, end=0.95, alpha=1) # lowest extreme of MAP
plot(p, cex=0.2, label.offset=5, no.margin=T, tip.col=cc) 
# tiplabels(pch = 21, bg = cc, cex = 0.8, adj = 4)
dev.off()

### continuous trait reconstruction ! ! ! TIMEWARN ! ! ! ~1-2 min per biovar
x   <- setNames(tr[,'bio1_q50'], rownames(tr)) # central value of MAT (thermal centroid)
(rng <- range(x))
cc  <- ecole::colvec(x, begin=0.1, end=0.95, alpha=1)
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

### rates of trait evolution: 
?ratebytree
ratebytree(trees, x, type='continuous', ...)
# `trees`	=
#   an object of class "multiPhylo". If x consists of a list of different traits
#   to be compared, then trees could also be a simple set of duplicates of the
#   same tree, e.g., rep(tree,length(x)).

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

### human-readable names, for publication
nm <- c('Ann T','T diurnal rng','Isothermality','T seasnlty','Max T warm mo',
        'Min T cold mo','T ann rng','T wet qtr','T dry qtr','T warm qtr',
        'T cold qtr','Ann P','P wet mo','P dry mo','P seasnlty','P wet qtr',
        'P dry qtr','P warm qtr','P cold qtr')


####    END    ####