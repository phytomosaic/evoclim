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
j   <- c('species','decimalLatitude','decimalLongitude')
fnm <- list.files(pth, pattern='.csv', full.names=T)
xy  <- rbindlist(lapply(fnm, f, keepcols=j), fill=T) # ! ! TIMEWARN ! ! 3-4 min
setnames(xy, names(xy), tolower(names(xy))) # cleanup column names
setnames(xy, c('decimallatitude','decimallongitude'), c('lat','lon')) # simplify
xy$species <- ecole::clean_text(xy$species) # cleanup taxon names
dim(xy) # 878681 observations

### filter occurrences (duplicates, artificial coords, etc)
sum(duplicated(xy))         # 7638 duplicates flagged for removal
xy  <- xy[!duplicated(xy),] # remove duplicates
flg <- CoordinateCleaner::clean_coordinates(
  x = xy, lon='lon', lat='lat',
  tests = c('capitals','centroids','equal','gbif','institutions','zeros'),
  capitals_rad=250, centroids_rad=250, centroids_detail='country', inst_rad=100)
sum(flg$.summary == 0)      # 3663 were flagged for removal
dim(flg)[[1]] - sapply(flg[,grep('\\.', names(flg))], function(i) sum(i)) 
xy <- xy[flg$.summary,]     # exclude 3663 records flagged by ANY test
dim(xy)                     # 867380 observations
rm(f,flg,fnm,j,pth)         # cleanup environment

### remove species with fewer than 10 occurrences
frq <- table(xy$species)    # species frequencies
set_par(2) ; hist(frq, breaks=55) ; hist(log10(1+frq), breaks=55)
length(unique(xy$species))         # 2854 taxa total
sum(frq < 10)                      # 490 taxa are too infrequent to analyze
length(j <- names(frq[frq >= 10])) # 2364 taxa to keep
xy <- xy[xy$species %in% j]
dim(xy)                     # 865255 observations

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
png('./fig/fig_00_climatevars.png', wid=8.5, hei=8.5, units='in',
    bg='transparent', res=500)
set_par(1)
plot(r, col=viridis::inferno(99, begin=0.1, end=0.95), axes=F, box=F)
data(wrld_simpl,package='maptools'); plot(wrld_simpl,add=T); rm(wrld_simpl)
dev.off()

### convert coordinates to spatial object
`make_xy` <- function(xy, crs, ...) {          # reproject points
  coordinates(xy) <- ~lon+lat
  proj4string(xy) <- sp::CRS('+init=epsg:4326')    # WGS84 dec deg
  return(sp::spTransform(xy, crs))
}
xy <- make_xy(xy, crs=trg_prj)
rm(make_xy,pth,fnm,trg_prj)
save(xy, file='./data/xy.rda')  # save intermediate result to save time
load(file='./data/xy.rda',verbose=T)

### extract env values at query coordinates
vals <- raster::extract(r,xy)      # ! ! ! TIMEWARN ! ! ! ~10 min
d    <- cbind.data.frame(xy, vals) # values and coordinates in new object
rm(xy, vals, r)
d$optional <- NULL
d <- d[!is.na(d$X2_5_mean_mean80sto00s_bio1),] # remove 50 NA values
names(d) <- gsub('X2_5_mean_mean80sto00s_','', names(d)) # clean names
d[,paste0('bio',c(1,2,4:11))] <- 
  d[,paste0('bio',c(1,2,4:11))] * 0.10 # rescale temperature units to degrees C
dim(d)  # 865205 observations
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
tr <- lapply(d[,paste0('bio',1:19)], function(var) { # for each biovar
  do.call('rbind', 
          tapply(var, d$species, # for each species
                 function(i) {
                   i   <- na.omit(i)
                   out <- c(stats::quantile(i, probs=c(0.05,0.50,0.95), na.rm=T),
                            stats::mad(i, na.rm=T))
                   names(out) <- c('q05','q50','q95','mad')
                   out
                 }))})
str(tr,1) # 19 matrices, where rows = species and cols = niche positions/breadth
head(tr[[1]])

### read in the phylogeny (Nelsen et al 2020 PNAS)
p <- read.tree('./data_raw/Lecanoromycetes_ML_Tree_Timescaled')
# plot(p, cex=0.2, type='fan')
ptax <- p$tip.label               # taxa names
# pgen <- sub('\\_.*', '', ptax)    # genera names
length(ptax) # full phylogeny has 3373 taxa

### prune to keep only GBIF-valid species
p <- keep.tip(p, tip = na.omit(match(unique(d$spe), ptax)))

### prune to keep only taxa with sufficient support (>100 occurrences?)
has_support <- NA
p <- keep.tip(p, tip = na.omit(match(has_support, ptax)))

### phylogenetic signal of climatic niches (climate niche conservatism)
?phytools::phylosig
lam <- sapply(n, FUN=function(j) {
  names(j) <- rownames(tr)
  l1 <- phytools::phylosig(p, j, method='K', test=T, nsim=999)
  l2 <- phytools::phylosig(p, j, method='lambda', test=T, nsim=999)
  round(c(k=l1$K, pval_k=l1$P, lambda=l2$lambda, pval_l=l2$P), 4)})
lam

### visualize ALL climate niche values across phylogeny
phylo.heatmap(p, tr, 0.4, colors=viridis::inferno(99), standardize=T)

### estimating diversification rates: plot number of lineages vs time
###    http://www.phytools.org/Cordoba2017/ex/11/LTT-and-gamma.html
(obj <- ltt(p,plot=FALSE))
plot(obj, log='y', log.lineages=FALSE)
# line of expecation under a pure-birth prediction:
h <- max(nodeHeights(p))
x <- seq(0,h,by=h/100)
b <- (log(Ntip(p))-log(2))/h
lines(x,2*exp(b*x),col=2,lty=2,lwd=2)
# simulate 100 similar trees:
trees <- pbtree(b=b,n=Ntip(p),t=h,nsim=100,method='direct',quiet=T)
obj   <- ltt(trees, plot=FALSE)
plot(obj, col='grey', log.lineages=T,main='Observed vs simulated LTTs')
lines(c(0,h),log(c(2,Ntip(p))),lty=2,lwd=2,col=4)
ltt(p,add=TRUE,col=2,lwd=2) # overlay original tree
# alternatively:
ltt95(trees, log=T, plot=F)
ltt(p,add=TRUE,col=2,lwd=2,log.lineages=F) # overlay original tree

### color tips by climate tolerance value
# tiff('phy_FIA_only.tif', wid=4.5, hei=5, unit='in', compr='lzw+p', res=600)
layout(matrix(1:2, 1,2))
plot(p, cex=0.3, label.offset=5, no.margin=T)
plot(p, cex=0.3, label.offset=5, no.margin=T, tip.col=cc)
# tiplabels(pch = 21, bg = cc, cex = 0.8, adj = 4)
# dev.off()

### continuous trait reconstruction _across branches_
trm <- contMap(p, x=tr, lims=c(16.5,33.5), plot=FALSE)
trm$cols[1:length(trm$cols)] <- viridis::inferno(length(trm$cols))
### continuous trait reconstruction _at each node_
fa  <- fastAnc(p, tr)

# across branches:
plot(trm, fsize=0.4, lwd=2, outline=F)
# at each node:
plot(p, cex=0.3, label.offset=5, no.margin=T, tip.col=cc,edge.width=2)
nodelabels(pch=19, col = colvec(fa), cex = 1)
tiplabels(pch=19, col = colvec(tr), cex = 1)

### plot ancestral trait reconstruction
png('./fig/phy_traitreconstruction.tif', 
    wid=3.5, hei=7.5, unit='in', bg='transparent', res=700)
plot(trm, fsize=0.4, lwd=3, outline=F,
     leg.txt='Upper thermal tolerance (\u00b0C)')
nodelabels(pch=19, col = colvec(fa), cex = 0.9)
tiplabels(pch=19, col = colvec(tr,alpha=1), cex = 0.8)
dev.off()

### traitgram: trait (y) vs time (x)
phenogram(p, tr, fsize=0.6, spread.costs=c(1,0.5))
phenogram(p, tr, fsize=0.6, spread.costs=c(1,0), colors=trm$cols)

### phylomorphospace
names(n)
par(las=1, bty='L', pty='s')
X <- n[,c('p95_mwmt','p95_cmd')]
phylomorphospace(p, X, label='horizontal', fsize=0.7)
# colored by trait
names(cc) <- p$tip.label
ccc <- c(cc[p$tip.label], rep('black', p$Nnode))
names(ccc) <- 1:(length(p$tip) + p$Nnode)
phylomorphospace(p, X, label='horizontal', fsize=0.7,
                 control=list(col.node=ccc))

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