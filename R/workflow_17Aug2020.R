######################################################################
#
#  Workflow of analysis
#
#    Rob Smith, phytomosaic@gmail.com, 17 Aug 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
# devtools::install_github('phytomosaic/ecole')
require(ecole)      # for plotting and convenience functions
require(phytools)   # for phylogenetic comparative tasks
require(sp)         # for spatial tasks
require(raster)     # for spatial tasks
require(FNN)        # to assign env values from nearest neighbors, when NA
require(data.table) # for memory-safe handling of large files
require(CoordinateCleaner) # for obvious
### read GBIF individual species files, binding them together
#     (rows = occurrences, cols = species+lon+lat)
`f` <- function(x, keepcols=NULL, ...) {
  cat(round(file.info(x)$size/1024^2,1),'MB; ')
  # cat('Time elapsed:',
  #     system.time({
        out <- data.table::fread(x, header = T, dec = '.', fill = T, 
                                 data.table = T, select = keepcols,
                                 integer64='character', ...)
      # })[[3]], '\n')
  return(out)
}
pth <- './data_raw/gbif/data_raw_gbif/'
j   <- c('key','species','decimalLatitude','decimalLongitude','issues','year')
fnm <- list.files(pth, pattern='.csv', full.names=T)
xy  <- rbindlist(lapply(fnm, f, keepcols=j), fill = T)
setnames(xy, names(xy), tolower(names(xy))) # cleanup column names
setnames(xy, c('decimallatitude','decimallongitude'), c('lat','lon')) # simplify
dim(xy)



anyDuplicated(xy)

### filter occurrences (duplicates, ocean-dwellers, etc)
# TODO
sum(duplicated(xy))        # 7638 duplicates
xy <- xy[!duplicated(xy),] # remove duplicates
xy$issues %in% c('')




### bin occurrences to grid (better than thinning, Smith et al 2020 J Biogeogr)
# TODO

# ### obtain WorldClim climate data
# r <- raster::getData('worldclim',var='bio',res=10)


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

### read in 1980-2010 30-y climate normals (from MERRAclim)
pth <- './data_raw/merraclim/'
fnm <- list.files(pth, pattern ='\\.tif$', full.names=T)
r   <- stack(paste0(pth, fnm))
trg_prj   <- projection(r)                     # raster projection



### convert coordinates to spatial object
`make_xy` <- function(xy, crs, ...) {          # reproject points
  xy <- xy[!(is.na(xy$lon) | is.na(xy$lat)),]  # rm NAs
  coordinates(xy) <- ~lon+lat
  proj4string(xy) <- CRS('+init=epsg:4269')    # NAD83 dec deg (as for FIA)
  return(spTransform(xy, crs))
}
xy <- make_xy(xy, crs=trg_prj) # reproject Albers equal-area proj

### extract env values at query coordinates
# vals <- extract(r,xy) ### SLOW WAY! 
### FAST WAY! extract values at query coordinates, in CHUNKS
`extract_raster` <- function(r, xy, breaks=10, ...){
  if(class(r)!='RasterLayer') {
    stop('r must be `RasterLayer` object')
  }
  if(!inherits(xy,'SpatialPoints')) {
    stop('xy must be `SpatialPoints` object')
  }
  # divide extent into equal area chunks for speed
  `make_chunk` <- function(e, breaks, ...){
    m  <- matrix(NA, nrow=breaks, ncol=4,
                 dimnames=list(NULL,
                               c('xmin','xmax','ymin','ymax')))
    rng <- range(c(e[1],e[2]))  # x dimension
    b   <- seq(rng[1], rng[2], length=breaks+1)
    for(i in 1:(breaks)){
      m[i,1:2] <- c(b[i], b[i+1])
    }
    rng <- range(c(e[3],e[4]))  # y dimension
    b   <- seq(rng[1], rng[2], length=breaks+1)
    for(i in 1:(breaks)){
      m[i,3:4] <- c(b[i], b[i+1])
    }
    m <- as.matrix(merge(m[,1:2],m[,3:4]))
    m
  }
  ck     <- make_chunk(e=extent(r), breaks=breaks)
  nchunk <- dim(ck)[[1]]
  npts   <- dim(xy@coords)[[1]]
  p      <- matrix(NA, nrow=npts, ncol=nchunk)
  pb     <- pbCreate(nchunk, progress='text', style=3, ...)
  for(i in 1:nchunk) {
    pbStep(pb, i)
    p[,i] <- raster::extract(crop(r, extent(ck[i,])), xy,
                             progress='text')
  }
  # collect chunks into one vector 'tmp'
  nonNA <- !is.na(p)
  tmp   <- vector('numeric', npts)
  for(i in 1:npts){
    if( all(!nonNA[i,]) ){
      tmp[i] <- NA   # some points were not queryable
    } else {
      tmp[i] <- p[i, min(which(nonNA[i,]))]
    }
  }
  # assemble matrix
  p <- matrix(c(1:npts,tmp),nrow=npts,ncol=2,dimnames=list(NULL,c('pt','val')))
  pbClose(pb)
  return(p)
}
vals <- extract_raster(r, xy, breaks=10)      # env values

### collect values and coordinates in new object
d <- cbind.data.frame(xy, vals)

### assign nearest-neighbor value if is NA (mostly at water boundaries)
i <- is.na(d$MAT)
x <- FNN::get.knn(d[,.(lon,lat)], k=50)$nn.index  # index 50 neighbors
x <- matrix(d$MAT[x], nrow=NROW(x), ncol=NCOL(x))[i,]
j <- apply(x, 1, function(x) which.min(is.na(x))) # col index nearest
d$MAT[i] <- x[cbind(1:NROW(x), j)]    # assign nearest non-NA
rm(i,j,x) # cleanup

### for each species, get niche positions and breadth (as 'traits')
tr <- sapply(sort(unique(d$spe)),
             function(i) {
               x <- d[d$spe==i,'MAP'] # TODO: more climate variables
               x <- na.omit(x)
               out <- c(quantile(x, probs=c(0.05,0.50,0.95)), stats::IQR(x))
               names(out) <- c('q05','q50','q95','iqr')
             })
cc <- colvec(tr)


### read in the phylogeny (Nelsen et al 2020 PNAS)
p <- read.tree('./data_raw/Lecanoromycetes_ML_Tree_Timescaled')
# plot(p, cex=0.2)
# plot(p, cex=0.2, type='fan')
# pspe <- p$tip.label               # species names
# pgen <- sub('\\_.*', '', pspe)    # genera names

### prune to keep only GBIF-valid species
has_gbif <- unique(d$spe)
p <- keep.tip(p, tip = na.omit(match(has_gbif, pspe)))

### prune to keep only taxa with sufficient support (>100 occurrences?)
has_support <- NA
p <- keep.tip(p, tip = na.omit(match(has_support, pspe)))

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

####    END    ####