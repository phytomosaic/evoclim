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
load(file='./data/d.rda', verbose=T)

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
save(tr, file='./data/tr.rda')

####    END    ####