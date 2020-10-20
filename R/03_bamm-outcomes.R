######################################################################
#
#  BAMM outcomes: analyze tests for rateshifts in climatic niche evolution
#
#    Rob Smith, phytomosaic@gmail.com, 19 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################



### preamble
rm(list=ls())
require(ecole)             # for plotting and convenience functions
require(BAMMtools)         # for BAMM
require(coda)              # for MCMC diagnostics
require(phytools)          # for ancestral state reconstruction
system( 'bamm --version' ) # BAMM 2.5.0 (2015-11-01)


### setup dirs and filenames
setwd('/home/rob/prj/evoclim/bamm')   # change into working dir (R)
fnm_phy     <- 'lichen_phy.tre'       # phylogeny
fnm_pos_tra <- 'lichen_tra_pos.txt'   # traits (niche position)
fnm_brd_tra <- 'lichen_tra_brd.txt'   # traits (niche breadth)
fnm_pos_event <- 'pos_event_data.txt' # event data (niche position)
fnm_brd_event <- 'brd_event_data.txt' # event data (niche breadth) 
fnm_div_event <- 'div_event_data.txt' # event data (niche breadth) 
fnm_pos_mcmc  <- 'pos_mcmc_out.txt'   # mcmc data (niche position)
fnm_brd_mcmc  <- 'brd_mcmc_out.txt'   # mcmc data (niche breadth)
fnm_div_mcmc  <- 'div_mcmc_out.txt'   # mcmc data (niche breadth)


### load tree and traits
phy       <- read.tree(fnm_phy)      # load tree
png('../fig/fig_04_pca_histograms.png', wid=6.5, hei=3.25, units='in',
    bg='transparent', res=700)
set_par_mercury(2)
hist(read.table(fnm_pos_tra)[,2], breaks=77, main='', xlab='PC1 niche positions') 
text(c(-12, 0, 12), c(100,100,100),c('Polar','Temperate','Tropical'), cex=0.7)
hist(read.table(fnm_brd_tra)[,2], breaks=77, main='', xlab='PC1 niche breadths') 
text(c(-14, 0), c(150,150),c('Generalist','Specialist'), cex=0.7)
dev.off()

# ### do quick ancestral state reconstructions
# lapply(1:NCOL(tr), function(j) { # ! ! ! TIMEWARN ! ! ! ~1-2 min per trait
#   cat(paste0(j, ' of ', NCOL(tr), '... '))
#   jnm <- j
#   x   <- tr[,j]
#   lim <- c(floor(min(x, na.rm=T)),ceiling(max(x, na.rm=T)))
#   # continuous trait reconstruction across branches 
#   trm <- phytools::contMap(p, x, lims=lim, plot=F)
#   trm$cols[1:length(trm$cols)] <- viridis::inferno(length(trm$cols))
#   # save output
#   fnm <- paste0('./res/asr/asr_',sprintf("%02d", jnm),'__',
#                 dimnames(tr)[[2]][jnm],'.rda')
#   save(trm, file=fnm)
#   # plot output (fan)
#   fnm <- paste0('./fig/asr/fig_',sprintf("%02d", jnm),'_asr__',
#                 dimnames(tr)[[2]][jnm],'.png')
#   png(fnm, wid=8.5, hei=8.5, units='in', bg='transparent', res=500)
#   plot(trm, lwd=1, outline=F, leg.txt=nm[jnm], type='fan', fsize=c(0.1,1))
#   dev.off()
# })


### load BAMM event data
ed_pos <- getEventData(phy = phy,
                   eventdata = fnm_pos_event,
                   burnin    = 0.10,  # 10% per PNAS
                   nsamples  = 2000,  # 2000 per PNAS
                   verbose   = TRUE,
                   type      = 'trait')
ed_brd <- getEventData(phy = phy,
                   eventdata = fnm_brd_event,
                   burnin    = 0.10,  # 10% per PNAS
                   nsamples  = 2000,  # 2000 per PNAS
                   verbose   = TRUE,
                   type      = 'trait')
ed_div <- getEventData(phy = phy,
                   eventdata = fnm_div_event,
                   burnin    = 0.10,  # 10% per PNAS
                   nsamples  = 2000,  # 2000 per PNAS
                   verbose   = TRUE,
                   type      = 'diversification')


### convergence diagnostics
mc_pos <- read.csv(fnm_pos_mcmc, header=T)
tail(mc_pos)
burnstart <- floor(0.1 * nrow(mc_pos)) # discard first 10% of samples as burnin
postburn  <- mc_pos[burnstart:nrow(mc_pos), ]
coda::effectiveSize(postburn$N_shifts) # effective sample sizes of number of shift events in each sample
coda::effectiveSize(postburn$logLik) # effective sample sizes of the log-likelihood 
set_par_mercury(2)
plot(mc_pos$logLik ~ mc_pos$generation) # sketch whether run converged 
abline(v=mc_pos$generation[burnstart], col=2)
plot(mc_pos$logPrior ~ mc_pos$generation)
abline(v=mc_pos$generation[burnstart], col=2)

# compare prior and posterior
set_par_mercury(1)
plotPrior(mc_pos, expectedNumberOfShifts=34) 
plotPrior(mc_pos, expectedNumberOfShifts=34) 
plotPrior(mc_pos, expectedNumberOfShifts=34) 


### identify 'best' configurations (3 alternatives):
# ### alternative 1 --- maximize Bayes factors for all pairwise models
# shift_probs <- summary(ed_pos)     # shift probabilities
# pr_null     <- 1/NROW(mc_pos) # probability of null no-shift model
# bfmat <- computeBayesFactors(mc_pos, expectedNumberOfShifts=34, burnin=0.1)
# plot_heatmap(bfmat)
# persp(bfmat, theta=135)
# persp(bfmat[,1:11], theta=135)
# (sbf <- stepBF(bfmat, step.size=1)) # optimal number of shifts
# set_par_mercury(3) # plot shift probabilities and Bayes factors
# plot(shift_probs[,1], shift_probs[,2], type='b', ylab='Shift probabilities')
# plot(shift_probs[,1], shift_probs[,2] / pr_null, type='b',
#      ylab='Min Bayes factor, relative to null no-shift model')
# plot(shift_probs[,1], bfmat[1,], type='b',
#      ylab='Bayes factors, relative to 9-shift model')
# abline(v=30)
# # abline(v=as.numeric(rownames(bfmat)[sbf$bestModel]))
# ### alternative 2 --- configuration that maximizes marginal probability of 
# ###    rate shifts along individual branches
# msc   <- maximumShiftCredibility(ed_pos, maximize='product')
# msc_i <- subsetEventData(ed_pos, index = msc$sampleindex)
# msc_i$eventData # 31 distinct shifts in this configuration
# plot(msc_i, 
#      method = 'polar',
#      pal=viridis::inferno(99),
#      logcolor=TRUE,
#      # breaksmethod='quantile',  # breaksmethod='jenks'
#      legend=T)
# addBAMMshifts(msc_i, method = 'polar', cex=1, bg='gold')
### alternative 3 --- configuration that maximizes a posteriori probability
best <- getBestShiftConfiguration(ed_pos, expectedNumberOfShifts=34)
plot(best, 
     method = 'polar',
     pal=viridis::inferno(99),
     logcolor=TRUE,
     # breaksmethod='quantile',  # breaksmethod='jenks'
     legend=T)
addBAMMshifts(best, method = 'polar', cex=1, bg='gold')
NROW(best$eventData[[1]]) # 31 distinct shifts in this configuration
best$eventData # which nodes have distinct shifts
getShiftNodesFromIndex(best)
ape::axisPhylo()



### mean, model-averaged rate of climatic trait evolution (for climate niche positions)
png('../fig/fig_04_rate_niche_pos.png', wid=10.5, hei=10.0, units='in',
    bg='transparent', res=700)
x <- plot(ed_pos, 
          method = 'polar', 
          pal=viridis::inferno(99), 
          breaksmethod='quantile',  # breaksmethod='jenks'
          legend=T) 
dev.off()
# den <- x$colordens # rate densities


# phylorate histogram 
`rate_hist` <- function (phylorates, xlim=NULL, ...) {
  j <-  c("coords", "colorbreaks","palette", "colordens")
  stopifnot(identical(names(phylorates),j))
  x <- phylorates$colordens[,1]
  y <- phylorates$colordens[,2]
  if(is.null(xlim)){
    xlim <- c(min(x), max(x))
  } 
  plot(NA, xlim = xlim, ylim = c(0, max(y)), ylab='density', xlab='trait rates', ...)
  segments(x, y, x, 0, lend = 2, col = phylorates$colordens[, 3], lwd = 3)
  lines(x,y, col='grey')
}
set_par_mercury(1)
rate_hist(x, xlim=c(0,2.2))


### rate shift probabilities as branch lengths
mst <- marginalShiftProbsTree(ed_pos)   # marginal shift probs tree
cst <- cumulativeShiftProbsTree(ed_pos) # cumulative shift probs tree
set_par_mercury(2)
plot.phylo(mst, no.margin=TRUE, show.tip.label=FALSE, type='fan')
plot.phylo(cst, no.margin=TRUE, show.tip.label=FALSE, type='fan')
# ### plot subtrees while preserving the original rates (colorbreaks)
# identify(phy) # locate node(s) graphically on tree...
# i_node <- 4263
# plot(subtreeBAMM(ed_pos, node = i_node), lwd = 1, colorbreaks = x$colorbreaks,
#      pal=viridis::viridis(99))


### 95 % credible set of rate shift configurations sampled with BAMM 
cset <- credibleShiftSet(ed_pos, expectedNumberOfShifts=34, threshold=3)
summary(cset)
cset$number.distinct
plot.credibleshiftset(cset, method='polar')


### macroevolutionary cohort analysis:
#     pairwise prob that 2 species share a common macroevolutionary rate dynamic
cmat <- getCohortMatrix(ed_pos)
cohorts(cmat, ed_pos, use.plot.bammdata=TRUE)


### plot evolutionary rates through time (separately for three clades)
ecole::set_par_mercury(3)
yl <- c(0,2.7)
st <- max(branching.times(phy))
plotRateThroughTime(ed_pos, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
text(x=200, y= 0.2, label='All taxa', font=4, cex=2.0, pos=4)
plotRateThroughTime(ed_pos, intervalCol='blue', avgCol='blue', start.time=st, 
                    node=i_node, nodetype = 'include', ylim=yl)
text(x=200, y= 0.2, label='Clade A', font=4, cex=2.0, pos=4)
plotRateThroughTime(ed_pos, intervalCol='green', avgCol='green', start.time=st, 
                    node=i_node, nodetype = 'exclude', ylim=yl)
text(x=200, y= 0.2, label='Clade B', font=4, cex=2.0, pos=4)


### plot evolutionary rates through time (overplot three clades)
ecole::set_par_mercury(1)
plotRateThroughTime(ed_pos, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
plotRateThroughTime(ed_pos, intervalCol='blue', avgCol='blue', start.time=st, 
                    node=i_node, nodetype = 'include', add=T)
# plotRateThroughTime(ed_pos, intervalCol='green', avgCol='green', start.time=st, 
#                     node=i_node, nodetype = 'exclude', add=T)






### plot geologic periods on fan phylogeny
`plot_geologic_rings` <- function(phy, legend = TRUE, ...) {
  require(phytools)
  require(plotrix)
  plotTree(phy, ftype='off', ylim=c(-0.2*Ntip(phy),Ntip(phy)),
           xlim=c(max(nodeHeights(phy)),0),
           direction='leftwards')
  obj <- phytools::geo.legend()
  r   <- max(obj$leg[,1]) - obj$leg[,2]
  plotTree(phy, type='fan', fsize=0.7, ftype='off', color='transparent')
  for(i in 1:nrow(obj$leg)) {
    plotrix::draw.circle(0,0,
                         radius=r[i],
                         col=obj$colors[i],
                         border='transparent')
  }
  par(fg='transparent')
  plotTree(phy, type='fan', add=TRUE, fsize=0.7, lwd=0.3, ftype='off')
  if (legend) {
    par(fg='black')
    add.simmap.legend(colors=obj$colors[rownames(obj$leg)],
    prompt=FALSE, x=0.95*par()$usr[1], y=100+0.6*par()$usr[3])
  }
}
plot_geologic_rings(phy)
### add shift nodes
shift_nodes <- getShiftNodesFromIndex(ed_pos, index = msc$sampleindex)
nodelabels(node=shift_nodes, pch=21, col='gold', cex=0.5)











################################################

### STRAPP analyses
strapp <- traitDependentBAMM(ephy = ed_div,
                             traits = tra,
                             reps = 999,
                             # return.full = TRUE,
                             method = 'spearman')
str(strapp)                      # correlation and its p-value
hist(strapp$obs.corr, breaks=55) # distribution of correlations from permutations
abline(v=strapp$estimate, lwd=3)



####    END    ####