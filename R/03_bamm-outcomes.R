######################################################################
#
#  BAMM outcomes: analyze tests for rateshifts in climatic niche evolution
#
#    Rob Smith, phytomosaic@gmail.com, 19 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################

#

### preamble
rm(list=ls())
require(ecole)             # for plotting and convenience functions
require(BAMMtools)         # for BAMM
require(coda)              # for MCMC diagnostics
require(phytools)          # for ancestral state reconstruction
system( 'bamm --version' ) # BAMM 2.5.0 (2015-11-01)


### setup dirs and filenames
setwd('/home/rob/prj/evoclim/bamm')   # change into working dir (R)
fnm_phy       <- 'lichen_phy.tre'     # phylogeny
fnm_pos_tra   <- 'lichen_tra_pos.txt' # traits (niche position)
fnm_brd_tra   <- 'lichen_tra_brd.txt' # traits (niche breadth)
fnm_pos_event <- 'pos_event_data.txt' # event data (niche position)
fnm_brd_event <- 'brd_event_data.txt' # event data (niche breadth) 
fnm_div_event <- 'div_event_data.txt' # event data (niche breadth) 
fnm_pos_mcmc  <- 'pos_mcmc_out.txt'   # mcmc data (niche position)
fnm_brd_mcmc  <- 'brd_mcmc_out.txt'   # mcmc data (niche breadth)
fnm_div_mcmc  <- 'div_mcmc_out.txt'   # mcmc data (niche breadth)
es            <- 34                   # expected number of shifts (set a priori)

### load tree and traits
phy <- read.tree(fnm_phy)      # load tree
tr  <- cbind(read.table(fnm_pos_tra,row.names=1),read.table(fnm_brd_tra)[,2])
names(tr) <- c('pc_position','pc_breadth')
head(tr)

### PCA histograms
png('../fig/fig_04_pca_histograms.png', wid=6.5, hei=3.25, units='in',
    bg='transparent', res=700)
set_par_mercury(2)
hist(read.table(fnm_pos_tra)[,2], breaks=77, main='', xlab='PC1 niche positions') 
text(c(-12, 0, 12), c(100,100,100),c('Polar','Temperate','Tropical'), cex=0.7)
hist(read.table(fnm_brd_tra)[,2], breaks=77, main='', xlab='PC1 niche breadths') 
text(c(-14, 0), c(150,150),c('Generalist','Specialist'), cex=0.7)
dev.off()

### do quick ancestral state reconstructions for PCA scores
trm <- lapply(1:NCOL(tr), function(j) { # ! ! ! TIMEWARN ! ! ! ~1-2 min per trait
  cat(paste0(j, ' of ', NCOL(tr), '... '))
  x   <- tr[,j]         # grab trait of interest
  names(x) <- dimnames(tr)[[1]]  # assign species names
  o <- phytools::contMap(phy, x, plot=F)  # continuous trait reconstruction
  o$cols[1:length(o$cols)] <- viridis::inferno(length(o$cols))  # change colors
  o
})
png('../fig/fig_05_anc_state_reconstrxn.png', wid=10.5, hei=5.25, units='in',
    bg='transparent', res=700)
set_par_mercury(2)
plot(trm[[1]], lwd=0.5, outline=F, leg.txt='PC1 niche positions', 
     type='fan', fsize=c(0.1,1))
plot(trm[[2]], lwd=0.5, outline=F, leg.txt='PC1 niche breadths', 
     type='fan', fsize=c(0.1,1))
dev.off()
rm(trm,tr)

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
mc_brd <- read.csv(fnm_brd_mcmc, header=T)
mc_div <- read.csv(fnm_div_mcmc, header=T)
`convergence_diagnostics` <- function(mc_pos){
  burnstart <- floor(0.1 * nrow(mc_pos)) # discard first 10% of samples as burnin
  postburn  <- mc_pos[burnstart:nrow(mc_pos), ]
  # effective sample sizes of number of shift events in each sample:
  print(coda::effectiveSize(postburn$N_shifts))
  # effective sample sizes of the log-likelihood:
  print(coda::effectiveSize(postburn$logLik))
  # sketch whether run converged 
  set_par_mercury(2)
  plot(mc_pos$logLik ~ mc_pos$generation) 
  abline(v=mc_pos$generation[burnstart], col=2)
  plot(mc_pos$logPrior ~ mc_pos$generation)
  abline(v=mc_pos$generation[burnstart], col=2)
}
convergence_diagnostics(mc_pos)
convergence_diagnostics(mc_brd)
convergence_diagnostics(mc_div)


# compare prior and posterior
set_par_mercury(3)
plotPrior(mc_pos, expectedNumberOfShifts=es) 
plotPrior(mc_brd, expectedNumberOfShifts=es) 
plotPrior(mc_div, expectedNumberOfShifts=es)


### identify 'best' configurations (3 alternatives):
# ### alternative 1 --- maximize Bayes factors for all pairwise models
# shift_probs <- summary(ed_pos)     # shift probabilities
# pr_null     <- 1/NROW(mc_pos) # probability of null no-shift model
# bfmat <- computeBayesFactors(mc_pos, expectedNumberOfShifts=es, burnin=0.1)
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
#      legend=T)
# addBAMMshifts(msc_i, method = 'polar', cex=1, bg='gold')
### alternative 3 --- configuration that maximizes a posteriori probability
best_pos <- getBestShiftConfiguration(ed_pos, expectedNumberOfShifts=es)
best_brd <- getBestShiftConfiguration(ed_brd, expectedNumberOfShifts=es)
best_div <- getBestShiftConfiguration(ed_div, expectedNumberOfShifts=es)
NROW(best_pos$eventData[[1]]) # 31 distinct shifts in this configuration
NROW(best_brd$eventData[[1]]) # 24 distinct shifts in this configuration
NROW(best_div$eventData[[1]]) # 7 distinct shifts in this configuration???

### TODO trying to plot geological periods as circles...
# set_par_mercury(1)
# plot(best_pos, method = 'polar', pal=viridis::inferno(99), 
#      breaksmethod = 'linear', color.interval = c(0,3), legend=F)
# set_par_mercury(1)
# plot(ed_pos, method = 'polar', pal=viridis::inferno(99), 
#      breaksmethod = 'linear', color.interval = c(0,3), legend=F)
# addBAMMshifts(best_pos, method = 'polar', cex=1, bg='gold')
# best_pos$eventData # which nodes have distinct shifts
# sort(getShiftNodesFromIndex(best_pos,1)) # doesnt match line below
# best_pos$eventData[[1]]$node             # doesnt match line above...
# # ape::axisPhylo()
# # require(plotrix)
# # for(i in c(0,.4,.9, 1.1)) {
# #   plotrix::draw.circle(0,0,
# #                        radius=i,
# #                        col='#00000040',
# #                        border='transparent')
# # }
# points(0,0,cex=1.1, col='blue')
# u <- par('usr')
# uu <- u * 0.95
# max(branching.times(phy))
# max(nodeHeights(phy))
# points(0,0,cex=1.1, col='blue')
# points(uu[1], uu[4], cex=1.1, col='red')
# points(uu[2], uu[3], cex=1.1, col='gold')
# `circle` <- function(x,y,r=1,border="black",lty="solid",lwd=1,fill=NULL) {
#   xapu <- sin(seq(0,pi,length=50)-pi/2)
#   for (i in 1:length(x)) {
#     xv1 <- x[i] + xapu*r[i]
#     xv2 <- x[i] - xapu*r[i]
#     yv1 <- sqrt(pmax(0,r[i]^2-(xv1-x[i])^2))+y[i]
#     yv2 <- -sqrt(pmax(0,r[i]^2-(xv2-x[i])^2))+y[i]
#     yv  <- c(yv1,yv2)
#     xv  <- c(xv1,xv2)
#     polygon(xv,yv,border=border,col=fill,lty=lty,lwd=lwd)
#   }
# }
# circle(0,0,1)
# circle(0,0,u[1])

### mean model-averaged rate of climatic trait evolution (for climate niche positions)
png('../fig/fig_06_mean_rate_niche_and_div.png', wid=12.5, hei=4.25, units='in',
    bg='transparent', res=700)
set_par_mercury(3)
plot(ed_pos, method = 'polar', pal=viridis::inferno(99), 
     breaksmethod = 'linear', color.interval = c(0,3), legend=F)
addBAMMshifts(best_pos, method = 'polar', cex=0.8, bg='gold', par.reset = F)
title('Rates of evolution of\nclimate niche position', cex.main=0.7)
plot(ed_brd, method = 'polar', pal=viridis::inferno(99), 
     breaksmethod = 'linear', color.interval = c(0,3), legend=F)
addBAMMshifts(best_brd, method = 'polar', cex=0.8, bg='gold', par.reset = F)
title('Rates of evolution of\nclimate niche breadth', cex.main=0.7)
plot(ed_div, method = 'polar', pal=viridis::inferno(99),
     breaksmethod = 'linear', 
     # color.interval = c(0,3), 
     legend=F)
addBAMMshifts(best_div, method = 'polar', cex=1, bg='gold')
title('Rate of diversification')
dev.off()


# ### phylorate histograms
# `rate_hist` <- function (phylorates, xlim=NULL, ...) {
#   j <-  c("coords", "colorbreaks","palette", "colordens")
#   stopifnot(identical(names(phylorates),j))
#   x <- phylorates$colordens[,1]
#   y <- phylorates$colordens[,2]
#   if(is.null(xlim)){
#     xlim <- c(min(x), max(x))
#   } 
#   plot(NA, xlim = xlim, ylim = c(0, max(y)), ylab='density', xlab='trait rates', ...)
#   segments(x, y, x, 0, lend = 2, col = phylorates$colordens[, 3], lwd = 3)
#   lines(x,y, col='grey')
# }
# set_par_mercury(3)
# rate_hist(x, xlim=c(0,3.2), add=T)
# rate_hist(y, xlim=c(0,2.2), add=T)
# rate_hist(z, xlim=c(0,2.2), add=T)


# ### rate shift probabilities as branch lengths
# mst <- marginalShiftProbsTree(ed_pos)   # marginal shift probs tree
# cst <- cumulativeShiftProbsTree(ed_pos) # cumulative shift probs tree
# set_par_mercury(2)
# plot.phylo(mst, no.margin=TRUE, show.tip.label=FALSE, type='fan')
# plot.phylo(cst, no.margin=TRUE, show.tip.label=FALSE, type='fan')
# # ### plot subtrees while preserving the original rates (colorbreaks)
# # identify(phy) # locate node(s) graphically on tree...
# # i_node <- 4263
# # plot(subtreeBAMM(ed_pos, node = i_node), lwd = 1, colorbreaks = x$colorbreaks,
# #      pal=viridis::viridis(99))


# ### 95 % credible set of rate shift configurations sampled with BAMM 
# cset <- credibleShiftSet(ed_pos, expectedNumberOfShifts=es, threshold=3)
# summary(cset)
# cset$number.distinct
# plot.credibleshiftset(cset, method='polar')


### macroevolutionary cohort analysis:
#     pairwise prob that 2 species share a common macroevolutionary rate dynamic
cohorts(getCohortMatrix(ed_pos), ed_pos, lwd=0.1, use.plot.bammdata=FALSE)
cohorts(getCohortMatrix(ed_brd), ed_brd, lwd=0.1, use.plot.bammdata=FALSE)
cohorts(getCohortMatrix(ed_div), ed_div, lwd=0.1, use.plot.bammdata=FALSE)


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


### overplot............
png('../fig/fig_07_rates_vs_timez.png', wid=4.75, hei=4.25, units='in',
    bg='transparent', res=700)
yl <- c(0,1)
st <- max(branching.times(phy))
ecole::set_par_mercury(1)
plotRateThroughTime(ed_pos, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
plotRateThroughTime(ed_brd, intervalCol='blue', avgCol='blue', start.time=st, 
                    add=T)
plotRateThroughTime(ed_div, intervalCol='forestgreen', avgCol='forestgreen', 
                    ratetype='netdiv', start.time=st, 
                    add=T)
legend('topleft',legend=c('Climate position','Climate breadth','Net diversification'),
       fill=c('red','blue','forestgreen'), bty='n')
dev.off()


### TODO plot geologic periods on fan phylogeny
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