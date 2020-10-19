######################################################################
#
#  BAMM test for rateshifts in climatic niche evolution
#
#    Rob Smith, phytomosaic@gmail.com, 16 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
require(ecole)       # for plotting and convenience functions
require(BAMMtools)   # for BAMM
require(coda)        # for MCMC diagnostics

### prune tree and traits to exact matches, then save as bamm inputs
wd <- '/home/rob/prj/evoclim/bamm'       # temporary working dir
load(file='./data/tr.rda', verbose=T)
p  <- read.tree('./data/Lecanoromycetes_ML_Tree_Timescaled')
p  <- keep.tip(p, tip = rownames(tr))
tr <- tr[match(p$tip.label, rownames(tr)),]
identical(rownames(tr), p$tip.label)      # expect TRUE
system( paste0('cd ', wd) )      # change into working dir (system)
setwd(wd)                        # change into working dir (R)

### PCA for either niche positions or niche breadth
is_breadth <- grepl('mad', dimnames(tr)[[2]]) # split niche position vs breadth
`pc` <- function(x, ...) { # PCA
  p <- prcomp(x, retx=T, center=T, scale.=T)
  print(summary(p))
  V <- t(apply(p$rotation,1,function(a,b){a*b},p$sdev)) 
  colnames(p$x) <- tolower(colnames(p$x))
  return(list(V=round(V,3), x=p$x))
}
pa <- pc(tr[, !is_breadth]) # trait position PCA
ecole::plot_heatmap(pa$V[,1:5], xexp=2.5) 
# PC1 = 'tropicality/boreality' (high vals = high mean, low seasonality/range)
pb <- pc(tr[, is_breadth])  # trait breadth PCA
ecole::plot_heatmap(pb$V[,1:5], xexp=2.5) 
# PC1 = 'generalist/specialist' (high vals = broader tolerances)
tr_position <- data.frame(rownames(pa$x), pa$x[,1])   # tropicality/boreality
tr_breadth  <- data.frame(rownames(pb$x), pb$x[,1])   # generalist/specialist
# names(tr_breadth) <- names(tr_position) <- c('spp','scr')
write.tree(p, file = 'lichen_phy.tre')                # phylogeny in bamm wd
write.table(tr_position, file = 'lichen_tra_pos.txt', sep = '\t', 
            quote = F, row.names=F, col.names=F)  # trait positions PCA
write.table(tr_breadth, file = 'lichen_tra_brd.txt', sep = '\t', 
            quote = F, row.names=F, col.names=F)  # trait breadth PCA
rm(is_breadth, pc, pa, pb, tr_position, tr_breadth)

### setup dirs and filenames
fnm_phy <- 'lichen_phy.tre'      # phylogeny
fnm_tra <- 'lichen_tra_pos.txt'  # traits
fnm_ctl <- 'lichen_ctl.txt'      # control (to be done below)
system( 'ls -F' )                # check which files are in here
list.files()                     # check which files are in here

### load tree and traits
phy <- read.tree(fnm_phy)
# tra <- read.table(fnm_tra)
# hist(tra[,2], breaks=77) # distribution of PC1 scores for niche position
# text(c(-14, 0, 12), c(100,100,100),c('Polar','Temperate','Tropical'))


### generate *custom* control file for LICHENS
(priors <- setBAMMpriors(phy     = phy,
                         traits  = 'lichen_tra_pos.txt',
                         outfile = NULL,
                         Nmax    = 9999))
generateControlFile(file = fnm_ctl,
                    type = 'trait',
                    params = list(
                      treefile            = fnm_phy,
                      traitfile           = fnm_tra, 
                      numberOfGenerations = '1000000', # 65 million per PNAS
                      mcmcWriteFreq       = '1000',   # 10,000 per PNAS
                      eventDataWriteFreq  = '1000',   # 10,000 per PNAS
                      acceptanceResetFreq = '1000',   # 10,000 per PNAS
                      printFreq           = '10000',   # how often print to console
                      outName = 'pos', # filename prefix
                      # lambdaInitPrior  = as.numeric(priors['lambdaInitPrior']),
                      # lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
                      # muInitPrior      = as.numeric(priors['muInitPrior']),
                      # lambdaIsTimeVariablePrior   = '0', # per PNAS
                      expectedNumberOfShifts = '34',     # per PNAS
                      betaInitPrior  = as.numeric(priors['betaInitPrior']),
                      betaShiftPrior = as.numeric(priors['betaShiftPrior']),
                      betaIsTimeVariablePrior = 0, # per PNAS and Matt, keep time-constant
                      useObservedMinMaxAsTraitPriors = '1'
                    ))
rm(priors)
system( paste0('cat ', fnm_ctl) ) # verify 



#############################################################################
###################################################################
#########################################################
### ! ! ! RUN BAMM ANALYSIS ! ! !
system( 'bamm --version' )              # BAMM 2.5.0 (2015-11-01)
system( paste0('bamm -c ', fnm_ctl) )
#########################################################
###################################################################
#############################################################################


### convergence diagnostics
mcmcout <- read.csv('pos_mcmc_out.txt', header=T)
plot(mcmcout$logLik ~ mcmcout$generation) # sketch whether run converged 
plot(mcmcout$logPrior ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout)) # discard first 10% of samples as burnin
postburn  <- mcmcout[burnstart:nrow(mcmcout), ]
coda::effectiveSize(postburn$N_shifts) # effective sample sizes of number of shift events in each sample
coda::effectiveSize(postburn$logLik) # effective sample sizes of the log-likelihood   


### load BAMM event data 
ed <- getEventData(phy = phy,
                   eventdata = 'pos_event_data.txt',
                   burnin = 0.10,   # 10% per PNAS
                   nsamples = NULL, # 2000 per PNAS
                   verbose = TRUE,
                   type = 'trait')


### shift probabilities and Bayes Factors
shift_probs <- summary(ed)
pr_null     <- 1/NROW(mcmcout) # probability of null no-shift model
# compare prior and posterior
set_par_mercury(1)
plotPrior(mcmcout, expectedNumberOfShifts=34) 

# Bayes factors for ALL pairwise models (highest BF best):
bfmat <- computeBayesFactors(mcmcout, expectedNumberOfShifts=34, burnin=0.1)
plot_heatmap(bfmat)
persp(bfmat, theta=135)
persp(bfmat[,1:11], theta=135)
stepBF(bfmat, step.size = 1, inputType = 'matrix') # optimal number of shifts
# plot shift probabilities and Bayes factors
set_par_mercury(3)
plot(shift_probs[,1], shift_probs[,2], type='b', ylab='Shift probabilities')
plot(shift_probs[,1], shift_probs[,2] / pr_null, type='b',
     ylab='Min Bayes factor, relative to null no-shift model')
plot(shift_probs[,1], bfmat[,1], type='b',
     ylab='Bayes factors, relative to 9-shift model')




### plot rate of climatic trait evolution (for climate niche positions)
x <- plot(ed, pal=viridis::viridis(99), breaksmethod='quantile', legend=T)
# x <- plot(ed, pal=viridis::viridis(99), breaksmethod='jenks')


### plot subtrees while preserving the original rates (colorbreaks)
identify(phy) # locate node(s) graphically on tree...
i_node <- 4263
plot(subtreeBAMM(ed, node = i_node), lwd = 1, colorbreaks = x$colorbreaks,
     pal=viridis::viridis(99))
# ### `dtRates` if you want to change tau and rates for entire tree
# tau <- 0.002
# ed  <- dtRates(ed, tau = tau)
# ecole::set_par(3)
# x  <- plot(ed, lwd = 3, spex = 's')             # full tree at finer tau
# sed <- dtRates(subtreeBAMM(ed, node = 421), tau) # subtree at finer tau
# plot(sed, lwd = 3, colorbreaks = x$colorbreaks)
# sed <- dtRates(subtreeBAMM(ed, node = 287), tau) # subtree at finer tau
# plot(sed, lwd = 3, colorbreaks = x$colorbreaks)


### 95 % credible set of rate shift configurations sampled with BAMM 
cset <- credibleShiftSet(ed, expectedNumberOfShifts=1, threshold=3)
summary(cset)
plot.credibleshiftset(cset, lwd=2.5)


### macroevolutionary cohort analysis:
#     pairwise prob that 2 species share a common macroevolutionary rate dynamic
cmat <- getCohortMatrix(ed)
cohorts(cmat, ed, use.plot.bammdata=TRUE)


### plot evolutionary rates through time (separately for three clades)
ecole::set_par_mercury(3)
yl <- c(0,2.7)
st <- max(branching.times(phy))
plotRateThroughTime(ed, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
text(x=200, y= 0.2, label='All taxa', font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, intervalCol='blue', avgCol='blue', start.time=st, 
                    node=i_node, nodetype = 'include', ylim=yl)
text(x=200, y= 0.2, label='Clade A', font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, intervalCol='green', avgCol='green', start.time=st, 
                    node=i_node, nodetype = 'exclude', ylim=yl)
text(x=200, y= 0.2, label='Clade B', font=4, cex=2.0, pos=4)


### plot evolutionary rates through time (overplot three clades)
ecole::set_par_mercury(1)
plotRateThroughTime(ed, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
plotRateThroughTime(ed, intervalCol='blue', avgCol='blue', start.time=st, 
                    node=i_node, nodetype = 'include', add=T)
# plotRateThroughTime(ed, intervalCol='green', avgCol='green', start.time=st, 
#                     node=i_node, nodetype = 'exclude', add=T)





################################################

### STRAPP analyses
strapp <- traitDependentBAMM(ephy = ed,
                             traits = tra,
                             reps = 999,
                             # return.full = TRUE,
                             method = 'spearman')
str(strapp)                      # correlation and its p-value
hist(strapp$obs.corr, breaks=55) # distribution of correlations from permutations
abline(v=strapp$estimate, lwd=3)



####    END    ####
