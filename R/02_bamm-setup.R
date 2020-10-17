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
write.tree(p, file = 'lichen_phy.tre')                # phylogeny in bamm wd
write.table(tr_position, file = 'lichen_tra_pos.txt') # trait positions PCA
write.table(tr_breadth, file = 'lichen_tra_brd.txt')  # trait breadth PCA
rm(is_breadth, pc, pa, pb, tr_position, tr_breadth)

### setup dirs and filenames
fnm_phy <- 'lichen_phy.tre'      # phylogeny
fnm_tra <- 'lichen_tra_pos.txt'  # traits
fnm_ctl <- 'lichen_ctl.txt'      # control (to be done below)
system( 'ls -F' )                # check which files are in here
list.files()                     # check which files are in here

### load tree and traits
phy <- read.tree(fnm_phy)
tra <- read.table(fnm_tra)
head(tra)
phy      # 2577 tips
str(tra)

### generate *custom* control file for LICHENS
(priors <- setBAMMpriors(phy    = phy,
                        traits  = tra,
                        outfile = NULL)) # TODO names must match............. 
generateControlFile(file = fnm_ctl,
                    type = 'trait',
                    params = list(
                      treefile   = fnm_phy,
                      traitfile  = fnm_tra, 
                      globalSamplingFraction = '1',
                      numberOfGenerations = '100000',
                      lambdaInitPrior  = as.numeric(priors['lambdaInitPrior']),
                      lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
                      muInitPrior      = as.numeric(priors['muInitPrior']),
                      expectedNumberOfShifts = '34',     # per PNAS
                      lambdaIsTimeVariablePrior   = '0', # per PNAS
                    ))
rm(priors)
system( paste0('cat ', fnm_ctl) ) # verify 
# ### generate a custom control file for PRIMATES
# priors <- setBAMMpriors(phy     = phy, 
#                         traits  = fnm_tra, 
#                         outfile = NULL)
# generateControlFile(fnm_ctl,
#                     type = 'trait',     # for traits analysis
#                     params = list(
#                       treefile               = fnm_phy,
#                       traitfile              = fnm_tra, 
#                       overwrite = '0',                   # do NOT overwrite
#                       # globalSamplingFraction = '0.98',
#                       numberOfGenerations    = '1000000',
#                       # lambdaInitPrior        = '1.889',
#                       # lambdaShiftPrior       = '0.032',
#                       # muInitPrior            = '1.889',
#                       expectedNumberOfShifts      = '1',
#                       # lambdaIsTimeVariablePrior   = '0',
#                       betaInitPrior  = as.numeric(priors['betaInitPrior']),
#                       betaShiftPrior = as.numeric(priors['betaShiftPrior']),
#                       useObservedMinMaxAsTraitPriors = 
#                         as.numeric(priors['useObservedMinMaxAsTraitPriors'])
#                     ))
# # numberOfGenerations = 100000 # Number of generations to perform MCMC simulation
# # mcmcWriteFreq = 5000 # Frequency in which to write the MCMC output to a file 
# system( paste0('cat ', fnm_ctl) ) # verify 



#############################################################################
### ! ! ! RUN THE BAMM ANALYSIS ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
system( 'bamm --version' )              # BAMM 2.5.0 (2015-11-01)
system( paste0('bamm -c ', fnm_ctl) )   # system( 'bamm -c traitcontrol.txt' )
### ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
#############################################################################





### load event data after BAMM
ed <- getEventData(phy = phy,
                   eventdata = 'event_data.txt',
                   burnin = 0.10,
                   nsamples = NULL,
                   verbose = TRUE,
                   type = 'trait')

### plot rate of climatic trait evolution (climate positions)
x <- plot(ed, pal=viridis::viridis(99), breaksmethod='quantile') # 'jenks'
### plot subtrees while preserving the original rates (colorbreaks)
identify(phy) # locate node(s) graphically on tree...
plot(subtreeBAMM(ed, node = 386), lwd = 3, colorbreaks = x$colorbreaks)
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
i_node  <- 386
ecole::set_par_mercury(3)
yl <- c(0,0.3)
st <- max(branching.times(phy))
plotRateThroughTime(ed, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
text(x=30, y= 0.8, label='All taxa', font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, intervalCol='blue', avgCol='blue', start.time=st, 
                    node=i_node, nodetype = 'include', ylim=yl)
text(x=30, y= 0.8, label='Clade A', font=4, cex=2.0, pos=4)
plotRateThroughTime(ed, intervalCol='green', avgCol='green', start.time=st, 
                    node=i_node, nodetype = 'exclude', ylim=yl)
text(x=30, y= 0.8, label='Clade B', font=4, cex=2.0, pos=4)


### plot evolutionary rates through time (overplot three clades)
ecole::set_par_mercury(1)
yl <- c(0,0.3)
st <- max(branching.times(phy))
plotRateThroughTime(ed, intervalCol='red', avgCol='red', start.time=st, ylim=yl)
plotRateThroughTime(ed, intervalCol='blue', avgCol='blue', start.time=st, 
                    node=i_node, nodetype = 'include', add=T)
plotRateThroughTime(ed, intervalCol='green', avgCol='green', start.time=st, 
                    node=i_node, nodetype = 'exclude', add=T)





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
