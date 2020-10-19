######################################################################
#
#  BAMM setup: run test for rateshifts in climatic niche evolution
#
#    Rob Smith, phytomosaic@gmail.com, 16 Oct 2020
#
##      GNU General Public License, Version 3.0    ###################


### preamble
rm(list=ls())
require(ecole)             # for plotting and convenience functions
require(BAMMtools)         # for BAMM
system( 'bamm --version' ) # BAMM 2.5.0 (2015-11-01)


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
rm(is_breadth, pc, pa, pb, tr_position, tr_breadth, p, tr, wd)


### setup dirs and filenames
fnm_phy     <- 'lichen_phy.tre'      # phylogeny
fnm_pos_tra <- 'lichen_tra_pos.txt'  # traits (niche position)
fnm_brd_tra <- 'lichen_tra_brd.txt'  # traits (niche breadth)
fnm_pos_ctl <- 'pos_ctl.txt'         # control (to be done below)
fnm_brd_ctl <- 'brd_ctl.txt'         # control (to be done below)
fnm_div_ctl <- 'div_ctl.txt'         # control (to be done below)
system( 'ls -F' )                    # check which files are in here
list.files()                         # check which files are in here
samp_freq <- '10000'     # '5000'      # 10,000 per PNAS
n_gens    <- '65000000'  # '5000000'   # 65 million per PNAS
phy       <- read.tree(fnm_phy)      # load tree


### generate custom control file for lichen climate niche POSITIONS
(priors <- setBAMMpriors(phy     = phy,
                         traits  = fnm_pos_tra,
                         outfile = NULL,
                         Nmax    = 9999))
generateControlFile(file = fnm_pos_ctl,
                    type = 'trait',
                    params = list(
                      treefile            = fnm_phy,
                      traitfile           = fnm_pos_tra, 
                      numberOfGenerations = n_gens,      # 65 million per PNAS
                      mcmcWriteFreq       = samp_freq,   # 10,000 per PNAS
                      eventDataWriteFreq  = samp_freq,   # 10,000 per PNAS
                      acceptanceResetFreq = samp_freq,   # 10,000 per PNAS
                      printFreq           = '50000',     # print to console
                      outName = 'pos',                   # filename prefix
                      expectedNumberOfShifts = '34',     # per PNAS
                      betaInitPrior  = as.numeric(priors['betaInitPrior']),
                      betaShiftPrior = as.numeric(priors['betaShiftPrior']),
                      betaIsTimeVariablePrior = 0, # per PNAS and Matt, keep time-constant
                      useObservedMinMaxAsTraitPriors = '1'
                    ))

### generate custom control file for lichen climate niche BREADTHS
(priors <- setBAMMpriors(phy     = phy,
                         traits  = fnm_brd_tra,
                         outfile = NULL,
                         Nmax    = 9999))
generateControlFile(file = fnm_brd_ctl,
                    type = 'trait',
                    params = list(
                      treefile            = fnm_phy,
                      traitfile           = fnm_brd_tra, 
                      numberOfGenerations = n_gens,      # 65 million per PNAS
                      mcmcWriteFreq       = samp_freq,   # 10,000 per PNAS
                      eventDataWriteFreq  = samp_freq,   # 10,000 per PNAS
                      acceptanceResetFreq = samp_freq,   # 10,000 per PNAS
                      printFreq           = '50000',     # print to console
                      outName = 'brd',                   # filename prefix
                      expectedNumberOfShifts = '34',     # per PNAS
                      betaInitPrior  = as.numeric(priors['betaInitPrior']),
                      betaShiftPrior = as.numeric(priors['betaShiftPrior']),
                      betaIsTimeVariablePrior = 0, # per PNAS and Matt, keep time-constant
                      useObservedMinMaxAsTraitPriors = '1'
                    ))

### generate custom control file for lichen species DIVERSIFICATION
(priors <- setBAMMpriors(phy     = phy,
                         traits  = NULL,
                         outfile = NULL,
                         Nmax    = 9999))
generateControlFile(file = fnm_div_ctl,
                    type = 'diversification',
                    params = list(
                      treefile            = fnm_phy,
                      numberOfGenerations = n_gens,      # 65 million per PNAS
                      mcmcWriteFreq       = samp_freq,   # 10,000 per PNAS
                      eventDataWriteFreq  = samp_freq,   # 10,000 per PNAS
                      acceptanceResetFreq = samp_freq,   # 10,000 per PNAS
                      printFreq           = '50000',     # print to console
                      outName = 'div',                   # filename prefix
                      lambdaInitPrior  = as.numeric(priors['lambdaInitPrior']),
                      lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
                      muInitPrior      = as.numeric(priors['muInitPrior']),
                      lambdaIsTimeVariablePrior   = '0', # per PNAS
                      expectedNumberOfShifts = '34'      # per PNAS
                    ))
rm(priors)


#############################################################################
###################################################################
#########################################################
### ! ! ! RUN BAMM ANALYSIS ! ! !
cat(paste0('Began at: ', (timebegin <- Sys.time()),'\n\n'))
system( paste0('bamm -c ', fnm_pos_ctl) )
cat(paste0('Ended at: ', Sys.time()), 'after', format(Sys.time() - timebegin),'\n\n')

cat(paste0('Began at: ', (timebegin <- Sys.time()),'\n\n'))
system( paste0('bamm -c ', fnm_brd_ctl) )
cat(paste0('Ended at: ', Sys.time()), 'after', format(Sys.time() - timebegin),'\n\n')

cat(paste0('Began at: ', (timebegin <- Sys.time()),'\n\n'))
system( paste0('bamm -c ', fnm_div_ctl) )
cat(paste0('Ended at: ', Sys.time()), 'after', format(Sys.time() - timebegin),'\n\n')
#########################################################
###################################################################
#############################################################################



####    END    ####
