######################################################################
# Using TREEBASE to assemble a lichen phylogeny
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 24 Jan 2019
##  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
#      useful site:  http://www.phytools.org/Cordoba2017/

###   preamble
rm(list=ls())
# devtools::install_github('phytomosaic/vuln', upgrade=F)
require(ecole)
require(picante)
require(phytools)
require(vuln)

data(lichen_spe, lichen_id, lichen_mex, crosswalk)
spe <- lichen_spe   # FIA species
id  <- lichen_id    # FIA sites
mex <- lichen_mex   # CNALH occurrences
rm(lichen_spe, lichen_id, lichen_mex) # clean up

setwd('C:/Users/Rob/Documents/inactive_prj/2019_phy_lichen/')

### available taxa, from Matt Nelsen, 06 Feb 2019
h   <- read.table('fia_taxa_to_robert.txt', sep = '\t')
h   <- as.character(unlist(h))
h   <- h[h!='']
has <- h[2:35]
not <- h[37:58]
has <- unlist(strsplit(has, '\\s+'))
has <- has[has!='']
has <- has[!grepl('\\[', has)]
not <- unlist(strsplit(not, '\\s+'))
not <- not[not!='']
not <- not[!grepl('\\[', not)]
may <- c('Cetraria_canadensis','Cetraria_oakesiana','Cetraria_pallidula',
         'Cetraria_pinastri','Cetraria_platyphylla','Cetraria_viridis',
         'Cladonia_squamosa_v_subsquamosa','Crespoa_crozalsiana','Enchylium_conglomeratum',
         'Heterodermia_galactophylla', 'Heterodermia_tropica',
         'Hypogymnia_hultenii','Hypogymnia_lophyrea','Hypotrachyna_horrescens',
         'Hypotrachyna_minarum','Lobaria_anomala','Lobaria_anthraspis',
         'Lobaria_scrobiculata','Nephromopsis_aurescens','Nephromopsis_chlorophylla',
         'Nephromopsis_ciliaris','Nephromopsis_coralligera', 'Nephromopsis_fendleri',
         'Nephromopsis_merrillii','Nephromopsis_orbata','Nephromopsis_sphaerosporella',
         'Nephromopsis_subalpina','Nephromopsis_weberi','Parmoterma_cetratum',
         'Physciella_melanchra','Polycauliona_tenax','Pyxine_caesiopruinosa',
         'Scytinium_cellulosum','Scytinium_lichenoides','Scytinium_palmatum',
         'Scytinium_polycarpum','Scytinium_teretiusculum','Sticta_wrightii',
         'Xanthomendoza_montana','Xanthomendoza_ulophyllodes')
mayr <- c('Vulpicida_canadensis','Usnocetraria_oakesiana','Ahtiana_pallidula',
          'Vulpicida_pinastri','Tuckermannopsis_platyphylla','Vulpicida_viridis',
          'Cladonia_squamosa','Canoparmelia_crozalsiana','Collema_conglomeratum',
          'Anaptychia_galactophylla','Anaptychia_tropica',
          'Cavernularia_hultenii','Cavernularia_lophyrea','Parmelinopsis_horrescens',
          'Parmelinopsis_minarum','Pseudocyphellaria_anomala','Pseudocyphellaria_anthraspis',
          'Lobarina_scrobiculata','Ahtiana_aurescens','Tuckermannopsis_chlorophylla',
          'Tuckermannopsis_ciliaris','Tuckermanella_coralligera','Tuckermanella_fendleri',
          'Kaernefeltia_merrillii','Tuckermannopsis_orbata','Ahtiana_sphaerosporella',
          'Tuckermannopsis_subalpina','Tuckermanella_weberi','Parmotrema_cetratum',
          'Phaeophyscia_melanchra','Xanthoria_tenax','Pyxine_albovirens',
          'Leptogium_cellulosum','Leptogium_lichenoides','Leptogium_palmatum',
          'Leptogium_polycarpum','Leptogium_teretiusculum','Dendriscosticta_wrightii',
          'Xanthoria_montana','Xanthoria_ulophyllodes')
key <- data.frame(may, mayr, stringsAsFactors=F)
not <- not[!not %in% may]
has <- c(has,may)
names(spe) <- crosswalk$species[match(names(spe), crosswalk$abbrev)]
relfrq <- colSums((spe>0)*1) / dim(spe)[[1]]
d     <- data.frame(from_matt  = c(has,not),
                  tree_has = c(rep(1,length(has)),rep(0,length(not))))
d$relfrq <- round(relfrq[match(d$from_matt,names(relfrq))],5)
fia <- read.csv('fia_lichen_species.csv', stringsAsFactors=F)
fia <- cbind(fia, d[match(fia$species, d$from_matt),])
fia$toquery <- key$mayr[match(fia$species,key$may)]
fia$toquery[is.na(fia$toquery)] <- fia$species[is.na(fia$toquery)]
fia$tree_has<-((fia$tree_has + as.numeric(fia$from_matt %in% may))>0)*1
fia <- fia[,c('abbr','species','toquery','tree_has','relfrq')]
names(fia)[2:3] <- c('old_query','new_query')
# write.csv(fia, 'fia_taxa_to_matt_06Feb2019.csv', row.names=F)
rm(h, has, not, may, relfrq, d)

# ### calculate vulnerability
# m_mwmt <- bingrid(mex, field='mwmt', nr=155, nc=200,
#                   xmn=-177.9, xmx=-52.63, ymn=15.37, ymx=82.4)
# m_cmd  <- bingrid(mex, field='cmd', nr=155, nc=200,
#                   xmn=-177.9, xmx=-52.63, ymn=15.37, ymx=82.4)
# v1 <- vuln(spe, y=id$mwmt, ybin=m_mwmt)
# v2 <- vuln(spe, y=id$cmd,  ybin=m_cmd)
# id <- cbind(id, v1=v1$t1, v2=v1$t2, v3=v1$t3)
# id <- cbind(id, vv1=v2$t1, vv2=v2$t2, vv3=v2$t3)
# rm(v1,v2,m_mwmt,m_cmd)

### manual corrections to species names
crosswalk$species[crosswalk$species=='Parmoterma_cetratum'] <-
     'Parmotrema_cetratum'
crosswalk$species[crosswalk$species=='Rimelia_subisidiosa'] <-
     'Parmotrema_subisidiosum'
crosswalk$species[crosswalk$species=='Canomaculina_subtinctoria'] <-
     'Parmotrema_subtinctorium'
crosswalk$species[crosswalk$species=='Dendriscocaulon_intricatulum'] <-
     'Ricasolia_quercizans'
crosswalk[crosswalk$abbrev == 'Brypik','species'] <- 'Bryoria_pikei'
species<- crosswalk$species[match(names(spe), crosswalk$species)]
genera <- sub('\\_.*', '', species)
names(spe) <- species
length(unique(genera)) # 73 genera

### read in the phylogeny
p    <- read.nexus('S22758.nex') ### Kraichak et al 2018
# p <- read.nexus('S15652.nex') ### Miadlikowska et al 2014
p <- read.tree('Lecanoromycetes_ML_Tree_Timescaled') # Nelsen et al 2020
plot(p, cex=0.4)
pspe <- p$tip.label
pgen <- sub('\\_.*', '', pspe)
unique(genera[!genera %in% pgen]) # lacks 23 of 73

### prune to keep only FIA species (or genera)
p  <- keep.tip(p, tip = na.omit(match(species, pspe)))

# ### collapse each genus to polytomy at its MRCA
# `collapse_genus` <- function(p, ...){
#      pspe <- p$tip.label
#      pgen <- sub('\\_.*', '', p$tip.label)
#      s <- aggregate(pspe ~ pgen, cbind(pspe, pgen), c)[[2]]
#      for (i in seq_along(s)){
#           tipgenus  <- pspe[s[[i]]]
#           if (!length(tipgenus) %in% c(1,2)){
#                node <- getMRCA(p, tip=tipgenus)
#                p    <- collapse.to.star(p, node)
#           }
#      }
#      return(p)
# }
# p <- collapse_genus(p)

# ### TODO add species from chklst to a tree, based on genus-matching
# `add_species` <- function(p = p, splist = species, ...) {
#      ### function to get ancestor of TIPS
#      `ancestor` <- function (x, node = NULL, ...) { ###
#           parents <- x$edge[, 1]
#           child   <- x$edge[, 2]
#           pvec    <- integer(max(x$edge))
#           pvec[child] <- parents
#           if(is.null(node)){
#                node <- child[child %in% 1:Ntip(x)]
#           }
#           return(pvec[node])
#      }
#      minlen <- min(p$edge.length[p$edge.length != 0])
#      pspe   <- p$tip.label
#      pgen   <- sub('\\_.*', '', p$tip.label)
#      addspe <- splist[!splist %in% pspe]
#      addgen <- sub('\\_.*', '', addspe)
#      has    <- addgen %in% pgen
#      addspe <- addspe[has]
#      addgen <- addgen[has]
#      ag     <- aggregate(pgen, list(pgen), length) # count per genus
#      ct     <- ag[match(pgen, ag[,1]),2]
#      x      <- data.frame(addgen, addspe, ct=ct[match(addgen,pgen)],
#                           stringsAsFactors=FALSE)
#      x$brlen <- ifelse(x$ct==1, minlen, 0)  # if monotypic
#      x       <- x[order(x$addgen, x$addspe),]
#      dimnames(x)[[1]] <- 1:dim(x)[[1]]
#      x$brlen <- (!duplicated(x$addgen)) * x$brlen
#      pp <- p # initialize second tree
#      for (j in 1:dim(x)[[1]]){  # for each addspe j
#           ####   for debugging   ####
#           if (j %in% c(157,165,227,236)) { next }
#           cat(j,'\n')
#           # recalculate parent node at every iteration
#           o <- pp$edge[pp$edge[,2] <= length(pp$tip.label), 2]
#           d <- data.frame(parentnode = ancestor(pp),
#                           pgen = sub('\\_.*', '', pp$tip.label[o]))
#           x$parentnode <- d$parentnode[match(x$addgen, d$pgen)]
#           pp <- bind.tip(tree      = pp,
#                          tip.label = x$addspe[j],
#                          where     = x$parentnode[j],
#                          position  = x$brlen[j])
#      }
#      return(pp)
# }
# p <- add_species(p, species)

# prune to keep only FIA species
sk <- seq_along(p$tip.label) * (p$tip.label %in% species)
p  <- keep.tip(p, tip = sk[sk != 0])  ;  rm(sk)

### phylogenetic signal of climatic niches (climate niche conservatism?)
n <- read.csv('~/_prj/4_vuln/pkg/vuln/fig/S2_nichesummary.csv',
              stringsAsFactors = FALSE)
n$Species[n$Abbreviation == 'Brypik'] <- 'Bryoria_pikei'
n$Species[n$Abbreviation == 'Pmosub'] <- 'Parmotrema_subtinctorium'
dimnames(n)[[1]] <- n$Species
n  <- n[,-c(1:3)]
n  <- n[na.omit(match(p$tip.label, dimnames(n)[[1]])),]
all.equal(p$tip.label, dimnames(n)[[1]])
tr <- n$p95_mwmt
names(tr) <- dimnames(n)[[1]]
cc <- colvec(tr)

?phytools::phylosig
### phylogenetic signal (climate niche conservatism?)
lam <- sapply(n, FUN=function(j){
     names(j) <- rownames(n)
     l1 <- phytools::phylosig(p,j,method='K',test=T,nsim=999)
     l2 <- phytools::phylosig(p,j,method='lambda',test=T,nsim=999)
     round(c(k=l1$K, pval_k=l1$P, lambda=l2$lambda, pval_l=l2$P),4)})
lam

### visualize ALL climate niche values across phylogeny
phylo.heatmap(p, nn, 0.4, colors=viridis::inferno(99), standardize=T)

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

### continuous trait reconstruction _acros s branches_
trm <- contMap(p, x=tr, lims=c(16.5,33.5), plot=FALSE)
trm$cols[1:length(trm$cols)] <- viridis::inferno(length(trm$cols))
### continuous trait reconstruction _at each node_
fa <- fastAnc(p, tr)

# across branches:
plot(trm, fsize=0.4, lwd=2, outline=F)
# at each node:
plot(p, cex=0.3, label.offset=5, no.margin=T, tip.col=cc,edge.width=2)
nodelabels(pch=19, col = colvec(fa), cex = 1)
tiplabels(pch=19, col = colvec(tr), cex = 1)

### Ancestral trait reconstruction
tiff('C:/Users/Rob/Desktop/phy_traitreconstruction.tif',
     wid=3.5, hei=7.5, unit='in',
     compr='lzw+p', res=700)
# png('C:/Users/Rob/Desktop/phy_traitreconstruction.png',
#           wid=3.5, hei=7.5, unit='in', bg='transparent', res=3000)
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
# as 3D
phylomorphospace3d(p, X=cbind(X, n[,c('iqr_mwmt')]),
                   control=list(spin=F,box=F,lwd=2), method='dynamic')

# ### community phylogenetic diversity ! ! ! ! TIMEWARN ! ! ! !
# id <- cbind(id, psd(spe, phy, compute.var=FALSE))
# head(id)
# print(mmap(id, 'v1', colorscale='inferno'))
# print(mmap(id, 'PSV', colorscale='inferno'))  # 0 = more related, 1 = unrelated
# print(mmap(id, 'PSC', colorscale='inferno'))  # 0 = more related, 1 = unrelated
# print(mmap(id, 'PSR', colorscale='inferno'))
# print(mmap(id, 'SR', colorscale='inferno'))
# set_par(4)
# with(id, plot(SR, v1, pch=19, col='#00000020'))
# with(id, plot(PSR, v1, pch=19, col='#00000020'))
# with(id, plot(PSV, v1, pch=19, col='#00000020'))
# with(id, plot(PSC, v1, pch=19, col='#00000020'))
# set_par(4)
# with(id, plot(SR, v2, pch=19, col='#00000020'))
# with(id, plot(PSR, v2, pch=19, col='#00000020'))
# with(id, plot(PSV, v2, pch=19, col='#00000020'))
# with(id, plot(PSC, v2, pch=19, col='#00000020'))
# set_par(4)
# with(id, plot(SR, v3, pch=19, col='#00000020'))
# with(id, plot(PSR, v3, pch=19, col='#00000020'))
# with(id, plot(PSV, v3, pch=19, col='#00000020'))
# with(id, plot(PSC, v3, pch=19, col='#00000020'))

###   END   ##########################################################