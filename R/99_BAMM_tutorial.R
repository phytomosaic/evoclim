
# File includes:
#	Code examples from manuscript
#	Code used to produce figures 1, 2, 3

##########################################



require(BAMMtools)
data(whales)
data(events.whales)
data(prior.whales)

# create bammdata object:
x <- getEventData(whales, events.whales, burnin=0.1)

summary(x) # summarize posterior

# mean phylorate plot:
plot.bammdata(x, lwd=2)
# see other plotting and color options: ?plot.bammdata

pset <- getBranchShiftPriors(whales, prior.whales)

# here is the credible set of shift configurations
css <- credibleShiftSet(x, threshold=pset) 

summary(css)

plot.credibleshiftset(css, plotmax=6, pal='temperature')

# or just call:
plot(css, plotmax=6, pal='temperature')

# Rate-through-time analysis
plotRateThroughTime(x)
#numerous plotting options at ?plotRateThroughTime

# Plotting for dolphins only (node 140)
plotRateThroughTime(x, node = 140) 

# plot while EXCLUDING dolphins
plotRateThroughTime(x, node=140, nodetype='exclude')

# clade specific rates:
# here, for dolphins
rts <- getCladeRates(x, node = 140)
mean(rts$lambda)
quantile(rts$lambda, c(0.05, 0.95))

# get non-dolphin rates
nd <- getCladeRates(x, node=140, nodetype = 'exclude')
mean(nd$lambda)
quantile(nd$lambda, c(0.05, 0.95))

#### Model selection stuff
pprobs <- summary(x) # posterior probabilities of models

# Bayes factors:
data(mcmc.whales, prior.whales)
bmat <- computeBayesFactors(mcmc.whales, prior.whales, burnin=0.1)

#####################################################################
#####################################################################

# FIgure 1. 6 Credible set of shift configurations
# Figure 2. Rate-through-time analysis.
# Figure 3. Rate histograms

library(BAMMtools)
data(whales, events.whales, prior.whales)

 
pshifts <- getBranchShiftPriors(whales, prior.whales)

css <- credibleShiftSet(ed, threshold=pshifts)

quartz.options(height=8, width=4)
zz <- plot.bammdata(x, lwd=1.5, legend=F, pal='temperature')

cc <- zz$colordens


########## Color histogram for rates:

quartz.options(height=4, width=4)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0, 0.76), ylim=c(0, 13))

for (i in 2:nrow(cc)){
	xv <- cc$kde.x[c(i-1,i-1,i,i)]
	yv <- c(0, cc$kde.y[c(i,i)], 0)
	polygon(xv, yv, col=cc$coldens[i], border=cc$coldens[i])
}
axis(1, at=seq(-0.2, 0.8, by=0.2), cex.axis=1.2)
axis(2, at=seq(-3, 15, by=2), cex.axis=1.2, labels=F)
mtext(side=1, text='Speciation rate', line=2.7, cex=1.4)
mtext(side=2, text='Frequency', line=1.5, cex=1.4)

########## Color bar to interpret rates:

quartz.options(height=5, width=5)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0,4), ylim=c(0, 0.25))

for (i in 2:nrow(cc)){
	
	xv <- c(0,0,1,1);
	yv <- cc$kde.x[c(i-1,i-1,i,i)]
 
	polygon(xv, yv, col=cc$coldens[i], border=cc$coldens[i])
}


dev.off(); # close graphics device when done.

################ Figure 1

quartz.options(height=7, width=10)
plot.credibleshiftset(css, pal = 'temperature', plotmax=6, lwd=1.5)

dev.off(); # Close graphics device when finished...
 
 
############### Figure 2
 
# node 140 dolphins + orca
# node 141, core dolphins

rtt <- getRateThroughTimeMatrix(x)
rtt.dolphin <- getRateThroughTimeMatrix(x, node=140, nodetype = 'include')
rtt.other <- getRateThroughTimeMatrix(x, node=140, nodetype = 'exclude')


xlim <- c(max(branching.times(whales)),0)
 
quartz.options(height=7, width=10)
plot.new()
par(oma=c(1,1,1,1))
par(mfrow=c(2,3))
 
 
plotRateThroughTime(rtt, intervalCol='red', avgCol='red', ylim=c(0,1), xlim=xlim, axis.labels=T);
mtext(side=3, line=-1, "All whales", cex=1.1, font=3);

plotRateThroughTime(rtt.dolphin, intervalCol='blue', avgCol='blue', ylim=c(0, 1), xlim=xlim);
mtext(side=3, line=-1, "Dolphins", cex=1.1, font=3);

plotRateThroughTime(rtt.other, intervalCol='darkgreen', avgCol='darkgreen', ylim=c(0, 1), xlim=xlim);
mtext(side=3, line=-1, "Non-dolphins", cex=1.1, font=3);

plotRateThroughTime(rtt, intervals=c(0.05,0.95), intervalCol='gray70', avgCol='black', ylim=c(0,1), xlim=xlim, opacity=1);
mtext(side=3, line=-1, "All whales", cex=1.1, font=3);

plotRateThroughTime(rtt.dolphin, intervals=c(0.05,0.95), intervalCol='gray70', avgCol='black', ylim=c(0,1), xlim=xlim, opacity=1);
mtext(side=3, line=-1, "Dolphins", cex=1.1, font=3);


plotRateThroughTime(rtt.other, intervals=c(0.05,0.95), intervalCol='gray70', avgCol='black', ylim=c(0,1), xlim=xlim, opacity=1);
mtext(side=3, line=-1, "Non-dolphins", cex=1.1, font=3);


dev.off(); # Close graphics device when finished...
 

##################### Figure 3
################### Frequency distribution

# Figure 3 is a composite figure. First, make the frequency distribution, 
#	then the Bayes factor colormatrix.

# Here is part 1 - the frequency distribution.

wid <- 0.4
quartz.options(height=6, width=4)
plot.new()
par(oma=c(1,1,1,1))
par(mar=c(5,5,1,1))
plot.window(xlim=c(-1,4.5), ylim=c(0,0.8))

for (i in 1:nrow(pprobs)){
	k <- pprobs$shifts[i];
	xv <- c(k-wid, k-wid, k+wid, k+wid)
	yv <- c(0, pprobs$prob[c(i,i)], 0)
	polygon(xv, yv, col='gray50', border='black')
}
axis(1, at=c(-2,0:6))
axis(2, at=seq(-0.2, 0.8, by=0.2), las=1)
mtext(side=1, text='Number of shifts', line =3, cex=1.4)
mtext(side=2, text='Probability', line =3, cex=1.4)

dev.off(); # close device


# Here is part 2 - the Bayes factor color matrix:


library(colorspace);
library(gplots);

# Some useful functions

getColorKey <- function(x, ncols = 16, units=2){
	
	if (ncols %% 2 != 0){
		stop("expecting even number of colors");
	}
	
	index <- (1:(ncols-1)) - (ncols/2);
	xv <- index * log(units);
	
	colset <- rich.colors(n=ncols);
	colmat <- matrix('', nrow=nrow(x), ncol=ncol(x));	
		
	for (i in 1:nrow(x)){
		for (j in 1:ncol(x)){
			if (!is.na(x[i,j])){
				if (x[i,j] <= xv[1]){
					colmat[i,j] <- colset[1];
				}else if (x[i,j] > xv[length(xv)]){
					colmat[i,j] <- colset[length(colset)];
				}else{
					for (k in 1:(length(xv)-1)){
						if (x[i, j] > xv[k] & x[i,j] <= xv[k+1]){
							colmat[i,j] <- colset[k];
						}
					}
				}

			}
		}
	}	
	
	return(colmat);
}

getColorBar <- function(ncols = 16, units=2){
	
	if (ncols %% 2 != 0){
		stop("expecting even number of colors");
	}
	
	index <- (1:(ncols-1)) - (ncols/2);
	xv <- index * log(units);
 
	colset <- rich.colors(n=ncols);
 	cx <- c(xv, max(xv)+log(units));
 
	return(data.frame(cuts=cx, cols=colset, stringsAsFactors=F));
}


bfmat <- computeBayesFactors(mcmc.whales, prior.whales, burnin=0.1, modelset=0:4)
 
diag(bfmat) <- rep(1, length=5) 

quartz.options(height=6, width=6)
plot.new()
par(oma=c(1,1,1,1))
par(mar=c(5,5,1,1))
plot.window(xlim=c(0,5), ylim=c(0,5), asp=1)

nc <- 22
bb <- 1.5

cmat <- getColorKey(log(bfmat), ncols=nc, units=bb)

quartz.options(height=7, width=10, dpi=72);
ll <- c(rep(1, 9), rep(2,3));
lmat <- matrix(ll, nrow=3, byrow=F);

plot.new();

layout(lmat);
plot.new();
par(mar=c(6,6,1,1));

plot.window(xlim=c(0,5), ylim=c(0,5), asp=1);
 

for (i in 1:nrow(bfmat)){
#for (i in 1:5){	
	for (j in 1:nrow(bfmat)){
		if (!is.na(bfmat[i,j])){
			xval <- as.numeric(rownames(bfmat)[i]);
			yval <- as.numeric(colnames(bfmat)[j])
			xco <- c(xval, xval, xval+1, xval+1);
			yco <- c(yval, yval+1, yval+1, yval);
			polygon(x=xco, y=yco, lwd=0.8, col=cmat[i, j], border='black')
		}
	}
}

axis(1, at=seq(-1, 5, by=1), cex.axis=1.2);
axis(2, at=seq(-1, 5, by=1), las=1, cex.axis=1.2);
mtext(side=1, text="Macroevolutionary regimes, numerator", cex=1.4, line=3.5);
mtext(side=2, text="Macroevolutionary regimes, denominator", cex=1.4, line=3.5);

# Add the color bar
cb <- getColorBar(ncols=nc, units=bb);

plot.new();
par(mar=c(7,1,1,7));

plot.window(xlim=c(0, 2.5), ylim=c(-8,8));

for (i in 2:nrow(cb)){
	xv <- c(0,0,1,1);
	yv <- c(cb$cuts[c(i-1,i,i,i-1)]);
	polygon(xv, yv, col=cb$cols[i]);
	
}

axis(4, at=seq(-5, 5, by=2), las=1, cex.axis=2,pos=1.5);
mtext(side=4, at=0, text="(Log) Bayes factor", cex=1.5)


dev.off(); #close graphics device when done.













