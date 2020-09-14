#!/usr/bin/env Rscript
# Experiment name
name = "detailed_crowding"
library("doParallel")
library("Rcpp")
library("dqrng")
setwd("/Users/jdraghi/Dropbox/adaptive_preferences/")
newDir = paste("/Users/jdraghi/Dropbox/adaptive_preferences/new_results/", name, sep="") 
system2("mkdir", newDir)
timeStamp = date()
temp = file.copy("/Users/jdraghi/Dropbox/adaptive_preferences/revised_parallel_model.R", to = paste0(name, "_", timeStamp, ".R"))

# Replicates
reps = 200
# Population size
N = 10000
# Non-neutral mutation rates
UA = 0.0005
UB = 0.0005
# Preference mutation rate 
mu = 0.0002
# Selection coefficients
sA = 0.01
sB = 0.01
# Starting ranks
kA = 70
kB = 70
# Recombination rate per block
rs = 0
# Frequency of environment 2
ps = 0.3

stepsizes = c(0.1, 0.25, 0.5, 1)

# Chance of dying per patch examined
c = 0.25
# Fraction of mutations that are beneficial
b = 1/20
# Numbers of loci per class
L = 100
# Crowding factor
crowding = TRUE
# Number of generations
maxT = 5e5
# Sampling interval
interval = 100
treatments = length(stepsizes)
totalRuns = treatments * reps
seeds = sample.int(2e9, totalRuns, replace=FALSE)
write.table(paste("id", "N", "UA", "UB", "mu", "sA", "sB", "kA", "kB", "r", "p", "c", "b", "L", "crowding", sep=" "), paste0(newDir, "/", "treatments.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(cbind(1:treatments, N, UA, UB, mu, sA, sB, kA, kB, rs, ps, c, b, L, crowding, stepsizes), paste0(newDir, "/", "treatments.txt"), append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)

sourceCpp("/Users/jdraghi/Dropbox/adaptive_preferences/pickParents.cpp")

hub = makeForkCluster(4)
registerDoParallel(hub)

temp = foreach(rep=1:totalRuns, .packages=("dqrng"), .noexport= c("pickParents", "pickParentsSoftSel", "pickParentsCrowding", "pickParentsFreePref")) %dopar%
{
	set.seed(seeds[rep])
	dqRNGkind("pcg64")
	dqset.seed(seeds[rep])
	
	treat = (rep - 1) %% treatments + 1
	filename = paste0(newDir, "/", name, "_", treat, ".txt")
	detailsFile = paste0(newDir, "/", name, "_details_", treat, ".txt")
	r = rs
	p = ps
	U = c(UA, UB)
	s = c(sA, sB)
	k = c(kA, kB)
	fits = matrix(0, ncol=2, nrow=L+1)
	for(i in 1:2)
	{
		fits[,i] = (1 + s[i]) ^ (L:0)
		fits[,i] = fits[,i] / fits[1,i]
	}
	stepSize = stepsizes[treat]

	pop = matrix(0, nrow=N, ncol=3)	
	newPop = pop
	pop[,1] = 0
	pop[,2] = k[1]
	pop[,3] = k[2]
	meanPref = rep(0, maxT / interval + 1)
	meanSpec = rep(0, maxT / interval + 1)
	meanAdapt = matrix(0, ncol=2, nrow=maxT / interval + 1)
	meanComp = matrix(0, ncol=2, nrow=maxT / interval + 1)
	parents = matrix(0, ncol=2, nrow=N)
	brks = seq(from = -0.05, to = 1.05, by = 0.1)
	
	for(t in 0:maxT)
	{	
		if(t %% interval == 0)
		{
			meanPref[t/interval+1] = mean(pop[,1])
			meanSpec[t/interval+1] = mean(pop[,1] >= 0.9)
			meanAdapt[t/interval+1,] = L - colMeans(pop[,-1]) + 1
			meanComp[t/interval+1,] = cbind(mean(fits[pop[,2],1]), mean(fits[pop[,3],2]))
			x = hist(pop[,1], breaks=brks, plot=FALSE)
			write.table(t(c(rep, t, x$counts)), detailsFile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
		}
		if(crowding == FALSE)
		{
			 #newPop = pickParents(N, pop, fits, p, c, r, dqrunif(5*N))
			 newPop = pickParentsFreePref(N, pop, fits, p, c, r, dqrunif(5*N))
		}
		if(crowding == TRUE) newPop = pickParentsCrowding(N, pop, fits, p, c, r, dqrunif(5*N))
		
		for(i in 1:2)
		{
			mutations = rpois(1, N * U[i])
			if(mutations > 0)
			{
				ids = sample(1:N, mutations, replace=TRUE)
				for(j in 1:mutations)
				{
					pBene = (newPop[cbind(ids[j], i+1)] - 1) / L * b
					pDel = 1 - (newPop[cbind(ids[j], i+1)] - 1) / L
					pNeutral = 1 - pBene - pDel
					newPop[cbind(ids[j], i+1)] = newPop[cbind(ids[j], i+1)] + sample(c(-1, 0, 1), 1, prob=c(pBene, pNeutral, pDel))
				}
			}
		}
		mutations = rpois(1, N * mu)
		if(mutations > 0)
		{
			ids = sample(1:N, mutations, replace=TRUE)
			for(i in 1:mutations)
			{
				#newPop[ids[i],1] = runif(1)
				newPop[ids[i],1] = newPop[ids[i],1] + sample(c(-1,1), 1) * stepSize
				if(newPop[ids[i],1] > 1) newPop[ids[i],1] = 1
				if(newPop[ids[i],1] < 0) newPop[ids[i],1] = 0
			}
		}
		pop = newPop
	}
	write.table(paste0("# ", seeds[rep]), filename, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table(cbind(rep, seq(from = 0, to = maxT, by = interval), meanPref, meanSpec, meanAdapt, meanComp), filename, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	
}

stopCluster(hub)




	