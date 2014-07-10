##################################################################################
################################ decide directory ################################
##################################################################################
#setwd("/Volumes/sim_back/cassava5")
setwd("~/Desktop/cassava5")


##################################################################################
######################### call script and package ################################
##################################################################################
source("kin.blup4.4.R")
source("breedingSchemeForCassavaGS140613_3.R")
library(rrBLUP)



####################################################################################
############### make base population based on specied information ##################
####################################################################################

#------------------------ make base population ------------------------------------
createBreedingBase <- function(seed=round(runif(1, 0, 1e9)), nMarkers, nQTL, effPopSize=100, propDomi, interactionMean, varEffects, corrEffects){
	print("Creating base data for breeding")
	set.seed(seed)
	
	# Parameters of the species (cassava)
	nChr <- 18
	lengthChr <- 110 # in cM
	minMAF <- 0.01
	nLoci <- nMarkers + nQTL * (interactionMean + 1) * 3     # "* 3" is for buck-up.
	
	# Coalescent simulation
	piecesPerM <- 10000   ## pieces per Morgan ##
	nPiecesPerChr <- lengthChr / 100 * piecesPerM
	recBTpieces <- 1 / piecesPerM
	
	# coalescent simulation
	# You have to multiply effPopSize by two because "genome" works on haploids
	coalSim <- getCoalescentSim(effPopSize=2*effPopSize, nMrkOrMut=nLoci, nChr=nChr, nPiecesPerChr=nPiecesPerChr, recBTpieces=recBTpieces, minMAF=minMAF, seed=seed)
	markers <- coalSim$markers
	map <- coalSim$map
	
	# Scramble ancestral state
	ancestralState <- rbinom(nLoci, 1, 0.5)
	markers[,ancestralState == 1] <- 1 - markers[,ancestralState == 1]
		
	# creating map data
	mapData <- makeMap(map, nLoci, nMarkers, nQTL, propDomi, interactionMean, varEffects, corrEffects)
	
	# Assemble the list
	return(list(mapData=mapData, founderHaps=markers))
}



#######################################################################
######################## GS breeding simulation #######################
#######################################################################

#----------------------- input -----------------------------------
nSim <- 50
propDomi <- 0
interactionMean <- 0
addSeedling <- F
nTraits <- 1
c12 <- 0.2
c23 <- 0.5
c31 <- 0.8
corrEffects <- matrix(c(1, c12, c31, c12, 1, c23, c31, c23, 1), 3, 3)[1:nTraits, 1:nTraits]

#------------------------- defined parameter ---------------------------
nPop0 <- 200
nPop1 <- nPop2 <- nPop3 <- nPop4 <- 600
nSel0 <- nSel1 <- nSel2 <- nSel3 <- 40
varSeedlings <- 36
varCET <- 16
varPYT <- 9
varYT <- 1
nMarkers <- 2500
nQTL <- 100
selectionWeights <- rep(1, nTraits)
varEffects <- rep(1, nTraits)
geneticVariance <- rep(1, nTraits)

#------------------------ simulation ----------------------------
ave <- matrix(NA, nSim * nTraits, 5)
top <- matrix(NA, nSim * nTraits, 5)
aveAll <- matrix(0, nSim, 5)
topAll <- matrix(0, nSim, 5)

if(addSeedling){
	pdf(paste("Fig_", nSim, "timesTrials_", nTraits, "traits_eachReplication.pdf", sep = ""))
}else{
	pdf(paste("Fig_", nSim, "timesTrials_", nTraits, "traits_eachReplication_withoutSeedlingsInTP.pdf", sep = ""))
}


for(sim in 1:nSim){
	print(paste("-----------------------------", sim, "-----------------------------------"))
	print(Sys.time())
	
	## create initial population ##
	baseData <- createBreedingBase(nMarkers = nMarkers, nQTL = nQTL, propDomi = propDomi, interactionMean = interactionMean, varEffects = varEffects, corrEffects = corrEffects)
	initData <- createInitPop(baseData = baseData, nPopInit = nPop0)
	initData <- adjustMap(initialPopulationData = initData, geneticVariance = geneticVariance, nTraits = nTraits)
	mapData <- initData$mapData
	genoC0 <- initData$geno
	print("Initial population")
	print(Sys.time())
		
	## C0 ##
	# training population
	for(trait in 1:nTraits){
		g <- calcGenotypicValue(geno = genoC0, mapData = mapData, trait = trait)
		y <- calcPhenotypicValue(gv = g, errorVar = varYT)
		phenoData <- data.frame(GID = 1:nPop0, g = g, y = y$pv, var = rep(y$var, nPop0))
		if(trait == 1){
			phenoDataAll <- phenoData
		}else{
			phenoDataAll <- cbind(phenoDataAll, phenoData)
		}
		ave[(trait - 1) * nSim + sim, 1] <- mean(g)
		top[(trait - 1) * nSim + sim, 1] <- max(g)
		if(trait == 1){
			gWeights <- g * selectionWeights[1]
		}else{
			gWeights <- gWeights + g * selectionWeights[trait]
		}
	} # trait
	aveAll[sim, 1] <- mean(gWeights) 
	topAll[sim, 1] <- max(gWeights)
	geno <- genoC0
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = geno, phenoData = phenoDataAll, nSelect = nSel0, selectPopGID = 1:nPop0, weights = selectionWeights)
	genoParents <- geno[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	# make next generation
	genoC1 <- randomMate(popSize = nPop1, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	print("C0")
	print(Sys.time())
		
	## C1 ##
	# phenotypic values of training population
	for(trait in 1:nTraits){
		g.add <- calcGenotypicValue(geno = genoC1, mapData = mapData, trait = trait)
	    if(addSeedling){
		    y.add <- calcPhenotypicValue(gv = g.add, errorVar = varSeedlings)
	    }else{
		    y.add <- list(pv = rep(NA, nPop1), var = varSeedlings)
	    }
		phenoData.add <- data.frame(GID = nPop0 + 1:nPop1, g = g.add, y = y.add$pv, var = rep(y.add$var, nPop1))
		if(trait == 1){
			phenoDataAll.add <- phenoData.add
		}else{
			phenoDataAll.add <- cbind(phenoDataAll.add, phenoData.add)
		}
		ave[(trait - 1) * nSim + sim, 2] <- mean(g.add)
		top[(trait - 1) * nSim + sim, 2] <- max(g.add)
		if(trait == 1){
			gWeights <- g.add * selectionWeights[1]
		}else{
			gWeights <- gWeights + g.add * selectionWeights[trait]
		}
	} # trait
	aveAll[sim, 2] <- mean(gWeights) 
	topAll[sim, 2] <- max(gWeights)
	phenoDataAll <- rbind(phenoDataAll, phenoDataAll.add)
	geno <- rbind(genoC0, genoC1)
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = geno, phenoData = phenoDataAll, nSelect = nSel1, selectPopGID = nPop0 + 1:nPop1, weights = selectionWeights)
	genoParents <- geno[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	genoRemC1 <- geno[-sort(c(parentsID * 2 - 1, parentsID * 2)), ][-(1:(nPop0 * 2)), ]
	# make next generation
	genoC2 <- randomMate(popSize = nPop2, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	print("C1")
	print(Sys.time())
	
	## C2 ##
	# phenotypic values of training population
	for(trait in 1:nTraits){
		g.add1 <- calcGenotypicValue(geno = genoRemC1, mapData = mapData, trait = trait)
		y.add1 <- calcPhenotypicValue(gv = g.add1, errorVar = varCET)
		phenoData.add1 <- data.frame(GID = (nPop0 + 1:nPop1)[-(parentsID - nPop0)], g = g.add1, y = y.add1$pv, var = rep(y.add1$var, nPop1 - nSel1))
		if(trait == 1){
			phenoDataAll.add1 <- phenoData.add1
		}else{
			phenoDataAll.add1 <- cbind(phenoDataAll.add1, phenoData.add1)
		}
		g.add2 <- calcGenotypicValue(geno = genoC2, mapData = mapData, trait = trait)
	    if(addSeedling){
		    y.add2 <- calcPhenotypicValue(gv = g.add2, errorVar = varSeedlings)
	    }else{
		    y.add2 <- list(pv = rep(NA, nPop2), var = varSeedlings)
	    }
		phenoData.add2 <- data.frame(GID = nPop0 + nPop1 + 1:nPop2, g = g.add2, y = y.add2$pv, var = rep(y.add2$var, nPop2))
		if(trait == 1){
			phenoDataAll.add2 <- phenoData.add2
		}else{
			phenoDataAll.add2 <- cbind(phenoDataAll.add2, phenoData.add2)
		}
		ave[(trait - 1) * nSim + sim, 3] <- mean(g.add2)
		top[(trait - 1) * nSim + sim, 3] <- max(g.add2)
		if(trait == 1){
			gWeights <- g.add2 * selectionWeights[1]
		}else{
			gWeights <- gWeights + g.add2 * selectionWeights[trait]
		}
	} # trait
	aveAll[sim, 3] <- mean(gWeights) 
	topAll[sim, 3] <- max(gWeights)
	phenoDataAll <- rbind(phenoDataAll, phenoDataAll.add1, phenoDataAll.add2)
	geno <- rbind(genoC0, genoC1, genoC2)
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = geno, phenoData = phenoDataAll, nSelect = nSel2, selectPopGID = nPop0 + nPop1 + 1:nPop2, weights = selectionWeights, reduce = T)
	genoParents <- geno[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	genoRemC2 <- geno[-sort(c(parentsID * 2 - 1, parentsID * 2)), ][-(1:((nPop0 + nPop1) * 2)), ]
	# make next generation
	genoC3 <- randomMate(popSize = nPop3, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	print("C2")
	print(Sys.time())
	
	## C3 ##
	# phenotypic values of training population
	for(trait in 1:nTraits){
		g.add1 <- calcGenotypicValue(geno = genoC1, mapData = mapData, trait = trait)
		y.add1 <- calcPhenotypicValue(gv = g.add1, errorVar = varPYT)
		phenoData.add1 <- data.frame(GID = nPop0 + 1:nPop1, g = g.add1, y = y.add1$pv, var = rep(y.add1$var, nPop1))
		if(trait == 1){
			phenoDataAll.add1 <- phenoData.add1
		}else{
			phenoDataAll.add1 <- cbind(phenoDataAll.add1, phenoData.add1)
		}
		g.add2 <- calcGenotypicValue(geno = genoRemC2, mapData = mapData, trait = trait)
		y.add2 <- calcPhenotypicValue(gv = g.add2, errorVar = varCET)
		phenoData.add2 <- data.frame(GID = (nPop0 + nPop1 + 1:nPop2)[-(parentsID - nPop0 - nPop1)], g = g.add2, y = y.add2$pv, var = rep(y.add2$var, nPop2 - nSel2))
		if(trait == 1){
			phenoDataAll.add2 <- phenoData.add2
		}else{
			phenoDataAll.add2 <- cbind(phenoDataAll.add2, phenoData.add2)
		}
		g.add3 <- calcGenotypicValue(geno = genoC3, mapData = mapData, trait = trait)
	    if(addSeedling){
		    y.add3 <- calcPhenotypicValue(gv = g.add3, errorVar = varSeedlings)
	    }else{
		    y.add3 <- list(pv = rep(NA, nPop3), var = varSeedlings)
	    }
		phenoData.add3 <- data.frame(GID = nPop0 + nPop1 + nPop2 + 1:nPop3, g = g.add3, y = y.add3$pv, var = rep(y.add3$var, nPop3))
		if(trait == 1){
			phenoDataAll.add3 <- phenoData.add3
		}else{
			phenoDataAll.add3 <- cbind(phenoDataAll.add3, phenoData.add3)
		}
		ave[(trait - 1) * nSim + sim, 4] <- mean(g.add3)
		top[(trait - 1) * nSim + sim, 4] <- max(g.add3)
		if(trait == 1){
			gWeights <- g.add3 * selectionWeights[1]
		}else{
			gWeights <- gWeights + g.add3 * selectionWeights[trait]
		}
	} # trait
	aveAll[sim, 4] <- mean(gWeights) 
	topAll[sim, 4] <- max(gWeights)
	phenoDataAll <- rbind(phenoDataAll, phenoDataAll.add1, phenoDataAll.add2, phenoDataAll.add3)
	geno <- rbind(genoC0, genoC1, genoC2, genoC3)
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = geno, phenoData = phenoDataAll, nSelect = nSel3, selectPopGID = nPop0 + nPop1 + nPop2 + 1:nPop3, weights = selectionWeights, reduce = T)
	genoParents <- geno[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	genoRemC3 <- geno[-sort(c(parentsID * 2 - 1, parentsID * 2)), ][-(1:((nPop0 + nPop1 + nPop2) * 2)), ]
	# make next generation
	genoC4 <- randomMate(popSize = nPop4, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	print("C3")
	print(Sys.time())
	
	## C4 ##
	# genotypic values of training population
	for(trait in 1:nTraits){
		g4 <- calcGenotypicValue(geno = genoC4, mapData = mapData, trait = trait)
		ave[(trait - 1) * nSim + sim, 5] <- mean(g4)
		top[(trait - 1) * nSim + sim, 5] <- max(g4)
		if(trait == 1){
			gWeights <- g4 * selectionWeights[1]
		}else{
			gWeights <- gWeights + g4 * selectionWeights[trait]
		}
	} # trait
	aveAll[sim, 5] <- mean(gWeights) 
	topAll[sim, 5] <- max(gWeights)
	
	## evaluate ##
	plot(aveAll[sim, ], type = "b", main = sim, ylim = c(min(aveAll[sim, ]), max(topAll[sim, ])))
	points(topAll[sim, ], type = "b", lty = 2)
	for(trait in 1:nTraits){
		points(ave[(trait - 1) * nSim + sim, ], type = "b", col = trait + 1)
		points(top[(trait - 1) * nSim + sim, ], type = "b", lty = 2, col = trait + 1)
	} # trait
	print("Evaluation")
	print(Sys.time())
	
} # sim

dev.off()


## draw pictures ##
if(addSeedling){
	pdf(paste("Fig_", nSim, "timesTrials_", nTraits, "traits.pdf", sep = ""))
}else{
	pdf(paste("Fig_", nSim, "timesTrials_", nTraits, "traits_withoutSeedlingsInTP.pdf", sep = ""))
}
if(nSim == 1){
  aveAllMu <- as.vector(aveAll)
  topAllMu <- as.vector(topAll)
}else{
  aveAllMu <- apply(aveAll, 2, mean)
  topAllMu <- apply(topAll, 2, mean)
}
a <- NULL
t <- NULL
for(trait in 1:nTraits){
  if(nSim == 1 | nSim == 0){
    a <- c(a, ave[trait, ])
    t <- c(t, top[trait, ])
  }else{
    a <- c(a, apply(ave[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    t <- c(t, apply(top[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
  }
} # trait
plot(aveAllMu, type = "b", lwd = 2.5, ylim = c(min(c(a, aveAllMu)), max(c(t, topAllMu))), main = "", xlab = "", ylab = "Genotypic value", axes = F)
points(topAllMu, type = "b", lwd = 2.5, lty = 2)
for(trait in 1:nTraits){
  if(nSim == 1 | nSim == 0){
    aveMu <- ave[trait, ]
    topMu <- top[trait, ]
  }else{
    aveMu <- apply(ave[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean)
    topMu <- apply(top[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean)
  }
  points(aveMu, type = "b", col = trait + 1)
  points(topMu, type = "b", col = trait + 1, lty = 2)
}
legend(x = 1, y = max(c(t, topAllMu)), legend = c("Mean", "Top"), lty = c(1, 2), pch = c(1, 1))
axis(1, labels = F)
mtext(paste("C", 0:4, sep = ""), side = 1, line = 1, at = 1:5)
axis(2)
box()

dev.off()

print("---------------------------------------- Summary ---------------------------------------------")
print(Sys.time())

