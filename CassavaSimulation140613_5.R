##################################################################################
################################ decide directory ################################
##################################################################################
#setwd("/Volumes/sim_back/cassava7")
setwd("~/Desktop/cassava7")


##################################################################################
######################### call script and package ################################
##################################################################################
source("kin.blup4.4.R")
source("breedingSchemeForCassavaGS140613_5.R")
library(rrBLUP)
library(ggplot2)
library(multicore)



####################################################################################
############### make base population and initial population ########################
####################################################################################

#------------------------ make base population of cassava ------------------------------------
createBreedingBase <- function(seed=round(runif(1, 0, 1e9)), nMarkers, nQTL, effPopSize=100, propDomi, interactionMean, varEffects, corrEffects){
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

#------------------------ make initial population of cassava ------------------------------------
createInitPopCassava <- function(nMarkers, nQTL, propDomi, interactionMean, varEffects, corrEffects, nPopInit, geneticVariance, nTraits){
	baseData <- createBreedingBase(nMarkers = nMarkers, nQTL = nQTL, propDomi = propDomi, interactionMean = interactionMean, varEffects = varEffects, corrEffects = corrEffects)
	initData <- createInitPop(baseData = baseData, nPopInit = nPopInit)
	initData <- adjustMap(initialPopulationData = initData, geneticVariance = geneticVariance, nTraits = nTraits)
	
	return(initData)
}



####################################################################################
############### make 4 selection cycles in cassava breeding ########################
####################################################################################
breedingCassava <- function(initData, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate, addSeedling){
	## read initial population and map data ##
	mapData <- initData$mapData
	genoC0 <- initData$geno
	
	## prepare box for results ##
	ave <- matrix(NA, nTraits, 5)
	top <- matrix(NA, nTraits, 5)
	gvar <- matrix(NA, nTraits, 5)
	aveAll <- rep(0, 5)
	topAll <- rep(0, 5)
	
	## C0 ##
	# genotype data
	nPop0miss <- (nPop0 * errorRate / 100) %/% 1
	genoC0true <- genoC0[1:(nPop0 * 2), ]
	if(errorRate > 0){
		genoC0false <- genoC0[c(1:((nPop0 - nPop0miss) * 2), (nPop0 * 2 + 1):((nPop0 + nPop0miss) * 2)), ]
	}else{
		genoC0false <- genoC0[1:(nPop0 * 2), ]
	}
	# training population
	for(trait in 1:nTraits){
		g <- calcGenotypicValue(geno = genoC0true, mapData = mapData, trait = trait)
		y <- calcPhenotypicValue(gv = g, errorVar = varYT)
		phenoData <- data.frame(GID = 1:nPop0, g = g, y = y$pv, var = rep(y$var, nPop0))
		if(trait == 1){
			phenoDataAll <- phenoData
		}else{
			phenoDataAll <- cbind(phenoDataAll, phenoData)
		}
		ave[trait, 1] <- mean(g)
		top[trait, 1] <- max(g)
		gvar[trait, 1] <- var(g)
		if(trait == 1){
			gWeights <- g * selectionWeights[1]
		}else{
			gWeights <- gWeights + g * selectionWeights[trait]
		}
	} # trait
	aveAll[1] <- mean(gWeights) 
	topAll[1] <- max(gWeights)
	genoT <- genoC0true
	genoF <- genoC0false
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = genoF, phenoData = phenoDataAll, nSelect = nSel0, selectPopGID = 1:nPop0, weights = selectionWeights)
	genoParents <- genoT[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	# make next generation
	genoC1 <- randomMate(popSize = nPop1 + (nPop1 * errorRate / 100) %/% 1, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	
	## C1 ##
	# genotype data
	nPop1miss <- (nPop1 * errorRate / 100) %/% 1
	genoC1true <- genoC1[1:(nPop1 * 2), ]
	if(errorRate > 0){
		genoC1false <- genoC1[c(1:((nPop1 - nPop1miss) * 2), (nPop1 * 2 + 1):((nPop1 + nPop1miss) * 2)), ]
	}else{
		genoC1false <- genoC1[1:(nPop1 * 2), ]
	}
	# phenotypic values of training population
	for(trait in 1:nTraits){
		g.add <- calcGenotypicValue(geno = genoC1true, mapData = mapData, trait = trait)
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
		ave[trait, 2] <- mean(g.add)
		top[trait, 2] <- max(g.add)
		gvar[trait, 2] <- var(g.add)
		if(trait == 1){
			gWeights <- g.add * selectionWeights[1]
		}else{
			gWeights <- gWeights + g.add * selectionWeights[trait]
		}
	} # trait
	aveAll[2] <- mean(gWeights) 
	topAll[2] <- max(gWeights)
	phenoDataAll <- rbind(phenoDataAll, phenoDataAll.add)
	genoT <- rbind(genoC0true, genoC1true)
	genoF <- rbind(genoC0false, genoC1false)
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = genoF, phenoData = phenoDataAll, nSelect = nSel1, selectPopGID = nPop0 + 1:nPop1, weights = selectionWeights)
	genoParents <- genoT[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	genoRemC1 <- genoT[-sort(c(parentsID * 2 - 1, parentsID * 2)), ][-(1:(nPop0 * 2)), ]
	# make next generation
	genoC2 <- randomMate(popSize = nPop2 + (nPop2 * errorRate / 100) %/% 1, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	
	## C2 ##
	# genotype data
	nPop2miss <- (nPop2 * errorRate / 100) %/% 1
	genoC2true <- genoC2[1:(nPop2 * 2), ]
	if(errorRate > 0){
		genoC2false <- genoC2[c(1:((nPop2 - nPop2miss) * 2), (nPop2 * 2 + 1):((nPop2 + nPop2miss) * 2)), ]
	}else{
		genoC2false <- genoC2[1:(nPop2 * 2), ]
	}
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
		g.add2 <- calcGenotypicValue(geno = genoC2true, mapData = mapData, trait = trait)
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
		ave[trait, 3] <- mean(g.add2)
		top[trait, 3] <- max(g.add2)
		gvar[trait, 3] <- var(g.add2)
		if(trait == 1){
			gWeights <- g.add2 * selectionWeights[1]
		}else{
			gWeights <- gWeights + g.add2 * selectionWeights[trait]
		}
	} # trait
	aveAll[3] <- mean(gWeights) 
	topAll[3] <- max(gWeights)
	phenoDataAll <- rbind(phenoDataAll, phenoDataAll.add1, phenoDataAll.add2)
	genoT <- rbind(genoC0true, genoC1true, genoC2true)
	genoF <- rbind(genoC0false, genoC1false, genoC2false)
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = genoF, phenoData = phenoDataAll, nSelect = nSel2, selectPopGID = nPop0 + nPop1 + 1:nPop2, weights = selectionWeights, reduce = T)
	genoParents <- genoT[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	genoRemC2 <- genoT[-sort(c(parentsID * 2 - 1, parentsID * 2)), ][-(1:((nPop0 + nPop1) * 2)), ]
	# make next generation
	genoC3 <- randomMate(popSize = nPop3 + (nPop3 * errorRate / 100) %/% 1, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	
	## C3 ##
	# genotype data
	nPop3miss <- (nPop3 * errorRate / 100) %/% 1
	genoC3true <- genoC3[1:(nPop3 * 2), ]
	if(errorRate > 0){
		genoC3false <- genoC3[c(1:((nPop3 - nPop3miss) * 2), (nPop3 * 2 + 1):((nPop3 + nPop3miss) * 2)), ]
	}else{
		genoC3false <- genoC3[1:(nPop3 * 2), ]
	}
	# phenotypic values of training population
	for(trait in 1:nTraits){
		g.add1 <- calcGenotypicValue(geno = genoC1true, mapData = mapData, trait = trait)
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
		g.add3 <- calcGenotypicValue(geno = genoC3true, mapData = mapData, trait = trait)
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
		ave[trait, 4] <- mean(g.add3)
		top[trait, 4] <- max(g.add3)
		gvar[trait, 4] <- var(g.add3)
		if(trait == 1){
			gWeights <- g.add3 * selectionWeights[1]
		}else{
			gWeights <- gWeights + g.add3 * selectionWeights[trait]
		}
	} # trait
	aveAll[4] <- mean(gWeights) 
	topAll[4] <- max(gWeights)
	phenoDataAll <- rbind(phenoDataAll, phenoDataAll.add1, phenoDataAll.add2, phenoDataAll.add3)
	genoT <- rbind(genoC0true, genoC1true, genoC2true, genoC3true)
	genoF <- rbind(genoC0false, genoC1false, genoC2false, genoC3false)
	# selection
	parentsID <- genomicSelection(mapData = mapData, geno = genoF, phenoData = phenoDataAll, nSelect = nSel3, selectPopGID = nPop0 + nPop1 + nPop2 + 1:nPop3, weights = selectionWeights, reduce = T)
	genoParents <- genoT[sort(c(parentsID * 2 - 1, parentsID * 2)), ]
	genoRemC3 <- genoT[-sort(c(parentsID * 2 - 1, parentsID * 2)), ][-(1:((nPop0 + nPop1 + nPop2) * 2)), ]
	# make next generation
	genoC4 <- randomMate(popSize = nPop4 + (nPop4 * errorRate / 100) %/% 1, geno = genoParents, pos = mapData$map$Pos)     # equal contribution is not assumed
	
	## C4 ##
	# genotype data
	genoC4true <- genoC4[1:(nPop4 * 2), ]
	# genotypic values of training population
	for(trait in 1:nTraits){
		g4 <- calcGenotypicValue(geno = genoC4true, mapData = mapData, trait = trait)
		ave[trait, 5] <- mean(g4)
		top[trait, 5] <- max(g4)
		gvar[trait, 5] <- var(g4)
		if(trait == 1){
			gWeights <- g4 * selectionWeights[1]
		}else{
			gWeights <- gWeights + g4 * selectionWeights[trait]
		}
	} # trait
	aveAll[5] <- mean(gWeights) 
	topAll[5] <- max(gWeights)
	
	## results of breeding ##
	return(list(ave = ave, top = top, gvar = gvar, aveAll = aveAll, topAll = topAll))
}




#######################################################################
######################## GS breeding simulation #######################
#######################################################################

#----------------------- input -----------------------------------
nSim <- 50
propDomi <- 0
interactionMean <- 0
nTraits <- 3
c12 <- 0
c23 <- 0
c31 <- 0
corrEffects <- matrix(c(1, c12, c31, c12, 1, c23, c31, c23, 1), 3, 3)[1:nTraits, 1:nTraits]
errorRate <- 20
nCore <- 10

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
## make initial population for breeding ##
initList <- list()
for(sim in 1:nSim){
	print(paste("Initial Population ------ ", sim, ": ", Sys.time(), sep = ""))
	initList[[sim]] <- createInitPopCassava(nMarkers = nMarkers, nQTL = nQTL, propDomi = propDomi, interactionMean = interactionMean, varEffects = varEffects, corrEffects = corrEffects, nPopInit = nPop0 + ((nPop0 * errorRate / 100) %/% 1), geneticVariance = geneticVariance, nTraits = nTraits)
} # sim

## breeding through 4 selection cycles ##
if(nCore == 0){
	print(paste("BREEDING ----------- ", Sys.time(), sep = ""))
	resSeedlings <- lapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate, addSeedling = T)
	print(paste("with Seedlings ------------ ", Sys.time(), sep = ""))
	resNoSeedlings <- lapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate, addSeedling = F)
	print(paste("without Seedlings ----------- ", Sys.time(), sep = ""))
	if(errorRate > 0){
		resSeedlings0 <- lapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate = 0, addSeedling = T)
		print(paste("with Seedlings 0 ------------ ", Sys.time(), sep = ""))
		resNoSeedlings0 <- lapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate = 0, addSeedling = F)
		print(paste("without Seedlings 0 ----------- ", Sys.time(), sep = ""))
	}
}else{
	print(paste("BREEDING ----------- ", Sys.time(), sep = ""))
	resSeedlings <- mclapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate, addSeedling = T, mc.cores = nCore)
	print(paste("with Seedlings ------------ ", Sys.time(), sep = ""))
	resNoSeedlings <- mclapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate, addSeedling = F, mc.cores = nCore)
	print(paste("without Seedlings ----------- ", Sys.time(), sep = ""))
	if(errorRate > 0){
		resSeedlings0 <- mclapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate = 0, addSeedling = T, mc.cores = nCore)
		print(paste("with Seedlings 0 ------------ ", Sys.time(), sep = ""))
		resNoSeedlings0 <- mclapply(initList, breedingCassava, nTraits, nPop0, nPop1, nPop2, nPop3, nPop4, nSel0, nSel1, nSel2, nSel3, varSeedlings, varCET, varPYT, varYT, selectionWeights, errorRate = 0, addSeedling = F, mc.cores = nCore)
		print(paste("without Seedlings 0 ----------- ", Sys.time(), sep = ""))
	}
}

## bind results from simulation replicates ##
aveS <- NULL
topS <- NULL
gvarS <- NULL
aveAllS <- NULL
topAllS <- NULL
aveN <- NULL
topN <- NULL
gvarN <- NULL
aveAllN <- NULL
topAllN <- NULL
if(errorRate > 0){
	aveS0 <- NULL
    topS0 <- NULL
    gvarS0 <- NULL
    aveAllS0 <- NULL
    topAllS0 <- NULL
    aveN0 <- NULL
    topN0 <- NULL
    gvarN0 <- NULL
    aveAllN0 <- NULL
    topAllN0 <- NULL
}
for(sim in 1:nSim){
	aveAllS <- rbind(aveAllS, resSeedlings[[sim]]$aveAll)
	topAllS <- rbind(topAllS, resSeedlings[[sim]]$topAll)
	aveAllN <- rbind(aveAllN, resNoSeedlings[[sim]]$aveAll)
	topAllN <- rbind(topAllN, resNoSeedlings[[sim]]$topAll)
	if(errorRate > 0){
	    aveAllS0 <- rbind(aveAllS0, resSeedlings0[[sim]]$aveAll)
	    topAllS0 <- rbind(topAllS0, resSeedlings0[[sim]]$topAll)
	    aveAllN0 <- rbind(aveAllN0, resNoSeedlings0[[sim]]$aveAll)
	    topAllN0 <- rbind(topAllN0, resNoSeedlings0[[sim]]$topAll)
	}
} # sim
for(trait in 1:nTraits){
	for(sim in 1:nSim){
		aveS <- rbind(aveS, resSeedlings[[sim]]$ave[trait, ])
		topS <- rbind(topS, resSeedlings[[sim]]$top[trait, ])
		gvarS <- rbind(gvarS, resSeedlings[[sim]]$gvar[trait, ])
		aveN <- rbind(aveN, resNoSeedlings[[sim]]$ave[trait, ])
		topN <- rbind(topN, resNoSeedlings[[sim]]$top[trait, ])
		gvarN <- rbind(gvarN, resNoSeedlings[[sim]]$gvar[trait, ])
		if(errorRate > 0){
		    aveS0 <- rbind(aveS0, resSeedlings0[[sim]]$ave[trait, ])
		    topS0 <- rbind(topS0, resSeedlings0[[sim]]$top[trait, ])
		    gvarS0 <- rbind(gvarS0, resSeedlings0[[sim]]$gvar[trait, ])
		    aveN0 <- rbind(aveN0, resNoSeedlings0[[sim]]$ave[trait, ])
		    topN0 <- rbind(topN0, resNoSeedlings0[[sim]]$top[trait, ])
		    gvarN0 <- rbind(gvarN0, resNoSeedlings0[[sim]]$gvar[trait, ])
		}
	} # sim
} # trait
print(paste("Bind Results ----------- ", Sys.time(), sep = ""))

## draw pictures ##
# genotyphic values
if(nSim == 1){
  aveAllMuS <- as.vector(aveAllS)
  topAllMuS <- as.vector(topAllS)
  aveAllMuN <- as.vector(aveAllN)
  topAllMuN <- as.vector(topAllN)
}else{
  aveAllMuS <- apply(aveAllS, 2, mean)
  topAllMuS <- apply(topAllS, 2, mean)
  aveAllMuN <- apply(aveAllN, 2, mean)
  topAllMuN <- apply(topAllN, 2, mean)
}
aveMuS <- NULL
topMuS <- NULL
aveMuN <- NULL
topMuN <- NULL
for(trait in 1:nTraits){
  if(nSim == 1){
    aveMuS <- c(aveMuS, aveS[trait, ])
    topMuS <- c(topMuS, topS[trait, ])
    aveMuN <- c(aveMuN, aveN[trait, ])
    topMuN <- c(topMuN, topN[trait, ])
  }else{
    aveMuS <- c(aveMuS, apply(aveS[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    topMuS <- c(topMuS, apply(topS[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    aveMuN <- c(aveMuN, apply(aveN[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    topMuN <- c(topMuN, apply(topN[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
  }
}
S <- c(aveAllMuS, aveMuS, topAllMuS, topMuS)
N <- c(aveAllMuN, aveMuN, topAllMuN, topMuN)
cycle <- rep(1:5, 2 + nTraits * 2)
trait <- rep(sort(rep(1:(nTraits + 1), 5)), 2)
line <- sort(rep(1:2, length(S) / 2))
value <- c(S, N)
cycle <- rep(cycle, 2)
trait <- rep(trait, 2)
linetype <- rep(line, 2)
shape <- sort(rep(1:2, length(value) / 2))
group <- sort(rep(1:((nTraits + 1) * 2 * 2), 5))
data <- data.frame(value = value, cycle = cycle, trait = trait, linetype = linetype, shape = shape, group = group, error = rep(1, length(value)))
if(errorRate > 0){
  if(nSim == 1){
    aveAllMuS <- as.vector(aveAllS0)
    topAllMuS <- as.vector(topAllS0)
    aveAllMuN <- as.vector(aveAllN0)
    topAllMuN <- as.vector(topAllN0)
  }else{
    aveAllMuS <- apply(aveAllS0, 2, mean)
    topAllMuS <- apply(topAllS0, 2, mean)
    aveAllMuN <- apply(aveAllN0, 2, mean)
    topAllMuN <- apply(topAllN0, 2, mean)
  }
  aveMuS <- NULL
  topMuS <- NULL
  aveMuN <- NULL
  topMuN <- NULL
  for(trait in 1:nTraits){
    if(nSim == 1){
      aveMuS <- c(aveMuS, aveS0[trait, ])
      topMuS <- c(topMuS, topS0[trait, ])
      aveMuN <- c(aveMuN, aveN0[trait, ])
      topMuN <- c(topMuN, topN0[trait, ])
    }else{
      aveMuS <- c(aveMuS, apply(aveS0[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
      topMuS <- c(topMuS, apply(topS0[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
      aveMuN <- c(aveMuN, apply(aveN0[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
      topMuN <- c(topMuN, apply(topN0[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    }
  }
  S <- c(aveAllMuS, aveMuS, topAllMuS, topMuS)
  N <- c(aveAllMuN, aveMuN, topAllMuN, topMuN)
  cycle <- rep(1:5, 2 + nTraits * 2)
  trait <- rep(sort(rep(1:(nTraits + 1), 5)), 2)
  line <- sort(rep(1:2, length(S) / 2))
  value <- c(S, N)
  cycle <- rep(cycle, 2)
  trait <- rep(trait, 2)
  linetype <- rep(line, 2)
  shape <- sort(rep(1:2, length(value) / 2))
  group <- sort(rep(1:((nTraits + 1) * 2 * 2), 5)) + (nTraits + 1) * 2 * 2
  data0 <- data.frame(value = value, cycle = cycle, trait = trait, linetype = linetype, shape = shape, group = group, error = rep(2, length(value)))
  data <- rbind(data, data0)
}
p1 <- ggplot(data = data, aes(x = cycle, y = value)) + geom_point(aes(colour = factor(trait), shape = factor(shape), size = factor(error))) + geom_line(aes(colour = factor(trait), linetype = factor(linetype), group = factor(group)))
p1 <- p1 + ylim(min(data$value), max(data$value))
p1 <- p1 + labs(title = "Genotypic value", x = "Cycle", y = "Genotypic value")
p1 <- p1 + scale_shape_discrete(name = "Seedlings in TP", breaks = factor(1:2), labels = c("YES", "NO"))
p1 <- p1 + scale_colour_discrete(name = "Trait", breaks = factor(1:(nTraits + 1)), labels = c("All", paste("Trait", 1:nTraits)))
p1 <- p1 + scale_linetype_discrete(name = "Average / Top", breaks = factor(c(1, 2)), labels = c("Average", "Top"))
if(errorRate > 0){
	p1 <- p1 + scale_size_manual(name = "Error", values = c(2, 3), labels = c(paste(errorRate, "%", sep = ""), "0%"))
}
p1 <- p1 + scale_x_continuous(breaks = 1:5, labels = paste("C", 0:4, sep = ""))
# genetic variance
gvarMuS <- NULL
gvarMuN <- NULL
for(trait in 1:nTraits){
  if(nSim == 1){
    gvarMuS <- c(gvarMuS, gvarS[trait, ])
    gvarMuN <- c(gvarMuN, gvarN[trait, ])
  }else{
    gvarMuS <- c(gvarMuS, apply(gvarS[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    gvarMuN <- c(gvarMuN, apply(gvarN[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
  }
} # trait
data <- data.frame(gvar = c(gvarMuS, gvarMuN), cycle = rep(rep(1:5, nTraits), 2), trait = rep(sort(rep(2:(nTraits + 1), 5)), 2), shape = sort(rep(1:2, length(gvarMuS))), group = sort(rep(1:(nTraits * 2), 5)), error = rep(1, length(c(gvarMuS, gvarMuN))))
if(errorRate > 0){
  gvarMuS <- NULL
  gvarMuN <- NULL
  for(trait in 1:nTraits){
    if(nSim == 1){
      gvarMuS <- c(gvarMuS, gvarS0[trait, ])
      gvarMuN <- c(gvarMuN, gvarN0[trait, ])
    }else{
      gvarMuS <- c(gvarMuS, apply(gvarS0[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
      gvarMuN <- c(gvarMuN, apply(gvarN0[((trait - 1) * nSim + 1):(trait * nSim), ], 2, mean))
    }
  } # trait
  data0 <- data.frame(gvar = c(gvarMuS, gvarMuN), cycle = rep(rep(1:5, nTraits), 2), trait = rep(sort(rep(2:(nTraits + 1), 5)), 2), shape = sort(rep(1:2, length(gvarMuS))), group = sort(rep(1:(nTraits * 2), 5)) + nTraits * 2, error = rep(2, length(c(gvarMuS, gvarMuN))))
  data <- rbind(data, data0)
}
p2 <- ggplot(data = data, aes(x = cycle, y = gvar)) + geom_point(aes(colour = factor(trait), shape = factor(shape), size = factor(error))) + geom_line(aes(colour = factor(trait), group = factor(group)))
p2 <- p2 + ylim(min(as.vector(as.matrix(data[, 1]))), max(as.vector(as.matrix(data[, 1]))))
p2 <- p2 + labs(title = "Genetic variance", x = "Cycle", y = "Genetic variance")
p2 <- p2 + scale_shape_discrete(name = "Seedlings in TP", breaks = factor(1:2), labels = c("YES", "NO"))
p2 <- p2 + scale_colour_discrete(name = "Trait", breaks = factor(2:(nTraits + 1)), labels = paste("Trait", 1:nTraits))
if(errorRate > 0){
	p2 <- p2 + scale_size_manual(name = "Error", values = c(2, 3), labels = c(paste(errorRate, "%", sep = ""), "0%"))
}
p2 <- p2 + scale_x_continuous(breaks = 1:5, labels = paste("C", 0:4, sep = ""))
pdf(paste("Fig_", nSim, "timesTrials_", nTraits, "traits_", errorRate, "error.pdf", sep = ""))
p1
p2
dev.off()


