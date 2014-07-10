############################################################
######################## crossing ##########################
############################################################

## make gamete ##
makeGamete <- function(geno, pos){
	diff <- diff(pos)
	rec <- (1 - exp(-2 * diff / 100)) / 2     # Haldane
	rec[rec < 0] <- 0.5
	rec <- c(0.5, rec)
	
	sample <- runif(length(rec))
	crossOver <- ((rec - sample) >= 0)
	
	selectHaplo <- cumsum(crossOver) %% 2
	selectHaplo <- selectHaplo + 1
	
	gamete <- geno[1, ]
	gamete[selectHaplo == 2] <- geno[2, selectHaplo == 2]
	return(gamete)     # vector
}

## make one progeny ##
makeProgeny <- function(genoMat, genoPat, pos){
	progeny <- rbind(makeGamete(genoMat, pos), makeGamete(genoPat, pos))
	return(progeny)     # matrix
}

## make progenies (combinations of parents are given) ##
makeProgenies <- function(parents, geno, pos){
	progenies <- matrix(NA, 2 * nrow(parents), ncol(geno))
	for(par in 1:nrow(parents)){
		genoMat <- geno[c(parents[par, 1] * 2 - 1, parents[par, 1] * 2), ]
		genoPat <- geno[c(parents[par, 2] * 2 - 1, parents[par, 2] * 2), ]
		progenies[c(par * 2 - 1, par * 2), ] <- makeProgeny(genoMat, genoPat, pos)
	} # par
	return(progenies)     # matrix
}

## random mating ##
randomMate <- function(popSize, geno, pos){
	parents <- t(sapply(rep(nrow(geno) / 2, popSize), sample, size=2))
	progenies <- makeProgenies(parents, geno, pos)
	return(progenies)
}

## random mating , then getting seeds from all parents ##
## The size of progeny population should not be smaller than that of parent population ##
randomMateAll <- function(popSize, geno, pos){
	parent1 <- rep(1:(nrow(geno) / 2), popSize %/% (nrow(geno) / 2))
	if(popSize %% (nrow(geno) / 2) != 0) parent1 <- c(parent1, 1:(popSize %% (nrow(geno) / 2)))
	parent2 <- rep(NA, popSize)
	for(pollen in 1:popSize){
		parent2[pollen] <- sample((1:(nrow(geno) / 2))[-parent1[pollen]], size = 1)
	} # pollen
	parents <- cbind(parent1, parent2)
	
	progenies <- makeProgenies(parents, geno, pos)
	return(progenies)
}
###################################################################################
################### make genotypic and phenotypic values ##########################
###################################################################################

#----------------- calculate genotypic values --------------------------
calcGenotypicValue <- function(geno, mapData, trait){
	nPop <- nrow(geno) / 2
	
	## genotypic value in one locus in one individual ##
	gv1pos <- function(geno1pos, actType, effect){
		geno1pos <- as.matrix(geno1pos, nrow = 2)
		coef <- ifelse(actType == 0, geno1pos[1, ] + geno1pos[2, ], geno1pos[1, ] * geno1pos[2, ])
		return(effect * prod(coef))
	}
	
	## genotypic value in one individual ##
	gv1ind <- function(genoVec, mapData, trait){
		nQTL <- max(mapData$effectID)
		geno1ind <- rbind(genoVec[1:(length(genoVec) / 2)], genoVec[(length(genoVec) / 2 + 1):length(genoVec)])
		
		power <- 0
		for(i in 1:nQTL){
			power <- power + gv1pos(geno = geno1ind[, mapData$effectivePos[mapData$effectID == i]], actType = mapData$actionType[mapData$effectID == i], effect = mapData$effects[i, trait])
		} # i
		return(power)
	}
	
	# genotypic value of all individuals
	genoVec <- cbind(geno[1:nPop * 2 - 1, ], geno[1:nPop * 2, ])
	gv <- apply(genoVec, 1, gv1ind, mapData = mapData, trait = trait)
	
	return(gv)
}


#---------------------- calculate phenotypic values -------------------------
# no correlation between environments!
calcPhenotypicValue <- function(gv, errorVar = NULL, H2 = NULL){
	if((!is.null(errorVar) & !is.null(H2)) | (is.null(errorVar) & is.null(H2))){
		stop("I cannot make phenotypic value!")
	}else{
		if(!is.null(errorVar)){
			pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
		}else{
			varG <- var(gv)
			errorVar <- varG * (1 - H2) / H2
			pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
		}
	}
	return(list(pv = pv, var = errorVar))
}


###################################################################################
############################ selection (GS: kin.blup4.4) ##########################
###################################################################################

## calculate score ##
makeScore <- function(geno){
	nPop <- nrow(geno) / 2
	
	#-------------------------- -1,0,1 ------------------------------
	score <- geno[1:nPop * 2 - 1, ] + geno[1:nPop * 2,]
	score <- score / 2     # Now, alleles are represented as -1, 1
	
	#------------------------- ID -------------------------------
	rownames(score) <- 1:nPop
	return(score)
}

## prediction ##
# phenoData includes all GID, phenotypic values, genotypic values and error variance.
genomicPrediction <- function(mapData, geno, phenoData, reduce = F){
	genoScore <- makeScore(geno)     # rownames are "1:nPop"
	genoScore <- genoScore[, mapData$markerPos]
	K <- A.mat(genoScore)
	model <- kin.blup4.4(data = phenoData, geno = "GID", pheno = "y", K = K, R = phenoData$var, reduce = reduce)
	pred <- model$g
	
	return(pred)
}

## selection ##
# selected candidate's GID should be shown in parameter.
genomicSelection <- function(mapData, geno, phenoData, reduce = F, nSelect, selectPopGID, weights = 1){
	nTraits <- length(weights)
	ncolPhenoData <- ncol(phenoData) / nTraits
	predAll <- rep(0, nrow(geno) / 2)
	for(trait in 1:nTraits){
		phenoData.trait <- phenoData[, ((trait - 1) * ncolPhenoData + 1):(trait * ncolPhenoData)]
		names.GID <- which(substr(names(phenoData.trait), 1, 3) == "GID")
		names.y <- which(substr(names(phenoData.trait), 1, 1) == "y")
		names.var <- which(substr(names(phenoData.trait), 1, 3) == "var")
		names(phenoData.trait)[c(names.GID, names.y, names.var)] <- c("GID", "y", "var")
		pred <- genomicPrediction(mapData = mapData, geno = geno, phenoData = phenoData.trait, reduce = reduce)
		predAll <- predAll + weights[trait] * pred
	} # trait
	orderGID <- order(predAll, decreasing = T)
	selectedGID <- NULL
	for(i in 1:length(orderGID)){
		if(sum(which(orderGID[i] == selectPopGID)) > 0) selectedGID <- c(selectedGID, orderGID[i])
		if(length(selectedGID) == nSelect) break
	} # i
	selectedGID <- sort(selectedGID)
	
	return(selectedGID)
}





####################################################################################################
#################### make positions and base population by "genome.exe" ############################
####################################################################################################

# -p 1 20 means one subpopulation out of which you extract 20 haplotypes
# -N 20 means Ne = 20
# a chr of 800 Mb is split into 20000 segments with a recombination frequency of 0.01 cM between them
# this corresponds therefore to a chr of 200 cM
# -s 1000 means put 1000 mutations on the tree per chromosome.
# Note that polymorphisms below a certain MAF will be
# filtered out below, so make sure you put enough on

# nPopsSamples is for the -p option; NULL will result in -p 1 effPopSize
# a single number, ps, will result in -p 1 ps
# a vector of two numbers, ps1 & ps2, will result in -p 2 ps1 ps2, etc.
# WARNING: not yet set up to return results for two subpopulations
# nMrkOrMut is for the number of markers or the mutation rate; 
# a single number, m, will result in -s m / nChr
# a vector of two numbers, l & u, will result in  -len l -mut u
# note that effPopSize can be a file name with a demographic scenario

## get loci and base population by coalescent (genome.exe) ##
getCoalescentSim <- function(nPopsSamples=NULL, effPopSize=100, nChr=7, nPiecesPerChr=15000, recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, seed=as.integer(Sys.time()), tree=0, extraOpt=NULL){
	# figure out where the executable is
	genomeExec <- "/Applications/genome/genome"
	
	# set up nPopsSamples to match effPopSize
	if (is.null(nPopsSamples)){
		if (mode(effPopSize) == "character"){ #must be a file name
			lastGen <- readLines(con=effPopSize, n=1)
			nPopsSamples <- as.integer(strsplit(lastGen, "[[:space:]]")[[1]])[-1]
		} else{
			nPopsSamples <- effPopSize
		}
	}
	# set up all the options for genome
	options <- c("-pop", paste(length(nPopsSamples), paste(nPopsSamples, collapse=" ")))
	options <- cbind(options, array(c("-N", effPopSize, "-c", nChr, "-pieces", nPiecesPerChr, "-rec", recBTpieces, "-maf", minMAF, "-seed", seed, "-tree", tree), c(2, 7)))
	if (length(nMrkOrMut) > 1){
		options <- cbind(options, c("-len", nMrkOrMut[1]), c("-mut", nMrkOrMut[2]))
	} else{
		options <- cbind(options, c("-s", round(8 * nMrkOrMut / nChr)))
	}
	outFile <- paste(as.integer(Sys.time()), runif(1), "GenomeOut.txt", sep="")
	systemCall <- paste(genomeExec, paste(options, collapse=" "), extraOpt, ">", outFile, sep=" ")
	#print(systemCall)
	system(systemCall)
	genOut <- readLines(outFile, n=-1)
	
	# pull out lines with Newick trees if asked
	if (tree){
		stripTreeLine <- function(treeLine){
			firstEq <- regexpr("=", treeLine) + 2
			return(substr(treeLine, firstEq, nchar(treeLine) - 1))
		}
		treeLines <- sapply(genOut[grep("The tree for fragment", genOut)], stripTreeLine)
	} else treeLines <- NULL
	
	# make a map of the SNP positions
	strtPos <- grep("SNP genetic position", genOut)
	map <- NULL
	for (chr in 1:nChr){
		# genome puts markers at both 0 and 1 so there is some inflation that needs to be deflated (sorry, strange calculation)
		posVec <- round(as.numeric(strsplit(genOut[strtPos[chr]+1], split=" ")[[1]]) * 100 * (nPiecesPerChr - 1) * recBTpieces, 6)
		if (length(posVec > 0)) map <- rbind(map, cbind(chr, posVec))
	}
	colnames(map) <- c("Chr", "Pos")
	
	# pull out the SNP themselves
	strtInd <- NULL
	for (pop in 1:length(nPopsSamples)){
		strtInd <- c(strtInd, grep(paste("POP", pop, ":", sep=""), genOut))
	}     ## 全部POP1なのに、まわす必要ある？ ##
	genOut <- substring(genOut[strtInd], 7)
	nSNP <- nchar(genOut[1])
	markers <- t(array(as.numeric(unlist(strsplit(genOut, split=""))), c(nSNP, sum(nPopsSamples))))
	## unlistで全部横並びになって、arrayで縦並びになるから、最後に転置 ##
	
	# remove markers with too low MAF
	freq <- colMeans(markers)
	maf <- abs((freq > 0.5) - freq)
	markers <- markers[, maf >= minMAF]
	map <- map[maf >= minMAF,]
	maf <- maf[maf >= minMAF]
	
	# if a fixed number of segregating sites, randomly pick nMrkOrMut markers
	if (length(nMrkOrMut) == 1){
		if (ncol(markers) < nMrkOrMut){
			print("Warning! Settings such that fewer markers simulated than demanded")
			print(paste("There were", ncol(markers), "markers"))
		} else{
			keepMrk <- sort(sample(ncol(markers), nMrkOrMut))
			markers <- markers[, keepMrk]
			map <- map[keepMrk, ]
			maf <- maf[keepMrk]
		}
	}
	nMrk <- ncol(markers)

	# find markers with exact same position and jitter them
	for (chr in 1:nChr){
		mrkThisChr <- map[,"Chr"] == chr
		uniqPos <- unique(map[mrkThisChr,"Pos"])
		for (pos in uniqPos){
			mrkAtPos <- which(mrkThisChr & map[,"Pos"] == pos)
			if (length(mrkAtPos) > 1) map[mrkAtPos,"Pos"] <- map[mrkAtPos,"Pos"] + sort(round(runif(length(mrkAtPos), 0, 100 * recBTpieces), 6))
		}
	}
	
	try(system(paste("rm", outFile)), silent = TRUE)
	
	map <- as.data.frame(map)
	return(list(markers=markers, map=map))
	
} #END getCoalescentSim


#######################################################################################
############################ make initial population ##################################
#######################################################################################

createInitPop <- function(baseData, nPopInit, seed=round(runif(1, 0, 1e9))){
	print("Creating an initial population for breeding")
	set.seed(seed)
	
	## double gamete ##
	doubleGametes <- function(gametes){
		genotypes <- matrix(NA, 2 * nrow(gametes), ncol(gametes))
		genotypes[1:nrow(gametes) * 2, ] <- gametes
		genotypes[1:nrow(gametes) * 2 - 1, ] <- gametes
	}
	
	# genotypes of base population
	# "genome" can create haplotypes, so we should create diploids for cassava.
	geno <- doubleGametes(baseData$founderHaps)
	
	# These genotypes are DHs and they need recombination.
	# First, we make F1. Then, we implement one cycle of random mating to produce recombinations.
	geno <- randomMate(popSize = nrow(geno) / 2, geno = geno, pos = baseData$mapData$map$Pos)
	geno <- randomMate(popSize = nPopInit, geno = geno, pos = baseData$mapData$map$Pos)
	
	# alleles should be represented as -1, 1
	# 0,1 -> -1,1
	geno <- geno * 2 - 1
	
	# return data
	mapData <- baseData$mapData
	return(list(mapData = mapData, geno = geno))
}


#######################################################################################
#################### make map information with QTL effects ############################
#######################################################################################

# varEffect shoul be a vector whose length is same as the number of traits.
# If there is one trait, corrEffects should be NULL.
# If there are two effects, corrEffects should be a vector whose length is one.
# If there are more than two effects, corrEffects should be matrix.
# corrEffects represents the upper triangular portion of correlation matrix.
#             1   c12 c13 c14
#             c12 1   c23 c24
#             c13 c23 1   c34
#             c14 c24 c34 1        like this.
# !!CAUTION!! The variance of effects in the first trait is not varEffects because of the weight, though those of the other traits are same as varEffects!
# All effective loci are pleiotropic.

## make marker loci and QTL (two traits; completely pleiotropic QTL) ##
makeMap <- function(map, nLoci, nMarkers, nQTL, propDomi = 0, interactionMean = 0, varEffects = 1, corrEffects = NULL){
	# select QTL and effect type
	nEffectiveLoci <- 1 + rpois(n = nQTL, lambda = interactionMean)     # main position + epistatic loci
	posEffectiveLoci <- sample(1:nLoci, sum(nEffectiveLoci))     # no replace
	actionType <- rbinom(sum(nEffectiveLoci), 1, propDomi)     # 0: additive, 1: dominance
	effectID <- NULL
	for(i in 1:nQTL){
		effectID <- c(effectID, rep(i, nEffectiveLoci[i]))
	} # i
	
	# make effects for the first trait
	if(length(varEffects) == 1){
		vEffect <- varEffects
	}else{
		vEffect <- varEffects[1]
	}
	effect <- rnorm(nQTL, 0, sqrt(varEffects))	
	
	# make effects for the other traits
	if((is.null(corrEffects)) | (!is.matrix(corrEffects))){
		effects <- as.matrix(effect, nQTL, 1)
		rownames(effects) <- 1:nQTL
		colnames(effects) <- "trait1"
	}else{
		varEffects[1] <- var(effect)
		cov <- corrEffects * sqrt(varEffects)     # multiply by row
		cov <- t(t(cov) * sqrt(varEffects))     # multiply by column
		if(sum(eigen(cov)$values < 0) == 0){
			chol.cov <- chol(cov)[, -1]
		}else{
			library(Matrix)
			cov.modified <- as.matrix(nearPD(cov)$mat)
			chol.cov <- chol(cov.modified)[, -1]
		}
		
		x1 <- scale(effect)
		x2 <- matrix(rnorm(nQTL * (length(varEffects) - 1)), nQTL, length(varEffects) - 1)
		xs <- cbind(x1, x2)
		effect2 <- xs %*% chol.cov
		
		effects <- cbind(effect, effect2)
		rownames(effects) <- 1:nQTL
		colnames(effects) <- paste("trait", 1:length(varEffects), sep = "")
	}
	
	# select markers
	if(nLoci - length(posEffectiveLoci) < nMarkers) print("Warning: Number of markers is not enough!")
	mrkPos <- sort(sample((1:nLoci)[-posEffectiveLoci], nMarkers, replace = F))
	
	# return data
	return(list(map = map, markerPos = mrkPos, effectID = effectID, effectivePos = posEffectiveLoci, actionType = actionType, effects = effects))
}



#################################################################################################
#################### remake map information by adjusting QTL effects ############################
#################################################################################################

adjustMap <- function(initialPopulationData, geneticVariance, nTraits){
	mapData <- initialPopulationData$mapData
	geno <- initialPopulationData$geno
	
	for(trait in 1:nTraits){
		g <- calcGenotypicValue(geno = geno, mapData = mapData, trait = trait)
		vg <- var(g)
		vg.ideal <- geneticVariance[trait]
		coef <- sqrt(vg.ideal / vg)
		mapData$effects[, trait] <- mapData$effects[, trait] * coef
	} # trait
	mapData$effects <- as.matrix(mapData$effects, ncol = nTraits)
	
	return(list(mapData = mapData, geno = geno))
}








