#### MATH 267 PROJECT ####
## Hongzhe Liu 
rm(list = ls())

# The main function, simulation the trials 10000 times and return a list contians the simulation information and average sample size of each trial.
simulation <- function(prior.dlt, true.dlt, NS = F, trtDlt = 0.2, delta = 0.05, r1 = 0.5, r2 = 0.95, nmax = 50, nmin = 10, a = 1.2, b = 1, c = 4) {

	#count the number of time that each combination were recommended as MTD during 10000 simulations	
	mtd.counter <- matrix(rep(0, times = 16), nrow = 4, byrow = T)

	total <- 0
	numofsim <- 10000
	for (k in 1:numofsim) {

		aij <- c * prior.dlt
		bij <- c * (1 - prior.dlt)

		## Step 1: find the first toxity dose or highest dose comb is reached
		i = 1
		j = 1
		smpSize = 0
		indic = 0
		repeat {

			# 0 no toxicity, 1 toxicity
			indic <- rbinom(1, 1, true.dlt[i, j])
			smpSize <- smpSize + 1

			# exit if reach the highest dose.
			if (i == 4 & j == 4) 
				break

			if (indic == 0) {
				dir <- rbinom(1, 1, 0.5)
				if (dir == 1 & i < 4) 
					i <- i + 1
				if (dir == 0 & j < 4) 
					j <- j + 1
			} else {
				# exit if the first toxicity occurs.
				break
			}
		}

		## Step 2: baysien algorithm
		repeat {

			if (smpSize > nmax) 
				break

			# update the aij,bij based on the comb that casused toxity patient
			idx <- 1:4
			if (indic == 1) {
				# toxity update aij
				aij[i:4, j:4] <- aij[idx >= i, idx >= j] + 1
			} else {
				# no toxity, updata bij
				bij[1:i, 1:j] <- bij[idx <= i, idx <= j] + 1
			}

			# stopping rules
			if (smpSize > nmin) {

				# rule 3
				if (pbeta(trtDlt + delta, aij[1, 1], bij[1, 1]) < (1 - r1)) {
					# No MTD is recommended
					i = -1
					j = -1
					break
				}

				#rule 4, if the trial reaches the highest dose, then we skip rule 4
				if (!(i == 4 & j == 4)) {

					gt.aij <- aij[idx >= i, idx >= j]
					gt.bij <- bij[idx >= i, idx >= j]
					Pr <- as.vector(1 - pbeta(trtDlt + delta, gt.aij, gt.bij))
					minPr <- min(Pr[2:length(Pr)])
					if (minPr > r2) 
						break
				}
			}

			# utility function to decide the amount of dose for next patient
			expUtl <- -(a + b) * (trtDlt * pbeta(trtDlt, aij, bij) - aij * pbeta(trtDlt, aij + 1, bij)/(aij + bij)) - b * (aij/(aij + bij) - trtDlt)

			#calculate the next dose comb based on the current data
			if (!NS) {
				#skipping method, which function returns the index of true logical object
				nextDose <- which(expUtl == max(expUtl), arr.ind = T)
			} else {
				#non-skipping method
				row <- ifelse(i == 4, i, i + 1)
				col <- ifelse(j == 4, j, j + 1)
				cand.dose <- c(expUtl[1:i, 1:j], expUtl[row, j], expUtl[i, col])
				nextDose <- which(expUtl == max(cand.dose), arr.ind = T)
			}

			i <- nextDose[1]
			j <- nextDose[2]

			# feed the new dose to next patient
			indic <- rbinom(1, 1, true.dlt[i, j])
			smpSize <- smpSize + 1
		}

		# i = -1, j = -1 represents that no MTD is recommended.
		if (i != -1 & j != -1) {
			mtd.counter[i, j] <- mtd.counter[i, j] + 1
		}

		# When the sample size reach 50, the counter smpSize is 51. I need to subtract 1 since the actual size is 50
		smpSize <- ifelse(smpSize == 51, smpSize - 1, smpSize)
		total <- total + smpSize
	}

	# return the info of MTD and average sample size
	averSize <- total/numofsim
	mtdInfo <- mtd.counter/numofsim
	simInfo <- list(mtdInfo = mtdInfo, averSize = averSize)
	invisible(simInfo)
}


# summary the simulation outputs
summary <- function(true.dlt, simInfo, tableType, trtDlt = 0.2) {

	mtdInfo <- simInfo$mtdInfo
	averSize <- simInfo$averSize

	# percentage of hitting the target in 10000 simulation
	error <- round(abs(true.dlt - trtDlt), digits = 2)
	atTarget <- mtdInfo[which(error == 0, arr.ind = T)]
	atTarget <- ifelse(length(atTarget) == 0, 0, atTarget)

	# percentage of deviating from the target among 1 to 5%
	ptsOf1to5 <- sum(mtdInfo[which(error >= 0.01 & error <= 0.05)])
	ptsOf1to5 <- ifelse(length(ptsOf1to5) == 0, 0, ptsOf1to5)

	# percentage of deviating from the target among 5 to 10%
	ptsOf5to10 <- sum(mtdInfo[which(error > 0.05 & error <= 0.1)])
	ptsOf5to10 <- ifelse(length(ptsOf5to10) == 0, 0, ptsOf5to10)

	# percentage of deviating from the target among 1% to 10%
	ptsOf1to10 <- ptsOf1to5 + ptsOf5to10

	# percentage of deviating from the target larger than 10% 
	ptsOfgt10 <- sum(mtdInfo[which(error > 0.1, arr.ind = T)])
	ptsOfgt10 <- ifelse(length(ptsOfgt10) == 0, 0, ptsOfgt10)

	# percentage of non-recommend
	nonRecom <- 1 - sum(mtdInfo)

	# summary, tableType is 3, output info in the form of table 3. Otherwise, output info by table 4
	if (tableType == 3) {
		summaryInfo <- round(c(atTarget, ptsOf1to10, ptsOfgt10, nonRecom, averSize/100), digits = 2) * 100
		names(summaryInfo) <- c("At target", " 1-10% of target", ">10% of target", "None recommended", "Average sample size")
	} else {
		summaryInfo <- round(c(atTarget, ptsOf1to5, ptsOf5to10, ptsOfgt10, nonRecom, averSize/100), digits = 2) * 100
		names(summaryInfo) <- c("At target", "1-5% of target", "5-10% of target", ">10% of target", "None recommended", "Average sample size")
	}

	return(summaryInfo)
}


############ DATA Preparation Start ################
scenA <- matrix(c(4,10,16,22,8,14,20,26,12,18,24,30,16,22,28,34), nrow = 4, byrow = T)/100	 #A
scenB <- matrix(c(2,5,8,11,4,7,10,13,6,9,12,15,8,11,14,17), nrow = 4, byrow = T)/100 		 #B
scenC <- matrix(c(10,25,40,55,20,35,50,65,30,45,60,75,40,55,70,85), nrow = 4, byrow = T)/100 #C

# LEF-
eij <- matrix(runif(16, 0.01, 0.05), ncol = 4)
preA.ls <- scenA - eij
preB.ls <- scenB - eij
preC.ls <- scenC - eij

preA.ls[preA.ls < 0.01] <- 0.01
preB.ls[preB.ls < 0.01] <- 0.01
preC.ls[preC.ls < 0.01] <- 0.01

# LEF+
preA.gt <- scenA + eij
preB.gt <- scenB + eij
preC.gt <- scenC + eij

# LEF*
#eij <- matrix(ifelse(rbinom(16,1,0.5)==1, runif(1,0.01,0.05), runif(1,-0.05,-0.01)), nrow=4, byrow = T)
eij <- matrix(runif(16, -0.05, 0.05), ncol = 4)
preA.rd <- scenA + eij
preB.rd <- scenB + eij
preC.rd <- scenC + eij

preA.rd[preA.rd < 0.01] <- 0.01
preB.rd[preB.rd < 0.01] <- 0.01
preC.rd[preC.rd < 0.01] <- 0.01
############ DATA Preparation END ################

############ Executation Start ################
# Scenario A 
true.dlt <- scenA
lfl.ns.a <- simulation(scenA, true.dlt, T)
lfl.s.a <- simulation(scenA, true.dlt)
lfl.ls.a <- simulation(preA.ls, true.dlt)
lfl.gt.a <- simulation(preA.gt, true.dlt)
lfl.rd.a <- simulation(preA.rd, true.dlt)

summary(scenA, lfl.ns.a, 3)
summary(scenA, lfl.s.a, 3)
summary(scenA, lfl.ls.a, 4)
summary(scenA, lfl.gt.a, 4)
summary(scenA, lfl.rd.a, 4)

# Scenario B
true.dlt <- scenB
lfl.ns.b <- simulation(scenB, true.dlt, T)
lfl.s.b <- simulation(scenB, true.dlt)
lfl.ls.b <- simulation(preB.ls, true.dlt)
lfl.gt.b <- simulation(preB.gt, true.dlt)
lfl.rd.b <- simulation(preB.rd, true.dlt)

summary(scenB, lfl.ns.b, 3)
summary(scenB, lfl.s.b, 3)
summary(scenB, lfl.ls.b, 4)
summary(scenB, lfl.gt.b, 4)
summary(scenB, lfl.rd.b, 4)

# Scenario C
true.dlt <- scenC
lfl.ns.c <- simulation(scenC, true.dlt, T)
lfl.s.c <- simulation(scenC, true.dlt)
lfl.ls.c <- simulation(preC.ls, true.dlt)
lfl.gt.c <- simulation(preC.gt, true.dlt)
lfl.rd.c <- simulation(preC.rd, true.dlt)

summary(scenC, lfl.ns.c, 3)
summary(scenC, lfl.s.c, 3)
summary(scenC, lfl.ls.c, 4)
summary(scenC, lfl.gt.c, 4)
summary(scenC, lfl.rd.c, 4)
############ Executation END ################
