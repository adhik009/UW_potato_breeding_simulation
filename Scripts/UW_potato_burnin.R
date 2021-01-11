rm(list=ls())

setwd("/Volumes/GoogleDrive/My Drive/UW_Madison/AlphaSimR/UW_potato_breeding_simulation/Scripts")
library(AlphaSimR)

## Create founder haplotypes
founderPop = runMacs2(nInd = 100, nChr = 12, segSites = 500, Ne = 100,
                      inbred = FALSE, ploidy = 4L)

## Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 10, mean = 0, var = 1)$setVarE(H2=0.4)

# add SNP chip
SP$addSnpChip(nSnpPerChr = 300)
SP$restrSegSites(overlap = TRUE, minSnpFreq = 0.05)

# double ploidy to create tetraploid population
diploid = newPop(founderPop)
Tetraploid = doubleGenome(diploid)

# Run 5 replications of the simulation
for(REP in 1:5){
  
  # setup empty shells for filling training population, parentpool and Parents_dropped 
  trainPop <- c()  
  ParentPool <- c()
  Parents_dropped <- c()
  breeding_efficiency <- c()
  
  ## Fill breeding program pipeline
  Parents = newPop(Tetraploid)
  
  # Year 4 - Stage 4
  #---------------------------------------------------------------------
  # 20 clones selected from Stage3 and planted for Stage 4 evaluation
  F1 = randCross(Parents, nCrosses=50, nProgeny=400) 
  Seedlings = F1
  Stage1 = setPheno(Seedlings, reps=1, varE=9)
  Stage2 = selectInd(Stage1, nInd=1000)
  Stage2 = setPheno(Stage2, reps=1, varE=2)
  Stage3 = selectInd(Stage2, nInd=200)
  Stage3 = setPheno(Stage3, reps=1, varE=1.15)
  Stage4 = selectInd(Stage3, nInd=20)
  Stage4 = setPheno(Stage4, reps=2, varE=1.3)
  
  # Year 3 - Stage 3
  #--------------------------------------------------------------------
  # 200 clones selected from Stage 2 and planted for Stage 3 evaluation
  F1 = randCross(Parents, nCrosses=50, nProgeny=400) 
  Seedlings = F1
  Stage1 = setPheno(Seedlings, reps=1, varE=9)
  Stage2 = selectInd(Stage1, nInd=1000)
  Stage2 = setPheno(Stage2, reps=1, varE=2)
  Stage3 = selectInd(Stage2, nInd=200)
  Stage3 = setPheno(Stage3, reps=1, varE=1.15)
  
  # Year 2 - Stage 2
  #---------------------------------------------------
  # 1000 clones selected from stage 1 evaluation and planted for stage 2 evaluation
  F1 = randCross(Parents, nCrosses=50, nProgeny=400)
  Seedlings = F1
  Stage1 = setPheno(Seedlings, reps=1, varE=9)
  Stage2 = selectInd(Stage1, nInd=1000)
  Stage2 = setPheno(Stage2, reps=1, varE=2)
  
  # Year 1 - Stage 1 
  #-----------------------------------------------------------------------
  #Seedling tubers planted and evaluated
  F1 = randCross(Parents, nCrosses=50, nProgeny=400)
  Seedlings = F1
  Stage1 = setPheno(Seedlings, reps=1, varE=9)
  
  # Year 0 
  #----------------------------------------------------------------------
  #Crossing and seedlings; seedling tubers harvested; no selection
  F1 = randCross(Parents, nCrosses=50, nProgeny=400)
  Seedlings = F1
  
  ## Create vectors for output (40 years total)
  meanStage1 = numeric(40)
  meanStage4 = numeric(40)
  varStage1 = numeric(40)
  
  
  ## Simulate 20 years of burn-in
  for(year in 1:20){
    
    # Year 4 - Stage 4
    Stage4 = selectInd(Stage3, nInd=20)
    Stage4 = setPheno(Stage4, reps=2, varE=1.3)
    
    # Year 3 - Stage 3
    Stage3 = selectInd(Stage2, nInd=200)
    Stage3 = setPheno(Stage3, reps=1, varE=1.15)
    
    # Year 2 - Stage 2
    Stage2 = selectInd(Stage1, nInd=2000)
    Stage2 = setPheno(Stage2, reps=1, varE=2)
    
    # Year 1 - Stage 1
    Stage1 = setPheno(Seedlings, reps=1, varE=9)
    
    # create parent pool
    # Parents from previous year and Stage4 of this year is appended in parentpool
    ParentPool <- append(Parents, Stage4) 
    New_parents = selectInd(ParentPool, nInd = 20, selectTop = TRUE)
    
    # Calculate number of old parents replaced in parentpool with new parents
    Parents_dropped[year] <- c(length(setdiff(Parents@id,New_parents@id)))
    
    # (Parents) until now is parents from previous year. It is updated with new parents below.
    # update parents for this year
    Parents <- New_parents
    
    # Year 0 (Crossing)
    F1 = randCross(New_parents, nCrosses=50, nProgeny=400)
    Seedlings = F1
    
    # Report mean and variance
    meanStage1[year] = meanG(Stage1)
    meanStage4[year] = meanG(Stage4)
    varStage1[year] = varG(Stage1)
    
    # Training Population for gBLUPs: populated over last 5 burn-in years
    if(year > 15){
      trainPop <- append(Stage3,trainPop)
    }
  }
  
  ## Save global environment
  save.image(paste0("/Volumes/GoogleDrive/My Drive/UW_Madison/AlphaSimR/UW_potato_breeding_simulation/Datasets/",
                    "BURNIN_PP",REP,".RData"))
  print("1 iteration complete !!")
}