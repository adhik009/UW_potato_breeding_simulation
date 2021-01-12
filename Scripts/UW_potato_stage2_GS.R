rm(list =ls())
setwd("C:/Users/ANeil/Google Drive/UW_Madison/AlphaSimR/UW Potato simulation/UW_potato_trial/GS_parents/Output")

library(AlphaSimR)
library(tibble)

# Using GS at stage 3

# Run 5 replications of the simulation
for(REP in 1:5) {
  load(paste0("BURNIN_PP",REP,".RData"))
  
  # Train GS model
  GS_accuracy <- c()
  
  ## Simulate years 21 through 40
  for(year in 21:40){
    
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
    #update training population
    trainPop <- append(trainPop[200:1000], Stage3)
    
    # train the GS model
    gs = RRBLUP(pop=trainPop, use = "pheno", useQtl = FALSE)
    
    # GS for parent selection at stage 2
    Stage2 = setEBV(Stage2, gs, value = "bv")
    Stage2_parents <- selectInd(Stage2, nInd = 20, use = "ebv", selectTop = TRUE)
    
    # Add the stage 2 parents to parentpool
    ParentPool <- append(Stage2_parents, Parents)
    
    # Use GS to select the best parents
    ParentPool = setEBV(ParentPool, gs, value = "bv")
    NewParents = selectInd(ParentPool, nInd = 20, use = "ebv", selectTop = TRUE)
    
    # Calculate number of old parents replaced in parentpool with new parents
    Parents_dropped[year] <- c(length(setdiff(Parents@id,New_parents@id)))
    
    # Update Parents
    Parents <- NewParents
    
    # Year 0 (Crossing)
    F1 = randCross(NewParents, nCrosses=50, nProgeny=400)
    Seedlings = F1
    
    # Report mean and variance
    meanStage1[year] = meanG(Stage1)
    meanStage4[year] = meanG(Stage4)
    varStage1[year] = varG(Stage1)
    GS_accuracy[year] = cor(gv(ParentPool), ebv(ParentPool))
  }
  
  ## Create tibble for final output
  # Mean centered on final burn-in year
  # Variance standardized to percentage of final burn-in year
  output = tibble(scenario=rep("Stage2_GS",40),
                  year=-19:20,
                  gain_stage1=meanStage1-meanStage1[20],
                  gain_stage4=meanStage4-meanStage4[20],
                  relVar=varStage1/varStage1[20]*100,
                  Parents_dropped = Parents_dropped,
                  GS_accuracy=GS_accuracy)
  saveRDS(output,paste0("Stage2_GS_",REP,".rds"))
  print("Iteration complete")
}