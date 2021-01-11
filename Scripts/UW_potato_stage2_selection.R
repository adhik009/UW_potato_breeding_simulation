rm(list =ls())
setwd("../Datasets")

library(AlphaSimR)
library(tibble)

# Changing crossing strategy to Stage2

# Run 5 replications of the simulation
for(REP in 1:5) {
  load(paste0("BURNIN_PP",REP,".RData"))
  
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
    Stage2_parents = selectInd(Stage2, nInd = 20, use = "pheno", selectTop = TRUE)
    
    # Year 1 - Stage 1
    Stage1 = setPheno(Seedlings, reps=1, varE=9)
    
    # create parent pool
    # Parents from previous year and Stage2 of this year is appended in parentpool
    ParentPool <- append(Parents, Stage2_parents) 
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
  }
  
  ## Create tibble for final output
  # Mean centered on final burn-in year
  # Variance standardized to percentage of final burn-in year
  output = tibble(scenario=rep("Stage2_sel",40),
                  year=-19:20,
                  gain_stage1=meanStage1-meanStage1[20],
                  gain_stage4=meanStage4-meanStage4[20],
                  relVar=varStage1/varStage1[20]*100,
                  Parents_dropped = Parents_dropped)
  print("Iteration complete !!")
  saveRDS(output,paste0("Stage2_sel_",REP,".rds"))
}
