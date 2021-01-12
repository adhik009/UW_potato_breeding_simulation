
rm(list =ls())
setwd("../Datasets")

library(dplyr) 
library(ggplot2) 

scenarios = c("Stage4_sel",
              "Stage3_sel",
              "Stage2_sel",
              "Stage1_sel",
              "Stage3_GS",
              "Stage2_GS")

# Read in output data
rawData = vector("list",5*length(scenarios))
i = 0L
for(SCENARIO in scenarios){
  for(REP in 1:5){
    i = i+1L
    FILE = paste0(SCENARIO,"_",REP,".rds")
    rawData[[i]] = readRDS(FILE)
  }
}
rawData = bind_rows(rawData)

# Genetic gain plot
pdf("../Plots/Genetic.gain.stage1.pdf")
gain = rawData %>%
  group_by(scenario,year) %>%
  summarise(
    MEAN = mean(gain_stage1),
    SE = sd(gain_stage1)/sqrt(5)
  )
ggplot(gain,aes(x=year,y=MEAN,color=scenario))+
  geom_ribbon(aes(x=year,ymin=MEAN-SE,ymax=MEAN+SE,
                  fill=scenario),alpha=0.2,linetype=0)+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()+
  scale_x_continuous("Year",limits=c(0,20))+
  scale_y_continuous("Genetic Gain", limits=c(-1,10))
dev.off()

pdf("../Plots/Genetic.gain.stage4.pdf")
gain = rawData %>%
  group_by(scenario,year) %>%
  summarise(
    MEAN = mean(gain_stage4),
    SE = sd(gain_stage4)/sqrt(5)
  )
ggplot(gain,aes(x=year,y=MEAN,color=scenario))+
  geom_ribbon(aes(x=year,ymin=MEAN-SE,ymax=MEAN+SE,
                  fill=scenario),alpha=0.2,linetype=0)+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()+
  scale_x_continuous("Year",limits=c(0,20))+
  scale_y_continuous("Genetic Gain", limits=c(-1,10))
dev.off()

# Genetic variance plot
pdf("../Plots/Genetic_variance_sel.stages.pdf")
relVar = rawData %>%
  group_by(scenario,year) %>%
  summarise(
    MEAN = mean(relVar),
    SE = sd(relVar)/sqrt(5)
  )
ggplot(relVar,aes(x=year,y=MEAN,color=scenario))+
  geom_ribbon(aes(x=year,ymin=MEAN-SE,ymax=MEAN+SE,
                  fill=scenario),alpha=0.2,linetype=0)+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()+
  scale_x_continuous("Year",limits=c(0,20))+
  scale_y_continuous("Relative Variance (%)",limits=c(0,200))
dev.off()

# Parent replacement plot
pdf("../Plots/Parent_replacement_parentpool.pdf")
PD = rawData %>%
  group_by(scenario,year) %>%
  summarise(
    MEAN = mean(Parents_dropped),
    SE = sd(Parents_dropped)/sqrt(5))

ggplot(PD,aes(x=year,y=MEAN,color=scenario))+
  geom_ribbon(aes(x=year,ymin=MEAN-SE,ymax=MEAN+SE,
                  fill=scenario),alpha=0.2,linetype=0)+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()+
  scale_x_continuous("Year",limits=c(1,20))+
  scale_y_continuous("Parent replacement", limits=c(1,20))
dev.off()

# Genomic selection accuracy
pdf("../Plots/GS_accuracy_plot.pdf")
GS_AC = rawData %>%
  group_by(scenario,year) %>%
  summarise(
    MEAN = mean(GS_accuracy),
    SE = sd(GS_accuracy)/sqrt(5)
  ) %>% filter(scenario %in% c("Stage2_GS","Stage3_GS")) %>% filter(year > 0)

ggplot(GS_AC,aes(x=year,y=MEAN,color=scenario))+
  geom_ribbon(aes(x=year,ymin=MEAN-SE,ymax=MEAN+SE,
                  fill=scenario),alpha=0.2,linetype=0)+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()+
  scale_x_continuous("Year",limits=c(1,20))+
  scale_y_continuous("GS Accuracy")
dev.off()

# Efficiency of converting genetic diversity to genetic gain
# Efficiency is calculated as decrease in genetic variance in (%) divided by genetic gain at year 20

# pdf("Efficienty.plot.SD_19.2020.pdf")
# Eff <- rawData %>% select(scenario,year,relVar,gain_stage4) %>% filter(year == 20) %>% 
#   group_by(scenario) %>%
#   summarise(
#     Variance = 100 - sd(relVar),
#     Gain = sd(gain_stage4)) 
# 
# Eff <- transform(Eff, Efficiency = Variance / Gain)
# 
# ggplot(Eff, aes(x = scenario, y = Efficiency, group =1)) + geom_line(size = 1, color = "blue") + geom_point(shape = 1, size = 2) +
#   labs(x = "Selection stage", y = "Efficiency") + theme_bw()
# dev.off()