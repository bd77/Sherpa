library(tidyverse)

setwd('D:/WORK/projects/1_urbIam/1_CODE_MATLAB/SHERPA/SHERPA-GITHUB/resAtlas/20171215_perSNAP_perPrec_aggAreas')

#read file with results
#sherpaAtlas is the case used in the Atlas
#sherpa120cls is the new case to be included in the Feb18 Sherpa
sherpaAtlas <- read_csv2("20171215_sherpa_chimere_PM25_perSNAP_perPrec_aggAreas.txt", 
                           skip = 5)
sherpaAtlas$relative_potential <- as.numeric(sherpaAtlas$relative_potential)
sherpaAtlas$snap <- as.factor(sherpaAtlas$snap)

sherpa120cls <- read_csv2("20171219_sherpa_chimere_noFW_120cls_PM25_perSNAP_perPrec_aggAreas.txt", 
                         skip = 5)
sherpa120cls$relative_potential <- as.numeric(sherpa120cls$relative_potential)
sherpa120cls$snap <- as.factor(sherpa120cls$snap)

#concatenate results
sherpaFinalResults <- rbind(sherpaAtlas, sherpa120cls)

#loop over cities
for (city in unique(sherpaFinalResults$target)) {
  print(paste0('processing ',city))
  #select useful columns, filter per city, group and summarize per source-snap, ggplot
  res1 <- sherpaFinalResults %>%
    select(model, target, source, snap, precursor, relative_potential) %>%
    filter(target==city, model=='sherpa_chimere') %>%
    group_by(source, snap) %>%
    summarise(relativePotential=sum(relative_potential))
    
  res2 <- sherpaFinalResults %>%
    select(model, target, source, snap, precursor, relative_potential) %>%
    filter(target==city, model=='sherpa_chimere_noFW_120cls') %>%
    group_by(source, snap) %>%
    summarise(relativePotential=sum(relative_potential))
  
  diff_res1_res2 <- res1
  diff_res1_res2$relativePotential <- (res1$relativePotential - res2$relativePotential)
  names(diff_res1_res2)[3]='difference between relative potentials'
  
  ggplot(data=diff_res1_res2, aes(fill=source, y=`difference between relative potentials`, x=snap)) + 
    geom_bar(stat='identity') +
    ylim(-10,10)
    labs(title=city) 

  ggsave(paste('./OUTPUT/diff_comparingSRR_atlas_vs_120cells_',city,'.png'))  
}

