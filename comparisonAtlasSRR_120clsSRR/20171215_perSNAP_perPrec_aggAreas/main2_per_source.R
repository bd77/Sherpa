library(tidyverse)

setwd('D:/WORK/projects/1_urbIam/1_CODE_MATLAB/SHERPA/SHERPA-GITHUB/Sherpa/comparisonAtlasSRR_120clsSRR/20171215_perSNAP_perPrec_aggAreas')
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
  sherpaFinalResults %>%
    select(model, target, source, snap, precursor, relative_potential) %>%
    filter(target==city) %>%
    group_by(model, source, snap) %>%
    summarise(relativePotential=sum(relative_potential)) %>%
    ggplot(aes(fill=snap, y=relativePotential, x=source)) + 
    geom_bar(stat='identity') +
    facet_grid( ~ model) +
    labs(title=city) +
    theme(axis.text=element_text(size=8, angle=90))
            
  ggsave(paste('./OUTPUT_per_source/comparingSRR_atlas_vs_120cells_',city,'.png'))  
}

