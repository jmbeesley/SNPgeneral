### 
### Generate annotated Manhattan plots

library(dplyr)
library(ggplot2)
library(readr)

# data input
homeDir            = "/Users/jonathanbeesley/GoogleDrive/work_informatics_projects/"
bcRegions          = read.table(paste(homeDir, "BCAC_FM_includesUCSCandEUGENE/data/BCAC_FM.regions.noHeader.bed", sep = ""), header=F)
eugeneKnown        = read.table(paste(homeDir, "BCAC_FM_includesUCSCandEUGENE/data/knownBCloci_EUGENEinfo", sep = ""), header=T)
eugeneSNPpositions = read.table(paste(homeDir, "BCAC_FM_includesUCSCandEUGENE/data/EUGENE.BREAST.snps.bed", sep = ""), header=F)
bcacFMdata         = read_delim(paste(homeDir, "BCAC_finemapping/bcac.oa.finemapping_added.rsIDs", sep = ""), delim= '\t', col_names=T)

#
eugeneKnown            = left_join(eugeneKnown, eugeneSNPpositions, by = c("EUGENE.eQTL" = "V4"))
eugeneKnown            = eugeneKnown %>% filter(EUGENE.analysis == "BREAST")
overallLoci            = eugeneKnown %>% select(locus, GWAS) %>% filter(GWAS == "breast.unconditional") %>% unique
ernegLoci              = eugeneKnown %>% select(locus, GWAS) %>% filter(GWAS == "breastERneg.unconditional") %>% unique
eugene.overall.regions = left_join(overallLoci, bcRegions, by = c("locus" = "V4"))
eugene.erneg.regions   = left_join(ernegLoci, bcRegions, by = c("locus" = "V4"))

################
### OVERALL RISK 

PlotManhattanOverall = function(i) {
 
assocLocus = eugene.overall.regions %>% filter(locus == i) %>% select(locus) %>% as.character
locusChr   = eugene.overall.regions %>% filter(locus == i) %>% select(V1)
minPos     = eugene.overall.regions %>% filter(locus == i) %>% select(V2) %>% as.numeric
maxPos     = eugene.overall.regions %>% filter(locus == i) %>% select(V3) %>% as.numeric
eugeneSNPs = eugeneKnown %>% filter(locus == i) %>% select(V3) %>% mutate(eugene = "EUGENE_SNP")
candidates = bcacFMdata %>% filter(locus == sub("BCAC_FM_", "", i), type != "excluded_SNP") %>% select(position, type)

assocData  = read.table(paste(homeDir, "BCAC_summaryStats/overall_position.space.pvalue_byChromosome/", locusChr[[1]], sep=""), header=F, )

graphData  = assocData %>% filter(V1 > minPos & V1 < maxPos)
graphData  = left_join(graphData, eugeneSNPs, by = c("V1" = "V3"))
graphData  = left_join(graphData, candidates, by = c("V1" = "position"))
graphData$type[is.na(graphData$type)] = "excluded_SNP"
graphData$eugene[is.na(graphData$eugene)] = "nonEUGENE_SNP"

# plot - SNP type in colour, EUGENE SNP is a triangle
ggplot(graphData) + 
  geom_point(aes(V1, -log10(V2), colour = type)) +
  geom_point(data = graphData %>% filter(eugene == "EUGENE_SNP"), aes(V1, -log10(V2)), shape = 2, size = 2, colour = 'red', stroke = 1) + 
  xlab(paste(locusChr[[1]]," position", sep='')) +
  ylab("-log10 P") +
  geom_hline(yintercept = -log10(5E-8), colour = 'red') + 
  scale_colour_manual(values = c("candidate_SNP"                = "blue", 
                                 "weaker_secondary_signals"     = "green", 
                                 "bcac_cimba_secondary_signals" = "violet",
                                 "cimba_secondary_signals"      = "orange",
                                 " "                            = "black"), 
                       na.value = "black") +
  ggtitle(assocLocus)

ggsave(paste(homeDir, "BCAC_FM_includesUCSCandEUGENE/manhattans/overall/",assocLocus, "_overallRisk", ".png", sep=""))

}

lapply(unique(eugene.overall.regions$locus), PlotManhattanOverall)

##############
### ERNEG RISK 

PlotManhattanERneg = function(i) {
 
assocLocus = eugene.erneg.regions %>% filter(locus == i)  %>% select(locus) %>% as.character
locusChr   = eugene.erneg.regions %>% filter(locus == i) %>% select(V1)
minPos     = eugene.erneg.regions %>% filter(locus == i) %>% select(V2) %>% as.numeric
maxPos     = eugene.erneg.regions %>% filter(locus == i) %>% select(V3) %>% as.numeric
eugeneSNPs = eugeneKnown %>% filter(locus == i) %>% select(V3) %>% mutate(eugene = "EUGENE_SNP")
candidates = bcacFMdata %>% filter(locus == sub("BCAC_FM_", "", i), type != "excluded_SNP") %>% select(position, type)

assocData  = read.table(paste(homeDir, "BCAC_summaryStats/erneg_position.space.pvalue_byChromosome/", locusChr[[1]], sep=""), header=F, )

graphData  = assocData %>% filter(V1 > minPos & V1 < maxPos)
graphData  = left_join(graphData, eugeneSNPs, by = c("V1" = "V3"))
graphData  = left_join(graphData, candidates, by = c("V1" = "position"))
graphData$type[is.na(graphData$type)] = "excluded_SNP"
graphData$eugene[is.na(graphData$eugene)] = "nonEUGENE_SNP"

# plot - SNP type in colour, EUGENE SNP is a triangle
ggplot(graphData) + 
  geom_point(aes(V1, -log10(V2), colour = type)) +
  geom_point(data = graphData %>% filter(eugene == "EUGENE_SNP"), aes(V1, -log10(V2)), shape = 2, size = 2, colour = 'red', stroke = 1) + 
  xlab(paste(locusChr[[1]]," position", sep='')) +
  ylab("-log10 P") +
  geom_hline(yintercept = -log10(5E-8), colour = 'red') + 
  scale_colour_manual(values = c("candidate_SNP"                = "blue", 
                                 "weaker_secondary_signals"     = "green", 
                                 "bcac_cimba_secondary_signals" = "violet",
                                 "cimba_secondary_signals"      = "orange",
                                 " "                            = "black"), 
                      na.value = "black") +
  ggtitle(assocLocus)

ggsave(paste(homeDir, "BCAC_FM_includesUCSCandEUGENE/manhattans/erneg/", assocLocus, "_ERnegRisk", ".png", sep=""))

}

lapply(unique(eugene.erneg.regions$locus), PlotManhattanERneg)





