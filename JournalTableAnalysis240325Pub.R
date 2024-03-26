##Combining Code for Data processing and Mapping##
library(miplicorn)
library(sf)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(tidyr)
library(scales)#only using with scale in radius (not nessary if not resizing pies)
#load required libraries for base maps
#library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
#library(ggrepel)
library(scatterpie)

#Analyze All Countries Separately because some HF have the same name and combining them makes it difficult to use dplyr functions

RW_ref <- "/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/RWAA2024/reference_AA_table.csv"
RW_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/RWAA2024/alternate_AA_table.csv"
RW_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/RWAA2024/coverage_AA_table.csv"
dataRW <- read_tbl_ref_alt_cov(RW_ref,RW_alt,RW_cov)
dataRW[dataRW<0]<- 0
manual_RW<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsRW.csv")
coordRW <-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesRWRepel.csv")
#coordRW <-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesRW.csv")
dataRW <- left_join(dataRW,left_join(manual_RW,coordRW))

DRC_ref <- "/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/DRCAAfiles/reference_AA_table.csv"
DRC_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/DRCAAfiles/alternate_AA_table.csv"
DRC_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/DRCAAfiles/coverage_AA_table.csv"
dataDRC <- read_tbl_ref_alt_cov(DRC_ref,DRC_alt,DRC_cov) 
dataDRC[dataDRC<0]<-0
manual_DRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsDRC.csv")
coordDRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesDRCRepel.csv")
#coordDRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesDRC.csv")
dataDRC <- left_join(dataDRC,left_join(manual_DRC,coordDRC))

UG_ref <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/UGAAfiles/reference_AA_table.csv"
UG_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/UGAAfiles/alternate_AA_table.csv"
UG_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/UGAAfiles/coverage_AA_table.csv"
dataUG <-read_tbl_ref_alt_cov(UG_ref,UG_alt,UG_cov)
dataUG[dataUG<0]<-0
manual_UG<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsUG.csv")
coordUG<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesUG.csv")
dataUG <- left_join(dataUG,left_join(manual_UG,coordUG))

TZ_ref <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/TZAAfiles/reference_AA_table.csv"
TZ_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/TZAAfiles/alternate_AA_table.csv"
TZ_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/TZAAfiles/coverage_AA_table.csv"
dataTZ <-read_tbl_ref_alt_cov(TZ_ref,TZ_alt,TZ_cov)
dataTZ[dataTZ<0]<-0
manual_TZ<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsTZ.csv")
coordTZ<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesTZ.csv")
dataTZ <- left_join(dataTZ,left_join(manual_TZ,coordTZ))

UMIthreshold <- function(data,n){
  #select significant UMIs/Coverage
  UMIs <- dplyr::filter(data, coverage > n | alt_umi_count > n)
  UMIs <- dplyr::group_by(UMIs,gene)
  #filter significant UMIs, calculate weighted frequency
  UMIs <- UMIs %>% mutate(totalUMI = ref_umi_count + alt_umi_count) %>% 
    mutate(frequency = alt_umi_count/totalUMI) %>% 
    mutate(weightedfreq = frequency*number) %>%
    drop_na(name)
  return (UMIs)
}
countryfreqforselectUMIoverN <- function(data,coord,select,n) {
  cutoff<-UMIthreshold(data,n) 
  #Finding Overall HF Frequencies
  sitegrouping <- cutoff %>% group_by(mutation_name, country, name) %>%add_count(Pool_Name)%>% reframe(number,totalUMI, weightedfreq,n)
  sitefrequency <- sitegrouping %>%group_by(mutation_name, country, name) %>%
    summarise_at(c("number", "totalUMI", "weightedfreq", "n"), ~ sum(., is.na(.), 0)) %>% 
    mutate(sitefreq = weightedfreq/number) %>%
    mutate(sitefreqWholeNum = sitefreq * 100) %>% rename(poolNumber = n)
  #select subset for managing graphing, combine with latitude longitude info, spread by the mutation to make multiple columns - innerjoin so that same HF names are not repeated
  forgraphing <- inner_join(sitefrequency, coord)
  #select the mutations for JID graphing
  forgraphing <- filter(forgraphing, grepl(select, mutation_name)) %>%
    drop_na(name) %>% 
    select(-totalUMI)%>% select(-weightedfreq) %>% select(-sitefreq) 
  #spread to export for viewing
  return (forgraphing)
}

admin10_rw <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/rwandaregions/rwa_adm1_2006_NISR_WGS1984_20181002.shp")
admin10_tz <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/Districts/Tanzania_District_wgs84.shp")
admin110_tz <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/Regions/Tanzania_Region_wgs84.shp")
admin110_cd <-st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/drcTerritory/cod_admbnda_adm2_rgc_20170711.shp")
admin110_ug <-st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/ugTerritories/Uganda_Districts-2020---136-wgs84.shp")
admin110_rw <-st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/rwandaregions/rwa_adm2_2006_NISR_WGS1984_20181002.shp")
natlpark <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/EastAfricaParks.kml")

k13 = "k13"
altDR = "dhps-Ala437Gly|dhps-Lys540Glu|dhps-Ala581Gly|dhfr-ts-Ile164Leu|dhfr-ts-Asn51Ile|dhfr-ts-Cys59Arg|dhfr-ts-Ser108Asn|mdr1-Asn86Tyr|mdr1-Tyr184Phe|mdr1-Asp1246Tyr|mdr1-Phe938Tyr|mdr1-Asn1042Asp|crt-Lys76Thr|crt-Arg371Ile|crt-Met74Ile|crt-Gln271Glu"

RWHFfreqs<-countryfreqforselectUMIoverN(dataRW,coordRW,altDR,10)
DRCHFfreqs <- countryfreqforselectUMIoverN(dataDRC,coordDRC,altDR,10)
TZHFfreqs <- countryfreqforselectUMIoverN(dataTZ, coordTZ,altDR,10)
UGHFfreqs <- countryfreqforselectUMIoverN(dataUG, coordUG,altDR,10)
DRmutations <- bind_rows(RWHFfreqs, TZHFfreqs, UGHFfreqs, DRCHFfreqs)

K13rw<-countryfreqforselectUMIoverN(dataRW,coordRW,k13,10)
K13tz<-countryfreqforselectUMIoverN(dataTZ,coordTZ,k13,10)
K13drc<-countryfreqforselectUMIoverN(dataDRC,coordDRC,k13,10)
K13ug<-countryfreqforselectUMIoverN(dataUG,coordUG,k13,10)
k13mutations <-bind_rows(K13drc,K13rw,K13tz,K13ug)

##MAKING JOURNAL TABLES####
#cleaning raw data tables
pubDataDRC <- dataDRC %>% select(!c(gene_id,gene,exonic_func,aa_change,targeted))
pubDataRW  <- dataRW %>% select(!c(gene_id,gene,exonic_func,aa_change,targeted)) %>% filter(!grepl("NTC|3d7|NTP|7g8",sample))
pubDataTZ <- dataTZ %>% select(!c(gene_id,gene,exonic_func,aa_change,targeted))
pubDataUG <- dataUG %>% select(!c(gene_id,gene,exonic_func,aa_change,targeted)) %>% filter(!grepl("NTC|3d7|NTP|7g8",sample))
write.csv(pubDataDRC,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/pubDataDRC.csv", row.names = TRUE)
write.csv(pubDataRW,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/pubDataRW.csv", row.names = TRUE)
write.csv(pubDataTZ,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/pubDataTZ.csv", row.names = TRUE)
write.csv(pubDataUG,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/pubDataUG.csv", row.names = TRUE)

#General Stats (Not in supplemental table)
totaln <- sum(manual_TZ$number, manual_DRC$number,manual_RW$number,manual_UG$number)
print(totaln)
npools <- sum(length(manual_TZ$number),length(manual_DRC$number),length(manual_RW$number),length(manual_UG$number))
avgpooln <- bind_rows(manual_TZ,manual_DRC,manual_RW,manual_UG) %>% mutate(avgerage = mean(number))
nsites <- sum(length(unique(manual_TZ$name)),length(unique(manual_DRC$name)),length(unique(manual_RW$name)),length(unique(manual_UG$name)))
countPools <- selectedMutations %>% group_by(mutation_name) %>%count(presence = (sitefreqWholeNum!=0))

#Health Center Averages
findAvgsk13 <- k13mutations %>% group_by(mutation_name, country,name) %>% summarise(avg= mean(sitefreqWholeNum))
samplesk13N <- k13mutations %>% group_by(mutation_name,country,name,poolNumber) %>% summarise(samples = sum(number))
mergednamesk13 <- left_join(findAvgsk13,samplesk13N)

findAvgsDR <- DRmutations %>% group_by(mutation_name, country,name) %>% summarise(avg= mean(sitefreqWholeNum))
samplesDRN <- DRmutations %>% group_by(mutation_name,country,name,poolNumber) %>% summarise(samples = sum(number))
mergednamesDR <- left_join(findAvgsDR,samplesDRN)

health_center_avgs <- bind_rows(mergednamesk13, mergednamesDR)
write.csv(health_center_avgs,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/AveragesByHealthCenters.csv", row.names = TRUE)

#Country HC Averages
findAvgsk13 <- k13mutations %>% group_by(mutation_name, country) %>% summarise(avg= mean(sitefreqWholeNum))
samplesk13N <- k13mutations %>% group_by(mutation_name,country) %>% summarise(samples = sum(number), pools = sum(poolNumber))
mergednamesk13 <- left_join(findAvgsk13,samplesk13N)

findAvgsDR <- DRmutations %>% group_by(mutation_name, country) %>% summarise(avg= mean(sitefreqWholeNum))
samplesDRN <- DRmutations %>% group_by(mutation_name,country) %>% summarise(samples = sum(number),pools=sum(poolNumber))
mergednamesDR <- left_join(findAvgsDR,samplesDRN)

country_avgs <- bind_rows(mergednamesk13, mergednamesDR)
write.csv(country_avgs,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/AveragesByCountry.csv", row.names = TRUE)

#Site Presence Counts by Country
countSitesk13 <- k13mutations %>% group_by(mutation_name,country) %>%count(presence = (sitefreqWholeNum!=0)) %>% rename(sites=n)
countSitesDR <- DRmutations %>% group_by(mutation_name,country)%>% count(presence = (sitefreqWholeNum!=0)) %>%rename(sites =n)
positiveHCcount <- bind_rows(countSitesk13, countSitesDR)
write.csv(positiveHCcount,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/FinalSuppXsl/PresenceCountSitesByCountry.csv", row.names = TRUE)

#Counting Pools with X (not used in supplemental)
RWgreater10 <- UMIthreshold(dataRW,10)
TZgreater10 <- UMIthreshold(dataTZ,10)
DRCgreater10 <-UMIthreshold(dataDRC,10)
UGgreater10 <- UMIthreshold(dataUG,10)
poolsseqd<-sum(length(unique(RWgreater10$sample)),length(unique(TZgreater10$sample)),length(unique(UGgreater10$sample)),length(unique(DRCgreater10$sample)))
alldataUMIfilter <- bind_rows(RWgreater10,TZgreater10,DRCgreater10,UGgreater10)
posk13Pools <- alldataUMIfilter%>% filter(grepl("k13",gene)) %>% group_by(mutation_name,country)%>% filter(alt_umi_count > 0)%>%summarise(posPools = length(unique(sample)))
write.csv(posk13Pools,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/2024Analysis01/NumberOfPositivePools.csv", row.names = TRUE)
posDRpools <- alldataUMIfilter%>% filter(grepl("dhps|mdr1|crt|dhfr",gene)) %>% group_by(mutation_name,country)%>% filter(alt_umi_count > 0)%>%summarise(posPools = length(unique(sample)))
write.csv(posDRpools,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/NumberOfPositiveDRPoolsByCountry.csv", row.names = TRUE)

#Calculate how many reads are involved in the output data (metric in discussion)
Repool_RW<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/repoolRW.csv")
Repool_DRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/repoolDRC.csv")
Repool_TZ<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/repoolTZ.csv")
Repool_UG<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/repoolUG.csv")
allRepool = bind_rows(Repool_DRC,bind_rows(Repool_RW, bind_rows(Repool_TZ,Repool_UG))) %>% filter(!grepl("3d7|NTC|NTP|7g8",Sample.ID)) %>%
  filter(Sample.ID%in%pools$sample)

##finding RW province averages (used in discussion)
K13rwprov <- K13rw %>% ungroup()%>% group_by(province,mutation_name) %>% summarize(provAvg = mean(sitefreqWholeNum))
K13drcprov <- K13drc %>% ungroup()%>% group_by(territory,mutation_name) %>% summarize(provAvg = mean(sitefreqWholeNum))

#permutation testing setup
allRawData <- bind_rows(dataRW, bind_rows(dataDRC,bind_rows(dataTZ,dataUG))) %>% filter(coverage > 10 | alt_umi_count >10)
pools <- as_tibble(allRawData$sample) %>% distinct()%>% filter(!grepl("3d7|NTC|NTP|7g8",value)) %>% rename(sample = value)
  
k13pools <- allRawData %>% mutate(WHO = if_else(grepl("k13-Arg561His|k13-Ala675Val|574Leu|469Tyr|Phe446Ile|458Tyr|476Ile|493His|539Thr|543Thr|553Leu|580Tyr|622Ile", mutation_name)>0, 2,if_else(grepl("469Phe|k13-Gly449Ala|441Leu|481Val|515Lys|527His|537Ile|538Val|568Gly", mutation_name)>0,1,0)))
valid <- filter(k13pools,WHO == 2) %>% mutate(validated = if_else(alt_umi_count >0, 1, 0)) %>% select(sample,validated) %>% group_by(sample) %>% summarise(whoValid = sum(validated))
candidate <- filter(k13pools,WHO==1) %>% mutate(cand449 = if_else(alt_umi_count > 0, if_else(grepl("449A",mutation_name),1,0),0)) %>% mutate(cand469 = if_else(alt_umi_count > 0, if_else(grepl("469P",mutation_name),1,0),0)) %>%
  mutate(cand441 = if_else(alt_umi_count > 0, if_else(grepl("441L",mutation_name),1,0),0)) %>%
  mutate(cand568 = if_else(alt_umi_count > 0, if_else(grepl("568G",mutation_name),1,0),0)) %>% group_by(sample)%>%
  summarise(cand= sum(cand449+cand469+cand441+cand568))
candidatebinary <- candidate %>%
  select(sample,cand449,cand469,cand441,cand568)%>%
  distinct() %>% group_by(sample)%>%
  summarise(count449 = sum(cand449), count469 =sum(cand469),count441 = sum(cand441),count568 = sum(cand568))

permutationdata <- left_join(pools, left_join(valid,left_join(candidatebinary,candidate))) %>% replace(is.na(.),0) %>% mutate(validbinary = if_else(whoValid >0,1,0)) %>% mutate(candbin = if_else(cand >0,1,0))
write_csv(permutationdata,"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/PermutationCalculations.csv")
