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
#RW####
RW_ref <- "/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/RWAA2024/reference_AA_table.csv"
RW_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/RWAA2024/alternate_AA_table.csv"
RW_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/RWAA2024/coverage_AA_table.csv"
dataRW <- read_tbl_ref_alt_cov(RW_ref,RW_alt,RW_cov)
dataRW[dataRW<0]<- 0
manual_RW<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsRW.csv")
coordRW <-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesRWRepel.csv")
#coordRW <-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesRW.csv")
dataRW <- left_join(dataRW,left_join(manual_RW,coordRW))

#DRC####
DRC_ref <- "/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/DRCAAfiles/reference_AA_table.csv"
DRC_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/DRCAAfiles/alternate_AA_table.csv"
DRC_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/DRCAAfiles/coverage_AA_table.csv"
dataDRC <- read_tbl_ref_alt_cov(DRC_ref,DRC_alt,DRC_cov) 
dataDRC[dataDRC<0]<-0
manual_DRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsDRC.csv")
coordDRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesDRCRepel.csv")
#coordDRC<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesDRC.csv")
dataDRC <- left_join(dataDRC,left_join(manual_DRC,coordDRC))

#UG####
UG_ref <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/UGAAfiles/reference_AA_table.csv"
UG_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/UGAAfiles/alternate_AA_table.csv"
UG_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/UGAAfiles/coverage_AA_table.csv"
dataUG <-read_tbl_ref_alt_cov(UG_ref,UG_alt,UG_cov)
dataUG[dataUG<0]<-0
manual_UG<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsUG.csv")
coordUG<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesUG.csv")
dataUG <- left_join(dataUG,left_join(manual_UG,coordUG))

#TZ####
TZ_ref <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/TZAAfiles/reference_AA_table.csv"
TZ_alt <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/TZAAfiles/alternate_AA_table.csv"
TZ_cov <-"/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/newTargsAll/TZAAfiles/coverage_AA_table.csv"
dataTZ <-read_tbl_ref_alt_cov(TZ_ref,TZ_alt,TZ_cov)
dataTZ[dataTZ<0]<-0
manual_TZ<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/manualDataPoolsTZ.csv")
coordTZ<-read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/coordinatesTZ.csv")
dataTZ <- left_join(dataTZ,left_join(manual_TZ,coordTZ))

#Functions for selecting data####
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
  sitegrouping <- cutoff %>% group_by(mutation_name, country, name) %>% reframe(number, totalUMI, weightedfreq)
  sitefrequency <- sitegrouping %>% group_by(mutation_name, country, name) %>%
    summarise_at(c("number", "totalUMI", "weightedfreq"), ~ sum(., is.na(.), 0)) %>% 
    mutate(sitefreq = weightedfreq/number) %>%
    mutate(sitefreqWholeNum = sitefreq * 100)
  #select subset for managing graphing, combine with latitude longitude info, spread by the mutation to make multiple columns - innerjoin so that same HF names are not repeated
  forgraphing <- inner_join(sitefrequency, coord)
  #select the mutations for JID graphing
  forgraphing <- filter(forgraphing, grepl(select, mutation_name)) %>%
    drop_na(name) %>% 
    select(-totalUMI)%>% select(-weightedfreq) %>% select(-sitefreq) 
  return (forgraphing)
  #to organize by district
  districtorg <- forgraphing %>% group_by(mutation_name, by)
  distrctFreq <- districtorg %>%
    summarise(distFreq=mean(sitefreqWholeNum))
  #return (distrctFreq)
}

#selecting focus data####
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

#Loading Country Boundaries####
#need to load admin110s first before you can run functions
#get admin shape file from naturalearth (admin10 is out of date)
rivers10 <- ne_download(scale = "large", type = 'rivers_lake_centerlines',
                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "large", type = 'lakes',
                       category = 'physical', returnclass = "sf")
# oceans10 <- ne_download(scale = "large", type = "coastline",
#                         category = 'physical', returnclass = "sf")
# sov110 <- ne_download(scale="large", type = "sovereignty",
#                       category = "cultural", returnclass = "sf")
# admin10 <- ne_download(scale="large", type = "admin_1_states_provinces_lines",
#                        category = "cultural", returnclass = "sf")

admin0_rw <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/rwandaregions/rwa_adm0_2006_NISR_WGS1984_20181002.shp")
admin10_rw <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/rwandaregions/rwa_adm1_2006_NISR_WGS1984_20181002.shp")
admin10_tz <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/Districts/Tanzania_District_wgs84.shp")
#admin110_tz <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/Regions/Tanzania_Region_wgs84.shp")
admin110_cd <-st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/drcTerritory/cod_admbnda_adm2_rgc_20170711.shp")
admin110_ug <-st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/ugTerritories/Uganda_Districts-2020---136-wgs84.shp")
admin110_rw <-st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/rwandaregions/rwa_adm2_2006_NISR_WGS1984_20181002.shp")
natlpark <- st_read("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/EastAfricaParks.kml")

pt.lim = data.frame(xlim = c(28.5, 32), ylim = c(1,-3))
pt.bbox <- st_bbox(c(xmin=pt.lim$xlim[1],
                     xmax=pt.lim$xlim[2],
                     ymin=pt.lim$ylim[1],
                     ymax=pt.lim$ylim[2]))
sf_use_s2(FALSE)
lakes10<- st_crop(lakes10, pt.bbox)
rivers10 <- st_crop(rivers10, pt.bbox)
admin110_cd <- st_crop(admin110_cd, pt.bbox)  
admin110_ug <- st_crop(admin110_ug, pt.bbox)
admin10_tz <- st_crop(admin10_tz, pt.bbox)

##PIE CHARTING####

#MAPPING base
{
  base <- ggplot()+
    #admin110RW off for faclab mal
    #Comment out sov110 because too many overlaps at DRC UG border and unneeded
    #geom_sf(data=sov110, color='grey80', size=20, alpha = 0.2, stroke=2.5) +
    geom_sf(data=admin10_tz, color="grey80", size= 0.6, alpha = 0.1) +
    geom_sf(data=admin110_cd, color="grey80", size= 0.6, alpha = 0.1) +
    geom_sf(data=admin110_ug, color="grey80", size= 0.6, alpha = 0.1) +
    #comment out national parks execept when graphing the health facilites and the co-occurances
    geom_sf(data=natlpark, color="palegreen", fill="green2", size=1,alpha =0.1) +
    geom_sf(data=rivers10, color="cyan4", size=0.5, alpha=0.5) +
    #geom_sf(data=admin110_rw, color="grey60", size= 0.6, alpha = 0.1) +
    geom_sf(data=admin10_rw, color="grey30", size= 0.6, alpha = 0.1) +
    geom_sf(data=lakes10, color="grey40", fill ="lightblue", size= 0.8) +
    geom_sf(data=admin0_rw, color= "black", size=0.6, linewidth =0.5, alpha=0.1)+
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.01, "in"), pad_y = unit(0.18, "in"),
                           style = north_arrow_fancy_orienteering)+
    annotate("text", x = 30, y = -1.1, label = "Uganda", 
             color="grey60", size=6 , fontface="italic") +
    annotate("text", x = 30.2, y = -2.6, label = "Burundi", 
             color="grey60", size=6 , fontface="italic") +
    annotate("text", x = 28.9, y = -1.5, label = "DRC", 
             color="grey60", size=6 , fontface="italic") +
    annotate("text", x = 30.8, y = -1.3, label = "Tanzania",
             color="grey60", size=6, fontface="italic")
}
{
  facilities <- base + 
    annotate("text", x = 29.7, y = -1.8, label = "Rwanda", 
             color="grey60", size=5 , fontface="italic") +
    geom_point(data = health_facilities, aes(x = lon, y = lat), size = 1.5,
             shape = 21, color= "grey40", fill = "cyan2", stroke = 1,
             alpha =0.8) +
    geom_text_repel(data = health_facilities, aes(x =lon, y= lat+.035, label = Health_Facilities),
                    point.size = NA,
                    size = 4,
                    color="grey20", fontface="bold",min.segment.length = 0.22, max.overlaps = 10) +
    #coord_sf(xlim = c(28.6, 31.6), ylim = c(0.3,-2.8), expand = TRUE) +
    #UG 2024 map coord_sf(xlim = c(29, 31.8), ylim = c(0,-2), expand = TRUE) +
    coord_sf(xlim = c(28.9, 30.85), ylim = c(-1.11,-2.81), expand = TRUE) +
    theme_void()+
    theme(legend.title=element_blank())+
    theme(legend.text=element_text(size=12))
  
  facilities
  ggsave("/media/nwernsma/ShareVolume/GEM_Rwanda_Mapping/facilities_mapped.svg", facilities, dpi=600, width=10, height=8.8)
}

#kelchpie Calculation
{
kelchpie <-filter(k13mutations, grepl("Arg561His|k13-Ala675Val|Gly449Ala|Pro441Leu|Cys469Phe", mutation_name))%>% 
  select(-c(territory,province,district,region,number))%>%
  spread(mutation_name, sitefreqWholeNum)%>%
  rename("R561H"="k13-Arg561His") %>%
  rename("A675V"="k13-Ala675Val")%>% 
  rename("G449A"="k13-Gly449Ala") %>%
  rename("P441L"="k13-Pro441Leu") %>%
  rename("C469F"="k13-Cys469Phe") %>%
  mutate(totalMT=R561H+A675V+G449A+P441L+C469F)%>%
  mutate(WT=100-totalMT)

#scale radius so that the higher R561H pies are layered on top
kelchpie$radius <- (rescale(kelchpie$R561H, to = c(3.01,3.0))/60)

#had issues with geom_sf so I rewrite and rename to avoid
write.csv(kelchpie,"/home/nwernsma/Documents/kelchpiewithRad.csv")
newkelchcsv<- read.csv("/home/nwernsma/Documents/kelchpiewithRad.csv")
}
#combined k13 pie simple
{
  #adjust the tanzania label in base MAP above when charted expanded points**
  # annotate("text", x = 31.1, y = -1.9, label = "Tanzania", 
  #          color="grey60", size=6 , fontface="italic")
  
  colorsline<- c("grey80","firebrick","cyan3","darkgoldenrod1","mediumpurple1","deeppink2")
  colorsline<- setNames(colorsline,c("WT","R561H","A675V","G449A","P441L","C469F"))
  
  multipie <- base +
    geom_scatterpie(aes(x=lon, y=lat, r=radius), 
                    data = newkelchcsv,
                    cols = c("R561H","A675V","G449A","P441L","C469F","WT"), 
                    sorted_by_radius = TRUE,
                    color = "grey20",
                    alpha=.8)+
    scale_fill_manual(values = colorsline)+
    #expanded coordinate set
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    # theme(legend.position = "none")
  #swap in once to get legend elements
    theme(legend.title=element_blank())+
    theme(legend.text=element_text(size=11),
          legend.position = c(.95,0),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,2,2,2),
          legend.box.background = element_rect(fill="white"))
  
  multipie
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240423kelchCombinedPiesExpandedRepelBBOX.pdf",multipie, dpi=600, width=11.5, height=8.8)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240423kelchCombinedPiesExpandedRepelBBOX.jpg",multipie, dpi=600, width=11.5, height=8.8)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240423kelchCombinedPiesExpandedRepelBBOX.svg",multipie, dpi=600, width=11.5, height=8.8)
}
{#DRC border
  
  multipie <-  ggplot() +
    geom_sf(data=admin10_tz, color="grey80", size= 0.6, alpha = 0.1) +
    geom_sf(data=admin110_cd, color="grey80", size= 0.6, alpha = 0.1) +
    geom_sf(data=admin110_ug, color="grey80", size= 0.6, alpha = 0.1) +
    geom_sf(data=rivers10, color="cyan4", size=0.5, alpha=0.5) +
    geom_sf(data=admin110_rw, color="grey60", size= 0.6, alpha = 0.1) +
    geom_sf(data=admin10_rw, color="grey60", size= 0.6, alpha = 0.1) +
    geom_sf(data=lakes10, color="grey40", fill ="lightblue", size= 0.8) +
    geom_sf(data=admin0_rw, color= "black", size=0.6, linewidth =.6,alpha=0.1)+
    #annotate("text", x = 29.2, y = -1.4, label = "DRC", 
    #         color="grey60", size=6 , fontface="italic") +
    #annotate("text", x = 29.45, y = -1.75, label = "Rwanda", 
    #         color="grey60", size=6 , fontface="italic") +
    geom_scatterpie(aes(x=lon, y=lat, r=radius), 
                    data = newkelchcsv,
                    cols = c("R561H","A675V","G449A","P441L","C469F","WT"), 
                    sorted_by_radius = TRUE,
                    color = "grey20",
                    alpha=.8)+
    scale_fill_manual(values = colorsline)+
    #DRC border
    coord_sf(xlim = c(28.97, 29.50), ylim = c(-1.45,-2.1), expand = TRUE) +
    theme_void()+
    theme(legend.position = "none")

  multipie
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240423kelchDRCBordSM.pdf",multipie, dpi=600, width=2, height=3)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240423kelchDRCBordSM.jpg",multipie, dpi=600, width=2, height=3)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240423kelchDRCBordSMallFresh.svg",multipie, dpi=600, width=2, height=3)
}

#justR561H
{R561H <- base +
  geom_scatterpie(aes(x=lon, y=lat, r = (radius/50)),
                  data = newkelchcsv,
                  cols = c("R561H", "WT"),
                  sorted_by_radius= TRUE,
                  color = "grey20",
                  legend_name = "Pfkelch13",
                  alpha=.8)+
  scale_fill_manual(values = c("firebrick","grey80")) +
  coord_sf(xlim = c(28.7, 31), ylim = c(-1,-2.9), expand = TRUE) +
  #geom_scatterpie_legend((figure2_k13561$radius/50)[1:18], x=28.8, y=-1.5, labeller = function(x) (x*100)) +
  #geom_scatterpie_legend((figure2_k13561$radius/50), x=28.8, y=-1.5)+
  theme_void()+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12))

ggsave("/home/nwernsma/Documents/230419RWMaps/561Hpies.png", R561H, dpi=600, width=7, height=5)
#Swap in use a set to small export and use to make inset
R561Hsmall <- R561H + 
  coord_sf(xlim = c(29.9, 30.4), ylim = c(-1.6,-2.3), expand = TRUE)
ggsave("/home/nwernsma/Documents/230419RWMaps/561Hpiesinset.png", R561Hsmall, dpi=600, width=7, height=5)
}
#675 & 561H border transmission
{#original HC locations
health_facilities <- read.csv("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/NWYMappingDemoFiles/dataRW23F_fac.csv")
#combined kelch pie chart 
multipie <- base +
  annotate("text", x = 29.2, y = -1.4, label = "DRC", 
           color="grey60", size=6 , fontface="italic") +
  annotate("text", x = 29.45, y = -1.75, label = "Rwanda", 
           color="grey60", size=6 , fontface="italic") +
  geom_sf(data=admin10_rw, color= "black", size=5, alpha=0.1)+
  geom_scatterpie(aes(x=lon, y=lat, r=(radius)), 
                  data = newkelchcsv,
                  cols = c("R561H","A675V","WT"), 
                  #sorted_by_radius = TRUE,
                  color = "grey20",
                  alpha=.8)+
  # geom_point(data = health_facilities, aes(x = lon, y = lat), size = 1.5,
  #            shape = 21, color= "grey40", fill = "black", stroke = 1,
  #            alpha =0.8) +
  scale_fill_manual(values =  c("firebrick","cyan3","grey80"))+
  #central
  #coord_sf(xlim = c(28.8, 31.3), ylim = c(-1.05,-2.85), expand = TRUE) +
  #DRC border
  coord_sf(xlim = c(28.97, 29.55), ylim = c(-1.35,-2.1), expand = TRUE) +
  theme_void()+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=11), 
        legend.position = c(.1,.9),
        #legend.justification = c("left","top"), 
        legend.margin = margin(1,2,2,2),
        legend.box.background = element_rect(fill="white"))

multipie
ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/561+675border.jpg",multipie, dpi=600, width=11.5, height=8.8)
ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/561+675border.pdf",multipie, dpi=600, width=11.5, height=8.8)
}

#Co-Ocurrance map
{k13occurance <- k13mutations %>% mutate(WHO = if_else(grepl("k13-Arg561His|k13-Ala675Val|574Leu|469Tyr|Phe446Ile|458Tyr|476Ile|493His|539Thr|543Thr|553Leu|580Tyr|622Ile", mutation_name)>0, 2,if_else(grepl("469Phe|k13-Gly449Ala|441Leu|481Val|515Lys|527His|537Ile|538Val|568Gly", mutation_name)>0,1,0)))
candidate<- filter(k13occurance,WHO==1) %>% filter(sitefreqWholeNum >0)%>% group_by(country,name)%>% count(name) %>% rename(cand = n) 
validate <- filter(k13occurance,WHO==2) %>% filter(sitefreqWholeNum >0)%>% group_by(country, name)%>% count(name) %>% rename(valid=n)
no_sites <- filter(k13occurance, WHO == 0)%>% ungroup()%>%select(-mutation_name) %>% group_by(country,name)%>% filter(!name %in% validate$name)%>% filter(!name %in% candidate$name)%>% count(name) %>% select(-n)
k13Cooccur <- right_join(candidate,validate) %>% bind_rows(no_sites)
k13Cooccur[is.na(k13Cooccur)] <-0
k13Cooccur <- k13Cooccur %>% mutate(mutants = if_else((cand >0 & valid >0),"Validated&Candidate",if_else((cand == 0 & valid > 0),"Validated Only","None")))
k13Cooccur <- left_join(k13Cooccur,full_join(coordRW,full_join(coordDRC, full_join(coordTZ,coordUG)))) %>% mutate(total = valid+cand) %>% arrange(mutants)
k13Cooccur$mutants <- factor(k13Cooccur$mutants, levels= c("Validated&Candidate","Validated Only","None"))
}
#VennDiagramCoOccur
#mapping CoOccur
{
  shapes = c(15,0,1)
  CoOccur <- base +
    geom_point(data = k13Cooccur, aes(x = lon, y = lat, shape = mutants), size = 2.5,
               stroke = 1, color = "grey20",
               alpha =0.8) +
    scale_shape_manual(values=shapes, name="WHO mutation")+
    #scale_color_manual(values = colors)+
    #coord_sf(xlim = c(28.6, 31.6), ylim = c(0.3,-2.8), expand = TRUE) +
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    #UG 2024 map coord_sf(xlim = c(29, 31.8), ylim = c(0,-2), expand = TRUE) +
    theme_void()+
    # theme(legend.title=element_blank())+
    # theme(legend.text=element_text(size=12))
    theme(legend.text=element_text(size=11),
          legend.title.align = 0.5,
          legend.position = c(.98,0.02),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,1,1,1),
          legend.box.background = element_rect(fill="white"))
  
  CoOccur
  
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/coOccur.pdf", CoOccur, dpi=600, width=6.34, height=7.64)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/coOccur.jpg", CoOccur, dpi=600, width=6.34, height=7.64)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/coOccurJIDImage.svg", CoOccur, dpi=600, width=6.34, height=7.64)
}

#mdr1 pie calculations
{
  mdr1pie <-filter(DRmutations, grepl("184",mutation_name)) %>% select(-number) %>% spread(mutation_name,sitefreqWholeNum) %>%rename("Y184F"="mdr1-Tyr184Phe")%>%
    mutate(WT=100-Y184F)

  mdr1pie$radius <- (rescale(mdr1pie$Y184F, to = c(3.01,3.0))/60)
}
#mdr1 184
{
  mdr184 <- base +
    geom_scatterpie(aes(x=lon, y=lat, r = radius),
                    data = mdr1pie,
                    cols = c("Y184F", "WT"),
                    sorted_by_radius= TRUE,
                    color = "grey20",
                    legend_name = "Pfkelch13",
                    alpha=.8)+
    scale_fill_manual(values = c("grey80","lightgoldenrod2"), name="mdr1") +
    #central
    #coord_sf(xlim = c(28.8, 31), ylim = c(-1.08,-2.85), expand = TRUE) +
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    #theme(legend.position = "none")
    theme(legend.title=element_blank())+
    theme(legend.text=element_text(size=11),
          legend.position = c(.98,0.02),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,2,2,2),
          legend.box.background = element_rect(fill="white"))
  
  mdr184
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/JID_Review_Edits/mdr1_184Pie240820.svg", mdr184, dpi=600, width=11.5, height=8.8)
}

#mdr1 N86Y 
{
  mdr1_86pie <- filter(DRmutations,grepl("86",mutation_name)) %>% select(-number) %>% spread(mutation_name, sitefreqWholeNum) %>% rename("N86Y" = "mdr1-Asn86Tyr") %>%
    mutate(WT= 100-N86Y)
  mdr1_86pie$radius <- (rescale(mdr1_86pie$N86Y,to=c(3.01,3.0))/60)
}
#mdr1 86
{
  mdr1_86 <- base +
    geom_scatterpie(aes(x=lon, y=lat, r = radius),
                    data = mdr1_86pie,
                    cols = c("N86Y", "WT"),
                    sorted_by_radius= TRUE,
                    color = "grey20",
                    legend_name = "Pfkelch13",
                    alpha=.8)+
    scale_fill_manual(values = c("grey80","dodgerblue"), name="mdr1") +
    #central
    #coord_sf(xlim = c(28.8, 31), ylim = c(-1.08,-2.85), expand = TRUE) +
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    #theme(legend.position = "none")
    theme(legend.title=element_blank())+
    theme(legend.text=element_text(size=11),
          legend.position = c(.98,0.02),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,2,2,2),
          legend.box.background = element_rect(fill="white"))
  
  mdr1_86
  #ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/mdr184Expandedpie240327.pdf", mdr184, dpi=600, width=11.5, height=8.8)
  #ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/mdr184Expandedpie240327.jpg", mdr184, dpi=600, width=11.5, height=8.8)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/JID_Review_Edits/mdr1_86Expandedpie240814.svg", mdr1_86, dpi=600, width=11.5, height=8.8)
}

#mdr1 1246
{
  mdr1_1246pie <- filter(DRmutations,grepl("1246",mutation_name)) %>% select(-number) %>% spread(mutation_name, sitefreqWholeNum) %>% rename("D1246Y" = "mdr1-Asp1246Tyr") %>%
    mutate(WT= 100-D1246Y)
  mdr1_1246pie$radius <- (rescale(mdr1_1246pie$D1246Y,to=c(3.01,3.0))/60)
}
#mdr1 1246
{
  mdr1246 <- base +
    geom_scatterpie(aes(x=lon, y=lat, r = radius),
                    data = mdr1_1246pie,
                    cols = c("D1246Y", "WT"),
                    sorted_by_radius= TRUE,
                    color = "grey20",
                    legend_name = "Pfkelch13",
                    alpha=.8)+
    scale_fill_manual(values = c("grey80","palevioletred1"), name="mdr1") +
    #central
    #coord_sf(xlim = c(28.8, 31), ylim = c(-1.08,-2.85), expand = TRUE) +
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    #theme(legend.position = "none")
    theme(legend.title=element_blank())+
    theme(legend.text=element_text(size=11),
          legend.position = c(.98,0.02),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,2,2,2),
          legend.box.background = element_rect(fill="white"))
  
  mdr1246
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/JID_Review_Edits/mdr1_1246Pie240820.svg", mdr1246, dpi=600, width=11.5, height=8.8)
}

#crt pie calculations
{
  crtpie <- filter(DRmutations, grepl("76",mutation_name)) %>%select(-number) %>% spread(mutation_name,sitefreqWholeNum) %>% rename(L76T = "crt-Lys76Thr")%>%
    mutate(WT = 100-L76T)
  crtpie$radius <-(rescale(crtpie$L76T, to = c(3.01,3.0))/60)
}
#crt 76
{
  crt76 <- base +
    geom_scatterpie(aes(x=lon, y=lat, r = radius),
                    data = crtpie,
                    cols = c("L76T", "WT"),
                    sorted_by_radius= TRUE,
                    color = "grey20",
                    legend_name = "Pfkelch13",
                    alpha=.8)+
    scale_fill_manual(values = c("grey80","olivedrab2")) +
    #central
    #coord_sf(xlim = c(28.8, 31), ylim = c(-1.08,-2.85), expand = TRUE) +
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    #theme(legend.position = "none")
    theme(legend.title=element_blank())+
    theme(legend.text=element_text(size=11),
        legend.position = c(.98,0.02),
        legend.justification = c("right", "bottom"),
        legend.margin = margin(1,2,2,2),
        legend.box.background = element_rect(fill="white"))

  crt76
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/crt76Expandedpie240322.pdf", crt76, dpi=600, width=11.5, height=8.8)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/crt76Expandedpie240322.jpg", crt76, dpi=600, width=11.5, height=8.8)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/crt76Expandedpie240430.svg", crt76, dpi=600, width=11.5, height=8.8)
}

#dhpspie
{##Working on dhps multi mutant
  #Thinking:
  # if pietotal < 100, then charting A581G, real 540, WT437, wedge437
  # else: then charting, WT437, ratioed 581 aND 540
  # 
  # 
  # 581 freq =  437,540&581haplotype
  # 540 freq- 581 freq = 437&540 haplotype
  # 437 freq - 540freq = 437 haplotype (if below 0 = 0 )
  # 
  # scale the 437% to the 581 and 540% 
  #   
  #   
  #   
  #   581 haplotype = 581 frequency
  # 540 haplotype = 540 measured minus 581 frequency
  # 537 haplotype =  537 freq - 540 frequency
  # 
  dhps437 <- filter(DRmutations, grepl("dhps-Ala437Gly",mutation_name)) %>% select(-number) %>% spread(mutation_name, sitefreqWholeNum)%>% rename("A437G"="dhps-Ala437Gly")
  dhps540 <- filter(DRmutations, grepl("dhps-Lys540Glu",mutation_name)) %>% select(-number)%>% spread(mutation_name, sitefreqWholeNum)%>% rename("K540E"="dhps-Lys540Glu")
  dhps581 <- filter(DRmutations, grepl("dhps-Ala581Gly",mutation_name)) %>% spread(mutation_name, sitefreqWholeNum)%>% rename("A581G"="dhps-Ala581Gly")
  
  dhps540 <- dhps540 %>% select(name,K540E)
  dhps437 <- dhps437 %>% select(name,A437G)
  combineddhpsPie <- left_join(dhps581,left_join(dhps540,dhps437)) %>% mutate(WT437 = 100-A437G) %>% mutate(MTratio = (K540E/A581G)+1)%>% mutate(chart581 = A437G/MTratio) %>% mutate(chart540 = A437G - chart581)%>%mutate(real540 = K540E - A581G) %>% mutate(totalMT = real540+A581G)%>% mutate(pietotal=WT437+totalMT) %>% mutate(wedge437 = 100-pietotal)
  combineddhpsPie$wedge437[combineddhpsPie$wedge437<0] <- 0
  combineddhpsPie <- select(combineddhpsPie,-c(pietotal,A437G,totalMT,MTratio,K540E)) %>% rename(WT=WT437) %>% rename(A437G = wedge437)
  
  #combineddhpsPie <- mutate(combineddhpsPie, pietotal = ifelse(wedge437 > 0,wedge437+A581G+WT437+real540,wedge437+WT437+chart581+chart540))
  
  allthree <- filter(combineddhpsPie, A437G >0) %>% select(-c(chart581,chart540)) %>% rename(K540E = real540)
  justtwo <- filter(combineddhpsPie,A437G==0) %>% select(-c(real540,A581G)) %>% rename(K540E = chart540)%>% rename(A581G =chart581)
  
  finaldhpsdata <- bind_rows(allthree,justtwo) 
  avg581 <- finaldhpsdata %>% filter(grepl("Rwanda",country)) %>% summarise(mean(A581G)) #in discussion
  finaldhpsdata$radius <- 0.0505
}
#dhps combined pie chart
{
  colorsline<- c("springgreen3","slateblue1","chocolate1","grey80")
  colorsline<- setNames(colorsline,c("A581G","K540E","A437G","WT"))
  
  multipie <- base +
    geom_scatterpie(aes(x=lon, y=lat, r=radius), 
                    data = finaldhpsdata,
                    cols = c("A581G","K540E","A437G","WT"), 
                    sorted_by_radius = TRUE,
                    color = "grey20",
                    alpha=.8)+
    scale_fill_manual(values = colorsline)+
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    theme(legend.position = "none")
  #swap in once to get legend elements
  # theme(legend.title=element_blank())+
  # theme(legend.text=element_text(size=11),
  #     legend.position = c(.95,0.05),
  #     legend.justification = c("right", "bottom"),
  #     legend.margin = margin(1,2,2,2),
  #     legend.box.background = element_rect(fill="white"))
  multipie
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/JID_Review_Edits/dhpsTriple240819.svg", multipie, dpi=600, width=11.5, height=8.8)
}

#dhfr pie 
{
  dhfr108 <- filter(DRmutations, grepl("108",mutation_name)) %>% select(-number) %>% spread(mutation_name, sitefreqWholeNum) %>%rename(S108N="dhfr-ts-Ser108Asn")
  dhfr51 <- filter(DRmutations, grepl("51",mutation_name)) %>% select(-number)%>% spread(mutation_name, sitefreqWholeNum)%>%rename(N51I="dhfr-ts-Asn51Ile")
  dhfr59 <- filter(DRmutations, grepl("59",mutation_name)) %>% select(-number)%>% spread(mutation_name, sitefreqWholeNum)%>%rename(C59R="dhfr-ts-Cys59Arg")
  dhfr164 <- filter(DRmutations, grepl("164",mutation_name)) %>% select(-number) %>%spread(mutation_name, sitefreqWholeNum)%>%rename(I164L="dhfr-ts-Ile164Leu")
  
  dhfr108 <- dhfr108 %>% select(name,S108N)
  dhfr51 <- dhfr51 %>% select(name,N51I)
  dhfr59 <- dhfr59 %>% select(name,C59R)
  combineddhfrPie <- left_join(dhfr164,left_join(dhfr59,left_join(dhfr51,dhfr108))) %>% mutate(WT108 = 100-S108N) %>% mutate(MTratio = (C59R/I164L)+1)%>% mutate(chart164 = S108N/MTratio) %>% mutate(chart59 = S108N - chart164)%>%mutate(real59 = C59R - I164L) %>% mutate(totalMT = real59+I164L)%>% mutate(pietotal=WT108+totalMT) %>% mutate(wedge108 = 100-pietotal)
  combineddhfrPiealt <- left_join(dhfr164,left_join(dhfr59,left_join(dhfr51,dhfr108))) %>% mutate(WT51 = 100-N51I) %>% 
    mutate(MTratio = (C59R/I164L)+1)%>% 
    mutate(chart164 = N51I/MTratio) %>% 
    mutate(chart59 = N51I - chart164)%>%
    mutate(real59 = C59R - I164L) %>% 
    mutate(totalMT = real59+I164L)%>% 
    mutate(pietotal=WT51+totalMT) %>% 
    mutate(wedge51 = 100-pietotal)
  combineddhfrPiealt$wedge51[combineddhfrPiealt$wedge51<0] <- 0
  combineddhfrPiealt <- select(combineddhfrPiealt,-c(pietotal,N51I,totalMT,MTratio,C59R)) %>% rename(WT=WT51) %>% rename(N51I = wedge51)
  allthree <- filter(combineddhfrPiealt, N51I >0) %>% select(-c(chart59,chart164)) %>% rename(C59R = real59)
  justtwo <- filter(combineddhfrPiealt,N51I==0) %>% select(-c(real59,I164L)) %>% rename(C59R = chart59)%>% rename(I164L =chart164)
  
  
  #Triplemutant finding only not GREAT
  tripledhfrPiealt <- left_join(dhfr164,left_join(dhfr59,left_join(dhfr51,dhfr108))) %>% mutate(WT108 = 100-S108N) %>% 
    mutate(MTratio = (N51I/C59R)+1)%>% 
    mutate(chart59 = S108N/MTratio) %>% 
    mutate(chart51 = S108N - chart59)%>%
    mutate(real51 = N51I - C59R) %>% 
    mutate(totalMT = real51+C59R)%>% 
    mutate(pietotal=WT108+totalMT) %>% 
    mutate(wedge108 = 100-pietotal)
  tripledhfrPiealt$wedge108[tripledhfrPiealt$wedge108<0] <- 0
  tripledhfrPiealt <- select(tripledhfrPiealt,-c(pietotal,S108N,totalMT,MTratio,N51I)) %>% rename(WT=WT108) %>% rename(S108N = wedge108)
  allthree <- filter(tripledhfrPiealt, S108N >0) %>% select(-c(chart51,chart59)) %>% rename(N51I = real51)
  justtwo <- filter(tripledhfrPiealt,S108N==0) %>% select(-c(real51,C59R)) %>% rename(N51I =chart51)%>% rename(C59R = chart59)
  DHFRforquad <- bind_rows(allthree,justtwo)
  #combineddhpsPie <- mutate(combineddhpsPie, pietotal = ifelse(wedge437 > 0,wedge437+A581G+WT437+real540,wedge437+WT437+chart581+chart540))
  
  #consolidating to a triplicate then finding ratio of triple to quad mutant
  tripledhfrPiealt <- left_join(dhfr164,left_join(dhfr59,left_join(dhfr51,dhfr108))) %>% mutate(triple = C59R+S108N+N51I) %>% mutate(tripleMTavg = (triple/300)*100) %>%
    mutate(WT=100-tripleMTavg) %>%
    mutate(quadRatio = I164L/tripleMTavg) %>%
    mutate(quadscaled = (quadRatio*tripleMTavg)/100) %>%
    mutate(chart164 = quadscaled*tripleMTavg) %>%
    mutate(chartTriple = tripleMTavg-chart164) %>%
    mutate(totalpie = chart164+WT+chartTriple) %>%
    select(-c(I164L,C59R,S108N,N51I,triple,tripleMTavg,quadRatio,quadscaled)) %>% rename(I164L = chart164) %>% rename(C59Rtrip =chartTriple)
  finaldhfrdata <- tripledhfrPiealt
  
  
  finaldhfrdata <- bind_rows(allthree,justtwo)
  finaldhfrdata$radius <- 0.0505
}
#Imputed dhfr164 as 164, 59 and 51
{
  colorsline<- c("brown2","forestgreen","gold1","grey80")
  colorsline<- setNames(colorsline,c("I164L","C59R","N51I","WT"))
  
  multipie <- base +
    geom_scatterpie(aes(x=lon, y=lat, r=radius), 
                    data = finaldhfrdata,
                    cols = c("I164L","C59R","N51I","WT"), 
                    sorted_by_radius = TRUE,
                    color = "grey20",
                    alpha=.8)+
    scale_fill_manual(values = colorsline)+
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    #theme(legend.position = "none")
    #swap in once to get legend elements
    theme(legend.title=element_blank())+
    #labs(fill="dhfr")+
    theme(legend.text=element_text(size=11),
          legend.title.align = 0.5,
          legend.position = c(.95,0.05),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,2,2,2),
          legend.box.background = element_rect(fill="white"))
  multipie
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240228dhfrPiesRepel.pdf",multipie, dpi=600, width=11.52, height=8.802)
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/24FebManuscript/240228dhfrPiesRepel.jpg",multipie, dpi=600, width=11.5, height=8.8)
}
#Imputed dhfr164 as tripleMT AVG
{
  colorsline<- c("brown2","forestgreen","grey80")
  colorsline<- setNames(colorsline,c("I164L","C59Rtrip","WT"))
  
  multipie <- base +
    geom_scatterpie(aes(x=lon, y=lat, r=radius), 
                    data = finaldhfrdata,
                    cols = c("I164L","C59Rtrip","WT"), 
                    sorted_by_radius = TRUE,
                    color = "grey20",
                    alpha=.8)+
    scale_fill_manual(values = colorsline)+
    #expanded
    coord_sf(xlim = c(28.9, 31.3), ylim = c(0.1,-2.8), expand = TRUE) +
    theme_void()+
    #theme(legend.position = "none")
    #swap in once to get legend elements
    theme(legend.title=element_blank())+
    #labs(fill="dhfr")+
    theme(legend.text=element_text(size=11),
          legend.title.align = 0.5,
          legend.position = c(.95,0.05),
          legend.justification = c("right", "bottom"),
          legend.margin = margin(1,2,2,2),
          legend.box.background = element_rect(fill="white"))
  multipie
  ggsave("/home/nwernsma/Documents/FinalArtR22RWPoolsAnalysis/JID_Review_Edits/dhfrTripleQuadNoInkscape240809.svg", multipie, dpi=600, width=11.5, height=8.8)
}

