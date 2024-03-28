# packages
library(readr)
library(dplyr)
library(tidyr)
library(soiltexture)
library(ggplot2)
library(stringr)
library(knitr)
library(glue)
library(terra)
library(futile.logger)

# read data: layers
flog.info("Reading data...")
dlay <- read_csv("/lustre/arpa/bonafeg/data/geo/SoilType/data/HWSD_2.0/HWSD2_DB/csv/HWSD2_LAYERS.csv")

# prepare SMU data, managing categories not defined by salt-silt-clay percentage
flog.info("Arranging data...")
dlay %>%
  distinct(
    HWSD2_SMU_ID, SEQUENCE, SHARE, WRB2,
    LAYER, TOPDEP, BOTDEP, SAND, SILT, CLAY, ORG_CARBON) %>%
  mutate(
    NoTextureCode=ifelse(SILT<0,-SILT,0),
    ShareTexture=ifelse(SILT<0,NA,SHARE),
    ShareOC=ifelse(ORG_CARBON<0,NA,SHARE),
    SILT=ifelse(SILT<0,NA,SILT),
    SAND=ifelse(SAND<0,NA,SAND),
    CLAY=ifelse(CLAY<0,NA,CLAY),
    ORG_CARBON=ifelse(ORG_CARBON<0,NA,ORG_CARBON)
  ) %>%
  ungroup() -> dsmu

# check organic carbon
flog.info("Checking organic carbon...")
org_thr <- 25
ggplot(dsmu%>%
         filter(SEQUENCE==1,LAYER%in%c("D1"),SHARE==100,ORG_CARBON>0), 
       aes(x=WRB2=="HS", y=ORG_CARBON))+
  geom_boxplot(varwidth=T, outlier.size = 0.5)+
  scale_y_continuous(sec.axis = sec_axis(~.,breaks=org_thr))+
  xlab("is the SMU classified as Histosol?")+
  ylab("organic carbon content (% weight)")+
  theme_bw()+
  geom_hline(yintercept = org_thr, linetype="dashed")+
  ggtitle("HWSD 2.0", subtitle = "0-20 cm depth") -> p
ggsave(p, filename = "hwsd2_orgcarbon_thr.pdf", height=4, width=4)

# summarize SMU data: texture
flog.info("Summarizing soil data defined by texture...")
dsmu %>%
  group_by(
    HWSD2_SMU_ID, LAYER, TOPDEP, BOTDEP
    ) %>%
  summarize(
    SAND=sum(SAND*ShareTexture,na.rm=T),
    SILT=sum(SILT*ShareTexture,na.rm=T),
    CLAY=sum(CLAY*ShareTexture,na.rm=T),
    ORG_CARBON=sum(ORG_CARBON*ShareOC,na.rm=T),
    TotShareTexture=sum(ShareTexture,na.rm=T),
    TotShareOC=sum(ShareOC,na.rm=T),
    SAND=SAND/TotShareTexture,
    SILT=SILT/TotShareTexture,
    CLAY=CLAY/TotShareTexture, 
    ORG_CARBON=ORG_CARBON/TotShareOC, 
    .groups = "drop"
  ) -> dsmu_tx

# summarize SMU data: no-texture categories
flog.info("Summarizing soil data not defined by texture...")
dsmu %>%
  filter(NoTextureCode>0) %>%
  group_by(
    HWSD2_SMU_ID, LAYER, TOPDEP, BOTDEP, 
    NoTextureCode
  ) %>%
  summarize(
    ShareNoTexture=sum(SHARE,na.rm=T), 
    .groups = "drop"
  ) -> dsmu_notx

# soil layers correspondances
flog.info("Splitting layers...")
#dest_br <- c(0,10,40,100,200)
dest_br <- c(0,30,200)
dest_br <- sort(dest_br)
orig_br <- sort(unique(c(dsmu$TOPDEP,dsmu$BOTDEP)))
spli_br <- sort(unique(c(dest_br, orig_br)))
ns <- length(spli_br)
spli_mid <- (spli_br[-1]+spli_br[-ns])/2
s2o <- cut(spli_mid, orig_br, labels = F)
s2d <- cut(spli_mid, dest_br, labels = F)
tibble(
  dest_bot=dest_br[-1][s2d],
  dest_top=dest_br[-ns][s2d],
  orig_bot=orig_br[-1][s2o],
  orig_top=orig_br[-ns][s2o],
  spli_bot=spli_br[-1],
  spli_top=spli_br[-ns]) %>%
  mutate(
    dest_thick=dest_bot-dest_top,
    orig_thick=orig_bot-orig_top,
    spli_thick=spli_bot-spli_top,
    dest_fract=spli_thick/dest_thick,
    orig_fract=spli_thick/orig_thick) -> slc

# split soil layers: no-texture categories
full_join(slc, dsmu_notx, 
          by=join_by(orig_top == TOPDEP, orig_bot == BOTDEP),
          relationship="many-to-many") -> dsmu_notx

# split soil layers: texture data
full_join(slc, dsmu_tx, 
          by=join_by(orig_top == TOPDEP, orig_bot == BOTDEP),
          relationship="many-to-many") -> dsmu_tx

# aggregate soil layers to the destination layers (no-texture categories)
flog.info("Aggregating layers...")
dsmu_notx  %>%
  group_by(HWSD2_SMU_ID,dest_bot,dest_top,
           NoTextureCode) %>% 
  summarise(ShareNoTexture=sum(ShareNoTexture*dest_fract), .groups = "drop") %>%
  group_by(HWSD2_SMU_ID,dest_bot,dest_top) %>% 
  arrange(desc(ShareNoTexture)) %>%
  mutate(NoTextureCodeSequence=paste(NoTextureCode,collapse=","),
         TotShareNoTexture=sum(ShareNoTexture,na.rm=T)) %>%
  slice(1) %>%
  ungroup() %>%
  rename(ShareNoTexture_1st=ShareNoTexture,
         NoTextureCode_1st=NoTextureCode)-> dsmu_notx

# aggregate soil layers to the destination layers (texture data and OC)
dsmu_tx  %>%
  group_by(HWSD2_SMU_ID,dest_bot,dest_top) %>% 
  summarise(
    TotShareTexture=sum(TotShareTexture*dest_fract), 
    TotShareOC=sum(TotShareOC*dest_fract), 
    SAND=sum(SAND*dest_fract), 
    SILT=sum(SILT*dest_fract), 
    CLAY=sum(CLAY*dest_fract), 
    ORG_CARBON=sum(ORG_CARBON*dest_fract), 
    .groups = "drop") %>% 
  filter(TotShareTexture>0) -> dsmu_tx

# merge texture and no-texture
full_join(dsmu_tx, dsmu_notx) %>%
  mutate(
    TotShareTexture=ifelse(is.na(TotShareTexture),0,TotShareTexture),
    TotShareNoTexture=ifelse(is.na(TotShareNoTexture),0,TotShareNoTexture),
    NoTextureDescr_1st=case_match(
      NoTextureCode_1st,
      1 ~ "Water bodies",
      2 ~ "Land ice and glaciers",
      3 ~ "Rock outcrops",
      4 ~ "Dunes/shifting sands",
      5 ~ "Salt flats",
      6 ~ "Rocky sublayers",
      7 ~ "Rocky sublayers",
      8 ~ "No data",
      9 ~ "Other")) -> dsmu

# remap to 12 soil categories defined by texture (USDA)
# install.packages("https://cran.r-project.org/src/contrib/Archive/soiltexture/soiltexture_1.5.1.tar.gz")
flog.info("Remapping to WRF categories...")
css_ok <- !is.na(dsmu$SILT+dsmu$SAND+dsmu$CLAY)
dd <- dsmu[css_ok,]
UsdaClass <- TT.points.in.classes(tri.data = as.data.frame(dd), 
                                  class.sys   = "USDA.TT", collapse=",", PiC.type = "t")
dsmu$UsdaClass <- NA
dsmu$UsdaClass[css_ok] <- UsdaClass

# remap USDA classes to Noah classes
stc <- read_csv("data/SoilTypeCategories.csv")
dsmu %>%
  mutate(UsdaClass=str_split_i(UsdaClass,",",1)) %>%
  left_join(
    stc %>% 
      filter(!is.na(UsdaClass)) %>% 
      dplyr::select(-NoahClassName)) -> dsmu

# assign to "Organic material" class used by Noah
OrgM_code <- stc$NoahClass[stc$NoahClassName=="Organic material"]
dsmu %>%
  mutate(NoahClass=ifelse(
    ORG_CARBON>=org_thr, 
    OrgM_code, 
    NoahClass)) -> dsmu

# assign to the other Noah special classes
notx_thr <- 50
Watr_code <- stc$NoahClass[stc$NoahClassName=="Water"]
Rock_code <- stc$NoahClass[stc$NoahClassName=="Bedrock"]
Othr_code <- stc$NoahClass[stc$NoahClassName=="Other"]
dsmu %>%
  mutate(
    NoahClass=case_when(
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 1) ~ Watr_code, # Water bodies         
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 2) ~ NoahClass, # Land ice and glaciers
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 3) ~ Rock_code, # Rock outcrops        
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 4) ~ NoahClass, # Dunes/shifting sands 
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 5) ~ NoahClass, # Salt flats           
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 6) ~ Rock_code, # Rocky sublayers      
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 7) ~ Rock_code, # Rocky sublayers      
      (TotShareNoTexture >= notx_thr) & (NoTextureCode_1st == 9) ~ NoahClass, # Other 
      .default = NoahClass
    ),
  NoahClass=ifelse(is.na(NoahClass), Othr_code, NoahClass)) %>%
  left_join(
    stc %>% 
      dplyr::select(-UsdaClass))-> dsmu

# synthesis
dsmu %>%
  group_by(NoahClass, NoahClassName) %>%
  summarize(Perc=signif(n()/nrow(dsmu)*100, 2)) %>%
  arrange(desc(Perc)) %>%
  mutate(Perc=sprintf("%g",Perc)) %>%
  View()

# read data: raster
flog.info("Reading raster...")
dras <- rast("/lustre/arpa/bonafeg/data/geo/SoilType/data/HWSD_2.0/HWSD2_RASTER/HWSD2.bil")

# area of interest, tiling
#bb_aoi <- c(minlon=10.4, maxlon=16.3, minlat=44.6, maxlat=48.1)
bb_aoi <- c(minlon=-180, maxlon=180, minlat=-90, maxlat=90)
tile_step <- c(90,90)
seq_lon <- seq(bb_aoi[1], bb_aoi[2], by=tile_step[1])
seq_lat <- seq(bb_aoi[3], bb_aoi[4], by=tile_step[2])
nx <- length(seq_lon)-1
ny <- length(seq_lat)-1

# join SMU data to raster and write data to TIFF tiles
flog.info("Joining data to raster...")
dsmu %>% distinct(dest_bot, dest_top) -> dest_lay
nlay <- nrow(dest_lay)
filename <- list("out/SoilType_depth{sprintf('%03i',dtop)}to{sprintf('%03i',dbot)}cm",
                 "_lon{sprintf('%+04g',minx)}to{sprintf('%+04g',maxx)}deg",
                 "_lat{sprintf('%+04g',miny)}to{sprintf('%+04g',maxy)}deg",
                 ".tif")
for (ix in 1:nx) {
  for (iy in 1:ny) {
    minx <- seq_lon[ix]
    maxx <- seq_lon[ix+1]
    miny <- seq_lat[iy]
    maxy <- seq_lat[iy+1]
    dras_aoi <- crop(dras, ext(minx,maxx,miny,maxy))
    for (il in 1:nlay){
      dtop <- dest_lay$dest_top[il]
      dbot <- dest_lay$dest_bot[il]
      dd <- dsmu %>% filter(dest_bot==dbot, dest_top==dtop)
      dsmu_aoi <- dras_aoi
      values(dsmu_aoi) <- dd$NoahClass[match(values(dras_aoi), dd$HWSD2_SMU_ID)]
      # save to GeoTIFF
      fout <- glue(glue_collapse(filename))
      flog.info(glue("Writing file {fout}..."))
      writeRaster(dsmu_aoi, glue(fout), overwrite=T)
    }
  }
}

# merge TIFFs in one single TIFF for each layer
for (il in 1:nlay){
  dtop <- dest_lay$dest_top[il]
  dbot <- dest_lay$dest_bot[il]
  rl <- Sys.glob(glue(filename[[1]],"_*_*.tif"))
  flog.info("Creating raster collection...")
  s <- sprc(rl)
  flog.info("Merging rasters...")
  m <- merge(s)
  fout <- glue(glue_collapse(filename[c(1,4)]))
  flog.info(glue("Writing file {fout}..."))
  writeRaster(m, glue(fout), overwrite=T)
}
