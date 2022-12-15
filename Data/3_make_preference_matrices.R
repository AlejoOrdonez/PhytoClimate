library(bigmemory)
library(doParallel)
library(rasterVis)
library(raster)
library(rgdal)
library(pals)


### Set projection info of equal area grid: World Eckert IV
projectinfo.ea <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"




##############
##############                  Section 1:   GROWTH FORM PREFERENCE MATRIX FOR TODAY
##############
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Load growth form data
traits <- readRDS("/home/shared/PROJECTS/biome_change/working/traits_final.rds")


### Prevent bigmemory package from printing warning message about down-casting from double to integer
options(bigmemory.typecast.warning = FALSE)


### TSS threshold
tss.th <- 0.7


### Exclude species with TSS < tss.th
traits <- traits[which(traits$tss >= tss.th), ]


### Exclude epiphytes (the don't have access to soil water, hence their models are not meaningful)
traits <- droplevels(traits[which(traits$ft != "epiphyte"), ])


### Split herbs in annulas, geophytes and herbs
traits$ft[ which(traits$ft == "H" & traits$lifeform == "geophyte") ] <- "HGeo"
traits$ft[ which(traits$ft == "H" & traits$lifeform == "therophyte")] <- "HThero"


### Growth forms in dataset 
gfs <- sort(unique(traits$ft))


### Time steps
timeid <- seq(-200, 20, 10)


# ### Make growth form preference matrix
# registerDoParallel(5)
# foreach(pro = 1:length(timeid)) %dopar% {
# 
#   projdir <- paste("/home/shared/PROJECTS/paleo_phytoclimates/range_projections/tid_", timeid[pro], sep = "")
#   projdir <- gsub("-", "neg", projdir)
#   f <- list.files(projdir, full.names = TRUE)
#   a <- readRDS(f[1])
#   pref.matrix <- matrix(nrow = length(a), ncol = length(gfs),
#                         dimnames = list(NULL, gfs))
#   for(f in 1:length(gfs)){
#     ### Make presence-absence dataframe
#     # Get model objects of all species of growthform f
#     gf.ss.fits <- paste(projdir, "/", rownames(traits)[which(traits$ft == gfs[f])], ".rds", sep = "")
#     # Set up big.matrix object. This is a presence-absence matrix for the species of growth form f. It does not have colnames. These would be gf.ss
#     pa.df <- big.matrix(nrow = length(a),
#                         ncol = length(gf.ss.fits),
#                         type = "integer", init = 0,
#                         shared = FALSE)
#     for(i in 1:length(gf.ss.fits)) {
#       # print(paste(f, i, sep="-"))
#       b <- readRDS(gf.ss.fits[i])
#       pa.df[, i] <- b
#     }
#   
#     ### Growth form preference for each cell
#     pref <- c()
#     for(i in 1:nrow(pa.df)){
#       pref[i] <- sum(pa.df[i,])/ncol(pa.df)
#     }
#   
#     ### Store in growth form preference matrix
#     pref.matrix[, f] <- pref
#   }
#   
#   ### Save prefmatrix
#   saveRDS(pref.matrix, paste("/home/shared/PROJECTS/paleo_phytoclimates/prefmatrices/tid_", timeid[pro], ".rds", sep = ""))
# }
# stopImplicitCluster()


### Make raster template
ch <- raster(list.files("/home/shared/GIS_DATA/climate_data/chelsa_v12/ambient/temp", full=T)[1])
cr <- raster(crs=crs(ch), ext=extent(ch), resolution=0.5, vals=NULL)
cr.ea <- projectRaster(cr, crs = projectinfo.ea)
load(paste("/home/shared/PROJECTS/paleo_phytoclimates/grid050/tid_", timeid[1], "_grid050.Rdata", sep=""))


### Make array of preference matrices
library(abind)
f <- list.files("/home/shared/PROJECTS/paleo_phytoclimates/prefmatrices")
timeid.new <- gsub("tid_", "", f)
timeid.new <- sort(as.numeric(gsub(".rds", "", timeid.new)))
pml <- list() # preference matrix list
for(pro in 1:length(timeid.new)){
  pml[[pro]] <- readRDS(paste("/home/shared/PROJECTS/paleo_phytoclimates/prefmatrices/tid_", timeid.new[pro], ".rds", sep = ""))
}


### Pm for 2100
pm85highres <- readRDS("/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange/gf_mean_rcp85_biomechange.rds")
grid025 <- readRDS("/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange/grid025_biomechange.rds")
r <- raster("/home/shared/PROJECTS/biome_change/custom_grids/coarse_res/CHELSA/ambient/tavg_1.tif")
for(i in 1:ncol(pm85highres)){
  print(i)
  myr <- rasterize(grid025[, 1:2], r, field = pm85highres[, i])
  resample(myr, cr.ea, filename = paste("/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange/", dimnames(pm85highres)[[2]][i], "_mean_rcp85.tif", sep = ""), overwrite = T)
}
plot(raster(list.files("/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange", pattern = "_mean_rcp85.tif", full.names = T)[14]))


### Make pm from gf stack and add to pml
rcpsuitstack <- stack(list.files("/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange", pattern = "_mean_rcp85.tif", full.names = T))
r <- rasterize(grid050[,1:2], cr.ea, field = pml[[1]][, 1])
rcpsuitstack <- mask(rcpsuitstack, r)
rcpsuitstack.pm <- getValues(rcpsuitstack)
rcpsuitstack.pm <- rcpsuitstack.pm[complete.cases(rcpsuitstack.pm),]
dimnames(rcpsuitstack.pm)[[2]] <- gsub("_mean_rcp85", "", dimnames(rcpsuitstack.pm)[[2]])
pml[[length(pml)+1]] <- rcpsuitstack.pm[, match(dimnames(rcpsuitstack.pm)[[2]], dimnames(pml[[1]])[[2]])]


### pml to array
pma <- abind(pml, along = 3)



### Data for Alejo
saveRDS(pml, "/home/shared/PROJECTS/paleo_phytoclimates/gf_suitab_matrices.rds")




######################### get row id for select locations

### Get coordinates of cities from different biomes
cities <- data.frame(Cityname = c("Atqasuk", "Rovaniemi", "Metz", "Augsburg", "Petropawl", "Foggia", "Hermannsburg", "Shangrao", "Lusaka", "Iquitos"),
                     CountryCode = c("US", "FI", "FR", "DE", "KZ", "IT", "AU", "CN", "ZM", "PE"),
                     lon = NA,
                     lat = NA)
library(RJSONIO)
for(i in 1:nrow(cities)){
  (mycity <- cities$Cityname[i])
  mycity <- gsub(" ", "", mycity)
  (mycntr <- cities$CountryCode[i])
  url <- paste("http://nominatim.openstreetmap.org/search?city=",
               mycity,
               "&countrycodes=",
               mycntr,
               "&limit=9&format=json",
               sep="")
  osm <- fromJSON(url)
  cities$lon[i] <- as.numeric(osm[[1]]$lon)
  cities$lat[i] <- as.numeric(osm[[1]]$lat)
}
cities <- SpatialPointsDataFrame(coords = cities[, c("lon", "lat")], data = cities, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
cities <- spTransform(cities, projectinfo.ea)
cities@data$lon <- coordinates(cities)[,1]
cities@data$lat <- coordinates(cities)[,2]



### Dataframe of mean annual temp in grid cells since 22kb
load(paste("/home/shared/PROJECTS/paleo_phytoclimates/grid050/tid_", timeid.new[1], "_grid050.Rdata", sep=""))
matm <- data.frame(matrix(nrow = nrow(grid050), ncol = length(timeid.new)+1,
                          dimnames = list(1:nrow(grid050), 1:(length(timeid.new)+1))))
matm[, 1] <- as.numeric(apply(grid050[,grep("tmean_",names(grid050))], 1, mean))
for(i in 2:length(timeid.new)){
  load(paste("/home/shared/PROJECTS/paleo_phytoclimates/grid050/tid_", timeid.new[i], "_grid050.Rdata", sep=""))
  matm[, i] <- as.numeric(apply(grid050[,grep("tmean_", names(grid050))], 1, mean))
}


### Add the data for 2070/2100
# grid025.list <- readRDS("/home/shared/PROJECTS/biome_change/working/grid025list_swapped.rds")
# grid025.list <- grid025.list[6:11] # the rcp85 lists
# 
# matm70 <- data.frame(matrix(nrow = nrow(grid025.list[[1]]), ncol = length(grid025.list), dimnames = list(1:nrow(grid025.list[[1]]), 1:length(grid025.list))))
# matm70[, 1] <- as.numeric(apply(grid025[, grep("tmean_", names(grid025))], 1, mean))
# for(i in 2:length(grid025.list)){
#   matm70[, i] <- as.numeric(apply(grid025.list[[i]][, grep("tmean_", names(grid025.list[[i]]))], 1, mean))
# }
# r <- raster("/home/shared/PROJECTS/biome_change/custom_grids/coarse_res/CHELSA/ambient/tavg_1.tif")
# mymat <- as.numeric(apply(matm70/10, 1, mean, na.rm = T))
# myr <- rasterize(grid025.list[[1]][, 1:2], r, field = mymat)
# resample(myr, cr.ea, filename = "/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange/mat_2070_rcp85.tif", overwrite = T)
myr.ea <- raster("/home/shared/PROJECTS/paleo_phytoclimates/datafrombiomechange/mat_2070_rcp85.tif")
matm[, ncol(matm)] <- extract(myr.ea, grid050[, 1:2])


### City data
mycity <- "Augsburg"
mycntr <- as.character(cities$CountryCode[which(cities$Cityname == mycity)])
my_gf <- "TD_cold"
ids_landcells <- cellFromXY(cr.ea, grid050[, 1:2])
id_mycity <- cellFromXY(cr.ea, cities[cities$Cityname == mycity,])
row_of_mycity <- which(ids_landcells == id_mycity)


my_city_through_time <- pma[row_of_mycity,,] # returns matrix with suitability for each GF at each time step at city coordinates



mytemp <- as.numeric(matm[row_of_mycity, ])
mytemp[16] <- mean(c(mytemp[15], mytemp[17]))
plot(kbp, mytemp,
     xlab = "kyr BP", ylab = "Mean annual\ntemperature [°C]", type = "b",
     pch = 21, bg = c(rep("blue", length(kbp)-1), "red"))
grid()




new_res <- 250
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_tdcold_rovaniemi.png",
    units = "in", width = 15, height = 18, pointsize = pntsz, res = new_res)
par(mfrow = c(2, 1), mar = c(4, 5, 0, 2) + 0.1)
kbp <- c(rev(seq(0,22,0.5))*-1, 0.1)
plot(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TD_cold" ), ],
     xlab = "", ylab = "Climatic suitability", type = "b",
     pch = 21, bg = c(rep("white", length(kbp)-1), "red"),
     main = "")
legend(pch = 21, bg = "white", legend = "Cold-decid. trees", bty = "n", "topleft")
plot(kbp, matm[row_of_mycity, ],
     xlab = "kyr BP", ylab = "Mean annual\ntemperature [°C]", type = "b",
     pch = 21, bg = c(rep("blue", length(kbp)-1), "red"))
abline(h=0, lty=2)
dev.off()


new_res <- 250
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_tdcold_c4_rovaniemi.png",
    units = "in", width = 15, height = 18, pointsize = pntsz, res = new_res)
par(mfrow = c(2, 1), mar = c(4, 5, 0, 2) + 0.1)
kbp <- c(rev(seq(0,22,0.5))*-1, 0.1)
plot(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TD_cold" ), ],
     xlab = "", ylab = "Climatic suitability", type = "b",
     pch = 21, bg = c(rep("white", length(kbp)-1), "red"),
     main = "")
points(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "G_C4" ), ],
       pch = 21, bg = c(rep("lightblue", length(kbp)-1), "red"), type = "b")
legend(pch = 21, pt.bg = c("white", "lightblue"), 
       legend = c("Cold-decid. trees", "C4 grasses"),  bty = "n", "topleft")
plot(kbp, matm[row_of_mycity, ],
     xlab = "kyr BP", ylab = "Mean annual\ntemperature [°C]", type = "b",
     pch = 21, bg = c(rep("blue", length(kbp)-1), "red"))
abline(h=0, lty=2)
dev.off()


new_res <- 250
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_tdcold_c4_te_rovaniemi.png",
    units = "in", width = 15, height = 18, pointsize = pntsz, res = new_res)
par(mfrow = c(2, 1), mar = c(4, 5, 0, 2) + 0.1)
kbp <- c(rev(seq(0,22,0.5))*-1, 0.1)
plot(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TD_cold" ), ],
     xlab = "", ylab = "Climatic suitability", type = "b",
     pch = 21, bg = c(rep("white", length(kbp)-1), "red"),
     main = "")
points(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "G_C4" ), ],
       pch = 21, bg = c(rep("lightblue", length(kbp)-1), "red"), type = "b")
points(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TE" ), ],
       pch = 21, bg = c(rep("darkgreen", length(kbp)-1), "red"), type = "b")
legend(pch = 21, pt.bg = c("white", "lightblue", "darkgreen"), 
       legend = c("Cold-decid. trees", "C4 grasses", "Evergreen trees"),  bty = "n", "topleft")
plot(kbp, matm[row_of_mycity, ],
     xlab = "kyr BP", ylab = "Mean annual\ntemperature [°C]", type = "b",
     pch = 21, bg = c(rep("blue", length(kbp)-1), "red"))
abline(h=0, lty=2)
dev.off()



new_res <- 250
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_tdcold_c4_te_tn_rovaniemi.png",
    units = "in", width = 15, height = 18, pointsize = pntsz, res = new_res)
par(mfrow = c(2, 1), mar = c(4, 5, 0, 2) + 0.1)
kbp <- c(rev(seq(0,22,0.5))*-1, 0.1)
plot(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TD_cold" ), ],
     xlab = "", ylab = "Climatic suitability", type = "b",
     pch = 21, bg = c(rep("white", length(kbp)-1), "red"),
     main = "")
points(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "G_C4" ), ],
       pch = 21, bg = c(rep("lightblue", length(kbp)-1), "red"), type = "b")
points(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TE" ), ],
       pch = 21, bg = c(rep("darkgreen", length(kbp)-1), "red"), type = "b")
points(kbp,  my_city_through_time[which( dimnames(my_city_through_time)[[1]] == "TN" ), ],
       pch = 21, bg = c(rep("orange", length(kbp)-1), "red"), type = "b")
legend(pch = 21, pt.bg = c("white", "lightblue", "darkgreen", "orange"), 
       legend = c("Cold-decid. trees", "C4 grasses", "Evergreen trees", "Needleaf trees"),  bty = "n", "topleft")
plot(kbp, matm[row_of_mycity, ],
     xlab = "kyr BP", ylab = "Mean annual\ntemperature [°C]", type = "b",
     pch = 21, bg = c(rep("blue", length(kbp)-1), "red"))
abline(h=0, lty=2)
dev.off()





mycex <- 1.5
suitvals <- my_city_through_time[, 45]

new_res <- 150
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/gf_profile_lusaka.png",
    units = "in", width = 10, height = 12, pointsize = pntsz, res = new_res)
par(mar = c(4,1.5,0,0))
barplot(suitvals[-c(1, 5, 6, 7, 8, 9)],
        main = "", ylab = "",
        xlab = "Climatic suitability",
        names.arg = "", # names(suitvals[-c(1, 5, 6, 7, 8, 9)]),
        col = c('#fdbf6f','#ff7f00', "#fb9a99", '#e31a1c', '#a6cee3','#b2df8a','#33a02c','#1f78b4'),
        horiz = TRUE,
        xlim = c(0, 0.32),
        cex.lab = mycex,
        cex.axis = mycex,
        space = 0)
dev.off()



local_forcing_global <- c()
future_forcing_global <- c()
dryasholocalib <- c()
for(i in 1:nrow(pma[,,1])){
  if(i %in% seq(1, 40000, 1000)) (print(i))
  my_cell_through_time <- pma[i,,] # returns matrix with suitability for each GF at each time step in cell i
  # Compute forcing between each time step
  dm <- as.matrix(dist(t(my_cell_through_time), diag = T, upper = 2))
  dryastoholo <- dm[35,21]
  forcings <- c()
  for(t1 in 2:ncol(dm)){
    t0 = t1-1
    forcings[t1-1] <- dm[t0, t1]
  }
  # Compare future to past forcings
  future_forcing <- forcings[length(forcings)]
  future_forcing_global[i] <- future_forcing
  qf <- quantile(forcings)
  if(future_forcing >= qf[1]) {local_forcing <- 1}
  if(future_forcing > qf[2]) {local_forcing <- 2}
  if(future_forcing > qf[3]) {local_forcing <- 3}
  if(future_forcing > qf[4]) {local_forcing <- 4}
  if(future_forcing == qf[5]) {local_forcing <- 5}
  local_forcing_global[i] <- local_forcing
  dryasholocalib[i] <- future_forcing/dryastoholo
}



library(rasterVis)
png("/home/shared/PROJECTS/paleo_phytoclimates/past_to_future_forcing.png",
    units = "in", width = 24, height = 12, pointsize = pntsz, res = new_res)
forcing.r <- rasterize(grid050[,1:2], cr.ea, field = local_forcing_global)
forcing.r <- as.factor(forcing.r)
rat <- levels(forcing.r)[[1]]
rat[["forcing"]] <- c("weak", "below average", "above average", "strong", "unprecedented")
levels(forcing.r) <- rat
cols <- c('#2c7bb6','#abd9e9', '#ffffbf', '#fdae61', '#d7191c')
levelplot(forcing.r, col.regions = cols,
          labels = rat$forcing,
          main = "",
          colorkey=list(labels=list(at=1:5, labels=rat[["forcing"]])))
dev.off()


dryasholocalib[which(!is.finite(dryasholocalib))] <- -999
dryasholocalib[which(is.na(dryasholocalib))] <- -999
dryasholocalib.r <- rasterize(grid050[,1:2], cr.ea, field = dryasholocalib)
dryasholocalib.r[dryasholocalib.r == -999] <- NA
dryasholocalib.r[dryasholocalib.r>20] <- NA

dryasholocalib.r[dryasholocalib.r>10] <- 100
dryasholocalib.r[dryasholocalib.r>5 & dryasholocalib.r <= 10] <- 90
dryasholocalib.r[dryasholocalib.r>2 & dryasholocalib.r <= 5] <- 80
dryasholocalib.r[dryasholocalib.r>1 & dryasholocalib.r <= 2] <- 70
dryasholocalib.r[dryasholocalib.r>0.5 & dryasholocalib.r <= 1] <- 60
dryasholocalib.r[dryasholocalib.r>0.25 & dryasholocalib.r <= 0.5] <- 50
dryasholocalib.r[dryasholocalib.r>=0 & dryasholocalib.r <= 0.25] <- 40
plot(dryasholocalib.r)

dryasholocalib.r <- as.factor(dryasholocalib.r)
rat <- levels(dryasholocalib.r)[[1]]
rat[["level"]] <- c("0-0.25", ">0.25-0.5", ">0.5-1", ">1-2", ">2-5", ">5-10", ">10-20")
levels(dryasholocalib.r) <- rat

myTheme <- rasterTheme(rev(pals::parula(n=255)))
maxofcolorkey <- ceiling(max(maxValue(dryasholocalib.r)))
cols <- c('#2b8cbe', '#f6eff7','#fdd49e','#fdbb84','#fc8d59','#e34a33','#b30000')
pntsz <- 250
png("/home/shared/PROJECTS/paleo_phytoclimates/dryas_to_future_forcing.png",
    units = "in", width = 24, height = 12, pointsize = pntsz, res = new_res)
levelplot(dryasholocalib.r, col.regions = cols,
          labels = rat$forcing,
          main = "",
          colorkey=list(labels=list(at=1:7, labels=rat[["level"]], cex = 3), width = 4, height = 0.5),
          xlab = NULL, ylab = NULL,
          scales = list(draw=FALSE)) +
  latticeExtra::layer(sp.polygons(cntr, lwd = 1) )
dev.off()



### City data
mycity <- "Rovaniemi"
mycntr <- as.character(cities$CountryCode[which(cities$Cityname == mycity)])
ids_landcells <- cellFromXY(cr.ea, grid050[, 1:2])
id_mycity <- cellFromXY(cr.ea, cities[cities$Cityname == mycity,])
row_of_mycity <- which(ids_landcells == id_mycity)




my_cell_through_time <- pma[row_of_mycity,,] # returns matrix with suitability for each GF at each time step in cell i
# Compute forcing between each time step
dm <- as.matrix(dist(t(my_cell_through_time), diag = T, upper = 2))
forcings <- c()
for(t1 in 2:ncol(dm)){
  t0 = t1-1
  forcings[t1-1] <- dm[t0, t1]
}
plot(forcings)



###










### GF stack for Cold-deciduous trees
gf = which(dimnames(pma[,,1])[[2]]=="TD_cold")
gf.stack <- rasterize(grid050[,1:2], r, field = pma[,gf,][,1])
for(i in 2:dim(pma)[3]){
  print(i)
  newr <- rasterize(grid050[,1:2], r, field = pma[,gf,][,i])
  gf.stack <- stack(gf.stack, newr)
}
plot(gf.stack[[c(1, 23, 35, 45)]], main = kbp[c(1, 23, 35, 45)]) 





### Figure preferences
land.lwd <- 0.1 # line width of the land polygon
ef <- 0.05 # Expansion factor for raster plotting extent
aoi <- c(-14571548 -14571548*ef, 16655656 +16655656*ef, -6622462 -6622462*ef, 8380642 +8380642*ef) # Are of interest to be shown on the maps
myTheme <- rasterTheme(rev(pals::parula(n=255)))
maxofcolorkey <- ceiling(max(maxValue(gf.stack[[c(1, 23, 35, 45)]]))*100)/100
lab.cex <- 3
axis.cex <- 3


### Country shape
cntr <- readOGR("/home/shared/GIS_DATA/maps_third_party/NaturalEarth/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp")
cntr <- spTransform(cntr, CRSobj = crs(r))
cntr <- disaggregate(cntr)
# # The next two lines are just for checking in QGIS which polygons need to be removed from cntr
# writeOGR(land, "/home/shared/PROJECTS/biome_change", "ne_110m_land_clipped", driver = "ESRI Shapefile")
# writeOGR(cntr, "/home/shared/PROJECTS/biome_change", "ne_110m_admin_0_countries_singlepart", driver = "ESRI Shapefile")
### Remove polygons for which me make no predictions (same as above but have different row IDs)
cntr <- cntr[-c(3,264:271,99), ]


### Shift cntr shapefiles to match raster
cntr <- shift(cntr, dx = res(r)[1]/2, dy = -res(r)[2]/2)

timesteps <- c(21, 26, 33, 38, 45, 46)
gf.stack.sub <- gf.stack[[timesteps]]
panelnames <- c(paste(gsub("-", "", kbp[timesteps]), "k-BP", sep = "")[1:4], c("1950", "2100"))


new_res <- 250
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_stack_TD.png",
    units = "in", width = 36, height = 16, pointsize = pntsz, res = new_res)
par(mar = c(0,0,0,0))
levelplot(crop(gf.stack.sub, aoi), mar = F,
          par.settings = myTheme,
          names.attr = panelnames,
          par.strip.text = list(cex = lab.cex, lines = 1.5), # edits the panel titles
          between = list(x = 0, y = 1.5), # space between panels
          at = seq(0, maxofcolorkey, 0.01),
          colorkey = list(#at = seq(0, 1, 0.01), # where colors change
            labels = list(at = seq(0, maxofcolorkey, 0.1),
                          cex = axis.cex), # where labels are located,
            cex = lab.cex,
            width = 2, height = 0.5,
            space = "bottom"),
          xlab = NULL, ylab = NULL,
          scales = list(draw=FALSE)) +
  latticeExtra::layer(sp.polygons(cntr, lwd = land.lwd) )
dev.off()


png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_stack_TD_rovaniemi.png",
    units = "in", width = 36, height = 16, pointsize = pntsz, res = new_res)
par(mar = c(0,0,0,0))
levelplot(crop(gf.stack.sub, aoi), mar = F,
          par.settings = myTheme,
          names.attr = panelnames,
          par.strip.text = list(cex = lab.cex, lines = 1.5), # edits the panel titles
          between = list(x = 0, y = 1.5), # space between panels
          at = seq(0, maxofcolorkey, 0.01),
          colorkey = list(#at = seq(0, 1, 0.01), # where colors change
            labels = list(at = seq(0, maxofcolorkey, 0.1),
                          cex = axis.cex), # where labels are located,
            cex = lab.cex,
            width = 2, height = 0.5,
            space = "bottom"),
          xlab = NULL, ylab = NULL,
          scales = list(draw=FALSE)) +
  latticeExtra::layer(sp.polygons(cntr, lwd = land.lwd) ) +
  latticeExtra::layer(sp.points(cities, lwd = land.lwd) ) 
dev.off()





add_shp=function(){plot(cntr, bg="transparent", add=TRUE, lwd = land.lwd)}
plot(crop(gf.stack.sub, aoi), addfun = add_shp)

animate(gf.stack[[-46]], pause=0.5, main = paste(gsub("-", "", kbp), "k-BP", sep = "")[-46], z = c(0,  max(maxValue(gf.stack))), maxpixels=50000, n=1)

spe <- readRDS("/home/shared/PROJECTS/paleo_phytoclimates/species.list.rds")
spe$ent1




### Atmospheric CO2 in last 800k years from Lüthi, D., Le Floch, M., Bereiter, B. et al. High-resolution carbon dioxide concentration record 650,000–800,000 years before present. Nature 453, 379–382 (2008). https://doi.org/10.1038/nature06949
library(openxlsx)
co2_luethi <- read.xlsx("/home/shared/PROJECTS/paleo_phytoclimates/co2_quarternary.xlsx", sheet = 3)
co2_luethi <- co2_luethi[co2_luethi$year_bp < 22000,]
co2_luethi$year_bp <- co2_luethi$year_bp/1000
newdf <- data.frame(year_bp = c(-150/1000, 0), co2_ppm = c(677, 310) )
co2_luethi <- rbind(newdf, co2_luethi)


globmean <- apply(matm, 2, mean)
deltaGSMT <- globmean - mean(matm[,45])
new_res <- 400
pntsz <- 12*new_res/72
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_mat.png",
    units = "in", width = 36, height = 16, pointsize = pntsz, res = new_res)
par(mar = c(4, 5, 0, 2) + 0.1)
plot(deltaGSMT ~ kbp, type = "n",
     xlab = "Age (kyr BP)", ylab = expression(paste(Delta, "Global mean\ntemperature (°C)")))
grid()
abline(h = 0, lty = 2, lwd = 2)
points(deltaGSMT[1:45] ~ kbp[1:45], type = "l", lwd = 2, col = "blue", pch = 21, bg = "blue")
points(deltaGSMT[1:45] ~ kbp[1:45], type = "b", lwd = 2, col = "blue", pch = 21, bg = "blue")
points(deltaGSMT[45:46] ~ kbp[45:46], type = "b", lty = 2, lwd = 2, col = "blue", pch = 21, bg = "blue")
points(deltaGSMT[46] ~ kbp[46], pch = 21, bg = "red")
dev.off()


co2 <- co2_luethi$co2_ppm
year_bp <- co2_luethi$year_bp*-1
png("/home/shared/PROJECTS/paleo_phytoclimates/timeseries_co2.png",
    units = "in", width = 36, height = 16, pointsize = pntsz, res = new_res)
par(mar = c(4, 5, 0, 2) + 0.1)
plot(co2 ~ year_bp, type = "n",
     xlab = "Age (kyr BP)", ylab = "CO2 (ppm)")
grid()
points(co2[-1] ~ year_bp[-1], type = "l", lwd = 6, col = "black", pch = 21, bg = "black")
points(co2[1:2] ~ year_bp[1:2], type = "l", lty = 2, lwd = 6, col = "black", pch = 21, bg = "black")
points(co2[1] ~ year_bp[1], pch = 21, bg = "red")
dev.off()


