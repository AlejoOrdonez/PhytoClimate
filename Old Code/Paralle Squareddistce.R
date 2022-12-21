rm(list=ls());gc()
library(snowfall)
library(terra)
library(analogue)
#setwd("~/Dropbox/Other Papers in progress/[Conradi] Phytoclimate  Velocity")
## Load the Data
### Load growth form data [21kaBP to present on 500Yrs intervals] last one is 2070 second to last is 1950
pml <- readRDS("./Data/gf_suitab_matrices.rds")
### Coordinates of the grid cells
xy <- readRDS("./Data/cell_coordinates.rds")
### Raster template
cr.ea <- rast("./Data/raster_template.tif")
### Matrix of names
NamesDtFrm <- data.frame(Acro1 = c("TE", "TDdry", "TDcold", "TN", "ShE", "ShDdry", "ShDcold","H","Geo", "Thero", "GC3", "GC4", "Suc", "Clim"),
                         Acro2 = c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C"),
                         Name = c("Evergreen trees", "Drought-deciduous trees", "Cold-deciduous trees", "Needleleaf trees", "Evergreen shrubs", "Drought-deciduous shrubs", "Cold-deciduous shrubs","Herbs", "Geophytes", "Therophytes", "C3 grasses", "C4 grasses", "Succulents", "Climbers"))

Maps21kaBPlist <- lapply(dimnames(pml[[1]])[[2]],
                         function(i){
                           rasterize(x = as.matrix(xy), # points as a matrix
                                     y = cr.ea, #  template raster
                                     values = pml[[1]][, i] # Values to aggregate
                           )
                         })
# Turn the list into a multy-layer SpatRaster
Maps21kaBP <- do.call("c",Maps21kaBPlist)
rm(Maps21kaBPlist)
# Add the GF names to each lager
names(Maps21kaBP) <- NamesDtFrm$Name[match(dimnames(pml[[1]])[[2]],NamesDtFrm$Acro2)]

## Make suitability raster for all growth forms for the present time step [1950] - 45th slot in the data list 
Maps0kaBPlist <- lapply(dimnames(pml[[45]])[[2]],
                        function(i){
                          rasterize(x = as.matrix(xy), # points as a matrix
                                    y = cr.ea, #  template raster
                                    values = pml[[45]][, i] # Values to aggregate
                          )
                        })
# Turn the list into a multy-layer SpatRaster
Maps0kaBP <- do.call("c",Maps0kaBPlist)
rm(Maps0kaBPlist);gc()
# Add the GF names to each lager
names(Maps0kaBP) <- NamesDtFrm$Name[match(dimnames(pml[[1]])[[2]],NamesDtFrm$Acro2)]
# Estimate the novelty of each cell in Maps0kaBP [the present conditions measured as 1950 climates]
rm(pml);rm(xy);rm(cr.ea);gc()
## Generate a data.frame of cells with data for present conditions
Maps0kaBPTble <- values(Maps0kaBP,
                        dataframe = TRUE,
                        na.rm = TRUE)
#rm(Maps0kaBP)
## Generate a data.frame of cells with data for past conditions
Maps21kaBPTble<- values(Maps21kaBP,
                        dataframe = TRUE,
                        na.rm = TRUE)
#rm(Maps21kaBP);gc()
CordDistRast <- rast("./Data/CordD_min.tif")
ChiDistRast <- rast("./Data/ChiD_min.tif")
#############
#Smp <- row.names(Maps0kaBPTble)[sample(1:dim(Maps0kaBPTble)[1],1000)]
#Maps0kaBPTble <- Maps0kaBPTble[Smp,]
############
a <- Sys.time()
sfInit(cpus=10,parallel=TRUE)
## Export packages
sfLibrary(analogue)
sfExport("Maps0kaBPTble")
sfExport("Maps21kaBPTble")

# Calculate the Squared Chord distance for each current cell in Maps0kaBP to all past cells in Maps21kaBP
CordDistSumList <- sfLapply(1:dim(Maps0kaBPTble)[1],#(x<-1)
                            function(x, DistMet = "chord"){
                              DistCalc <- analogue::distance(x = Maps0kaBPTble[x,],
                                                             y = Maps21kaBPTble,
                                                             method = DistMet,
                                                             double.zero = TRUE)
                              
                              min(DistCalc,na.rm=T)
                            })
CordDistSum <- data.frame(Dist = do.call("c",CordDistSumList),
                          CellID = as.numeric(rownames(Maps0kaBPTble)))
write.csv(CordDistSum,"./Data/CordDist_min.csv")

# Turn the Min Chord distance per cell into a raster
#CordDistSum <- read.csv("./Data/CordDist_min.csv") # Load data
CordDistRast <- rast(Maps0kaBP[[1]]) # create the empty raster
values(CordDistRast) <- NA # fill with empty data
names(CordDistRast) <- "minCordDist" # change layer name
values(CordDistRast)[CellWithData] <- CordDistSum$Dist # add the ChordD.min data
# Save the 
writeRaster(CordDistRast,
            "./Data/CordD_min.tif",
            overwrite=TRUE)

Sys.time()-a
sfStop()
a <- Sys.time()
sfInit(cpus=10,parallel=TRUE)
## Export packages
sfLibrary(analogue)
sfExport("Maps0kaBPTble")
sfExport("Maps21kaBPTble")

chiSqurDistSumList <- sfLapply(1:dim(Maps0kaBPTble)[1],#(x<-1)
                               function(x, DistMet = "chi.square"){
                                 DistCalc <- analogue::distance(x = Maps0kaBPTble[x,],
                                                                y = Maps21kaBPTble,
                                                                method = DistMet,
                                                                double.zero = TRUE)
                                 
                                 min(DistCalc,na.rm=T)
                               })
chiSqurDistSum <- data.frame(Dist = do.call("c",chiSqurDistSumList),
                             CellID = as.numeric(rownames(Maps0kaBPTble)))
write.csv(chiSqurDistSum,"./Data/ChiSqur_min.csv")


sfStop()