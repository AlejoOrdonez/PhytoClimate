rm(list=ls());gc()
library(terra)
require(sp)
## Load the functions
source("./CODE/R-Code/VelocityFnc.R")

#setwd("/Users/alejandroordonez/Dropbox/Other Papers in progress/[Conradi] Phytoclimate  Velocity/PhytoClimate")
## Load the Data
### Load growth form data [22kaBP to present on 500Yrs intervals] last one is 2070 second to last is 1950
pml <- readRDS("./Data/RawData/LGM/gf_suitab_matrices.rds")
### Coordinates of the grid cells
xy <- readRDS("./Data/RawData/LGM/cell_coordinates.rds")
### Raster template
cr.ea <- rast("./Data/RawData/LGM/raster_template.tif")

### Matrix of names
NamesDtFrm <- data.frame(Acro1 = c("TE", "TDdry", "TDcold", "TN", "ShE", "ShDdry", "ShDcold","H","Geo", "Thero", "GC3", "GC4", "Suc", "Clim"),
                         Acro2 = c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C"),
                         Name = c("Evergreen trees", "Drought-deciduous trees", "Cold-deciduous trees", "Needleleaf trees", "Evergreen shrubs", "Drought-deciduous shrubs", "Cold-deciduous shrubs","Herbs", "Geophytes", "Therophytes", "C3 grasses", "C4 grasses", "Succulents", "Climbers"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## Per period
Dates <- list(All = seq(-22.0,1,0.5), # All time periods
              GS2toNow = seq(-15.0,0,0.5), # Crop to 15kaBPO 
              Holocene = seq(-11.5,0,0.5), # Holocene
              GS1 = seq(-13, -11.5,0.5), # The same as YD
              GI1 = seq(-14.5,-13,0.5), # The same as BA
              GS2 = seq(-21.0,-14.5,0.5), # PostLGM
              WarmBA = seq(-15.5,-13.5,0.5), # BA onset warming
              WarmPostYD = seq(-12.5,-9.5,0.5) # BA onset warming
              )

MethodsUse <- c("AbsDif", "Anomaly1", "Anomaly2", "ARM", "lm")

# I will loop over all Methods 
for(MethodUse in MethodsUse){
#(MethodUse <- "AbsDif")
# I will loop over all periods 
  for(TimePer in names(Dates)){#(TimePer<- names(Dates)[2])
    TimePerUse <- Dates[[TimePer]] # Extract the dates
    TimePerUse <- which(seq(-21,0,0.5)%in%TimePerUse) # Turn dates into points in the chronology 
    for(GF.Use in NamesDtFrm$Acro2){#(GF.Use <- NamesDtFrm$Acro2[1])
## **Zero**: Generate a SpatRaster with the suitability rasteres
      a<-Sys.time()
      MapsList <- lapply(TimePerUse,
                         function(i){
                           rasterize(x = as.matrix(xy), # points as a matrix
                                     y = cr.ea, #  template raster
                                     values = pml[[i]][,GF.Use] # Values to aggregate
                                     )
                           })
# make a multi-band SpatRaster
      MapsPer.GF <- do.call("c",MapsList)
## **First**: Generate a SpatRaster that estimates the temporal trend for each cell using the app function from terra
# IMPORTANT: here the approach can be changed - for ease I use a Median of all 500Yrs differences.
      a<-Sys.time()
      TempHetRast <- app(MapsPer.GF,
                         fun = function(i, ff) ff(i,method="Anomaly2"),#
                         ff=TempGradFnc,
                         cores = 10 # the function is run automatically in parallel in 10 cores
                         )
      TempHetRast <- mask(TempHetRast, # Mask based on the input raster
                          MapsPer.GF[[1]])
# Save the Raster file
      writeRaster(TempHetRast,
                  paste0("./Data/LGM/Velocity/TempChng_",MethodUse,"/",
                         TimePer,"_",GF.Use,"_TempChng.tif"),
                  overwrite=TRUE)
      Sys.time()-a
## **Second**: Estimate the spatial gradients magnitude using a using the maximum average technique [Burrough & McDonnell 1998].
      a<-Sys.time()
      if(!paste0(TimePer,"_",GF.Use,"_SpatHet.tif")%in%dir("./Data/LGM/Velocity/SpatHet/")){
        SpaceHetRast <- SpatHetFnc(MapsPer.GF)
# Save the Raster file
        writeRaster(SpaceHetRast,
                    paste0("./Data/LGM/Velocity/SpatHet/",TimePer,"_",GF.Use,"_SpatHet.tif"),
                    overwrite=TRUE)
      } else {
        SpaceHetRast <- rast(paste0("./Data/LGM/Velocity/SpatHet/",TimePer,"_",GF.Use,"_SpatHet.tif"))
      }
      print(Sys.time()-a)
## **Third**: Estimate the spatial gradients direction using a using the maximum average technique [Burrough & McDonnell 1998].
##            Here 90 degrees is poleward direction, so that 0 degrees is East in the north and West in the south
      a<-Sys.time()
      if(!paste0(TimePer,"_",GF.Use,"_Bearing.tif")%in%dir("./Data/LGM/Velocity/Bearing/")){
        BearingRast <- BearingFnc(MapsPer.GF)
# Save the Raster file
        writeRaster(BearingRast,
                    paste0("./Data/LGM/Velocity/Bearing/",TimePer,"_",GF.Use,"_Bearing.tif"),
                    overwrite=TRUE)
      } else {
        BearingRast <- rast(paste0("./Data/LGM/Velocity/Bearing/",TimePer,"_",GF.Use,"_Bearing.tif"))
      }
      print(Sys.time()-a)
## **Forth**: Estimate the velocity magnitude (i.e., speed) as the ratio between the temporal and spatial gradient.
      Velocity <- VelocityFnc (TempHetRast,
                               SpaceHetRast)

      # Save the Raster file
      writeRaster(Velocity,
                  paste0("./Data/LGM/Velocity/Velocity_Anomaly2/",TimePer,"_",GF.Use,"_Velocity.tif"),
                  overwrite=TRUE)
      Sys.time()-a
      rm(list=c("Velocity","BearingRast","SpaceHetRast", "TempHetRast","MapsPer.GF","MapsList"))
    }
  }
}
#plot(log10(Velocity))

## **Fifth**: Estimate the displacement magnitude (i.e., speed).
for(MethodUse in MethodsUse){#(MethodUse <- "AbsDif")
  for(TimePer in names(Dates)){#(TimePer<- names(Dates)[2])
    if(length(dir(paste0("./Data/LGM/Velocity/Velocity_",MethodUse,"/"),
                  pattern=TimePer))>0){
# Load the Velocity vectors
      VelocityByGF <- lapply(NamesDtFrm$Acro2,
                             function(GF.Use){#(GF.Use <- NamesDtFrm$Acro2[1])
                               rast(paste0("./Data/LGM/Velocity/Velocity_",
                                           MethodUse,
                                           "/",TimePer,"_",
                                           GF.Use,
                                           "_Velocity.tif"))
                             })
# Merge all values 
      VelocityByGF <- do.call("c",VelocityByGF)
      names(VelocityByGF) <- NamesDtFrm$Name
# Estimate the Displacement <-  geometric mean of velocities
      Displacement <- 10^mean(log10(VelocityByGF))
# Save the Displacement raster
      writeRaster(Displacement,
                  filename = paste0("./Data/LGM/Displacement/Displacement_",
                                    MethodUse,
                                    "/",
                                    TimePer,
                                    "_Displacement.tif"),
                  overwrite = TRUE)
    }
  }
}
## **Sxith**: Estimate the Divergence of velocity vectors (i.e., sd).
for(TimePer in names(Dates)){#(TimePer<- names(Dates)[2])
# Load the Velocity vectors
    BearingByGF <- lapply(NamesDtFrm$Acro2,
                           function(GF.Use){#(GF.Use <- NamesDtFrm$Acro2[1])
                             rast(paste0("./Data/LGM/Velocity/Bearing/",
                                         TimePer,"_",GF.Use,
                                         "_Bearing.tif"))
                           })
# Merge all values 
    BearingByGF <- do.call("c",BearingByGF)
    names(BearingByGF) <- NamesDtFrm$Name
# Estimate the Displacement <-  geometric mean of velocities
    Divergence <- app(BearingByGF,sd,na.rm=T)
# Save the Displacement raster
    writeRaster(Divergence,
                filename = paste0("./Data/LGM/Divergence/",
                                  TimePer,
                                  "_Divergence.tif"),
                overwrite = TRUE)
}    
