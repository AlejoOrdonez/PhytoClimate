#Create The GEOTIFFS
rm(list=ls());gc()
require(terra)
## LGM
setwd("/Users/au467796/Library/CloudStorage/Dropbox/Other Papers in progress/[Conradi] Phytoclimate  Velocity/PhytoClimate")
## Load the Data
### Load growth form data [21kaBP to present on 500Yrs intervals] last one is 2070 second to last is 1950
pml <- readRDS("./Data/RawData/LGM/gf_suitab_matrices.rds")
### Coordinates of the grid cells
xy <- readRDS("./Data/RawData/LGM/cell_coordinates.rds")
### Raster template
cr.ea <- rast("./Data/RawData/LGM/raster_template.tif")
### Matrix of names
NamesDtFrm <- data.frame(Acro1 = c("TE", "TDdry", "TDcold", "TN", "ShE", "ShDdry", "ShDcold","H","Geo", "Thero", "GC3", "GC4", "Suc", "Clim"),
                         Acro2 = c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C"),
                         Name = c("Evergreen trees", "Drought-deciduous trees", "Cold-deciduous trees", "Needleleaf trees", "Evergreen shrubs", "Drought-deciduous shrubs", "Cold-deciduous shrubs","Herbs", "Geophytes", "Therophytes", "C3 grasses", "C4 grasses", "Succulents", "Climbers"))

for (GF.Use in NamesDtFrm$Acro2[1]){
## **Zero**: Generate a SpatRaster with the suitability rasteres
  MapsList <- lapply(1:44,
                     function(i){
                       rasterize(x = as.matrix(xy), # points as a matrix
                                 y = cr.ea, #  template raster
                                 values = pml[[i]][,GF.Use] # Values to aggregate
                                 )
                       })
# make a multi-band SpatRaster
  Mapsout <- do.call("c",MapsList)
  writeRaster(Mapsout,
              paste0("./Data/LGM/LatePleistocene/suitability_",
                     GF.Use,
                     "LatPleist.tif"),
              overwrite=TRUE)
print(GF.Use)
}
