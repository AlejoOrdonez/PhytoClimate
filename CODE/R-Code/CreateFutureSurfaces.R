library(raster)

setwd("/Users/au467796/Library/CloudStorage/Dropbox/Other Papers in progress/[Conradi] Phytoclimate  Velocity/PhytoClimate")
setwd("/Users/alejandroordonez/Dropbox/Other Papers in progress/[Conradi] Phytoclimate  Velocity/PhytoClimate2")

### Grid with coordinates of grid cells (needed for the rasterization below)
scenarios <- c("rcp26_CCSM4", "rcp26_CNRM-CM5", "rcp26_FGOALS-g2",
               "rcp26_MIROC-ESM", "rcp26_MPI-ESM-LR", "rcp85_CCSM4",
               "rcp85_CNRM-CM5", "rcp85_FGOALS-g2", "rcp85_MIROC-ESM",
               "rcp85_MPI-ESM-LR")
grid025.list <- readRDS("./Data/RawData/Future/grid025list_swapped.rds")
names(grid025.list) <- c("today", scenarios)
unlist(lapply(grid025.list, nrow))


### Raster template "high res" (circa 25 km x 25 km)
r.high <- raster("./Data/RawData/Future/tavg_1.tif")


### Raster template "low res" (circa 50 km x 50 km)
r.low <- raster("./Data/RawData/Future/coarse_raster_template.tif")


### Preference matrix for present-day climatology (mean of 1979-2013)
### (contains suitability values for the growth forms)
pm <- readRDS("./Data/RawData/Future/pref_matrix_today.rds")


### Make suitability raster for each growth form and write it to disk
for(i in 1:ncol(pm)){
  myr <- rasterize(grid025.list[[1]][, 1:2], r.high, field = pm[, i])
  resample(myr, r.low,
           filename = paste("./Data/Future climate/present-day/suitability_", dimnames(pm)[[2]][i], "_ambient.tif", sep = ""),
           overwrite = TRUE)
}


### Do the same for preference matrices for future climatologies
### (mean of 2061-2080; henceforth 2070) (takes 5 min)
f <- dir("./Data/RawData/Future/",pattern = "pref_matrix_2070", full.names = TRUE)
for (i in 1:length(f)) {
  print(i)
  pm <- readRDS(f[i])
  scen <- gsub("\\.rds", "", gsub("pref_matrix_2070_", "", basename(f[i])))
  gcm <- strsplit(scen, "_")[[1]][[1]]
  rcp <- strsplit(scen, "_")[[1]][[2]]
  # add cells to pm if incomplete
  if(nrow(pm) != nrow(grid025.list[[1]])){
    new.matrix <- matrix(NA, nrow = nrow(grid025.list[[1]]), ncol = ncol(pm), dimnames = dimnames(pm))
    grid.id <- intersect(grep(gcm, names(grid025.list)), grep(rcp, names(grid025.list)))
    scores <- which(!is.na(grid025.list[[grid.id]][, "tmax_1"])) 
    for(r in 1:length(scores)){
      new.matrix[scores[r], ] <- pm[r, ]
    }
    pm <- new.matrix
  }
  # rasterize columns of pref matrix
  for(k in 1:ncol(pm)){
    myvals <- pm[, k]
    if(sum(is.na(myvals)) > 0) {
      myvals[ which(is.na(myvals)) ] <- -1
      myr <- rasterize(grid025.list[[1]][, 1:2], r.high, field = myvals)
      myr[Which(myr == -1, cells = T) ] <- NA
      resample(myr, r.low,
               filename = paste("./Data/Future climate/",
                                str_to_upper(strsplit(scen,"_")[[1]][2]),
                                "/All/suitability_", dimnames(pm)[[2]][k], "_", scen, ".tif", sep = ""),
               overwrite = TRUE)
    } else {
      myr <- rasterize(grid025.list[[1]][, 1:2], r.high, field = myvals)
      resample(myr, r.low,
               filename = paste("./Data/Future climate/",
                                str_to_upper(strsplit(scen,"_")[[1]][2]),
                                "/All/suitability_", dimnames(pm)[[2]][k], "_", scen, ".tif", sep = ""), overwrite = TRUE)
    }
  }
}


### Ensemble MEDIANS for RCP 2.6 and 8.5, respectively, do this (takes 2 min)
pm <- readRDS("./Data/RawData/Future/pref_matrix_today.rds")
rcps <- c("rcp26", "rcp85")
for(r in 1:length(rcps)){
  f <- list.files(paste0("./Data/Future climate/",str_to_upper(strsplit(scen,"_")[[1]][2]),"/All"),
                  pattern = rcps[r],
                  full.names = TRUE)
  for(i in 1:ncol(pm)){
    gf <- dimnames(pm)[[2]][i]
    gfs <- stack(f[grep(paste("_", gf, "_", sep = ""), f)])
    gfs.vals <- extract(gfs, grid025.list[[1]][,1:2])
    gfs.median <- apply(gfs.vals, 1, median, na.rm = TRUE)
    if(sum(is.na(gfs.median)) > 0) {
      gfs.median[ which(is.na(gfs.median)) ] <- -1
      myr <- rasterize(grid025.list[[1]][, 1:2], r.high, field = gfs.median)
      myr[Which(myr == -1, cells = T) ] <- NA
      resample(myr, r.low,
               filename = paste("./Data/Future climate/",
                                str_to_upper(strsplit(scen,"_")[[1]][2]),
                                "suitability_", dimnames(pm)[[2]][i], "_MEDIAN_", rcps[r], ".tif", sep = ""), overwrite = TRUE)
    } else {
      myr <- rasterize(grid025.list[[1]][, 1:2], r.high, field = gfs.median)
      resample(myr, r.low, filename = paste("./Data/Future climate/",
                                            str_to_upper(strsplit(scen,"_")[[1]][2]),
                                            "suitability_",
                                            dimnames(pm)[[2]][i], "_MEDIAN_", rcps[r], ".tif", sep = ""), overwrite = TRUE)
    }
  }
}


