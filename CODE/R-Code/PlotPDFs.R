rm(list=ls());gc()
library(terra)
require(maptools)
Ice <- rast("./Data/IceCover21kaBP.tif") 
## Countries shapefile
data(wrld_simpl)
wrld_simpl <- spTransform(x = wrld_simpl,
                          CRSobj = CRS("+proj=eck4"))

pdf("./PDF Figures/LGM/TempHet.pdf",width = 10, height = 6)
#Plot Temp Heterogenities
for(MethodUse in c("AbsDif", "Anomaly1", "Anomaly2", "ARM", "lm")){#(MethodUse <- "AbsDif")
  # I will loop over all periods 
  for(TimePer in c("All", "Holocene", "GS1", "GI1", "GS2")){#(TimePer<- "All")
    for (GF.Use in c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C")){#(GF.Use<- "TE")
      if(paste0(TimePer, "_", GF.Use,"_TempChng.tif")%in%dir(paste0("./Data/Velocity/TempChng_",MethodUse))){
        RastPlot <- rast(paste0("./Data/Velocity/TempChng_",MethodUse, "/",
                                TimePer, "_", GF.Use,"_TempChng.tif"))
        plot(RastPlot,
             main = paste0("Temp Heterogenity\n",
                           "[",MethodUse,"-",TimePer,"-",GF.Use,"]"),
             col= rev(hcl.colors(100,"RdYlBu")) # set the colors
             )
        plot(wrld_simpl,add=T,lwd=1.2)
        plot(Ice, 
             col = gray(0.3,alpha = 0.6),
             legend = F,
             add = T)
      }
    }
  }
}
dev.off()

pdf("./PDF Figures/LGM/Velocity.pdf",width = 10, height = 6)
#Plot Velocity
for(MethodUse in c("AbsDif", "Anomaly1", "Anomaly2", "ARM", "lm")){#(MethodUse <- "AbsDif")
  # I will loop over all periods 
  for(TimePer in c("All", "Holocene", "GS1", "GI1", "GS2")){#(TimePer<- "All")
    for (GF.Use in c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C")){#(GF.Use<- "TE")
      if(paste0(TimePer, "_", GF.Use,"_Velocity.tif")%in%dir(paste0("./Data/Velocity/Velocity_",MethodUse))){
        RastPlot <- rast(paste0("./Data/Velocity/Velocity_",MethodUse, "/",
                                TimePer, "_", GF.Use,"_Velocity.tif"))
        plot(log10(RastPlot),
             main = paste0("Velocity\n",
                           "[",MethodUse,"-",TimePer,"-",GF.Use,"]"),
             col= rev(hcl.colors(100,"RdYlBu")) # set the colors
        )
        plot(wrld_simpl,add=T,lwd=1.2)
        plot(Ice, 
             col = gray(0.3,alpha = 0.6),
             legend = F,
             add = T)
        }
      }
  }
}
dev.off()

pdf("./PDF Figures/LGM/Displacement.pdf",width = 10, height = 6)
#Plot Displacement
for(MethodUse in c("AbsDif", "Anomaly1", "Anomaly2", "ARM", "lm")){#(MethodUse <- "Anomaly1")
  # I will loop over all periods 
  for(TimePer in c("All", "Holocene", "GS1", "GI1", "GS2")){#(TimePer<- "All")
    if(paste0(TimePer, "_Displacement.tif")%in%dir(paste0("./Data/Displacement/Displacement_",MethodUse))){
      RastPlot <- rast(paste0("./Data/Displacement/Displacement_",MethodUse, "/",
                              TimePer, "_Displacement.tif"))
      plot(log10(RastPlot),
           main = paste0("Displacement\n",
                         "[",MethodUse,"-",TimePer,"]"),
           col= rev(hcl.colors(100,"RdYlBu")) # set the colors
      )
      plot(wrld_simpl,add=T,lwd=1.2)
      plot(Ice, 
           col = gray(0.3,alpha = 0.6),
           legend = F,
           add = T)
    }
  }
}
dev.off()

pdf("./PDF Figures/LGM/Bearing.pdf",width = 10, height = 6)
#Plot Bearings
#Plot Temp Heterogenities
# I will loop over all periods 
for(TimePer in c("All", "Holocene", "GS1", "GI1", "GS2")){#(TimePer<- "All")
  for (GF.Use in c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C")){#(GF.Use<- "TE")
    if(paste0(TimePer, "_", GF.Use,"_Bearing.tif")%in%dir("Data/Velocity/Bearing/")){
      RastPlot <- rast(paste0("./Data/Velocity/Bearing/",
                                TimePer, "_", GF.Use,"_Bearing.tif"))
      plot(RastPlot,
           main = paste0("Temp Bearing\n",
                         "[",TimePer,"-",GF.Use,"]"),
           col= rev(hcl.colors(100,"RdYlBu")) # set the colors
      )
      plot(wrld_simpl,add=T,lwd=1.2)
      plot(Ice, 
           col = gray(0.3,alpha = 0.6),
           legend = F,
           add = T)
    }
  }
}
dev.off()


pdf("./PDF Figures/LGM/Divergence.pdf",width = 10, height = 6)
#Plot Divergence
for(TimePer in c("All", "Holocene", "GS1", "GI1", "GS2")){#(TimePer<- "All")
  if(paste0(TimePer, "_Divergence.tif")%in%dir("./Data/Divergence/")){
    RastPlot <- rast(paste0("./Data/Divergence/",
                            TimePer, "_Divergence.tif"))
    plot(RastPlot,
         main = paste0("Divergence\n",
                       "[",TimePer,"]"),
         col= rev(hcl.colors(100,"RdYlBu")) # set the colors
    )
    plot(wrld_simpl,add=T,lwd=1.2)
    plot(Ice, 
         col = gray(0.3,alpha = 0.6),
         legend = F,
         add = T)
  }
}
dev.off()


##############################################################################
rm(list=ls());gc()
library(terra)
require(maptools)
## Countries shapefile
data(wrld_simpl)
wrld_simpl <- spTransform(x = wrld_simpl,
                          CRSobj = CRS("+proj=eck4"))


#Plot Temp Heterogeneity
for(TimePer in c("RCP26", "RCP85")){#(TimePer<- "RCP26")
  pdf(paste0("./PDF Figures/",TimePer,"/TempHet.pdf"),width = 10, height = 6)
  RastPlot <- rast(paste0("./Data/Future climate/Velocity/TempChng/",
                          TimePer,
                          "_AllGF_TempChng.tif"))
  for (GF.Use in c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C")){#(GF.Use<- "TE")
    plot(RastPlot[[GF.Use]],
         main = paste0("Temp Heterogenity\n",
                       "[",TimePer,"-",GF.Use,"]"),
         col= rev(hcl.colors(100,"RdYlBu")) # set the colors
         )
    plot(wrld_simpl,add=T,lwd=1.2)
    }
    dev.off()
}

#Plot Temp Velocity
for(TimePer in c("RCP26", "RCP85")){#(TimePer<- "RCP26")
  pdf(paste0("./PDF Figures/",TimePer,"/Velocity.pdf"),width = 10, height = 6)
  RastPlot <- rast(paste0("./Data/Future climate/Velocity/Velocity/",
                          TimePer,
                          "_AllGF_Velocity.tif"))
  for (GF.Use in c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C")){#(GF.Use<- "TE")
    plot(log10(RastPlot[[GF.Use]]),
         main = paste0("Velocity\n",
                       "[",TimePer,"-",GF.Use,"]"),
         col= rev(hcl.colors(100,"RdYlBu")) # set the colors
         )
    plot(wrld_simpl,add=T,lwd=1.2)
    }
dev.off()
}

#Plot Temp Displacement
for(TimePer in c("RCP26", "RCP85")){#(TimePer<- "RCP26")
  pdf(paste0("./PDF Figures/",TimePer,"/Displacement.pdf"),width = 10, height = 6)
  RastPlot <- rast(paste0("./Data/Future climate/Displacement/",
                          TimePer,
                          "_AllGF_Displacement.tif"))
  plot(log10(RastPlot),
         main = paste0("Displacement\n",
                       "[",TimePer,"]"),
         col= rev(hcl.colors(100,"RdYlBu")) # set the colors
    )
    plot(wrld_simpl,add=T,lwd=1.2)
  dev.off()
}


pdf("./PDF Figures/LGM/Bearing.pdf",width = 10, height = 6)

#Plot Bearings
for(TimePer in c("RCP26", "RCP85")){#(TimePer<- "RCP26")
  pdf(paste0("./PDF Figures/",TimePer,"/Bearing.pdf"),width = 10, height = 6)
  RastPlot <- rast(paste0("./Data/Future climate/Velocity/Bearing/",
                          TimePer,
                          "_AllGF_Bearing.tif"))
  for (GF.Use in c("TE", "TD_dry", "TD_cold", "TN", "ShrE", "ShrD_dry", "ShrD_cold", "H","HGeo", "HThero", "G_C3", "G_C4", "Suc", "C")){#(GF.Use<- "TE")
    plot(RastPlot[[GF.Use]],
         main = paste0("Temp Bearing\n",
                       "[",TimePer,"-",GF.Use,"]"),
         col= rev(hcl.colors(100,"RdYlBu")) # set the colors
    )
    plot(wrld_simpl,add=T,lwd=1.2)
  }
  dev.off()
}


#Plot Temp Displacement
for(TimePer in c("RCP26", "RCP85")){#(TimePer<- "RCP26")
  pdf(paste0("./PDF Figures/",TimePer,"/Divergence.pdf"),width = 10, height = 6)
  RastPlot <- rast(paste0("./Data/Future climate/Divergence/",
                          TimePer,
                          "_AllGF_Divergence.tif"))
  plot(RastPlot,
       main = paste0("Divergence\n",
                     "[",TimePer,"]"),
       col= rev(hcl.colors(100,"RdYlBu")) # set the colors
  )
  plot(wrld_simpl,add=T,lwd=1.2)
  dev.off()
}



