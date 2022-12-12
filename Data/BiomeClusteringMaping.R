dbiome <- raster("biome_mclust_nodapc_18.tif")


### Define colors
### Green to yellow series
c17 <- "#33562c" # Amazon
c07 <- "#2d872c" # Arc
c18 <- "#33a02c" # Seasonal
c06 <- "#9ce633" # Cerrado
c13 <- "#ddf21e" # N-India/Zambezi
c12 <- "#f2f2a8" # Caatinga


### Orange series
c09 <- "#fec44f" # Sahel
c14 <- "#fe9929" # Karroo
c04 <- "#fc4e2a" # Maghreb desert
c03 <- "#dc0f06" # Sahara


### Purple series
c01 <- "#9e13de" # cool subtropical to temperate
c08 <- "#e077bd" # humid subtropical


### Blue series
c16 <- "#0f0c61" # tundra
c11 <- "#251ef2" # boreal
c15 <- "#258cf2" # hemi-boreal
c02 <- "#bdbdbd" # temperate
c02 <- "#9e9ac8"


### Brown series
c10 <- "#cb7d11" # Great Plains
c05 <- "#cba573" # Kazakhstan


### Put in vector
biomecols <- c(c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16, c17, c18)


### appealing order of zones in legend
nice.numbering <- c(17,7,18,6,13,12, # green series
                    8,1, # purple series
                    9,14,4,3, # orange series
                    2,15,11,16, # blue series inverse
                    10,5) # brown series


### Zone map
par(mar = c(0,0,0,0))
plot(dbiome, legend = FALSE, col = biomecols, axes = FALSE, box = FALSE)
legend(x= -14000000, y=-5, ncol= 3,
       bty="n", legend = c(1:maxValue(dbiome)),
       fill = biomecols[nice.numbering],
       title = expression(bold("Zone")), title.adj=0.45,
       x.intersp = 0.2)