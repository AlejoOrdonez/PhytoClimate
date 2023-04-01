rm(list=ls());gc()
require(tidyverse)
require(analogue)
require(mgcv)
require(maptools)
require(terra)
data(wrld_simpl)
NorthAmeric <- vect(wrld_simpl)
NorthAmeric2 <- terra::crop(NorthAmeric,
                     ext(c(range(PollenTRees$long),range(PollenTRees$lat))))
# Load the data
PollenTRees <- read_csv("./Data/Pollen data/all_woody_counts_final_fullTS.csv")

#Names of the variables
names(PollenTRees)

# used functions
## Turn counts/proportions into presence absence matrix
PresAbsFnc <- function(x){ifelse(x==0,0,1)}
# Remove a given proportion
Lesthan5Per <- function(x){ifelse(x<0.005,0,x)}

#=====================================================================================
#### Time coverage of the cores

# Count the number of samples per core 
PollenTRees %>% 
  count(dataset.id) %>%  # Count the number of times a core is been used
  summarise(n.Min=min(n), #Min number of times
            nMax=max(n) #Max number of times
            )
### NOTE: All cores have 21 time steps - Any of these is empty?


# Count the number of samples per core removing those that do not have pollen counts
PollenTRees %>% 
  rowwise() %>% # Make summaries rowwise
  mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains
  filter(N.pollCnt>0) %>% 
  count(dataset.id) %>%  # Count the number of times a core is been used
  ungroup() %>%
  summarise(n.Min=min(n), #Min number of times
            nMax=max(n) #Max number of times
  )
## NOTE: The range of time bins changes between cores with some having as low as only 4 points in time


# Check the time span (number of years) for each core.
corePerTim <- PollenTRees %>%
  rowwise() %>% # Make summaries row-wise
  mutate( N.pollCnt = sum(Alnus:Acacia)) %>% #Summ the pollen grains
  select(dataset.id,N.pollCnt) %>% #Only use The Id and counts
  filter(N.pollCnt!=0) %>% #Remove time/core combinations with no counts
  count(dataset.id) %>% #Number of times per core
  ungroup() %>% #Remove the per row summaries
  count(n,
        name = "StPrTim")  # Count the number of cores 
## Plot the distribution of time bins per core
ggplot(corePerTim) +
  aes(y=StPrTim,
      x=n) +
  geom_col()
### NOTE: Most cores have 13 t0 14 time steps


# Define the start end time for each core
StrEndCore <- PollenTRees %>%
  rowwise() %>% # Make summaries row-wise
  mutate( N.pollCnt = sum(Alnus:Acacia)) %>% #Summ the Pollen grains
  select(dataset.id,N.pollCnt,age_group) %>% #Only use The Id and counts
  filter(N.pollCnt!=0) %>% #Remove time/core combinations with no counts
  ungroup() %>% #Remove the per row summaries
  group_by(dataset.id) %>% # Grup per IDE
  summarise(MinT = min(age_group), # get the Summary per Core
            MaxT = max(age_group))

table(StrEndCore$MinT) # most cores end at 1ka BP
table(StrEndCore$MaxT) # most cores start between at 13 and 14ka BP

# A simple visualization of the time extend for each core
df <- tibble(StrEnd = c(StrEndCore$MaxT,StrEndCore$MinT),
             Core = as.factor(c(1:dim(StrEndCore)[1],1:dim(StrEndCore)[1])))
ggplot(df,aes(x=Core,y=StrEnd)) + 
  geom_boxplot(coef=0,
               lower = StrEndCore$MinT, # set the low value of the bar to the min.
               upper = StrEndCore$MaxT # set the upper value of the bar to the max.
  )
### NOTE: Based on this I will focus on changes between 16kaBP to 1kaBP.

## Conclusion: I should only focus on 16kaBP to now
#=====================================================================================

#=====================================================================================
## Evaluating the taxonomic coverage of the Pollen data

# Summary of Pollen counts core/time? 
PollenTRees %>% 
  filter(age_group < 16) %>% # Remove samples older tha 16kaBP
  rowwise() %>% # Make summaries rowwise
  mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains
  filter(N.pollCnt>0) %>% # Filter cores with no spp data
  select(N.pollCnt) %>% #only use the N.pollCnt sumamry
  ungroup() %>%  # Remote the riw wise
  summarise(Min = min(N.pollCnt),
            lowQ = quantile(N.pollCnt,0.25),
            Median = quantile(N.pollCnt,0.5),
            highQ = quantile(N.pollCnt,0.75),
            Lessthan10 = sum(N.pollCnt<10),
            Lessthan50 = sum(N.pollCnt<50),
            Lessthan100 = sum(N.pollCnt<100))
## NOTE: There are a cores with only one pollen grains. With half of the cores having a median of 691 pollen grain


# Assessing number of time bins with low pollen counts
SumPercore <- PollenTRees %>% 
              filter(age_group < 16) %>% # Remove samples older tha 16kaBP
              rowwise() %>% # Make summaries rowwise
              mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains
              filter(N.pollCnt>0) %>% 
              ungroup() %>%
              group_by(dataset.id) %>%  
              summarise(Lessthan10 = sum(N.pollCnt<10),
                        Lessthan50 = sum(N.pollCnt<50),
                        Lessthan100 = sum(N.pollCnt<100))  

#Number of cores/times with less than 10 pollen grains
sum(SumPercore$Lessthan10>0)
#Number of cores/times with less than 50 pollen grains
sum(SumPercore$Lessthan50>0)



# Total number of taxa per core/time
PresAbsFnc <- function(x){ifelse(x==0,0,1)}
SppSPerCore <- PollenTRees %>% 
                filter(age_group < 16) %>% # Remove samples older tha 16kaBP
                mutate(across(Alnus:Acacia, # turn the obs into presence absence
                              PresAbsFnc)) %>% #count of non zero values 
                rowwise() %>% 
                mutate( N.Spp = rowSums(across(Alnus:Acacia))) %>%
                filter(N.Spp>0)

# plot the number of spp per core/time
ggplot(SppSPerCore, aes(x=N.Spp))+
  geom_density()
## NOTE: For the median number of sp is 12

# Summary of richness per core/time
PollS.Summ <- SppSPerCore %>%
              filter(age_group < 16) %>% # Remove samples older than 16kaBP
              group_by(dataset.id) %>% # Make summaries by Core IDE
              summarise(AvSppRich = mean(N.Spp),
                        MinSppRich = min(N.Spp),
                        MaxSppRich = max(N.Spp),
                        lowQ.S = quantile(N.Spp,0.25),
                        Median.S = quantile(N.Spp,0.5),
                        highQ.S = quantile(N.Spp,0.75))

#Plot the values
ggplot(PollS.Summ, aes(x=MinSppRich))+
  geom_density()

ggplot(PollS.Summ, aes(x=MaxSppRich))+
  geom_density()

## NOTE: For most cores, the minimum number of spp at any given time bin is ~5, but could be as low as 1. On the top end, cores could have as much as ~10 or ~20 spp.
##       However, all these are estimates based on the idea that one pollen grain means a presence.


## Richness based on proportions
Lesthan5Per <- function(x){ifelse(x<0.005,0,x)}
PercentSppDbs <- PollenTRees %>% 
                  filter(age_group < 16) %>% # Remove samples older tha 16kaBP
                  rowwise() %>% # Make summaries rowwise
                  mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains
                  filter(N.pollCnt>0) %>% # remove sites with no spp data
                  mutate(across(Alnus:Acacia, # Estimate relative abundances
                                ~./N.pollCnt)) %>%
                  mutate(across(Alnus:Acacia, #make low proportion Zero
                                Lesthan5Per)) %>%
                  mutate(across(Alnus:Acacia, # turn the obs into presence absence
                                PresAbsFnc)) %>% #count of non zero values 
                  rowwise() %>% 
                  mutate( N.Spp = rowSums(across(Alnus:Acacia))) %>% # estimate the spp richness
                  filter(N.Spp>0) %>% # remove cells with no spp
                  select(N.Spp)
                  
ggplot(SppSPerCore, aes(x=N.Spp))+
  geom_density()

## Conclusion: Remove spp with low pollen counts (less that 5% of the original count), and use these filtered dataset for estimates of similarity/beta diversity
#=====================================================================================

# Build a Spp proportion data-set
PercentSppDbs <- PollenTRees %>% 
  filter(age_group < 16) %>% # Remove samples older than 16kaBP
  rowwise() %>% # Make summaries rowwise
  mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains
  filter(N.pollCnt>0) %>% # remove sites with no spp data
  mutate(across(Alnus:Acacia, # Estimate relative abundances
                ~./N.pollCnt)) %>%
  mutate(across(Alnus:Acacia, #make low proportion Zero
                Lesthan5Per)) %>%
  mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains to make sure all all sacle to a 100%
  mutate(across(Alnus:Acacia, # Estimate relative abundances after removing the low values
                ~./N.pollCnt)) %>%
  mutate( N.pollCnt = rowSums(across(Alnus:Acacia))) %>% #Summ the pollen grains
  filter(N.pollCnt>0) # remove sites with no spp data
## NOTE: PercentSppDbs object has re-scales all the pollen count to relative abundances for contrast  

# Filter all sites with less than a 14kaYrs coverage (i.e., less tna 14 time points)
PercentSppDbs <- PercentSppDbs[PercentSppDbs$dataset.id%in%names(table(PercentSppDbs$dataset.id))[which(table(PercentSppDbs$dataset.id)>14)],]

 
# Visualize the cumulative dissimilarity 

# Estimate the off diagonal of the distance matrix 
PollenTReesList <- lapply(sort(unique(PercentSppDbs$dataset.id)),
                          function(x){temp <- PollenTRees[PollenTRees$dataset.id==x,] # Get the full time series of a given core
                                      temp <- temp[rowSums(temp[,-c(1:6)])>0,] # Remove sites with no spp data
                                      temp <- temp[temp$age_group<=16,] # remove sites older than 16kaBO
                                      return(temp)
                            })

# estimate the commutative dissimilarity 
DistUse <- "SQeuclidean"#"SQchord" "SQchi.square", "SQeuclidean","bray"

DistList <- lapply(PollenTReesList,
                   function(x){
                     xSpp <- x[,-c(1:6)]
                     xSpp <- (xSpp[ , colSums(xSpp)>0]/rowSums(xSpp)) # filter spp with No data and scale the metrics as %of the total pollen count
                     Distcores <- analogue::distance(xSpp,
                                                     method = DistUse)
                     dimnames(Distcores) <- list(x$age_group,x$age_group)
                     tibble(Core = unique(x$dataset.id),
                            Age = -(rev(x$age_group)),
                            CummDis = c(0,cumsum(sapply(dim(Distcores)[1]:2,
                                                  function(x){Distcores[x-1,x]}))))
                                       })

plot.new()
plot.window(xlim = c(-22,0),
            ylim = c(0,1.5))#round(max(sapply(DistList[-c(22,56)],
                          #              function(x){max(x$CummDis)})))))
axis(1);axis(2)
mtext(DistUse,side=3,outer=T,
      xpd=NA,
      cex = 2, font = 2,
      line=-5)
for(i in c(1:length(DistList))[-c(22,56)]){
x <- DistList[[i]]  
points(x = x$Age,
     y = x$CummDis,
     type = "b")
}



# Pot one core at a time
for(i in c(1:length(DistList))[-c(22,56)]){#i<-1
  x <- DistList[[i]]  # load the core info
  
  
  
  #Plot the Core Accumulative dissimilarity
  par(fig=c(0,1,0,1),new=F,
      mar=c(3.1, 3.1, 2.1, 2.1))
  plot(x = x$Age,
       y = x$CummDis,
       type = "b",
       pch=19,
       xlim = c(-22,0),
       ylim = c(0,1.5),
       main = paste0("Core ID - ",x$Core[1]))
  GAM.Mod <- gam(CummDis~s(Age),data=x)
  lines(predict(GAM.Mod)~x$Age,
        lwd=2,
        col="red")
  legend("topright",
         paste0("EDF=", round(summary(GAM.Mod)$edf),1))
  plot.window(xlim = c(0,1),
              ylim = c(0,1))
  par(fig=c(0,0.5,0.5,1),
      new=T)
  plot(NorthAmeric2,
       axes=F)
  box()
  points(PollenTRees[which(PollenTRees$dataset.id==x$Core[1])[1],c("long","lat")],
         pch=19)
  # Load the displacement Info
  #abs(min(x$Age)) # start age
  #Disp <- 
}

# sort out Sites with the max cumulative disimialrity
order(sapply(DistList,
             function(x){max(x[,'CummDis'])}))

# Max dissimilarity
PollenTReesList[[56]] # Long -130.096; Lat = 57.59498 --> Kinaskan lake (Yukon)

PollenTReesList[[22]] # Long = -144.6583 Lat = 63.94 --> Healy Lake (Alaska)



All_Displacement.tif

# Suitability extract
a <- rast("./Data/LGM/LatePleistocene/suitability_CLatPleist.tif")


b <-project(a,"+proj=longlat +datum=WGS84")

d <- extract(b,PollenTRees[which(PollenTRees$dataset.id==x$Core[1])[1],c("long","lat")])





d[seq(-21,0.5,by=0.5)%in%x$Age]
