rm(list=ls());gc()
require(raster)
set.seed(180)
####---------------------------------------------------------------------------------

#TEMPERATURE
#### Define the Rate of change (change in units per time interval) for each cell
RtChngTemp<-sort(runif(100,0,0.3))
#### Generate an example temperate data series for 10-interval periods
TempGrad.SN<-raster(matrix(rep(1:10,each=10),nrow=10,byrow=T)+matrix(runif(100,-0.5,0.5),nrow=10))
TempGrad.SN<-stack(TempGrad.SN,raster(matrix(TempGrad.SN[]+RtChngTemp+ runif(100,-0.025,0.025),nrow=10,byrow=T)))
for (i in 3:10){TempGrad.SN<-stack(TempGrad.SN,raster(matrix(TempGrad.SN[[c(i-1)]][]+RtChngTemp+ runif(100,-0.025,0.025),nrow=10,byrow=T)))}
# dev.new(width=8,height=6)
# plot(TempGrad.SN,main=paste("Temperature -",1:10),col=colfunc(100),cex.main=2, axes=F)
####---------------------------------------------------------------------------------
####---------------------------------------------------------------------------------
#PRECIPITATION
#### Define the Rate of change (change in units per time interval) for each cell
RtChngPrec<-c(sort(runif(50,-0.5,0.5),T),sort(runif(50,-0.5,0.5),T))
RtChngPrec<-c(sort(c(runif(25,-0.3,0),runif(25,0,0.3)),T),sort(c(runif(25,-0.3,0),runif(25,0,0.3)),T))
#### Generate an example temperate data series for 10-interval periods
#PresGrad.SN <-raster(matrix(rep(10:1,each=10),nrow=10,byrow=T)+matrix(runif(100,-0.5,0.5),nrow=10))
PresGrad.SN <-raster(matrix(rep(c(10:1,1:10),each=5),nrow=10,byrow=T)+matrix(runif(100,-0.5,0.5),nrow=10))
PresGrad.SN <-stack(PresGrad.SN,raster(matrix(PresGrad.SN[]+ RtChngPrec + runif(100,-0.025,0.025),nrow=10,byrow=T)))
for (i in 3:10){PresGrad.SN <-stack(PresGrad.SN,raster(matrix(PresGrad.SN[[c(i-1)]][]+ RtChngPrec + runif(100,-0.025,0.025),nrow=10,byrow=T)))}
# dev.new(width=8,height=6)
# plot(PresGrad.SN,main=paste("Precipitation -",1:10),col= rev(colfunc(100)),cex.main=2, axes=F)
####---------------------------------------------------------------------------------
####---------------------------------------------------------------------------------
#### Define the Spatial Heterogeneity and directionality of the historical gradient 

#### TEMPERATURE: Estimating Spatial Heterogeneity-based on initial year [i.e., 1901] conditions 
### The direction of the gradients goes as follows: positive values indicate a northward vector (north<south) or westwards vectors (west<east) while negative indicates the opposite. The idea is a high to low directionality  
    SpcHet.EastWest <- focal(TempGrad.SN[[1]], w=matrix(1, nrow=3, ncol=3), fun=function(x){mean(c((x[2]-x[1]),(x[3]-x[2]),(x[5]-x[4]),(x[6]-x[5]),(x[8]-x[7]),(x[9]-x[8])),na.rm=T)}, pad=TRUE, padValue=NA)
    SpcHet.NorthSouth <- focal(TempGrad.SN[[1]], w=matrix(1, nrow=3, ncol=3), fun=function(x){mean(c((x[4]-x[1]),(x[7]-x[4]),(x[5]-x[2]),(x[8]-x[5]),(x[6]-x[3]),(x[9]-x[6])),na.rm=T)}, pad=TRUE, padValue=NA)
	SpcHetTemp <- sqrt((SpcHet.NorthSouth^2)+(SpcHet.EastWest^2))/mean(res(TempGrad.SN))

#### TEMPERATURE: Estimating the bearing of the velocity vector based on initial conditions
    Bearing.Temp <- atan2(x=SpcHet.EastWest, y= SpcHet.NorthSouth)*(180/pi)
    Bearing.Temp[Bearing.Temp[]<0] <- (360+ Bearing.Temp[Bearing.Temp[]<0])

#### PRECIPITATION: Estimating Spatial Heterogeneity-based on initial year [i.e., 1901] conditions 
### The direction of the gradients goes as follows: positive values indicate a northward vector (north<south) or westwards vectors (west<east) while negative indicates the opposite. The idea is a high to low directionality  
    SpcHet.EastWest <- focal(PresGrad.SN[[1]], w=matrix(1, nrow=3, ncol=3), fun=function(x){mean(c((x[2]-x[1]),(x[3]-x[2]),(x[5]-x[4]),(x[6]-x[5]),(x[8]-x[7]),(x[8]-x[9])),na.rm=T)}, pad=TRUE, padValue=NA)
    SpcHet.NorthSouth <- focal(PresGrad.SN[[1]], w=matrix(1, nrow=3, ncol=3), fun=function(x){mean(c((x[4]-x[1]),(x[7]-x[4]),(x[5]-x[2]),(x[8]-x[5]),(x[6]-x[3]),(x[9]-x[6])),na.rm=T)}, pad=TRUE, padValue=NA)
	SpcHetPrec <- sqrt((SpcHet.NorthSouth^2)+(SpcHet.EastWest^2))/mean(res(PresGrad.SN))

#### PRECIPITATION: Estimating the bearing of the velocity vector based on initial conditions
    Bearing.Prec <- atan2(x=SpcHet.EastWest, y= SpcHet.NorthSouth)*(180/pi)
    Bearing.Prec[Bearing.Prec[]<0] <- (360+ Bearing.Prec[Bearing.Prec[]<0])


#### Local Spatial Gradient
	Diverg<-abs(Bearing.Temp-Bearing.Prec)
	Diverg<-calc(Diverg,fun=function(x){ifelse(x<=180,x,360-x)})


#### TEMPERATURE: Determine the end point of a vector representing the Bearing 
VecPosTemp<-t(sapply(Bearing.Temp[],function(x){
			x1<-ifelse(x<90,x,ifelse(x<180,180-x,ifelse(x<270,270-x,ifelse(x<360,360-x))))
			out<-c(cos(x1*(pi/180)), sin(x1*(pi/180)))							
			out<-out*			c(ifelse(x<90,1,ifelse(x<180,-1,ifelse(x<270,-1,ifelse(x<360,1)))),ifelse(x<90,1,ifelse(x<180,1,ifelse(x<270,-1,ifelse(x<360,-1)))))
#This quadrant requires the inversion of the X and Y sides
			if(x>180 & x<270){out<-out[c(2,1)]}
			return(t(out))}))

####PRECIPITATION: Determine the end point of a vector representing the Bearing 
	VecPosPrec<-t(sapply(Bearing.Prec[],function(x){
				x1<-ifelse(x<90,x,ifelse(x<180,180-x,ifelse(x<270,270-x,ifelse(x<360,360-x,x))))
				out<-c(cos(x1*(pi/180)), sin(x1*(pi/180)))							
				out<-out*			c(ifelse(x<90,1,ifelse(x<180,-1,ifelse(x<270,-1,ifelse(x<360,1,1)))),ifelse(x<90,1,ifelse(x<180,1,ifelse(x<270,-1,ifelse(x<360,-1,1)))))
#This quadrant requires the inversion of the X and Y sides
				if(x>180 & x<270){out<-out[c(2,1)]}
				return(t(out))}))

### TEMPERATURE: Time integrated bearing The direction of change would change if the temporal tendency is to a reduction in the variable of interest
TempChng.Coef.OLS.fnc<-function(x){if(sum(is.na(x))==length(x)){out<-NA}
								   else{out<-as.numeric(coef(lm(x~I(1:length(x)),na.action=na.omit))[2])}
								   return(out)
								   }
ChngRastTmp <-calc(TempGrad.SN,fun= TempChng.Coef.OLS.fnc)

TempChng.Sig.OLS.fnc<-function(x){if(sum(is.na(x))==length(x)){out<-NA}
								  else{out<-as.numeric(coef(lm(x~I(1:length(x)),na.action=na.omit))[2])}
								  return(out)
								  }
PvalChngRastTmp <-calc(TempGrad.SN,fun= TempChng.Sig.OLS.fnc)

### PRECIPITATION: Time integrated bearing The direction of change would change if the temporal tendency is to a reduction in the variable of interest
PrecChng.Coef.OLS.fnc<-function(x){if(sum(is.na(x))==length(x)){out<-NA}
								   else{out<-as.numeric(coef(lm(x~I(1:length(x)),na.action=na.omit))[2])}
								   return(out)
								   }
ChngRastPrec <-calc(PresGrad.SN,fun= PrecChng.Coef.OLS.fnc)

PrecChng.Sig.OLS.fnc<-function(x){if(sum(is.na(x))==length(x)){out<-NA}
								  else{out<-as.numeric(coef(lm(x~I(1:length(x)),na.action=na.omit))[2])}
								  return(out)
								  }
PvalChngRastPrec <-calc(PresGrad.SN,fun= PrecChng.Sig.OLS.fnc)

### TEMPERATURE: Defining rotation of the bearing based on the temporal change
	Bearing.Temp.Chng<-Bearing.Temp
	Bearing.Temp.Chng[(ChngRastTmp[]<0) & (Bearing.Temp[]<=180)]<-Bearing.Temp.Chng[(ChngRastTmp[]<0) & (Bearing.Temp[]<=180)]+180
	Bearing.Temp.Chng[(ChngRastTmp[]<0) & (Bearing.Temp[]>180)]<-Bearing.Temp.Chng[(ChngRastTmp[]<0) & (Bearing.Temp[]>180)]-180
### PRECIPITATION:  Defining rotation of the bearing based on the temporal change
	Bearing.Prec.Chng<-Bearing.Prec
	Bearing.Prec.Chng[(ChngRastPrec[]<0) & (Bearing.Prec[]<=180)]<-Bearing.Prec.Chng[(ChngRastPrec[]<0) & (Bearing.Prec[]<=180)]+180
	Bearing.Prec.Chng[(ChngRastPrec[]<0) & (Bearing.Prec[]>180)]<-Bearing.Prec.Chng[(ChngRastPrec[]<0) & (Bearing.Prec[]>180)]-180


### TEMPERATURE: Determine the end point of a vector representing the Bearing 
	VecPosTimINtegTemp<-t(sapply(Bearing.Temp.Chng[],function(x){
						x1<-ifelse(x<90,x,ifelse(x<180,180-x,ifelse(x<270,270-x,ifelse(x<360,360-x,x))))
						out<-c(cos(x1*(pi/180)), sin(x1*(pi/180)))							
						out<-out*			c(ifelse(x<90,1,ifelse(x<180,-1,ifelse(x<270,-1,ifelse(x<360,1,1)))),ifelse(x<90,1,ifelse(x<180,1,ifelse(x<270,-1,ifelse(x<360,-1,1)))))
    #This quadrant requires the inversion of the X and Y sides
						if(x>180 & x<270){out<-out[c(2,1)]}
						return(t(out))}))

### PRECIPITATION:  Determine the end point of a vector representing the Bearing 
	VecPosTimINtegPrec<-t(sapply(Bearing.Prec.Chng[],function(x){
						x1<-ifelse(x<90,x,ifelse(x<180,180-x,ifelse(x<270,270-x,ifelse(x<360,360-x,x))))
						out<-c(cos(x1*(pi/180)), sin(x1*(pi/180)))							
						out<-out*			c(ifelse(x<90,1,ifelse(x<180,-1,ifelse(x<270,-1,ifelse(x<360,1,1)))),ifelse(x<90,1,ifelse(x<180,1,ifelse(x<270,-1,ifelse(x<360,-1,1)))))
    #This quadrant requires the inversion of the X and Y sides
						if(x>180 & x<270){out<-out[c(2,1)]}
						return(t(out))}))

#Temporal SPATIAL GRADIENT
	TempDiverg<-abs(Bearing.Temp.Chng-Bearing.Prec.Chng)
	TempDiverg <-calc(TempDiverg,fun=function(x){ifelse(x<=180,x,360-x)})




#dev.new(width=18,height=6)
#png("~/Desktop/Panel1.png",width=18,height=6,units="in",res=250,bg="transparent")
#par(mfrow=c(1,3),oma=c(0.5,0.5,0.5,1))
#TEMPERATURE
colfunc <- colorRampPalette(c("#053061","#92c5de","#f7f7f7","#f4a582","#67001f"))
png("~/Desktop/SuppMat_S1/Panel_1.png",width=6,height=6,units="in",res=250,bg="transparent")
plot(TempGrad.SN[[1]],main="Temperature",col=colfunc(100),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=c(0,15), legend.only=TRUE, legend.lab="Co", legend.line=2, legend.width = 1.5,col=colfunc(100))
dev.off()
#PRECIPITATION
png("~/Desktop/SuppMat_S1/Panel_2.png",width=6,height=6,units="in",res=250,bg="transparent")
plot(PresGrad.SN[[1]],main="Precipitation",col=rev(colfunc(100)),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=range(PresGrad.SN[[1]][]), col=rev(colfunc(100)),legend.only=TRUE, legend.lab="mm per year", legend.line=2, legend.width = 1.5)
dev.off()

#Local gradinet
png("~/Desktop/SuppMat_S1/Panel_3.png",width=6,height=6,units="in",res=250,bg="transparent")
plot(Diverg,main="Spatial angle",col=(colfunc(100)),zlim=c(0,180),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=c(0,180), col=(colfunc(100)),legend.only=TRUE, legend.lab="Degrees", legend.line=2.2, legend.width = 1.5)
#### TEMPERATURE: plot the bearing vector
arrows(x0= coordinates(TempGrad.SN)[,1],
		y0= coordinates(TempGrad.SN)[,2],
		x1=c(coordinates(TempGrad.SN)[,1]+c(VecPosTemp[,1]*0.025)),
		y1=c(coordinates(TempGrad.SN)[,2]+c(VecPosTemp[,2]*0.025)),
		length=0.025,
		angle=60,
		code=2,col="black",lwd=2)
		
### PRECIPITATION: plot the bearing vector
	arrows(x0= coordinates(PresGrad.SN)[,1],
			y0= coordinates(PresGrad.SN)[,2],
			x1=c(coordinates(PresGrad.SN)[,1]+c(VecPosPrec[,1]*0.025)),
			y1=c(coordinates(PresGrad.SN)[,2]+c(VecPosPrec[,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="dark grey",lwd=2)
dev.off()

### TEMPERATURE: Plot the Velocity of change - including the directionality of change
png("~/Desktop/SuppMat_S1/Panel_4.png",width=6,height=6,units="in",res=250,bg="transparent")
plot((ChngRastTmp/SpcHetTemp),main="Temperature - Magnitude",col= colfunc(100),zlim=c(-max(abs(ChngRastTmp/SpcHetTemp)[]),max(abs(ChngRastTmp/SpcHetTemp)[])),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=c(-max(abs(ChngRastTmp/SpcHetTemp)[]),max(abs(ChngRastTmp/SpcHetTemp)[])), col=(colfunc(100)),legend.only=TRUE, legend.lab="km * year", legend.line=2.5, legend.width = 1.5)

### TEMPERATURE: plot the bearing vector No change
	arrows(x0= coordinates(TempGrad.SN)[ChngRastTmp[]>=0,1],
			y0= coordinates(TempGrad.SN)[ChngRastTmp[]>=0,2],
			x1=c(coordinates(TempGrad.SN)[ChngRastTmp[]>=0,1]+c(VecPosTimINtegTemp[ChngRastTmp[]>=0,1]*0.025)),
			y1=c(coordinates(TempGrad.SN)[ChngRastTmp[]>=0,2]+c(VecPosTimINtegTemp[ChngRastTmp[]>=0,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="black",lwd=2)

### TEMPERATURE: plot the bearing vector Change
	arrows(x0= coordinates(TempGrad.SN)[ChngRastTmp[]<0,1],
			y0= coordinates(TempGrad.SN)[ChngRastTmp[]<0,2],
			x1=c(coordinates(TempGrad.SN)[ChngRastTmp[]<0,1]+c(VecPosTimINtegTemp[ChngRastTmp[]<0,1]*0.025)),
			y1=c(coordinates(TempGrad.SN)[ChngRastTmp[]<0,2]+c(VecPosTimINtegTemp[ChngRastTmp[]<0,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="purple",lwd=2)
legend("bottomleft",legend=c("Increase","Decrease"),col=c("black","purple"),lty=1,xpd=NA,inset=c(0,-0.125),ncol=2,lwd=2)
dev.off()

### PRECIPITATION: Plot the Velocity of change - including the directionality of change
png("~/Desktop/SuppMat_S1/Panel_5.png",width=6,height=6,units="in",res=250,bg="transparent")
plot((ChngRastPrec/SpcHetPrec),main="Precipitation - Magnitude ",col= rev(colfunc(100)),zlim=c(-max(abs(ChngRastPrec/SpcHetPrec)[]),max(abs(ChngRastPrec/SpcHetPrec)[])),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=c(-max(abs(ChngRastPrec/SpcHetPrec)[]),max(abs(ChngRastPrec/SpcHetPrec)[])), col= rev(colfunc(100)),legend.only=TRUE, legend.lab="km * year", legend.line=2.5, legend.width = 1.5)


### PRECIPITATION: plot the bearing vector No change
	arrows(x0= coordinates(PresGrad.SN)[ChngRastPrec[]>=0,1],
			y0= coordinates(PresGrad.SN)[ChngRastPrec[]>=0,2],
			x1=c(coordinates(PresGrad.SN)[ChngRastPrec[]>=0,1]+c(VecPosTimINtegPrec[ChngRastPrec[]>=0,1]*0.025)),
			y1=c(coordinates(PresGrad.SN)[ChngRastPrec[]>=0,2]+c(VecPosTimINtegPrec[ChngRastPrec[]>=0,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="black",lwd=2)

### PRECIPITATION: plot the bearing vector change
	arrows(x0= coordinates(PresGrad.SN)[ChngRastPrec[]<0,1],
			y0= coordinates(PresGrad.SN)[ChngRastPrec[]<0,2],
			x1=c(coordinates(PresGrad.SN)[ChngRastPrec[]<0,1]+c(VecPosTimINtegPrec[ChngRastPrec[]<0,1]*0.025)),
			y1=c(coordinates(PresGrad.SN)[ChngRastPrec[]<0,2]+c(VecPosTimINtegPrec[ChngRastPrec[]<0,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="purple",lwd=2)
legend("bottomleft",legend=c("Increase","Decrease"),col=c("black","purple"),lty=1,xpd=NA,inset=c(0,-0.125),ncol=2,lwd=2)
dev.off()


#Temporal Factor of change
png("~/Desktop/SuppMat_S1/Panel_6.png",width=6,height=6,units="in",res=250,bg="transparent")
DirChng<-raster(ChngRastTmp)
DirChng[]<-ifelse(((ChngRastTmp[]>0)&(ChngRastPrec[]>0)),0,ifelse(((ChngRastTmp[]<0)&(ChngRastPrec[]<0)),0,180))
plot(DirChng,main="Temporal angle",col= colfunc(100),zlim=c(0,180),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=c(0,180), col=(colfunc(100)),legend.only=TRUE, legend.lab="Degrees", legend.line=2.2, legend.width = 1.5)
dev.off()


### Plot the Mean displacment (Average speed) and the time integrated direction of the velocity vectors
png("~/Desktop/SuppMat_S1/Panel_7.png",width=6,height=6,units="in",res=250,bg="transparent")
colfunc <- colorRampPalette(c("#f7f7f7","#f4a582","#67001f"))
plot(mean(stack(abs(ChngRastTmp/SpcHetTemp),abs(ChngRastPrec/SpcHetPrec))),main="Displacement",col=(colfunc(100)),cex.main=2, axes=F,legend=F, zlim=c(0,0.06))
fields::image.plot( zlim=c(0,0.06), col=(colfunc(100)),legend.only=TRUE, legend.lab="km * year", legend.line=2.5, legend.width = 1.5)
dev.off()

### Plot the Divergence between gradients
png("~/Desktop/SuppMat_S1/Panel_8.png",width=6,height=6,units="in",res=250,bg="transparent")
colfunc <- colorRampPalette(c("#053061","#92c5de","#f7f7f7","#f4a582","#67001f"))
plot(TempDiverg,main="Temporal angle",col= colfunc(100),zlim=c(0,180),cex.main=2, axes=F,legend=F)
fields::image.plot( zlim=c(0,180), col=(colfunc(100)),legend.only=TRUE, legend.lab="Degrees", legend.line=2.2, legend.width = 1.5)
### TEMPERATURE:  plot the bearing vector No change
	arrows(x0= coordinates(TempGrad.SN)[,1],
			y0= coordinates(TempGrad.SN)[,2],
			x1=c(coordinates(TempGrad.SN)[,1]+c(VecPosTimINtegTemp[,1]*0.025)),
			y1=c(coordinates(TempGrad.SN)[,2]+c(VecPosTimINtegTemp[,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="black",lwd=2)

### PRECIPITATION: Plot the bearing vector No change
	arrows(x0= coordinates(PresGrad.SN)[,1],
			y0= coordinates(PresGrad.SN)[,2],
			x1=c(coordinates(PresGrad.SN)[,1]+c(VecPosTimINtegPrec[,1]*0.025)),
			y1=c(coordinates(PresGrad.SN)[,2]+c(VecPosTimINtegPrec[,2]*0.025)),
			length=0.025,
			angle=60,
			code=2,col="dark grey",lwd=2)
legend("bottomleft",legend=c("Temperature","Precipitation"),col=c("black","dark grey"),lty=1,xpd=NA,inset=c(0,-0.1),ncol=2,lwd=2)
dev.off()


# plot(Diverg-TempDiverg,main="Divergence",col=(colfunc(100)),zlim=c(-180,180),cex.main=2, axes=F,legend=F)
# fields::image.plot( zlim=c(-180,180),legend.only=TRUE, legend.lab="Degrees", legend.line=2, legend.width = 1.5,col=(colfunc(100)))

# #### TEMPERATURE: plot the bearing local vector
	# arrows(x0= coordinates(TempGrad.SN)[,1],
			# y0= coordinates(TempGrad.SN)[,2],
			# x1=c(coordinates(TempGrad.SN)[,1]+c(VecPosTimINtegTemp[,1]*0.025)),
			# y1=c(coordinates(TempGrad.SN)[,2]+c(VecPosTimINtegTemp[,2]*0.025)),
			# length=0.025,
			# angle=60,
			# code=2,col="black",lwd=2)
# #### PRECIPITATION: plot the Local vector precipitation
	# arrows(x0= coordinates(PresGrad.SN)[,1]+0.005,
			# y0= coordinates(PresGrad.SN)[,2],
			# x1=c(coordinates(PresGrad.SN)[,1]+c(VecPosTimINtegPrec[,1]*0.025)),
			# y1=c(coordinates(PresGrad.SN)[,2]+c(VecPosTimINtegPrec[,2]*0.025)),
			# length=0.025,
			# angle=60,
			# code=2,col="black",lwd=2)

# #### TEMPERATURE: plot the bearing Temporal vector
# arrows(x0= coordinates(TempGrad.SN)[,1],
		# y0= coordinates(TempGrad.SN)[,2],
		# x1=c(coordinates(TempGrad.SN)[,1]-c(VecPosTemp[,1]*0.025)),
		# y1=c(coordinates(TempGrad.SN)[,2]-c(VecPosTemp[,2]*0.025)),
		# length=0.025,
		# angle=60,
		# code=2,col="dark grey",lwd=2)
		
# ### PRECIPITATION: plot the Temporal vector
	# arrows(x0= coordinates(PresGrad.SN)[,1],
			# y0= coordinates(PresGrad.SN)[,2],
			# x1=c(coordinates(PresGrad.SN)[,1]-c(VecPosPrec[,1]*0.025)),
			# y1=c(coordinates(PresGrad.SN)[,2]-c(VecPosPrec[,2]*0.025)),
			# length=0.025,
			# angle=60,
			# code=2,col="dark grey",lwd=2)

# legend("bottomleft",legend=c("Spatial vectors","Temporal vectors"),col=c("black","dark grey"),lty=1,xpd=NA,inset=c(0,-0.1),ncol=2,lwd=2)
# dev.off()



