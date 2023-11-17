# R-Code-STATS-
This is the R code for executing ordinary kriging using the Swiss rainfall data
#Exploratory Spatial Data Analysis
library(geoR)
points(sic.100, borders=sic.borders,col="green")
points(sic.367, borders=sic.borders,col="red",add=TRUE)
plot.geodata(sic.100,bor=sic.borders)


# Variogram Analysis
library(geoR)
library(fields)
# Assuming you are trying to use variog() function with some arguments
vario.b <- variog(sic.100, option = c("bin", "cloud", "smooth"), bin.cloud = TRUE)
vario.c <- variog(sic.100, op="cloud")
bplot.xy(vario.c$u,vario.c$v, breaks=vario.b$u,col="grey80", lwd=2,cex=0.1,outline=FALSE)
 
library(geoR)
vario.ex <- variog(sic.100, bin.cloud=TRUE)
plot(vario.ex)
vario4<-variog4(sic.100)
plot(vario4,same=FALSE)

#Semi-variogram modeling
library(geoR)
vario.ex<- variog(sic.100,option="bin")
vario.sphe<-(variofit(vario.ex,cov.model= "spher", ini.cov.pars=c(15000,200)))
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp =c (2,1,0))

plot(vario.ex,main="Spherical")
lines.variomodel(cov.model="sphe",cov.pars=c(15000,100), nug=0,max.dist=350)

plot(vario.ex,main="Exponential")
lines.variomodel(cov.model="exp",cov.pars=c(15000,100), nug=0,max.dist=350)

plot(vario.ex,main="Exponential with nugget effect")
lines.variomodel(cov.model="exp",cov.pars=c(10000,100), nug=5000,max.dist=350)

plot(vario.ex,main="Matern")
lines.variomodel(cov.model="matern",cov.pars=c(10000,100), nug=0,max.dist=350,kappa=0.5)


#Ordinary Kriging
#Spherical variogram
library(geoR)
vario.ex<- variog(sic.100, bin.cloud=TRUE)
plot(vario.ex,main="")
lines.variomodel(cov.model="spher",cov.pars=c(15000,50), nug=0,max.dist=300)

#raw data
library(geoR)
pred.grid <- expand.grid(seq(0,350, l=51),seq (0,220, l=51))
rgb.palette <- colorRampPalette(c("blue", "lightblue",
                                  "orange", "red"),space = "rgb")
kc<- krige.conv(sic.100, loc = pred.grid,
                krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
image(kc, loc = pred.grid,col =rgb.palette(20) ,xlab="Coord X",
      ylab="Coord Y",borders=sic.borders,main="Estimation")
krige.var <- kc$krige.var
image(kc, krige.var,loc = pred.grid,col=rgb.palette(20),
      xlab="Coord X",ylab="Coord Y",borders=sic.borders,
      main="Kriging variance")

library(geoR)
kc1<- krige.conv(sic.100, loc = sic.100$coords,
                 krige=krige.control(cov.model="spherical",cov.pars=c(16000,47)))
kc2<- krige.conv(sic.100, loc = sic.367$coords,
                 krige=krige.control(cov.model="spherical",cov.pars=c(16000,47)))
plot(sic.100$data,kc1$predict,xlab="Observed",ylab="Estimated",
     main="Control sample")
abline(a=0,b=1,col="red")
plot(sic.367$data,kc2$predict,xlab="Observed",ylab="Estimated",
     main="Control")
abline(a=0,b=1,col="red")

## CALCULATING RMSE
# Calculate RMSE for kriging estimates on the observed data
rmse_sample <- sqrt(mean((sic.100$data - kc1$predict)^2))
rmse_control <- sqrt(mean((sic.367$data - kc2$predict)^2))

# Print RMSE
cat("RMSE (Sample):", rmse_sample, "\n")
cat("RMSE (Control sample):", rmse_control, "\n")


#Box-Cox Transformation 
library(geoR)
plot.geodata(sic.100,bor=sic.borders,lambda=0.5)

library(geoR)
vario.ext<- variog(sic.100,option="bin",lambda=0.5)
plot(vario.ext)
lines.variomodel(cov.m = "matern",cov.pars =c (105, 36), nug = 6.9,
                 max.dist = 300,kappa = 1, lty = 1)

##Kriging
library(geoR)
kct<- krige.conv(sic.100, loc = pred.grid,
                 krige=krige.control(cov.model="matern",cov.pars=c(105, 36),
                                     kappa=1,nugget=6.9,lambda=0.5))
pred.grid <- expand.grid(seq(0,350, l=51),seq (0,220, l=51))
rgb.palette <- colorRampPalette(c("blue", "lightblue",
                                  "orange", "red"),space = "rgb")
image(kct, loc = pred.grid,col =rgb.palette(20) , xlab="Coord X",
      ylab="Coord Y",borders=sic.borders,main="Estimation")


image(kct, krige.var,loc = pred.grid,col =rgb.palette(20) ,
      xlab="Coord X",ylab="Coord Y",borders=sic.borders,
      main="Kriging variance")


legend_colors <- rgb.palette(20)
legend_breaks <- seq(min(kct$krige.var), max(kct$predict), length.out = length(legend_colors) + 1)
image.plot(x = NULL, y = NULL, z = matrix(legend_breaks, nrow = 1), col = rgb.palette(20),
           axes = FALSE, legend.only = TRUE, horizontal = FALSE, legend.width = 1, legend.shrink = 0.8,
           legend.mar = 5, labels("Legend"))
           



library(geoR)
kct1<- krige.conv(sic.100, loc = sic.100$coords,
                  krige=krige.control(cov.model="matern",cov.pars=c(16000,47),
                                      kappa=1,nugget=6.9,lambda=0.5))
kct2<- krige.conv(sic.100, loc = sic.367$coords,
                  krige=krige.control(cov.model="matern",cov.pars=c(16000,47),
                                      kappa=1,nugget=6.9,lambda=0.5))
plot(sic.100$data,kct1$predict,xlab="Observed",ylab="Estimated",
     main="Sample")
abline(a=0,b=1,col="red")
plot(sic.367$data,kct2$predict,xlab="Observed",ylab="Estimated",
     main="Control sample")
abline(a=0,b=1,col="red")

# Calculate RMSE for kriging estimates on the observed data
rmse_sample <- sqrt(mean((sic.100$data - kct1$predict)^2))
rmse_control <- sqrt(mean((sic.367$data - kct2$predict)^2))

# Print RMSE
cat("RMSE (Sample):", rmse_sample, "\n")
cat("RMSE (Control sample):", rmse_control, "\n")
