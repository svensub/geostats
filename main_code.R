rm(list = ls())

################################################################################
#                         HANDLING SPATIAL DATA                                #
################################################################################

# step 1: set working directory

setwd("")

# step 2: download relevant libraries

library(PrevMap)
library(sf)
library(tmap)
library(raster)
library(ggplot2)
library(dplyr)

# step 3: load the data

# prevalence data

tz <- read.csv("Remo_TANZANIA.csv")

# spatial data

tz.adm0 <- st_read("TZA_adm/geoBoundaries-TZA-ADM0.shp")
tz.adm0 <- st_transform(tz.adm0,crs=3857)

# covariates data, keep only tanzania

all_cov <- readRDS("covariates_10km.rds")
levels(all_cov$ADMIN0)

cov <- subset(all_cov, ADMIN0== "Tanzania (Mainland)"|ADMIN0== "Tanzania (Zanzibar)") 

cov$logprecip <- log(cov$AnnualPrecip)
cov$logPopDens2020adj <- log(cov$PopDens2020adj)

# step 4: create prevalence

tz$prev <- tz$y/tz$n

# step 5: create e.logit

tz$e.logit <- log((tz$y+0.5)/(tz$n-tz$y+0.5))

# step 6: convert to utm_x/utm_y + scaling

tz.sf <- st_as_sf(tz,coords=c("lat","long"),crs=4326)
tz.sf <- st_transform(tz.sf,crs=3857)
tz$web_x <- st_coordinates(tz.sf)[,1]
tz$web_y <- st_coordinates(tz.sf)[,2]

# step 7: prevalence map

map0 <- tm_shape(tz.adm0) + 
  tm_borders(lwd=3) 
map0

map0+tm_shape(tz.sf)+tm_dots(size=0.5)

map.with.points <- map0+tm_shape(tz.sf)+
  tm_bubbles("prev", col = "prev", 
             border.col = "black",
             style="fixed", 
             breaks=seq(0,0.4,0.05),
             palette="-RdYlBu",
             title.size="Prevalence", 
             scale = 1,
             title.col="") 

map.with.points+
  tm_compass(type="8star", 
             position = c("right","top"))+
  tm_scale_bar(breaks = c(0,100,200),text.size=1,
               position=c("center","bottom")) 

tmap_mode("view")

map.with.points
tmap_mode("plot")

# step 7: include waterways

tz.wl <- st_read("TZA_wat/hotosm_tza_waterways_lines.shp")
tz.wl <- st_transform(tz.wl,crs = 3857)

map.with.points+tm_shape(tz.wl)+
  tm_lines(col="blue", palette="dodgerblue3", 
           title.col="Waterways")

# step 8: create grid
tz.grid <- st_make_grid(tz.adm0,
                        cellsize = 10000,
                        what="centers")

# step 9: subset to include only grid locations that fall inside Tanzania
tz.inout <- st_intersects(tz.grid,
                          tz.adm0,
                          sparse = FALSE)
tz.grid <- tz.grid[tz.inout]

# step 10: create raster of log precipitation

logprecip <- rasterFromXYZ(cbind(cov[,c("x","y","logprecip")]), 
                           crs=CRS("+init=epsg:3857"))

tz.logprecip <- mask(logprecip, as(tz.adm0,"Spatial"))

tz.logprecip <- extract(logprecip,
                        st_coordinates(tz.grid)) # predictions for all of tz

tz.logprecip.data <- extract(logprecip,
                             st_coordinates(tz.sf))# logprecip for specified loci

length(which(is.na(tz.logprecip)))

length(which(is.na(tz.logprecip.data)))

ind.na <- which(is.na(tz.logprecip))
tz.grid <- tz.grid[-ind.na] # remove NAs
tz.logprecip <- tz.logprecip[-ind.na]

ind.na <- which(is.na(tz.logprecip.data))
tz.logprecip.data <- tz.logprecip.data[-ind.na]

plot(tz.logprecip)
plot(tz.logprecip.data)

# step 11: plot with all key attributes

tm_shape(logprecip)+
  tm_raster(title="Annual precipitation (log)")+
  tm_shape(tz.adm0)+
  tm_borders(lwd=2)+
  tm_shape(tz.sf)+
  tm_bubbles("prev", col = "prev", 
             border.col = "black",
             style="fixed", 
             breaks=seq(0,0.4,0.05),
             palette="-RdYlBu",
             title.size="Prevalence", 
             scale = 1,
             title.col="")+
  tm_scale_bar(breaks = c(0,100,200),text.size=1,
               position=c("left","bottom")) +
  tm_layout(legend.outside = TRUE) +
  tm_compass(type="8star", 
             size = 2,
             position = c("right","top"))

# step 12: extract values to tz to plot

tz$logprecip <- extract(logprecip,tz.sf)

# step 13: scatter plots of empirical logit against log precipitation

tz <- na.omit(tz)

plot(tz$logprecip,tz$prev,
     xlab="log(Annual Precipitation)", ylab="Prevalence")

plot(tz$logprecip,tz$e.logit,
     xlab="log(Annual Precipitation)",
     ylab="Empirical logit")     

plot.lpr <- ggplot(tz, aes(x = logprecip, 
                           y = e.logit)) + geom_point() +
  labs(x="log Annual Precipitation", y="Empirical logit") +
  stat_smooth(method = "gam", formula = y ~ s(x), se=FALSE) +
  stat_smooth(method = "lm", formula = y ~ x ,
              col="red",lty="dashed",se=FALSE) 

plot.lpr

# step 14: assume model is binomial model, with log (p/1-p) = beta0 + beta1*d_i
# + beta2max{d_i-c,0}, where d_i is log precip and c is the value of precip (log) 
# which corresponds to the change in slope of a linear spline. use glm to find 
# which values of c 

glmx1 <- glm(cbind(y,n-y) ~ logprecip + 
               I((logprecip-7)*(logprecip>7)),
             data=tz, family = binomial)

glmx2 <- glm(cbind(y,n-y) ~ logprecip + 
               I((logprecip-7.2)*(logprecip>7.2)),
             data=tz, family = binomial)

glmx3 <- glm(cbind(y,n-y) ~ logprecip + 
               I((logprecip-7.4)*(logprecip>7.4)),
             data=tz, family = binomial)

logLik(glmx1)
logLik(glmx2)
logLik(glmx3)

# step 15: create plot with chosen model + a quadratic model

plot.lpr + stat_smooth(method = "lm", formula = y ~ x + I((x-7.2)*(x>7.2)),
                       col="green",lty="dashed",se=FALSE)

################################################################################
#                                 ANALYSIS                                     #
################################################################################

# model 1: linear model

lm.fit <- lm(e.logit ~ logprecip, data=tz)
predictors.tz <- data.frame(logprecip=tz.logprecip)
pred.lm <- 1/(1+exp(-predict(lm.fit,newdata=predictors.tz)))
tz.pred.lm <- rasterFromXYZ(cbind(st_coordinates(tz.grid)/1000,pred.lm))

plot(tz.pred.lm)

# prior to building linear geostats model, plot point map and variogram

point.map <- point.map(tz,~logprecip,coords=~I(web_x/1000)+I(web_y/1000),
                       pt.divide="quintiles")

## point.map + tm_shape(tz.adm0) + tm_borders(lwd=3)

vari <- variogram(tz, e.logit~logprecip,
                  coords=~I(web_x/1000)+I(web_y/1000),
                  uvec=seq(10,150,length=15))

plot(vari,type="b")

# model 2: linear geostatistical model

## first we fit linear geostatistical model to the empirical logit

spat.corr.diagnostic(e.logit~1,
                     data=tz,
                     coords=~I(web_x/1000)+I(web_y/1000),
                     likelihood = "Gaussian",
                     lse.variogram = TRUE)

sigma2.guess <- 5.306937e+05
phi.guess <- 1.026837e+08
tau2.guess <- 1.769597e+00 

fit.mle <- linear.model.MLE(e.logit ~ 1,
                            coords=~I(web_x/1000)+I(web_y/1000),
                            kappa=0.5,
                            start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                            data=tz, method="nlminb")

summary(fit.mle,log.cov.pars=FALSE)

## we then predict nodule prevalence across Tanzania + show the exceedance 
## probability for a 20% threshold.
pred.mle.lm <- spatial.pred.linear.MLE(fit.mle,grid.pred = st_coordinates(tz.grid)/1000,
                                       standard.errors = TRUE,
                                       scale.predictions = "prevalence",
                                       n.sim.prev = 1000,
                                       thresholds = 0.2, scale.thresholds = "prevalence")

plot(pred.mle.lm,"prevalence","predictions")
plot(pred.mle.lm,"prevalence","standard.errors")
plot(pred.mle.lm,summary="exceedance.prob")

## then we introduce log-transformed annual precipitation as a linear predictor

spat.corr.diagnostic(e.logit~logprecip,
                     data=tz,
                     coords=~I(web_x/1000)+I(web_y/1000),
                     likelihood = "Gaussian",
                     lse.variogram = TRUE)

sigma2.guess <- 1.947061
phi.guess <- 77.983029
tau2.guess <- 1.237006

fit.mle.lpr <- linear.model.MLE(e.logit~logprecip,
                                coords=~I(web_x/1000)+I(web_y/1000),
                                kappa=0.5,
                                start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                                data=tz, method="nlminb")

summary(fit.mle.lpr,log.cov.pars=FALSE)
summary(fit.mle.lpr,log.cov.pars=TRUE)

exp(5.51459)
3*exp(5.51459)

pred.mle.lm.lpr <- 
  spatial.pred.linear.MLE(fit.mle.lpr,grid.pred = st_coordinates(tz.grid)/1000,
                          predictors = predictors.tz,
                          scale.predictions = "prevalence",
                          thresholds = 0.2,n.sim.prev = 1000,
                          scale.thresholds = "prevalence")

plot(pred.mle.lm.lpr,"prevalence","predictions")
plot(pred.mle.lm.lpr,summary="exceedance.prob")

# model 3: Binomial statistical model

## here we consider a binomial model for prevalence, with log-transformed annual
## precipitation as a linear predictor

c.mcmc <- control.mcmc.MCML(n.sim = 10000,
                            burnin=2000,
                            thin=8)

glm.fit <- glm(cbind(y, n-y) ~ logprecip, data=tz, family = binomial)
par0 <- c(coef(glm.fit), 2.3, 248, 1.2) 
fit.bin.lpr <- binomial.logistic.MCML(y ~ logprecip,
                                      units.m = ~ n,
                                      coords=~I(web_x/1000)+I(web_y/1000),
                                      kappa=0.5, control.mcmc = c.mcmc,
                                      par0=par0, 
                                      start.cov.pars = c(248, 1.2), 
                                      data=tz, method="nlminb")
summary(fit.bin.lpr)

pred.mle.bin.lpr <- 
  spatial.pred.binomial.MCML(fit.bin.lpr, grid.pred = st_coordinates(tz.grid)/1000,
                             control.mcmc = c.mcmc,
                             predictors = predictors.tz,
                             scale.predictions = "prevalence",
                             thresholds = 0.2,
                             scale.thresholds = "prevalence")

par(mfrow=c(1,2))
plot(pred.mle.lm.lpr$prevalence$predictions,
     pred.mle.bin.lpr$prevalence$predictions,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

plot(pred.mle.lm.lpr$exceedance.prob,
     pred.mle.bin.lpr$exceedance.prob,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

library("splancs")

par(mfrow = c(1, 2))
plot(pred.mle.bin.lpr, type = "prevalence", summary = "predictions",
     zlim = c(0, 0.45),
     main = "Prevalence predictions")
contour(pred.mle.bin.lpr, type = "prevalence", summary = "predictions",
        levels = c(0.05, 0.1, 0.2, 0.3), add = TRUE)

plot(pred.mle.bin.lpr, summary = "exceedance.prob",
     zlim = c(0,1),
     main = "Exceedance probabilities")
contour(pred.mle.bin.lpr, summary = "exceedance.prob",
        levels = c(0.25, 0.75), add = TRUE)


variog.diagnostic.glgm(
  fit.bin.lpr,
  n.sim = 1000,
  uvec = NULL,
  plot.results = TRUE,
  which.test = "variogram"
)

exp(5.90618)
3*exp(5.90618)