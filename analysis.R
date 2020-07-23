# -------------------------------------------------------------------------- #
# This script is meant to reproduce the analysis in the manuscript:          #
#                                                                            # 
# Current flood-rich period exceptional compared to past 500 years in Europe #
#                                                                            #
#   by Bloeschl et al.                                                       #
# -------------------------------------------------------------------------- #

# Summary
# PART 1: INTERPOLATION
# PART 2: FLOOD PERIOD SELECTION AND RANKING
# PART 3: TEMPERATURE AND SEASONALITY ANALYSES


# --------------------------------------------------------------------------------- #
# PART 1: INTERPOLATION
#  i.e., from event intensities i_f, to annual intensities i_a, to interpolated intensities i_i

# Load the data
sites <- read.csv('sites.csv')
i_f <- read.csv('i_f.csv')
bias <- read.csv('bias.csv')

# Event intensities
# I consider only class 2 and 3 floods for the interpolation
i_f$numIndex <- c(NA, 2, 3, NA, NA, NA)[as.numeric(as.factor(i_f$Index))]

# Annual intensities
# Flood-intensities of multiple flood events per year aggregated via euclidean norm
i_a <- matrix(NA, nrow=517, ncol=103)
rootsumsquares <- function(x) {sqrt(sum(x^2))}
 rownames(i_a) <- seq(1500, 2016)
 colnames(i_a) <- unique(i_f$Code)
for (j in 1:ncol(i_a)) {
 series <- i_f[i_f$Code == colnames(i_a)[j], ]
 anni <- series$Year
 cehflood <- series$numIndex
 dummy <- tapply(cehflood, anni, rootsumsquares)
 i_a[rownames(i_a) %in% as.character(anni), j] <- dummy
}

# Project the sites in a km space and create the grid for the interpolation:
library(sp)
coor <- array(cbind(sign((sites$lon4 == 'O') - 0.5)*(sites$lon1 + sites$lon2/60 + sites$lon3/3600),
                    sites$lat1 + sites$lat2/60 + sites$lat3/3600),
              dim=c(length(sites$code), 2), dimnames=list(sites$code, c('lon', 'lat')))
coordinates(sites) <- coor
 proj4string(sites) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
lonC = 7
latC = 51
sites_km103 <- spTransform(sites, CRS(paste('+proj=aeqd +lat_0=', latC, ' +lon_0=', lonC, ' +units=km', sep='')))
sites_km <- sites_km103[sites_km103$Interp == 1,]  # only the sites that are used for interpolation
library(raster)
library(rgeos)
bbbox <- as(extent(as.vector(t(matrix(c(-1450, -1600, 1150, 1500), ncol=2)))), 'SpatialPolygons')
 proj4string(bbbox) <- paste('+proj=aeqd +lat_0=', latC, ' +lon_0=', lonC, ' +units=km', sep='')
tr <- seq(1500, 2018, by=4)
kmperyear = 50
xr <- seq(bbox(bbbox)[1,1], bbox(bbbox)[1,2], length=length(tr)/2)/kmperyear
yr <- seq(bbox(bbbox)[2,1], bbox(bbbox)[2,2], length=length(tr)/2)/kmperyear
griglia <- expand.grid(x=xr, y=yr, t=tr)
pixelsize <- diff(seq(bbox(bbbox)[1,1], bbox(bbbox)[1,2], length=length(tr)/2))[1] * diff(seq(bbox(bbbox)[2,1], bbox(bbbox)[2,2], length=length(tr)/2))[1]  # km^2
voxelsize <- pixelsize*diff(tr)[1]  # km^2*yr

# Interpolation
pow0 = 10
scaleU = 2
powU = 2
kappaW = c(0.2,1,1.5)
tapering = 20
smoothing = 10
# Note that the following loop takes 1 hour and 40 min on a INTEL COREi7 vPro 8th Gen processor
library(fields)
timestamp()
for (seed in 11:60) {
 cat('\n\n****************************************************************\n************************ seed n.', seed, '****************************\n****************************************************************\n\n')
 set.seed(seed)  # for repeatability
 timewindow=100  # one century
 i_a_add0 <- i_a
 for (j in 1:ncol(i_a_add0)) {
  # bias correction
  numIndex <- i_a_add0[,j]
  firstlastfloods <- c(1, length(numIndex))
   isfloodgr2 <- numIndex; isfloodgr2[is.na(isfloodgr2)] <- 0
   isfloodgr2[isfloodgr2 <= 1] <- 0; isfloodgr2[isfloodgr2 > 1] <- 1
  floodfrequency <- filter(isfloodgr2, filter=rep(1, timewindow))/timewindow
   # to remove the NA
   dummy <- range(which(!is.na(floodfrequency)))
   floodfrequency[1:dummy[1]] <- floodfrequency[dummy[1]]
   floodfrequency[dummy[2]:length(floodfrequency)] <- floodfrequency[dummy[2]]
    floodfrequency[floodfrequency < 1/timewindow] <- 1/timewindow  # to avoid 0s
  # sample from Bernoulli distribution with probability proportional to floodfrequency
  dummyi <- rep(1, length(numIndex)); dummyi[bias[,j+1] < 1.1] <- 0  # no no-flood where there is no info
  probab0 <- floodfrequency*dummyi;
  probab0 <- 1 - (1 - probab0)^pow0
   probab0[probab0 > 1] <- 1
  add0 <- as.logical(rbinom(n=length(numIndex), size=1, prob=probab0))
   dummy <- numIndex[add0]; dummy[is.na(dummy)] <- 0
  numIndex[add0] <- dummy
  # remove 0s outside the observational period
  if (firstlastfloods[1] > 1) {
   numIndex[1:(firstlastfloods[1] - 1)] <- NA
  }
  if (firstlastfloods[2] < length(numIndex)) {
   numIndex[(firstlastfloods[2] + 1):length(numIndex)] <- NA
  }
  i_a_add0[,j] <- numIndex
 }

 # prepare input for fastTps
 coor <- coordinates(sites_km)
 anni <- as.integer(rownames(i_a_add0))
 tabellone_km_add0 <- NULL
 for (j in 1:length(sites_km)) {
  codice <- as.character(sites_km$code[j])
  dummy <- data.frame(code=codice, x=coor[codice,1], y=coor[codice,2], yr=anni, flood=i_a_add0[,codice],
                      u=sites_km$Repres[j])
  tabellone_km_add0 <- rbind(tabellone_km_add0, dummy[!is.na(dummy$flood),])
 }
  coordinates(tabellone_km_add0) = ~x+y
  proj4string(tabellone_km_add0) <- paste('+proj=aeqd +lat_0=', latC, ' +lon_0=', lonC, ' +units=km', sep='')

 weights00 <- (tabellone_km_add0$u/scaleU)^powU
 weights00 <- weights00 * as.numeric(as.character(
                           cut(tabellone_km_add0$flood,
                               c(0,1.5,2.5,10), include.lowest=TRUE,
                               labels=kappaW)))
 # Thin plate spline interpolation for the 3D space-time cube:
 tps_intrp <- fastTps(x=cbind(coordinates(tabellone_km_add0)[,1]/kmperyear,
                              coordinates(tabellone_km_add0)[,2]/kmperyear,
                              tabellone_km_add0$yr),
                      Y=tabellone_km_add0$flood,
                      m=2, theta=tapering, lambda=smoothing,
                      weights=weights00)
 out_intrp <- predict(tps_intrp, griglia)  # this takes time and memory
 intrp_cube <- cbind(griglia, round(out_intrp, 4))
  names(intrp_cube)[4] <- 'floodIntensity'
 save(intrp_cube, file=paste('intrp_cubeSEED', seed, '.RData', sep=''), compress='xz')
 rm(intrp_cube)
 gc()  # to free memory
}
timestamp()

# Postprocessing
load('intrp_cubeSEED11.RData')
intrp_cube_50 <- intrp_cube
for (seed in 12:60) {
 load(paste('intrp_cubeSEED', seed, '.RData', sep=''))
 intrp_cube_50 <- cbind(intrp_cube_50, intrp_cube[,4])
}
names(intrp_cube_50)[-c(1:3)] <- paste('flInt', 11:60, sep='')
save('intrp_cube_50', file='intrp_cube_50.RData', compress='xz')  # this takes some time
# after saving on intrp_cube_50.RData, one can remove the 50 intrp_cubeSEEDXX.RData files

# Averaging of the seeds:
i_i_cube <- cbind(intrp_cube_50[,1:3], i_i=apply(intrp_cube_50[,-c(1:3)], 1, median))
# this is the space time interpolated flood intensity



# --------------------------------------------------------------------------------- #
# PART 2: FLOOD PERIOD SELECTION AND RANKING
# i.e., from interpolated intensities i_i to the 10 biggest flood periods

# Event identification
thr_q = quantile(i_i_cube$i_i, prob=.95)
thr_a = 10  # 1 pixel is 40.625km*48.4375km*4yr
minoverlap = 0.15
minfraction = 0.05

yrs <- unique(i_i_cube$t)
x <- unique(i_i_cube$x)*kmperyear
y <- unique(i_i_cube$y)*kmperyear
matrice <- array(i_i_cube$i_i, dim=c(length(x), length(y), length(yrs)), dimnames=list(round(x, 1), round(y, 1), yrs))
# Algorithm to connect neighbouring voxels from Haslinger and Bloeschl (2017, doi:10.1002/2017WR020797)
if (exists('df01')) rm(df01)
for (ii in 1:length(yrs)) {
 dat1 <- list()
 dat1$x <- x
 dat1$y <- y
 dat1$z <- matrice[,,ii]
 r <- raster(dat1)
 r[r < thr_q] <- NA
 cl <- clump(r)
 n_cl <- data.frame(freq(cl))
 ind <- which(is.na(n_cl$value) == FALSE & n_cl$count >= thr_a, arr.ind=TRUE)
 if (length(ind) >= 1) {
  for (i in 1:length(ind)) {
   if (i == 1) {
    idc <- Which(cl == ind[i], cells=TRUE)
    d <- data.frame(id=rep(i, length(idc)),
                    yr=rep(yrs[ii], length(idc)),
                    x=coordinates(r)[,1][idc],
                    y=coordinates(r)[,2][idc],
                    igp=idc,
                    value=r[idc])
   } else {
    idc = Which(cl == ind[i], cells=TRUE)
    d1 <- data.frame(id=rep(i, length(idc)),
                     yr=rep(yrs[ii], length(idc)),
                     x=coordinates(r)[,1][idc],
                     y=coordinates(r)[,2][idc],
                     igp=idc,
                     value=r[idc])
                     d<-rbind.data.frame(d, d1)
   }
  }
  if (!exists('df01')) { # first time step with detection
   df01 <- d
  } else {
   d$id <- d$id + max(df01$id)
   df01 <- rbind.data.frame(df01, d)
  }
 }
}
events3 <- df01

# tracking flood areas for identifying coherent flood-rich periods in space-time, 
# the method of Haslinger and Bloeschl (2017, doi:10.1002/2017WR020797) is employed
time <- yrs
for (t in 2:length(time)) {
 ind0 <- which(events3$yr == time[t-1])
 ind1 <- which(events3$yr == time[t])
 if (length(ind0) > 0 & length(ind1) > 0) {
  df0 <- events3[ind0,]
  df1 <- events3[ind1,]
  int <- intersect(df0$igp, df1$igp)
  ev0 <- unique(df0$id[df0$igp %in% int])
  ev1 <- unique(df1$id[df1$igp %in% int])
  if (length(ev0) > length(ev1)) {
   ev1 <- NULL
   for (i in 1:length(ev0)) {
    ev1[i] <- df1$id[df1$igp %in% intersect(df0$igp[df0$id==ev0[i]], df1$igp)][1]
   }
  } else {
   ev0<-NULL
   for (i in 1:length(ev1)) {
    ev0[i] <- df0$id[df0$igp %in% intersect(df1$igp[df1$id==ev1[i]], df0$igp)][1]
   }
  }

  for (i in 1:length(ev1)) {
   size0 <- length(df0$id[df0$id == ev0[i]])
   size1 <- length(df1$id[df1$id == ev1[i]])
   if (size0 >= size1) fraction <- size1/size0 else fraction <- size0/size1
   if (size0 >= size1) {
    overlap <- length(intersect(df0$igp[df0$id == ev0[i]], df1$igp[df1$id == ev1[i]]))/size1
    mean_of_larger <- mean(df0$value[df0$igp %in% df1$igp])
   } else {
    overlap <- length(intersect(df0$igp[df0$id == ev0[i]], df1$igp[df1$id == ev1[i]]))/size0
    mean_of_larger <- mean(df1$value[df1$igp %in% df0$igp])
   }
   if (is.na(overlap) == FALSE) {
    if ((fraction >= minfraction & overlap >= minoverlap)) {
     events3$id[events3$id == ev1[i]] <- ev0[i]
    }
   }
  }
 }
}

# compiling event table
events <- events3
un <- unique(events$id)
nt<-length(yrs)
nx<-r@ncols
ny<-r@nrows

for (i in 1:length(un)) {
 ind <- which(events$id == un[i])
 ev <- events[ind,]

 event_a <- data.frame(id=ev$id[1])   # id
 event_a$yr_st <- ev$yr[1]            # yr_start
 event_a$yr_en <- ev$yr[length(ind)]  # yr_end
 event_a$dur <- (event_a$yr_en - event_a$yr_st) + 1   # duration

 event_a$sumstrength <- round(sum(ev$value), 1)
 event_a$maxstrength <- round(max(ev$value), 3)
 event_a$meanstrength <- round(mean(ev$value), 3)
 event_a$maxarea <- length(unique(ev$igp))
 event_a$volume <- length(ev$value)

 if (i > 1) event_c <- rbind.data.frame(event_c, event_a) else event_c <- event_a
}
events <- event_c
events$scaledmeanstrength <- round((events$meanstrength - min(events$meanstrength))/(max(events$meanstrength) - min(events$meanstrength)), 3)
events$scaledvolume <- round((events$volume - min(events$volume))/(max(events$volume) - min(events$volume)), 3)
events$maxarea_Mkm2 <- round(events$maxarea*pixelsize*1e-6, 3)
events$volume_Mkm2yr <- round(events$volume*voxelsize*1e-6, 2)

# Selection and ranking of the largest events:
criterio <- events$scaledmeanstrength + events$scaledvolume
max10events <- events[order(criterio, decreasing=TRUE),][1:10,]
 max10events$rank <- 1:10
max10events3 <- merge(events3, max10events, all.y=TRUE)
max10events$name <- NA; max10events[order(max10events$yr_st),'name'] <- c('I','II','III','IV','Va','Vb','VI','VII','VIII','IX')
 print(max10events[order(max10events$yr_st),], row.names=FALSE)
# this should be consistent with Extended Data Table 2 in the paper



# --------------------------------------------------------------------------------- #
# PART 3: TEMPERATURE AND SEASONALITY ANALYSES

# 3.1) Temperature anomalies within and outside of flood-rich periods
 
# Load the temperature data
CEtempanom500yEXT <- read.csv('Delta_T.csv')

eventsXregion <- data.frame(c=c(T,T,F,T,T,F,F,T,F,T), s=c(F,F,T,F,T,F,T,F,F,T), w=c(T,T,F,T,T,F,T,F,F,T), n=c(F,F,F,F,F,T,F,F,T,F))
 rownames(eventsXregion) <- c('I','II','III','IV','Va','Vb','VI','VII','VIII','IX')
CEtempanom500y_periods <- data.frame(CEtempanom500yEXT, c=NA, s=NA, w=NA, n=NA)
for (reg in c('c','s','w','n')) {
 dummy <- max10events[order(max10events$id),][eventsXregion[,reg],]
 for (i in 1:nrow(dummy)) {
  if ((i == 1)&(dummy$yr_st[i] > 1500)) CEtempanom500y_periods[CEtempanom500y_periods$year < dummy$yr_st[i], reg] <- paste('STto', dummy$name[i], sep='')
  CEtempanom500y_periods[(CEtempanom500y_periods$year >= dummy$yr_st[i])&(CEtempanom500y_periods$year <= dummy$yr_en[i]), reg] <- dummy$name[i]
  if (i < nrow(dummy)) CEtempanom500y_periods[(CEtempanom500y_periods$year > dummy$yr_en[i])&(CEtempanom500y_periods$year < dummy$yr_st[i+1]), reg] <- 
                       paste(dummy$name[i], 'to', dummy$name[i+1], sep='')
  if (i == nrow(dummy)) CEtempanom500y_periods[CEtempanom500y_periods$year > dummy$yr_en[i], reg] <- paste(dummy$name[i], 'toEND', sep='')
 }
}
# central EU
floodperiodmeananomaly_c <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$c, mean)[c('I','II','IV','Va','VII','IX')]
floodperiodsdanomaly_c <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$c, sd)[c('I','II','IV','Va','VII','IX')]
floodperiodlength_c <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$c, length)[c('I','II','IV','Va','VII','IX')]
nofloodperiodmeananomaly_c <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$c, mean)[c('ItoII','IItoIV','IVtoVa','VatoVII','VIItoIX','IXtoEND')]
nofloodperiodsdanomaly_c <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$c, sd)[c('ItoII','IItoIV','IVtoVa','VatoVII','VIItoIX','IXtoEND')]
nofloodperiodlength_c <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$c, length)[c('ItoII','IItoIV','IVtoVa','VatoVII','VIItoIX','IXtoEND')]
 # the following values are used in Fig. 4a 
 print(round(floodperiodmeananomaly_c, 3))
 print(round(floodperiodmeananomaly_c - 1.645*floodperiodsdanomaly_c/sqrt(floodperiodlength_c), 3))  # lower confidence bound
 print(round(floodperiodmeananomaly_c + 1.645*floodperiodsdanomaly_c/sqrt(floodperiodlength_c), 3))  # upper confidence bound
 print(round(nofloodperiodmeananomaly_c, 3))
 print(round(nofloodperiodmeananomaly_c - 1.645*nofloodperiodsdanomaly_c/sqrt(nofloodperiodlength_c), 3))  # lower confidence bound
 print(round(nofloodperiodmeananomaly_c + 1.645*nofloodperiodsdanomaly_c/sqrt(nofloodperiodlength_c), 3))  # upper confidence bound
# southern EU
floodperiodmeananomaly_s <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$s, mean)[c('III','Va','VI','IX')]
floodperiodsdanomaly_s <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$s, sd)[c('III','Va','VI','IX')]
floodperiodlength_s <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$s, length)[c('III','Va','VI','IX')]
nofloodperiodmeananomaly_s <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$s, mean)[c('STtoIII','IIItoVa','VatoVI','VItoIX','IXtoEND')]
nofloodperiodsdanomaly_s <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$s, sd)[c('STtoIII','IIItoVa','VatoVI','VItoIX','IXtoEND')]
nofloodperiodlength_s <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$s, length)[c('STtoIII','IIItoVa','VatoVI','VItoIX','IXtoEND')]
 # the following values are used in Extended Data Fig. 4a 
 print(round(floodperiodmeananomaly_s, 3))
 print(round(floodperiodmeananomaly_s - 1.645*floodperiodsdanomaly_s/sqrt(floodperiodlength_s), 3))  # lower confidence bound
 print(round(floodperiodmeananomaly_s + 1.645*floodperiodsdanomaly_s/sqrt(floodperiodlength_s), 3))  # upper confidence bound
 print(round(nofloodperiodmeananomaly_s, 3))
 print(round(nofloodperiodmeananomaly_s - 1.645*nofloodperiodsdanomaly_s/sqrt(nofloodperiodlength_s), 3))  # lower confidence bound
 print(round(nofloodperiodmeananomaly_s + 1.645*nofloodperiodsdanomaly_s/sqrt(nofloodperiodlength_s), 3))  # upper confidence bound
# western EU
floodperiodmeananomaly_w <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$w, mean)[c('I','II','IV','Va','VI','IX')]
floodperiodsdanomaly_w <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$w, sd)[c('I','II','IV','Va','VI','IX')]
floodperiodlength_w <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$w, length)[c('I','II','IV','Va','VI','IX')]
nofloodperiodmeananomaly_w <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$w, mean)[c('ItoII','IItoIV','IVtoVa','VatoVI','VItoIX','IXtoEND')]
nofloodperiodsdanomaly_w <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$w, sd)[c('ItoII','IItoIV','IVtoVa','VatoVI','VItoIX','IXtoEND')]
nofloodperiodlength_w <- tapply(CEtempanom500y_periods$temp_anom, CEtempanom500y_periods$w, length)[c('ItoII','IItoIV','IVtoVa','VatoVI','VItoIX','IXtoEND')]
 # the following values are used in Fig. 4c 
 print(round(floodperiodmeananomaly_w, 3))
 print(round(floodperiodmeananomaly_w - 1.645*floodperiodsdanomaly_w/sqrt(floodperiodlength_w), 3))  # lower confidence bound
 print(round(floodperiodmeananomaly_w + 1.645*floodperiodsdanomaly_w/sqrt(floodperiodlength_w), 3))  # upper confidence bound
 print(round(nofloodperiodmeananomaly_w, 3))
 print(round(nofloodperiodmeananomaly_w - 1.645*nofloodperiodsdanomaly_w/sqrt(nofloodperiodlength_w), 3))  # lower confidence bound
 print(round(nofloodperiodmeananomaly_w + 1.645*nofloodperiodsdanomaly_w/sqrt(nofloodperiodlength_w), 3))  # upper confidence bound


# 3.2) Flood seasonality within and outside of flood periods

# for the seasonality, we use all events
allevents00 <- i_f
allevents00$numIndex <- c(1, 2, 3, NA, NA, NA)[as.numeric(as.factor(allevents00$Index))]
allevents00 <- allevents00[!is.na(allevents00$numIndex),]

dummy <- as.data.frame(sites); rownames(dummy) <- dummy$code
allevents00$regionEU <- dummy[as.character(allevents00$Code), 'regionEU']
 allevents00$period <- NA
for (reg in c('c','s','w','n')) {
 dummy <- max10events[order(max10events$id),][eventsXregion[,reg],]
 for (i in 1:nrow(dummy)) {
  allevents00[(allevents00$regionEU == reg)&(allevents00$Year >= dummy$yr_st[i])&(allevents00$Year <= dummy$yr_en[i]), 'period'] <- dummy$name[i]
 }
}
allevents00$Season[allevents00$Season == ''] <- 'undated'
allevents00$period[is.na(allevents00$period)] <- 'noFrich'

# Seasonality analysis
weights <- c(1,2,3)
# (weighted) proportion of main seasonality in the 3 main regions
seasonarray <- tapply(allevents00$numIndex, 
                      list(allevents00$period, allevents00$Season, allevents00$regionEU, allevents00$numIndex), 
                      length)[c('I','II','III','IV','Va','Vb','VI','VII','VIII','IX','noFrich'),c('spring','summer','autumn','winter'),c('c','s','w'),]
 seasonarray[is.na(seasonarray)] <- 0
 
seas00 <- data.frame(summerCeu=seasonarray[,'summer','c','1'] +
                               seasonarray[,'summer','c','2'] +
                               seasonarray[,'summer','c','3'],
                     autumnSeu=seasonarray[,'autumn','s','1'] +
                               seasonarray[,'autumn','s','2'] +
                               seasonarray[,'autumn','s','3'],
                     winterWeu=seasonarray[,'winter','w','1'] +
                               seasonarray[,'winter','w','2'] +
                               seasonarray[,'winter','w','3'])
seasWW <- data.frame(summerCeu=weights[1]*seasonarray[,'summer','c','1'] +
                               weights[2]*seasonarray[,'summer','c','2'] +
                               weights[3]*seasonarray[,'summer','c','3'],
                     autumnSeu=weights[1]*seasonarray[,'autumn','s','1'] +
                               weights[2]*seasonarray[,'autumn','s','2'] +
                               weights[3]*seasonarray[,'autumn','s','3'],
                     winterWeu=weights[1]*seasonarray[,'winter','w','1'] +
                               weights[2]*seasonarray[,'winter','w','2'] +
                               weights[3]*seasonarray[,'winter','w','3']) * apply(seasonarray, c(1,3), sum) / 
          apply(weights[1]*seasonarray[,,,'1'] + weights[2]*seasonarray[,,,'2'] + weights[3]*seasonarray[,,,'3'], c(1,3), sum)
 seasWW[is.nan(as.matrix(seasWW))] <- 0
 
prop_all_events <- round(seasWW / apply(seasonarray, c(1,3), sum), 2)
 prop_all_events[is.nan(as.matrix(prop_all_events))] <- NA

# values for Figure 5b and Extended Data Fig. 5bd:
seas00_3 <- rbind(seas00['noFrich',], pastFrich=apply(seas00[1:9,], 2, sum), seas00['IX',])
seasWW_3 <- rbind(seasWW['noFrich',], pastFrich=apply(seasWW[1:9,], 2, sum), seasWW['IX',])
n <- rbind(noFrich=apply(seasonarray, c(1,3), sum)['noFrich',], 
           pastFrich=apply(apply(seasonarray, c(1,3), sum)[1:9,], 2, sum), 
           IX=apply(seasonarray, c(1,3), sum)['IX',])
p_hat <- seasWW_3/n
 print(round(p_hat, 3))  # these are the proportions in Figure 5b and Extended Data Fig. 5bd.
# confidence bounds based on clt, see e.g. section 9.4. in introduction to probability and statistics using R (2010) by by GJ Kerns 
std_p_hat <- sqrt(p_hat*(1 - p_hat)/n)  # sqrt(p*(1 - p)/n)
 print(round(p_hat - 1.645*std_p_hat, 3))  # lower confidence bound
 print(round(p_hat + 1.645*std_p_hat, 3))  # upper confidence bound

