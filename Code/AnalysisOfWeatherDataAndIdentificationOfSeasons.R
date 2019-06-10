# weatherprocessing.R
# generate weather and season data files for:
#   The seasonality of diarrheal pathogens: A retrospective study of seven sites over three years. bioRxiv. 2019.
#   doi: 10.1101/541581
#
# Purpose:
#   1. Convert raw data downloaded from the Global Summary of the Day (GSOD) to temperature, rainfall, relative humidity, and specific humidity
#   2. Divide the year into three different kinds of "season" for each site
#      Note: The seasons are not labeled
#
# Input data should be in ../Data, and output files go to the current directory, which can be moved to ../Data if desired.

#####################
# Read GSOD station data
sites <- data.frame(countryname=c("The Gambia","Mali","Mozambique","Kenya","India","Bangladesh","Pakistan"),
                    startdate=c("2007-12-10","2007-12-17","2007-12-10","2008-01-31","2007-12-01","2007-12-02","2008-03-03"), # date of first clinic visit in the study
                    noaa.stationid=c(616870,612910,673410,637080,428090,419230,417800), # USAF/ID of noaa stations closest to the study sites
                    noaa.elevation=c(49.1,380.1,44.2,1208.5,4.9,7.3,30.5)) # weather station elevation in meters, from isd-history.csv
sites$startdate <- as.Date(sites$startdate, format="%Y-%m-%d")

# read data from weather stations (downloaded from gsod)
files <- list.files(path="../Data", pattern="^gsod-", full.names=T, recursive=FALSE)
ids <- NULL
for (x in files) {
    temp <- read.csv(x, header=T) # load file
    ids <- append(ids,temp$STN...[1])
}
files.gsod <- data.frame(filename=files, id=ids)
                        
# weather conversion functions
relativetospecifichumidity <- function(elevation, temperature, rh) { # temperatures should be in celsius, elevation in m. Function from James Tamerius
    P <- 101.3*exp(-elevation/8200) # approximate air pressure from elevation
    Es <- 0.622 * exp((17.502*temperature)/(temperature+240.97)) # saturation vapor pressure (kPa)
    E <- (0.01*rh)*Es # vapor pressure (kPA)
    W <- (0.622*E)/(P-E) # mixing ratio
    1000 * W/(1+W) # specific humidity (g/kg)
}# Low humidity is defined as less than 11–12 g/kg

dewpointtorelativehumidity <- function(dewpoint,temperature) { # temperatures should be in celsius
    beta <- (112 - (0.1 * temperature) + dewpoint)/(112 + (0.9 * temperature)) # from Bosen
    relative.humidity <- 100 * beta^8
}

vaporpressuretorelativehumidity <- function(vp, temperature) { # temperature in celsius, vp in hectopascals
    actualvaporpressure <- vp
    Es <- 0.622 * exp((17.502*temperature)/(temperature+240.97)) # saturation vapor pressure (kPa)
    saturationvaporpressure <- Es*10 # hectopascals
    relativehumidity <- 100*actualvaporpressure/saturationvaporpressure # in %
}

# chopflag - removes flags from numeric data from GSOD
chopflag <- function(s) { l <- nchar(s)
                          ch <- substr(s,l,l)
                          if ((ch>='A' && ch<='I') || ch=='*') {
                              as.numeric(substr(s,1,l-1))
                          } else {
                            as.numeric(s)
                          }}

# write daily weather data from GSOD sites to file
cat("country", "station id", "date", "rain (in)", "temperature (F)", "relative humidity (%)", "specific humidity (g/kg)", sep=",", file="weather-daily.csv", append=FALSE, fill=1000)
for (sitenum in 1:7) {
    data.gsod <- read.table(as.character(files.gsod$filename[match(sites$noaa.stationid[sitenum],files.gsod$id)]), sep=",", strip.white=TRUE, header=TRUE)
    gsod.date <- as.Date(as.character(data.gsod$YEARMODA), format="%Y%m%d")
    gsod.dayssincestart <- gsod.date-sites$startdate[sitenum]
    data.gsod$dewp <- data.gsod$DEWP
    data.gsod$dewp[data.gsod$dewp==9999.9] <- NA
    data.gsod$temp <- data.gsod$TEMP
    data.gsod$temp[data.gsod$temp==9999.9] <- NA
    data.gsod$prcp <- sapply(levels(data.gsod$PRCP)[data.gsod$PRCP], chopflag) # ignore the QC flags
    data.gsod$prcp[data.gsod$prcp==99.99] <- NA
    data.gsod$rh <- dewpointtorelativehumidity(dewpoint=(data.gsod$dewp-32)*5/9, temperature=(data.gsod$temp-32)*5/9)
    elevation <- sites$noaa.elevation[sitenum]
    data.gsod$sh <- relativetospecifichumidity(elevation=rep(elevation,nrow(data.gsod)), temperature=(data.gsod$temp-32)*5/9, rh=data.gsod$rh)

    for (day in -60:(365*3+60)) {
        cat(as.character(sites$countryname[sitenum]),
            as.character(sites$noaa.stationid[sitenum]),
            format(sites$startdate[sitenum]+day,"%m-%d-%Y"),
            data.gsod$prcp[gsod.dayssincestart==day],
            data.gsod$temp[gsod.dayssincestart==day],
            data.gsod$rh[gsod.dayssincestart==day],
            data.gsod$sh[gsod.dayssincestart==day],
            sep=",", file="weather-daily.csv", append=TRUE, fill=TRUE)
    }
}


##############
# identifying seasons (based on monthly gsod weather)
# read weather data
data.gsod.daily <- read.csv(file="../Data/weather-daily.csv", header=TRUE)
data.gsod.daily$DATE <- as.Date(data.gsod.daily$date, format="%m-%d-%Y")
gsod.months <- as.numeric(format(data.gsod.daily$DATE, "%m"))
gsod.years <- as.numeric(format(data.gsod.daily$DATE, "%Y"))

# generate matrices with weather data for each site
weather.monthly <- NULL
weather.monthly.mean0 <- NULL # 0 mean variance 1 version of weather data
for (countrynum in 1:7) {
    country <- levels(data.gsod.daily$country)[countrynum]
    temp <- matrix(data=NA, nrow=12*5, ncol=4)
    colnames(temp) <- c("rain","temp", "RH", "SH")
    temp[,1] <- sapply(1:(5*12), function(m) {sum(data.gsod.daily$rain..in.[data.gsod.daily$country==country & gsod.years==trunc((m-1)/12)+2007 & gsod.months==((m-1)%%12)+1], na.rm=TRUE)})# rain
    temp[,2] <- sapply(1:(5*12), function(m) {mean(data.gsod.daily$temperature..F.[data.gsod.daily$country==country & gsod.years==trunc((m-1)/12)+2007 & gsod.months==((m-1)%%12)+1], na.rm=TRUE)})# temp
    temp[,3] <- sapply(1:(5*12), function(m) {mean(data.gsod.daily$relative.humidity....[data.gsod.daily$country==country & gsod.years==trunc((m-1)/12)+2007 & gsod.months==((m-1)%%12)+1], na.rm=TRUE)})# RH
    temp[,4] <- sapply(1:(5*12), function(m) {mean(data.gsod.daily$specific.humidity..g.kg.[data.gsod.daily$country==country & gsod.years==trunc((m-1)/12)+2007 & gsod.months==((m-1)%%12)+1], na.rm=TRUE)})# SH
    weather.monthly[[countrynum]] <- temp
    temp[,1] <- (temp[,1]-mean(temp[1,],na.rm=TRUE))/sqrt(var(temp[,1]))
    temp[,2] <- (temp[,2]-mean(temp[2,],na.rm=TRUE))/sqrt(var(temp[,2],na.rm=TRUE))
    temp[,3] <- (temp[,3]-mean(temp[3,],na.rm=TRUE))/sqrt(var(temp[,3],na.rm=TRUE))
    temp[,4] <- (temp[,4]-mean(temp[4,],na.rm=TRUE))/sqrt(var(temp[,4],na.rm=TRUE))
    weather.monthly.mean0[[countrynum]] <- temp
}

postscript("weather-pca-propvariance.ps",
           width=8,height=4,horizontal=FALSE)
par(mfrow=c(2,4))
for (countrynum in 1:7) { 
    countryname <- levels(data.gsod.daily$country)[countrynum]
    numclusters <- 3
    v <- !is.na(rowSums(weather.monthly.mean0[[countrynum]]))
    result.pca <- prcomp(weather.monthly.mean0[[countrynum]][v,]) # PCA
    variances <- (result.pca$sdev)^2
    par(mar=c(4,3,3,1)+0.1)
    par(mgp=c(2,0.7,0))
    plot(cumsum(variances)/sum(variances), main=countryname, ylim=c(0.6,1.0), ylab="cumulative prop var", xlab="principal components", type="b") # variances
    print(cumsum(variances)/sum(variances))
}
dev.off()

# identify 3 seasons per site and write to file
# Note that kmeans is stochastic, so results can change
# Also note that the clusters/seasons are labeled as 1, 2, and 3, and human-readable labels need to be added by the user
cat(c("",paste(rep(1:12,length(temp)/12),rep(2007:2011,each=12),sep="-")), file="gems-seasons.csv", sep=",", fill=2000, append=FALSE)
for (countrynum in 1:7) { 
    countryname <- levels(data.gsod.daily$country)[countrynum]
    numclusters <- 3
    v <- !is.na(rowSums(weather.monthly.mean0[[countrynum]]))
    result.pca <- prcomp(weather.monthly.mean0[[countrynum]][v,]) # PCA
    result.kmeans <- kmeans(result.pca$x[,1:3], numclusters) # only use first three PCA components
    # average weather of and months covered by seasons
    for (i in sort(unique(result.kmeans$cluster))) {
        cat(countryname, " season ",i," avg monthly rainfall: ", round(mean(weather.monthly[[countrynum]][!is.na(temp) & temp==i,1], na.rm=TRUE),2), sep="", fill=TRUE)
        cat(countryname, " season ",i," avg temperature: ", round(mean(weather.monthly[[countrynum]][!is.na(temp) & temp==i,2],na.rm=TRUE),1), sep="", fill=TRUE)
        cat(countryname, " season ",i," avg relative humidity: ", round(mean(weather.monthly[[countrynum]][!is.na(temp) & temp==i,3],na.rm=TRUE),1), sep="", fill=TRUE)
        cat(countryname, " season ",i," avg specific humidity: ", round(mean(weather.monthly[[countrynum]][!is.na(temp) & temp==i,4]),1), sep="", fill=TRUE)
        cat("months in season",i,":", sort(unique(rep(1:12,5)[!is.na(temp) & temp==i])), fill=TRUE)
    }                                        
    # write newly identified seasons to file
    cat(c(levels(data.gsod.daily$country)[countrynum],temp), file="gems-seasons.csv", sep=",", fill=2000, append=TRUE)
}
