# analyses.R
# analyses for:
#   The seasonality of diarrheal pathogens: A retrospective study of seven sites over three years. bioRxiv. 2019.
#   doi: 10.1101/541581
#
# Purpose:
#   1. Summarize weather data for different seasons.
#   2. Compare fortnights with high vs low numbers of MSD cases
#   3. Compare numbers of cases in different seasons
#
# Input data files should be in ../Data.

################
# Read data

# read GEMS pathogen case count data
data.pathogens.fortnightly <- read.csv(file="../Data/estimated-positive-fortnightly.csv", header=TRUE)
data.pathogens.monthly <- read.csv(file="../Data/estimated-positive-monthly.csv", header=TRUE)

# read weather data
data.gsod.daily <- read.csv(file="../Data/weather-daily.csv", header=TRUE)
data.gsod.daily$DATE <- as.Date(data.gsod.daily$date, format="%m-%d-%Y")

# read season data
seasons.with.labels <- read.csv("../Data/gems-seasons.csv", header=TRUE, row.names=1)
temp <- rownames(seasons.with.labels)
seasons.with.labels <- data.frame(lapply(seasons.with.labels, as.character), stringsAsFactors=FALSE)
rownames(seasons.with.labels) <- temp

################
# summarize weather conditions for each season

# print mean and interquartile range of monthly weather conditions for each country-season pair
tempmonth <- as.numeric(format(data.gsod.daily$DATE, "%m"))
tempyear <- as.numeric(format(data.gsod.daily$DATE, "%Y"))
for (sitenum in 1:nrow(seasons.with.labels)) {
    countryname <- rownames(seasons.with.labels)[sitenum]
    for (season in unique(as.character(seasons.with.labels[sitenum,]))) {
        if (!is.na(season)) {
            temp <- (seasons.with.labels[sitenum,]==season)
            temp[is.na(temp)] <- FALSE
            m <- rep(1:12,5)[temp]
            y <- rep(2007:2011,each=12)[temp]
            temperatures <- sapply(1:length(m), function(i) {mean(data.gsod.daily$temperature..F.[data.gsod.daily$country==countryname & tempmonth==m[i] & tempyear==y[i]], na.rm=TRUE)})
            rainfall <- sapply(1:length(m), function(i) {sum(data.gsod.daily$rain..in.[data.gsod.daily$country==countryname & tempmonth==m[i] & tempyear==y[i]], na.rm=TRUE)})
            rh <- sapply(1:length(m), function(i) {mean(data.gsod.daily$relative.humidity....[data.gsod.daily$country==countryname & tempmonth==m[i] & tempyear==y[i]], na.rm=TRUE)})
            sh <- sapply(1:length(m), function(i) {mean(data.gsod.daily$specific.humidity..g.kg.[data.gsod.daily$country==countryname & tempmonth==m[i] & tempyear==y[i]], na.rm=TRUE)})
            cat(countryname, season, "rainfall (in)", round(mean(rainfall),2), round(quantile(rainfall,c(0.25,0.75)),2), fill=TRUE)
            cat(countryname, season, "temperature (F)", round(mean(temperatures),1), round(quantile(temperatures,c(0.25,0.75)),1), fill=TRUE)
            cat(countryname, season, "relative humidity (%)", round(mean(rh),1), round(quantile(rh,c(0.25,0.75)),1), fill=TRUE)
            cat(countryname, season, "specific humidity (g/kg)", round(mean(sh),1), round(quantile(sh,c(0.25,0.75)),1), fill=TRUE)
        }
    }
}

################
# fortnightly weather vs pathogen prevalence

# add fortnightly weather conditions to data.pathogens.fortnightly data frame
fortnightdates <- as.Date(data.pathogens.fortnightly$date,format="%m-%d-%Y")
gsoddates <- as.Date(data.gsod.daily$date,format="%m-%d-%Y")
data.pathogens.fortnightly$rain.in <- NA
data.pathogens.fortnightly$temp.f <- NA
data.pathogens.fortnightly$relative.humidity <- NA
data.pathogens.fortnightly$specific.humidity <- NA
for (i in 1:nrow(data.pathogens.fortnightly)) {
    data.pathogens.fortnightly$rain.in[i] <- sum(data.gsod.daily$rain..in.[data.gsod.daily$country==data.pathogens.fortnightly$country[i] & (as.numeric(gsoddates-fortnightdates[i]) %in% 0:13)], na.rm=TRUE)
    data.pathogens.fortnightly$temp.f[i] <- mean(data.gsod.daily$temperature..F.[data.gsod.daily$country==data.pathogens.fortnightly$country[i] & (as.numeric(gsoddates-fortnightdates[i]) %in% 0:13)], na.rm=TRUE)
    data.pathogens.fortnightly$relative.humidity[i] <- mean(data.gsod.daily$relative.humidity....[data.gsod.daily$country==data.pathogens.fortnightly$country[i] & (as.numeric(gsoddates-fortnightdates[i]) %in% 0:13)], na.rm=TRUE)
    data.pathogens.fortnightly$specific.humidity[i] <- mean(data.gsod.daily$specific.humidity..g.kg.[data.gsod.daily$country==data.pathogens.fortnightly$country[i] & (as.numeric(gsoddates-fortnightdates[i]) %in% 0:13)], na.rm=TRUE)
}

# compare weather between fortnights with high vs low numbers of cases
for (countryname in levels(data.pathogens.fortnightly$country)) {
    for (colnum in grep("eligible", names(data.pathogens.fortnightly))) {
        pos <- data.pathogens.fortnightly[,colnum] # fortnightly case data for all sites
        v <- data.pathogens.fortnightly$country==countryname # rows with data for one site
        if (sum(pos[v])>20) { # must have at least 20 cases
            thresh <- quantile(pos[v],c(0.25,0.75)) # top and bottom quartiles of cases
            if (wilcox.test(data.pathogens.fortnightly$rain.in[v & pos <=thresh[1]], data.pathogens.fortnightly$rain.in[v & pos>thresh[2]])$p.value<0.05) {
                hi <- (mean(data.pathogens.fortnightly$rain.in[v & pos <=thresh[1]],na.rm=TRUE)<mean(data.pathogens.fortnightly$rain.in[v & pos>thresh[2]],na.rm=TRUE))
                cat(countryname, names(data.pathogens.fortnightly)[colnum], ifelse(hi, "rainy", "not rainy"), "p<0.05", fill=TRUE)
            }
            if (wilcox.test(data.pathogens.fortnightly$temp.f[v & pos <=thresh[1]], data.pathogens.fortnightly$temp.f[v & pos>thresh[2]])$p.value<0.05) {
                hi <- (mean(data.pathogens.fortnightly$temp.f[v & pos <=thresh[1]],na.rm=TRUE)<mean(data.pathogens.fortnightly$temp.f[v & pos>thresh[2]],na.rm=TRUE))
                cat(countryname, names(data.pathogens.fortnightly)[colnum], ifelse(hi, "hot", "cold"), "p<0.05", fill=TRUE)
            }
            if (wilcox.test(data.pathogens.fortnightly$relative.humidity[v & pos <=thresh[1]], data.pathogens.fortnightly$relative.humidity[v & pos>thresh[2]])$p.value<0.05) {
                hi <- (mean(data.pathogens.fortnightly$relative.humidity[v & pos <=thresh[1]],na.rm=TRUE)<mean(data.pathogens.fortnightly$relative.humidity[v & pos>thresh[2]],na.rm=TRUE))
                cat(countryname, names(data.pathogens.fortnightly)[colnum], ifelse(hi,"RH humid","RH dry"), "p<0.05", fill=TRUE)
            }
            if (wilcox.test(data.pathogens.fortnightly$specific.humidity[v & pos <=thresh[1]], data.pathogens.fortnightly$specific.humidity[v & pos>thresh[2]])$p.value<0.05) {
                hi <- (mean(data.pathogens.fortnightly$specific.humidity[v & pos <=thresh[1]],na.rm=TRUE)<mean(data.pathogens.fortnightly$specific.humidity[v & pos>thresh[2]],na.rm=TRUE))
                cat(countryname, names(data.pathogens.fortnightly)[colnum], ifelse(hi,"SH humid","SH dry"), "p<0.05", fill=TRUE)
            }
        } else if (sum(pos[v])==0) { # no cases
            cat(countryname, names(data.pathogens.fortnightly)[colnum], "NOT ENOUGH CASES", fill=TRUE)
        }
    }
}


################
# seasons vs pathogen prevalence

library(dunn.test) # Dunn's test

season.month <- rep(1:12,5)
season.year <- rep(2007:2011,each=12)
season.month.year <- paste(rep(1:12,5),rep(2007:2011,each=12),sep="-")
pathogens.month.year <- paste(data.pathogens.monthly$month,data.pathogens.monthly$year,sep="-")

for (countryname in levels(data.pathogens.monthly$country)) {
    seasonindex <- match(countryname,rownames(seasons.with.labels)) # which row of season data corresponds to this countryname?
    for (pathogenname in levels(data.pathogens.monthly$pathogen)) {
        positives <- list()
        for (seasonnum in 1:3) {
            seasonname <- unique(unlist(seasons.with.labels[seasonindex,]))[seasonnum+1]
            v <- (data.pathogens.monthly$country==countryname & data.pathogens.monthly$pathogen==pathogenname)
            positives[[seasonname]] <- (data.pathogens.monthly$positive.ages.0.23 + data.pathogens.monthly$positive.ages.24.59)[(data.pathogens.monthly$country==countryname & data.pathogens.monthly$pathogen==pathogenname) & pathogens.month.year %in% season.month.year[seasons.with.labels[seasonindex,]==seasonname]]
        }
        results.dunn <- dunn.test(positives, kw=FALSE, table=FALSE,
                                  method = "Bonferroni") # we do a Bonferroni correction

        cat(countryname, pathogenname, "p values", results.dunn$P.adjusted[1], results.dunn$P.adjusted[2], results.dunn$P.adjusted[3], fill=TRUE)
        if (results.dunn$P.adjusted[1]<0.05) {
            cat(countryname, pathogenname, unique(unlist(seasons.with.labels[seasonindex,]))[c(1,2)+1], "p<0.05", fill=TRUE)
        }
        if (results.dunn$P.adjusted[2]<0.05) {
            cat(countryname, pathogenname, unique(unlist(seasons.with.labels[seasonindex,]))[c(1,3)+1], "p<0.05", fill=TRUE)
        }
        if (results.dunn$P.adjusted[3]<0.05) {
            cat(countryname, pathogenname, unique(unlist(seasons.with.labels[seasonindex,]))[c(2,3)+1], "p<0.05", fill=TRUE)
        }
    }
}
