Repository of code and data for "The seasonality of diarrheal pathogens: A retrospective study of seven sites over three years" (https://doi.org/10.1101/541581)

The "Data" directory has weather data and estimated case count data.
 * The seven files "gsod-*.csv" are daily weather data for the stations closest to the study sites. They were downloaded from: https://data.noaa.gov/dataset/global-surface-summary-of-the-day-gsod.
 * "weather-daily.csv" is the processed weather data for all seven stations, in which relative humidity estimates were added.
 * "gems-seasons.csv" contains the season labels for the seven sites for the time period of the study.
 * "estimated-positive-fortnightly.csv" has fortnightly counts of the number of children with moderate-to-severe diarrhea (MSD) seen at each site and estimates of the number of these children who were positive for different pathogens.
 * "estimated-positive-monthly.csv" has monthly counts of the number of children with MSD and positive for pathogens. Separate counts for children who were 0 to 23 months and 24 to 59 months are given.

The directory "Code" contains the MatLab script "AnalysisOfPeriodicityByDiseaseByCountry" required to detect significant periodicities in the study data and to generate Figure 3.  The Matlab script "PlottingPeriodogramForRotavirusFig1" generates the Lomb-Scargle periodogram for rotavirus seen in Figure 1 of the paper.  

The "Code" directory contains R scripts for processing weather data and identifying seasonality in "AnalysisOfWeatherDataAndIdentificationOfSeasons" and to detect associations between weather and pathogen prevalence in "AnalysisAssociationOfWeatherAndSeasonsToPrevalence".
  * AnalysisOfWeatherDataAndIdentificationOfSeasons.R contains code for generating "weather-daily.csv" from the raw "gsod-*.csv" files, generating a version of "gems-seasons.csv" from the weather data, and summarizing the weather conditions for each season.
  * AnalysisAssociationOfWeatherAndSeasonsToPrevalence.R contains code for analyzing the associations between weather or season and diarrheal case counts or pathogen prevalence.

The "Extra Figures" directory contains additional Figures not contained in the supporting material of the manuscript.
