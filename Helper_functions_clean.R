##Libraries=====================================================================
require(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(overlap)
library(maptools)
library(sp)
library(scales)
library(anytime)  
library(ggpubr)
library(survival) 
library(survminer)
library(sjPlot)


## List of functions ===========================================================

# controlInterval: This function sets a minimum time interval between successive
# camera trap images. Used to ensure the images are temporally independent. 
 
# overlapPlotCI: used to create customized overlap plots for the diel activity 
# of two species. This function is built using functions from the Overlap package.  

# overlapCI: Used to obtain the confidence interval for the overlap coefficient
# of the activity of two species. Built using functions form the overlap package. 

# overlapPlotST: This function creates overlap plots similar to those of 
# the OverlapPlotCI function, but uses the standard plot desgin of the overlap 
# package. 

# overapPlotGrid: This function creates a grid of overlap plots built using the 
# OverlapPlotCI fucntion 
 
# getTimes: For a given pair of species (sp1 and sp2) and a .csv file exported
# from Timelapse, this function creates a new dataframe with three new columns.
# This dataframe only includes rows for images of sp1. The Time_diff column 
# shows the time between the observation of sp1 and the next observation of sp1 
# or sp2 by the camera.The Next_sp column shows whether sp1 and sp2 was detected 
# after sp1. The Censor column shows whether the time interval between the two 
# observations is censored or not (1 indicates not censored, 0 indicates 
# censored). All rows where sp1 was detected after sp1 are censored. 
 
# mergeLayers: This function will merge the information from a dataframe with 
# spatial data and a dataframe obtained with getTimes or imported from 
# Timelapse. The name of the camera trapping site (ine the RelativePath column)
# is used as the key for the merge. 
 
# coxModel: function to run cox proportional hazard models using character 
# formulas and store the results and related plots in a list. 
 
# addFreqSp: This function is used to add columns with the weekly and monthly
# frequency of detentions for a species into a dataframe exported from Timelapse 
 
# addFrqHum: This function is used to add columns with the weekly and monthly
# frequency of detentions for humans into a dataframe exported from Timelapse 

# sample data: 
sdata1 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Research/Professor_Molnar/Canid_temporal_segregation_manuscript/Datasets/TimelapseDataTUWall_Nsp_over0_V7.csv")
sdata2 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Research/Professor_Molnar/Canid_temporal_segregation_manuscript/Datasets/AllCameras_SpatialData_500mbuffer.csv")

## Function to remove images within a time interval ============================

controlInterval <- function(data,lim, keep = c()){
  
## Notes =======================================================================  

  # data: dataframe obtained by importing .csv file generated with Timelapse
  # lim: length of the time interval (in minutes) following a single camera detection during which no other pictures will be retained. 
  # keep: add the names of species of interest (sp 1 and 2 in OverlapPlotIi and overlapCI) to exclude all other detection

  
  # Times will be added to the date column
  # Times will be converted to the POSIXct format 
  # Use the keep argument to exclude empty images ort those that have only humans 
  
## Test call ===================================================================

    # x <- controlInterval(sdata1, 10)
  # x <- controlInterval(sdata1, 15, keep = c("fox","dog", "coyote"))

## Function body ===============================================================  

  ##format date, convert to POSIXct if it is not already
  as.character(data$Species) -> data$Species
  if (class(data$DateTime)[1] == "character"){ 
    data$DateTime <- anytime(data$DateTime)
  }
  
  #remove unwanted pictures 
  if (length(keep) != 0){
    data <- filter(data, (data$Species %in% keep))
  }
  
  ##sort following camera and time
  data <- data %>% dplyr::arrange(RelativePath, DateTime)
  
  # calculate the time differences between images 
  data <- data %>% mutate(Time_diff = DateTime - lag(DateTime))
  
  #find the cameras of each preceding row
  data <- data %>% mutate(Prevcam = lag(RelativePath))
  
  ## convert the time difference from difftime to hours
  data$Time_diff <- as.numeric(data$Time_diff, units = "mins")
  
  ## set Time_diff to limit for the first observation of each new camera so it is not removed later
  data$Time_diff[data$RelativePath != data$Prevcam] <- lim
  
  # drop rows within time limit
  output_data <- data[!(data$Time_diff < lim),]
  
  # drop rows where the previous camera is from a different site 
  output_data <- output_data[(output_data$RelativePath == output_data$Prevcam),]
  
  # dropthe first row, which is not followed by any camera 
  output_data <- output_data[2:nrow(output_data),]
  
  #drop columns used in the function 
  drop <- c("Time_diff","Prevcam")
  output_data <- output_data[,!(names(output_data) %in% drop)]
  
  ##renumber rows
  rownames(output_data) <- NULL
  
  return(output_data)
}

## Function for overlap plots with confidence intervals ========================

overlapPlotCI <- function(sp1, sp2, data, grid = 100, nb = 1000, cam = c(),
                          type = 1, rug = TRUE, stat = TRUE, title= TRUE, legend = TRUE) {

## Notes =======================================================================  
  
  # sp1: first species in the pair
  # sp2: second species in the pair
  # data: dataframe obtained by importing .csv file generated with Timelapse
  # cam: vector with camera names that will be used for the plot. if empty, will use all cameras
  # grid: number of points for plotting the Von Mises density 
  # nb: number of bootstrap samples 
  # plot_type: can be type 1 or 2
  # rug: includes a rug on the plot when TRUE
  # stat: adds the value of the overlap to the plot when TRUE
  # title: adds a title to the plot when TRUE 
  # legend: adds a legend to the plot when TRUE

  # This function uses sun time. The location is set to a single point in 
  # Toronto, Ontario, Canada by default
  # The function will not work for cameras that captured less than 2 individuals of a given species
  
## Test call====================================================================
  
    # overlapPlotCI(sp1 = "fox", sp2 = "coyote", data = sdata1, grid = 100, nb = 1000, type = 1)
  # overlapPlotCI(sp1 = "fox", sp2 = "coyote", data = sdata1, grid = 100, nb = 1000, type = 2)
  # overlapPlotCI(sp1 = "dog", sp2 = "fox", data = sdata1, grid = 100, nb = 1000, type = 1, cam = c("TUW28"))
  
## Function body ===============================================================  
  
  locale = locale(date_format = "%Y-%m-%d",
                  time_format = "%H:%M:%S",
                  tz = "UTC")
  #filter cameras
  if (length(cam) != 0){
    data <- filter(data, data$RelativePath %in% cam)
  }
  
  ##format date, convert to POSIXct if it is not already
  as.character(data$Species) -> data$Species
  if (class(data$DateTime)[1] == "character"){ 
    data$DateTime <- anytime(data$DateTime)
  }
  
  
  ##Convert time to value between 0 and 2pi (radian)
  data$Time <- ((hour(data$DateTime)*60 + minute(data$DateTime))/1440) * 2*pi
  
  ##dates to a POSIXct object with the right time zone (GMT):
  dates_toronto <- data$DateTime
  
  ##create spatial point for the study area (latitude and longitude of Toronto, assume one point is accurate enoguh for whole study area)
  coords <- matrix(c(-79.3831843, 43.653226), nrow = 1)
  point <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
  
  ##Assign sun time to data$Time (replace clock time by sun time)  
  data$Time <- sunTime(data$Time, dates_toronto, point)
  
  ##create vectors of observation times for each species 
  sp1_ob <- data$Time[data$Species  == sp1] #Species A
  sp2_ob <- data$Time[data$Species  == sp2] #Species B
  
  ##Check if the number of detections is higher than 2 
  if (length(sp1_ob) < 2 || length(sp2_ob) < 2){
    return("Number of dections for sp1 or sp2 is less than 2")
  }
  
  ##create empty matrix to store bootstrapped values
  mat1 <- matrix(data = NA, nrow = length(sp1_ob), ncol = nb)
  mat2 <- matrix(data= NA, nrow = length(sp2_ob), ncol = nb)
  
  ##sample sun times with replacements and add the values to the matrix  
  for (i in 1:nb){ 
    mat1[1:length(sp1_ob), i] <- sample(sp1_ob, replace = TRUE)
    mat2[1:length(sp2_ob), i] <- sample(sp2_ob, replace = TRUE)
  }
  
  ##fit sun times to Von Mises distribution
  #create empty matrices for density 
  mat_densityA <- matrix(data = NA, nrow = grid, ncol = nb)
  mat_densityB <- matrix(data = NA, nrow = grid, ncol = nb)
  
  #fit density for all sun times in the matrices
  for (i in 1:nb){
    
    bwA <- getBandWidth(sp1_ob, kmax = 3)
    bwB <- getBandWidth(sp2_ob, kmax = 3)
    
    mat_densityA[1:grid, i] <- densityFit(mat1[1:nrow(mat1),i], seq(0,2*pi,length= grid), bwA) /(100/(2*pi))
    mat_densityB[1:grid, i] <- densityFit(mat2[1:nrow(mat2),i], seq(0,2*pi,length= grid), bwB) / (100/(2*pi))
    
  }
  
  ##matrix with 4 columns to contain the upper and lower bounds for the 95% confidence interval 
  confA =  matrix(data = NA, nrow = grid, ncol = 3)
  confB = matrix(data = NA, nrow = grid, ncol = 3)
  
  ##Get the mean and boundaries of the confidence interval for the density distributions 
  for (i in 1:grid){
    
    ##method 3 (based on input from Phillip) 
    sorted_pointsA <- sort(mat_densityA[i, 1:ncol(mat_densityA)])
    
    confA[i, 3] = mean(mat_densityA[i, 1:ncol(mat_densityA)]) #mean 
    confA[i, 1] = sorted_pointsA[round(0.025 * length(sorted_pointsA))] #2.5% quantile 
    confA[i, 2] = sorted_pointsA[round(0.975 * length(sorted_pointsA))] #97.5% quantile
    
    
    sorted_pointsB <- sort(mat_densityB[i, 1:ncol(mat_densityB)])
    
    confB[i, 3] = mean(mat_densityB[i, 1:ncol(mat_densityB)]) #mean 
    confB[i, 1] = sorted_pointsB[round(0.025 * length(sorted_pointsB))] #2.5% quantile 
    confB[i, 2] = sorted_pointsB[round(0.975 * length(sorted_pointsB))] #97.5% quantile
    
  }
  
  ##dataframe for the plots 
  df <- data.frame(estimateA = confA[1:nrow(confA),3], confA_up = confA[1:nrow(confA),1],  confA_down = confA[1:nrow(confA),2],
                   estimateB = confB[1:nrow(confA),3], confB_up = confB[1:nrow(confB),1], confB_down = confB[1:nrow(confB),2])
  
  if (rug == TRUE){
    ##dataframes for rugs
    df2 <- data[(data$Species  == sp1),] # rug for sp1
    df2$Time <- df2$Time*100/(2 *pi)
    
    df3 <- data[(data$Species  == sp2),]# rug for sp2 
    df3$Time <- df3$Time*100/(2 *pi)
  }
  
  if (stat == TRUE){
    ## calculate overlap value
    if ((length(sp1_ob) > 50) | (length(sp2_ob) > 50)){
      overlap <- overlapEst(sp1_ob, sp2_ob, type =  "Dhat4")  
    } else {
      overlap <- overlapEst(sp1_ob, sp2_ob, type =  "Dhat1")
    }
    
  }
  
  #legend names
  legend_sp1 = paste(sp1, "density")
  legend_sp2 = paste(sp2, "density")
  
  ## single lines plot with overlap (type 1)
  if (type == 1){
    
  plot = ggplot(data = df, mapping = aes(x = seq(0,100,length = grid))) +
      geom_area(aes(y = pmin(confA_down, confB_down)),
                alpha = 0.2) +
      geom_area(aes(y = pmin(estimateA, estimateB)),
                alpha = 0.2) +
      geom_area(aes(y = pmin(confA_up, confB_up)),
                alpha = 0.2) +
      geom_line(aes(y = estimateA, colour = "blue")) +
      geom_line(aes(y = confA_up), colour = "blue", linetype = "dashed") +
      geom_line(aes(y = confA_down), colour = "blue", linetype = "dashed") +
      geom_line(aes(y = estimateB, colour = "red")) +
      geom_line(aes(y = confB_up), colour = "red", linetype = "dashed") +
      geom_line(aes(y = confB_down), colour = "red", linetype = "dashed") +
      {if(rug == TRUE) geom_rug(data = df2, aes(x = Time), sides = "t", colour = "blue")} +
      {if(rug == TRUE) geom_rug(data = df3, aes(x = Time), sides = "b", colour = "red")} +
      scale_color_manual("", labels = c(legend_sp1, legend_sp2), values = c("blue", "red")) +
      {if (legend == FALSE) guides(fill = "none", color = "none", linetype = "none", shape = "none")} +
      scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = c("Midnight", 
                                    "Sunrise", "Noon", "Sunset", "Midnight"))+
      labs(x = "Time of the day",
           y = "Density")+
    {if(title == TRUE) labs(title = paste("Kernel density for", sp1, "and", sp2, "activity, with confidence interval"))} +
    {if(stat == TRUE) labs(subtitle = paste("Overlap estimate = ", round(overlap, digits = 2), 
                                            ", N", sp1, " = ", length(sp1_ob),", N", sp2," = ", length(sp2_ob), sep = ""))} +
    {if(stat == TRUE & length(cam) == 1) labs(subtitle = paste(cam[1],", overlap estimate = ", round(overlap, digits = 2), 
                                            ", N", sp1, " = ", length(sp1_ob),", N", sp2," = ", length(sp2_ob), sep = ""))} +
      theme_bw() +
      theme(text = element_text(size= 11)) #change all font sizes
    
  }

  ## multiple line plot (type 2)
  if (type == 2){
    
  df_A <- gather(as.data.frame(mat_densityA))
  df_A$time <- rep(seq(0,100,length = grid), 100)
  
  df_B <- gather(as.data.frame(mat_densityB))
  df_B$time <- rep(seq(0,100,length = grid), 100)
  
  plot = ggplot() + 
    geom_area(data = df, aes(y = pmin(estimateA, estimateB),
                            x = rep(seq(0,100,length = grid))), alpha = 0.3) +
    geom_line(data = df_A, aes(x = time, y = value, group = key, colour = legend_sp1))+
    geom_line(data = df_B, aes(x = time, y = value, group = key, colour = legend_sp2))+
    scale_color_manual("", labels = c(legend_sp2, legend_sp1), values = alpha(c("red", "blue"), 0.03))+
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = c("Midnight", 
                                                                    "Sunrise", "Noon", "Sunset", "Midnight"))+
    labs(x = "Time of the day",
         y = "Density",
         title = paste("Kernel density for", sp1, "and", sp2, "activity, with confidence interval"))+
    theme_bw()+
    theme(text = element_text(size= 11)) #change all font sizes
  }

  return(plot)
  
}

## Function to calculate the confidence interval limits for the overlap ========

overlapCI <- function(sp1, sp2, data, nb = 1000, ci= 0.95, cam = c(), samp = FALSE) {
  
## Notes =======================================================================  
  # sp1: first species in the pair
  # sp2: second species in the pair
  # data: dataframe obtained by importing .csv file generated with Timelapse
  # nb: number of boostrap samples 
  # ci: size of the confidence interval
  # cam: vector with camera names that will be used for the plot. If empty, will use all cameras.
  # samp: whether or not to return the estimated overlap for the samples. 
  
## Test call ===================================================================
  
   # overlapCI("fox", "coyote", sdata1, nb = 1000, ci= 0.95, cam = c("TUW28"))
  # overlapCI("fox", "coyote", sdata1, nb = 1000, ci= 0.95, cam = c("TUW28"), samp = TRUE)
  # overlapCI("fox", "coyote", sdata1, nb = 1000, ci= 0.95, cam = c("TUW28", "TUW27"))
  
  # This function uses sun time. The location is set to Toronto, Ontario, Canada by default

## Function body ===============================================================  
  
  locale = locale(date_format = "%Y-%m-%d",
                  time_format = "%H:%M:%S",
                  tz = "UTC")
  
  #filter cameras
  if (length(cam) != 0){
    data <- filter(data, data$RelativePath %in% cam)
  }
  
  ##format date, convert to POSIXct if it is not already
  data$Species <- as.character(data$Species)
  if (class(data$DateTime)[1] == "character"){ 
    data$DateTime <- anytime(data$DateTime)
  }
  
  ##Convert time to value between 0 and 2pi (radian)
  data$Time <- ((hour(data$DateTime)*60 + minute(data$DateTime))/1440) * 2*pi
  
  ##dates to a POSIXct object with the right time zone (GMT):
  dates_toronto <- data$DateTime
  
  ##create spatial point for the study area (latitude and longitude of Toronto, assume one point is accurate enoguh for whole study area)
  coords <- matrix(c(-79.3831843, 43.653226), nrow = 1)
  point <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
  
  ##Assign sun time to data$Time (replace clock time by sun time)  
  data$Time <- sunTime(data$Time, dates_toronto, point)
  
  
  ##create vectors of observation times for each species 
  sp1_ob <- data$Time[data$Species  == sp1] 
  sp2_ob <- data$Time[data$Species  == sp2] 
  
  ##Check if the number of detections is higher than 2 
  if (length(sp1_ob) < 2 || length(sp2_ob) < 2){
    return("Number of dections for sp1 or sp2 is less than 2")
  }
  
  ##Calculate bootstrap estimates of the overlap using (bootstrap() from the Overlap package)
  if ((length(sp1_ob) > 50) | (length(sp2_ob) > 50)){
    bootstat <- bootstrap(sp1_ob, sp2_ob, nb, type =  "Dhat4")  
  } else {
    bootstat <- bootstrap(sp1_ob, sp2_ob, nb, type =  "Dhat1")
  }
  
  ##overlap for the original dataset (overlapEst() from the Overlap package)
  if ((length(sp1_ob) > 50) | (length(sp2_ob) > 50)){
    overlap <- overlapEst(sp1_ob, sp2_ob, type =  "Dhat4")  
  } else {
    overlap <- overlapEst(sp1_ob, sp2_ob, type =  "Dhat1")
  }
  
  ##return the bootstrap interval
  if (samp == FALSE){
  return(bootCI(overlap, bootstat, ci))
  } else{
    return(overlap)
  }
}

## Function for regular overlap plots with sun times ===========================

overlapPlotST <- function(sp1, sp2, data, cam = "") {
  
  ## Notes =======================================================================  
  
  # sp1: first species in the pair
  # sp2: second species in the pair
  # data: dataframe obtained by importing .csv file generated with Timelapse
  # cam: string to specify the camera for which the data will be used in plotting 
  
  # This function uses sun time. The location is set to toronto, Ontario, Canada by default
  
  
  ## Test call====================================================================
  
   # overlapPlotST(sp1 = "fox", sp2 = "coyote", sdata1, cam = "TUW28") 
  # overlapPlotST(sp1 = "fox", sp2 = "coyote", sdata1, cam = "TUW27")
  
  ## Function body ===============================================================  
  
  locale = locale(date_format = "%Y-%m-%d",
                  time_format = "%H:%M:%S",
                  tz = "UTC")
  
  #filter to  camera
  data <- data[(data$RelativePath == cam),]
  
  
  ##format date, convert to POSIXct if it is not already
  data$Species <- as.character(data$Species)
  if (class(data$DateTime)[1] == "character"){ 
    data$DateTime <- anytime(data$DateTime)
  }
  
  ##Convert time to value between 0 and 2pi (radian)
  data$Time <- ((hour(data$Date)*60 + minute(data$Date))/1440) * 2*pi
  
  ##dates to a POSIXct object with the right time zone (GMT):
  dates_toronto <- data$DateTime
  
  ##create spatial point for the study area (latitude and longitude of Toronto, assume one point is accurate enough for whole study area)
  coords <- matrix(c(-79.3831843, 43.653226), nrow = 1)
  point <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs "))
  
  ##Assign sun time to data$Time (replace clock time by sun time)  
  data$Time <- sunTime(data$Time, dates_toronto, point)
  
  ##create vectors of observation times for each species 
  sp1_ob <- data$Time[data$Species  == sp1] #Species A
  sp2_ob <- data$Time[data$Species  == sp2] #Species B
  
  par(mfrow=c(1,1))
  
  plot = overlapPlot(sp1_ob, sp2_ob, rug = TRUE,
                     linetype = c(1, 2),
                     xscale = 360,
                     linecol = c("black", "blue"),
                     ylab = "Density",
                     xlab = "Sun position (degrees)",
                     main = paste(cam[1],", overlap estimate =", round(overlap, digits = 2), 
                                  ", Nsp1 =", length(sp1_ob),", Nsp2 =", length(sp2_ob)))
  
  legend("topright", legend = c(toString(sp1), toString(sp2)), lty=c(1, 2),
         col=c("black", "blue"),
         bg="white")
  
}


## Function to generate a grid overlap plot with sun times =====================

overlapPlotGrid <- function(sp1, sp2, data, grid, nb, cams = c()){

## Notes =======================================================================  

# sp1: first species in the pair
# sp2: second species in the pair
# data: dataframe obtained by importing .csv file generated with Timelapse
# grid: number of points for plotting the Von Mises density 
# nb: number of bootstrap samples 
# cams: specify the cameras to include in the plot grid 
  
# This function uses sun time. The location is set to Toronto, Ontario,
# Canada by default


## Test call====================================================================

  # overlapPlotGrid(sp1 = "fox", sp2 = "coyote", data = sdata1, grid = 100, nb = 500, cams = c("TUW11", "TUW1", "TUW4", "TUW09", "TUW09b"))
  # overlapPlotGrid(sp1 = "fox", sp2 = "coyote", data = sdata1, grid = 100, nb = 500, cams = c("TUW36", "TUW36b", "TUW37", "TUW37b", "TUW35a", "TUW35b", "TUW34", "TUW10", "TUW13", "TUW14"))
  # overlapPlotGrid(sp1 = "fox", sp2 = "coyote", data = sdata1, grid = 100, nb = 500, cams = c("TUW27", "TUW29b", "TUW28", "TUW25", "TUW24", "TUW23", "TUW19"))
  
## Function body ===============================================================  

  plot_list <- list()
  n <- 1 
  
  for (i in cams){
    
    plot = overlapPlotCI(sp1, sp2, data = x, grid, nb, type = 1, cam = c(i), title = FALSE, legend = FALSE)
    
  if (class(plot)[1] != "character"){ 
    
    plot_list[[n]] <- plot  
    n = n + 1  
    }
  }
  
  return(ggarrange(plotlist = plot_list))
  
}

## Function to generate times-to-encounter =====================================

getTimes <- function(sp1, sp2, data, interval = 0, rcens = NA, rtrunc = NA, excludecen = FALSE){
  
## Notes =======================================================================  

# sp1: first species in the pair
# sp2: second species in the pair
# interval: minimum interval (minutes) between detections of each species
# data: dataframe obtained by importing .csv file generated with Timelapse
# rcens: maximum time-to-encounter (hours) at which observations will be right censored
# rtrunc: maximum time-toe-encounter (hours) at which images will be excluded from the dataset 
# excludecen: if true, censored values are not included in the output dataset 

  
# This function returns a dataframe with a new column, the time difference 
# between each detection of sp2 after sp1. Only rows for detection with sp1 are 
# kept. 
  
# If a the sampling period after sp1 was detected and sp2 is not spotted, 
# then the time difference will be the time to the end and the observation will
# be right-censored (independent of rcens value).
  
# If sp1 is detected after sp1, then the observation is right-censored 
# (independent of rcens value).
  
## Test call====================================================================

# times <- getTimes("coyote", "fox", sdata1, 10, rcens  = 168)
# times <- getTimes("coyote", "fox", sdata1, 10, rtrunc = 168)
# times <- getTimes("dog", "fox", sdata1, 10, rtrunc = 10)
# times <- getTimes("dog", "fox", sdata1, 10, rtrunc = 168, excludecen = TRUE)
# times <- getTimes("dog", "fox", sdata1, 20, rtrunc = 168)
  
## Function body ===============================================================  

  
  data_org <- data
  
  data <- controlInterval(data_org, interval, keep = c(sp1, sp2))
  
  locale = locale(date_format = "%Y-%m-%d",
                  time_format = "%H:%M:%S",
                  tz = "UTC")
  
  ##sort data following time
  data <- data %>% dplyr::arrange(DateTime)
  
  times_df = data[0,]
  
  for (i in unique(data$RelativePath)){
    
    ##subset so only the two species of interest are included
    subdata <- filter(data, data$RelativePath == i)
    
    ##sort following time
    subdata <- subdata %>% dplyr::arrange(DateTime)
    
    subdata <- subdata %>% mutate(Time_diff = lead(DateTime) - DateTime)
    subdata <- subdata %>% mutate(Next_sp = lead(Species))
    
    ## convert from difftime to hours
    subdata$Time_diff <- as.numeric(subdata$Time_diff, units = "hours")
    
    
    ##getlast sampling time
    sub_org <- data_org[(data_org$RelativePath == i),]
    sub_org <- sub_org %>% dplyr::arrange(DateTime)
    t_last <- anytime(sub_org$DateTime[nrow(sub_org)])
    
    
    subdata$Time_diff[nrow(subdata)] <- (t_last - subdata$DateTime[nrow(subdata)])
    subdata$Next_sp[nrow(subdata)] <- "End"
    
    for (y in 1:nrow(subdata)){
      if ((subdata$Species[y] == subdata$Next_sp[y]) | (subdata$Next_sp[y] == "End")){
        subdata$Censor[y] <- 0
      } else {
        subdata$Censor[y] <- 1 
      }
      if ((is.na(rcens) == FALSE) & (subdata$Time_diff[y] > rcens)){
        subdata$Censor[y] <- 0
      }
    }
    
    if (excludecen == TRUE){
      subdata <- filter(subdata, subdata$Censor == 1)
    }
    
    times_df <- rbind(times_df, subdata)
  }
  
  ##Keep only observations when sp1 is the one that is observed 
  times_df <- times_df[(times_df$Species == sp1),]

  
  ## if a rcens value was entered, censor observations when the time interval before sp2 is higher than rcens
  if (is.na(rcens) == FALSE){
    for (i in 1:nrow(times_df)){
      if (times_df$Time_diff[i] > rcens){
        times_df$Time_diff[i] <- rcens 
      } 
      
    }
  }
  
  ## if a rtrunc value was entered, truncate the data above a certain maximum time point 
  if (is.na(rtrunc) == FALSE){
      times_df <- times_df[(times_df$Time_diff < rtrunc),] 
  }
  
  
  ## this is added to remove the time differences with the end of the sampling period that are smaller than the interval for the pictures. 
  times_df <- times_df[(times_df$Time_diff > (interval/60)),]
  
  return(times_df)
}

## Function to generate merge layers into timelapse export dataset=============

mergeLayers <- function(data1, data2){

## Notes =======================================================================  
   
  # data1: exported .csv file with the detections from TimeLapse
  # data2: .csv file with the covariates from the layers

## Test call====================================================================

  # data1 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Research/Professor_Molnar/Canid_temporal_segregation_manuscript/GitHubWildCanids/Datasets/TimelapseDataTUWall_Nsp_over0_V7.csv")
  # data2 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Research/Professor_Molnar/Canid_temporal_segregation_manuscript/GitHubWildCanids/Datasets/AllCameras_SpatialData_500mbuffer.csv")
  # mData <- mergeLayers(sdata1, sdata2) 

## Function body =============================================================== 

  #Change the name of the column that will be used to join the dataframes 
  data2 <- data2 %>% rename (RelativePath = Name)
  
  #merge the dataframes 
  output_data <- merge(data1, data2, col = "RelativePath")
  
  return(output_data)

}

## Function to run a cox proportional model=============

coxModel <- function(moddata, formula, m.null = FALSE){
  
  ## Notes =======================================================================  
  
  # moddata: Dataframe with the time-to-encounters 
  # formula: The formula for the cox proportional hazard model as a string 
  # m.null: set as TRUE if the model is null  
  
  ## Test call====================================================================
  
  # sdata3 <- getTimes("dog", "coyote", sdata1, 20, rtrunc = 48)
  # sdata3 <- data1[(data1$RelativePath %in% c("TUW19", "TUW26", "TUW28")),]
  # sdata3 <-  filter(data1, DateTime < anytime("2021-05-01 00:00:00"))
  # mData <- mergeLayers(sdata3, sdata2)
  
  # mData$Built_PA <- as.numeric(mData$Built_PA)
  # mData$Tree_PA <- as.numeric(mData$Tree_PA)
  
  #m_output <- coxModel(mData, "Surv(Time_diff, Censor) ~ Built_PA + Tree_PA")
  
  ## Function body =============================================================== 
  
  if (m.null == TRUE){
    
    #Kaplan Meier plot 
    km.model <- survfit(Surv(Time_diff, Censor) ~ RelativePath, type = "kaplan-meier", data = moddata)
    
    km.plot <- ggsurvplot(km.model, data = moddata, conf.int = TRUE)
    
    #cox proportional hazard model
    c.model <-  coxph(as.formula(formula), data = moddata)
    c.sum <- summary(c.model)
    
    #calculate the AIC for the model 
    c.AIC <- extractAIC(c.model)
    
    return(list("km.plot" = km.plot, "km.model" = km.model,
                "c.model" = c.model, "c.sum" = c.sum, "c.AIC" = c.AIC))
  } else {
  
    #Kaplan Meier plot 
    km.model <- survfit(Surv(Time_diff, Censor) ~ RelativePath, data = moddata, type = "kaplan-meier")
    
    km.plot <- ggsurvplot(km.model, data = moddata, conf.int = TRUE)
    
    #cox proportional hazard model
    c.model <-  coxph(as.formula(formula), data = moddata)
    
    c.sum <- summary(c.model)
    
    #linearity plot for the CPH model
    #c.linplot <- ggcoxfunctional(c.model, ylim = c(-3.5, 3.5), moddata)
    c.linplot <- 1
    
    #Schoenfeld individual test 
    c.proptest <- cox.zph(c.model)
    
    c.propplot <- ggcoxzph(cox.zph(c.model))
    
    #Deviance residuals
    c.devres <- ggcoxdiagnostics(c.model, type = "deviance",
                     linear.predictions = FALSE, ggtheme = theme_bw())
    
    #calculate the AIC for the model 
    c.AIC <- extractAIC(c.model)
    
    #model deviance
    c.dev <- anova(c.model)
    
    #return list of model information 
    results <- list("km.plot" = km.plot, "km.model" = km.model,
                    "c.model" = c.model, "c.sum" = c.sum, "c.linplot" = c.linplot,
                    "c.proptest" = c.proptest, "c.propplot" = c.propplot,
                    "c.devres" = c.devres, "c.AIC" = c.AIC, "c.dev" =  c.dev)  
    return(results)
  }
}


## Function to get species weekly and monthly frequency columns ================

addFreqSp <- function(data, sp){
  
  ## Notes =====================================================================
  
  # data: Dataframe with species detections from timelapse  (not the an output from getTimes)
  # sp: The name of the species for which a weekly and monthly frequency column 
  # will be added
  
  # This functions creates a weekly frequency and weekly column for the selected
  # species. A month and a weekly index columns are also added to the dataframe. 

  ## Test call====================================================================
  
  # outputdf <- addFreqSp(sdata1, c("fox", "coyote", "dog"))

  ## Function body =============================================================== 
  
  # Add an index column for weeks
  data$WeekIndex <- strftime(data$DateTime, format = "%Y-W%V")
  
  #Add an index column for months 
  data$MonthIndex <- strftime(data$DateTime, format = "%Y-%m")
  
  #calculate the frequencies of observations of each species per week 
  Weekdf <-  data %>% 
    group_by(Species, RelativePath, WeekIndex) %>%
    summarise(Nspecies = sum(Nspecies)) %>%
    #add rows with zeroes for when the species is not detected 
    complete(WeekIndex = data$WeekIndex, fill = list(Nspecies = 0))
  
  Weekdf <- data.frame(Weekdf)
  rename(Weekdf, WeekFreq = WeekIndex)
  
  #calculate the frequencies of observations of each species per month
  Monthdf <- data %>% 
    group_by(Species, RelativePath, MonthIndex) %>%
    summarise(Nspecies = sum(Nspecies)) %>%
    #add rows with zeroes for when the species is not detected 
    complete(MonthIndex = data$MonthIndex, fill = list(Nspecies = 0))
  
  Monthdf <- data.frame(Monthdf)
  rename(Monthdf, MonthFreq = MonthIndex)
  
  # for each species for which the weekly and monthly frequency was calculated,
  # merge columns with its monthly and weekly frequencies into the main 
  # dataset and rename them accordingly. 
  
  for (i in 1:length(sp)){
    
    #get all monthly and weekly frequencies and index for one species 
    Weekfreqs <- Weekdf[(Weekdf == sp[i]),] 
    Monthfreqs <- Monthdf[(Monthdf == sp[i]),]
    
    #merge the new columns to the main dataset using the weekly/monthly index 
    #column as the key
    temp1 <- merge(data, Weekfreqs, by = c("RelativePath", "WeekIndex"), all.x = TRUE)
    temp2 <- merge(data, Monthfreqs, by = c("RelativePath", "MonthIndex"), all.x = TRUE)
    
  
    #create and assign column names for the new column that contain the name of
    #the species whose frquency has been calculated  
    weekcolname <- paste("WeekFreq", sp[i], sep ="")
    monthcolname <- paste("MonthFreq", sp[i], sep ="")
    
    temp1 <- arrange(temp1, DateTime)
    temp2 <- arrange(temp2, DateTime)
    data <- arrange(data, DateTime)
    
    data[weekcolname] <- temp1$Nspecies.y
    data[monthcolname] <- temp2$Nspecies.y
    
  } 
  
  return(data) 
  
}

## Function to get human weekly and monthly frequency columns ==================

addFreqHum <- function(spdf, humdf, lim){
  
  ## Notes =======================================================================  
  
  # spdf: Dataframe with species detections from timelapse (not the output from getTimes) 
  # humdf: The name of the species for which a weekly and monthly frequency column 
  # will be added
  # lim: the independance interval between detection of indidividual humans  
  
  ## Test call====================================================================
  
  # Humans <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Research/Professor_Molnar/Canid_temporal_segregation_manuscript/Datasets/TimelapseData_all_imagesleaned.csv")
  # Species <- sdata1
  # outputdf <- addFreqHum(Species, Humans , 0.5)
  
  ## Function body =============================================================== 
  
  humdf <- humdf[(humdf$humans == "true"),]
  
  ##format date, convert to POSIXct if it is not already
  humdf$DateTime <- anytime(humdf$DateTime)
  
  ##sort following camera and time
  humdf<- humdf %>% dplyr::arrange(RelativePath, DateTime)
  
  # calculate the time differences between images 
  humdf <- humdf %>% mutate(Time_diff = DateTime - lag(DateTime))
  
  #find the cameras of each preceding row
  humdf <- humdf %>% mutate(Prevcam = lag(RelativePath))
  
  ## convert the time difference from difftime to hours
  humdf$Time_diff <- as.numeric(humdf$Time_diff, units = "mins")
  
  ##set Time_diff to limit for the first observation of each new camera so it is not removed later
  humdf$Time_diff[humdf$RelativePath != humdf$Prevcam] <- lim
  
  ## drop rows within time limit
  humdf2 <- humdf[!(humdf$Time_diff < lim),]
  humdf2 <- humdf2[(humdf2$RelativePath == humdf2$Prevcam),]
  humdf2 <- humdf2[2:nrow(humdf2),]
  
  #drop columns used in the function 
  drop <- c("Time_diff","Prevcam")
  humdf2 <- humdf2[,!(names(humdf2) %in% drop)]
  
  ##renumber rows
  rownames(humdf2) <- NULL
  
  # convert the human entries (which contains only "true") to 1 for summation
  humdf2$humans <- 1
  
  # Add indexes for months and weeks 
  humdf2$WeekIndex <- strftime(humdf2$DateTime, format = "%Y-W%V")
  humdf2$MonthIndex <- strftime(humdf2$DateTime, format = "%Y-%m")
  spdf$WeekIndex <- strftime(spdf$DateTime, format = "%Y-W%V")
  spdf$MonthIndex <- strftime(spdf$DateTime, format = "%Y-%m")
  
  # summarize the weekly frequency of human detections and store it in a dataframe 
  Weekdf <-  humdf2 %>% 
    group_by( RelativePath, WeekIndex) %>%
    summarise(humans = sum(humans)) %>%
    complete(WeekIndex = humdf2$WeekIndex, fill = list(humans = 0))
  
  Weekdf <- data.frame(Weekdf)
  rename(Weekdf, WeekFreq = WeekIndex)
  
  # summarize the monthly frequency of human detections and store it in a dataframe 
  Monthdf <- humdf2 %>% 
    group_by(RelativePath, MonthIndex) %>%
    summarise(humans = sum(humans)) %>%
    complete( MonthIndex = humdf2$MonthIndex, fill = list(humans = 0))
  
  Monthdf <- data.frame(Monthdf)
  rename(Monthdf, MonthFreq = MonthIndex)
  
  # merge the species detection dataframe with the ones containing the frequencies to make  temporary dataframes
  temp1 <- merge(spdf, Weekdf, by = c("RelativePath", "WeekIndex"), all.x = TRUE)
  temp2 <- merge(spdf, Monthdf, by = c("RelativePath", "MonthIndex"), all.x = TRUE)
  
  #sort the dataframes 
  temp1 <- arrange(temp1, DateTime)
  temp2 <- arrange(temp2, DateTime)
  spdf <- arrange(spdf, DateTime)
  
  #add the frequency columns to the original Species dataframe 
  spdf["WeekFreqhumans"] <- temp1$humans.y
  spdf["MonthFreqhumans"] <- temp2$humans.y
  
  return(spdf) 
  
}
