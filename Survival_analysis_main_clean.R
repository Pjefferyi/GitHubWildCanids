library(nlme)
library(sjPlot)
library(webshot)
library(splines)
library(MuMIn)
library(simPH)
library(car)
require(survival)
library(circular)


#====================== constructing the dataset ==============================
# Load data
spatialdata <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/Layers/AllCameras_AllBuffers_500B.csv")
humanobs <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/TimelapseData_all_images_v2.csv")
spobs <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/2022_full/TimelapseData_2021_Combined.csv")
#spobs <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/TimelapseDataTUWall_Nsp_over0_V7.csv")

# subset to get only sites and period of interest
spobs <- spobs[(spobs$RelativePath %in% c("TUW28", "TUW26", "TUW19")),]
spobs <-  filter(spobs, DateTime < anytime("2021-05-01 00:00:00"))
spobs <-  filter(spobs, DateTime > anytime("2020-09-30 23:59:59"))

# set species detection interval of confidence (in minutes)
spobs1 <- controlInterval(spobs, 2, keep = c("dog"))
spobs2 <- controlInterval(spobs, 30, keep = c("fox"))
spobs3 <- controlInterval(spobs, 30, keep = c("coyote"))
spobs <- rbind(spobs1, spobs2, spobs3)

# Add species frequency covariates
spobs <- addFreqSp(spobs, c("fox", "coyote", "dog"))
spobs <- addFreqHum(spobs, humanobs, 0.5) # 0.5 minutes is the independence interval for humans

# Add a column for time of day in radians (clocktime)
spobs$Time <-((hour(spobs$DateTime)*60 + minute(spobs$DateTime))/1440) * 2*pi

# Add a column for time of day in radians (Suntime)
dates_toronto <- spobs$DateTime
coords <- matrix(c(-79.3831843, 43.653226), nrow = 1)
point <- SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat + datum=WGS84"))
spobs$SunTime <- sunTime(spobs$Time, dates_toronto, point)

# Add a column for day and night
spobs$Period <- "night"
spobs[((spobs$SunTime > 0.5* pi) & (spobs$SunTime < 1.5*pi)),]$Period <- "day"

# Add spatial covariates
spobs <- mergeLayers(spobs, spatialdata)

# drop columns that are not useful
drop <- c("X.1","ImageQuality", "DeleteFlag", "melanistic", "mange",
          "SpeciesExtra", "UWIN", "Flag0", "Flag1", "speciesextraname",
          "checked", "flagged", "vehicles", "X", "OBJECTID..", "SHAPE..",
          "y", "x", "status", "ORIG_FID", "SHAPE_Area", "SHAPE_Length" )
spobs <- spobs[, !(names(spobs) %in% drop)]

write.csv(spobs, "C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/TimelapseDataTUWall_WaitTimeDataset.csv" )

#======================  Cox proportional hazard models ========================
# load dataset 
dataset <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/TimelapseDataTUWall_WaitTimeDataset.csv")
dataset2 <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/TimelapseDataTUWall_WaitTimeDataset2.csv")


#commands for model info:"km.plot", "km.model", "c.model", "c.sum", "c.linplot",
#"c.proptest", "c.propplot", "c.devres", "c.AIC", "c.dev"

#format:  modelname$[command]
#e.g. modfd1$c.model, modfd1$km.plot


# Dog after fox models

# import waiting times 
datafoxdog <- getTimes("fox", "dog", interval = 0, dataset, rtrunc = 168)

# construct models
modfdnull <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ 1", m.null = TRUE)
modfd1 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA")
modfd2 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Built_PA")
modfd3 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ WeekFreqdog")
modfd4 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ MonthFreqdog") # best frequency predictor 
modfd5 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ WeekFreqfox") # best frequency predictor 
modfd6 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ MonthFreqfox")
modfd7 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ WeekFreqcoyote")
modfd8 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ MonthFreqcoyote")
modfd9 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ WeekFreqfox + WeekFreqcoyote")
modfd10 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqdog")
modfd11 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqfox")
modfd12 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Built_PA +  WeekFreqfox")
modfd13 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqfox")
modfd14 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ WeekFreqdog + WeekFreqfox")
modfd15 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ WeekFreqdog+ WeekFreqcoyote")
modfd16 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqdog + WeekFreqfox +WeekFreqcoyote")
modfd17 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Built_PA + WeekFreqdog + WeekFreqfox + WeekFreqcoyote ")
modfd18 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqdog + WeekFreqfox + WeekFreqcoyote")
modfd19 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqdog+ WeekFreqfox")
modfd20 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Built_PA + WeekFreqdog + WeekFreqfox")
modfd21 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqdog+ WeekFreqfox")
modfd22 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqdog")
modfd23 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Built_PA + WeekFreqdog")
modfd24 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqcoyote")
modfd25 <- coxModel(datafoxdog, "Surv(Time_diff, Censor) ~ Built_PA +  WeekFreqcoyote")

# model results 
tab_model(list(modfd1$c.model, modfd2$c.model, modfd3$c.model, modfd4$c.model,
               modfd5$c.model, modfd6$c.model, modfd7$c.model, modfd8$c.model,
               modfd9$c.model, modfd10$c.model, modfd11$c.model, 
               modfd12$c.model, modfd13$c.model, modfd14$c.model,
               modfd15$c.model, modfd16$c.model, modfd17$c.model,
               modfd18$c.model, modfd19$c.model, modfd20$c.model,
               modfd21$c.model, modfd22$c.model,modfd23$c.model,
               modfd24$c.model, modfd25$c.model),
          show.aicc = TRUE, show.dev = TRUE, show.loglik = TRUE)  
#file = "C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/foxdogmodels.html")



# model comparison table
mod_tablefd <- model.sel(list(modfd1$c.model, modfd2$c.model, modfd3$c.model, 
                              modfd4$c.model, modfd5$c.model, modfd6$c.model, 
                              modfd7$c.model, modfd8$c.model, modfd9$c.model, 
                              modfd10$c.model, modfd11$c.model,modfd12$c.model,
                              modfd13$c.model, modfd14$c.model,modfd15$c.model,
                              modfd16$c.model, modfd17$c.model,modfd18$c.model,
                              modfd19$c.model, modfd20$c.model,modfd21$c.model,
                              modfd22$c.model, modfd23$c.model,modfd24$c.model,
                              modfd25$c.model, modfdnull$c.model))


# Dog after Coyote model
# import waiting times 
datacoydog <- getTimes("coyote", "dog", interval = 0, dataset, rtrunc = 168)

# construct models
modcdnull <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ 1", m.null = TRUE)
modcd1 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA")
modcd2 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Built_PA")
modcd3 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ WeekFreqdog") 
modcd4 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ MonthFreqdog") # best freq predictor 
modcd5 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ WeekFreqfox")
modcd6 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ MonthFreqfox")
modcd7 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ WeekFreqcoyote")
modcd8 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ MonthFreqcoyote")
modcd9 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Built_PA + WeekFreqdog")
modcd10 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ WeekFreqdog + WeekFreqcoyote")
modcd11 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Built_PA + WeekFreqdog + WeekFreqcoyote")
modcd12 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Built_PA + WeekFreqcoyote")
modcd13 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqdog")
modcd14 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqdog + WeekFreqcoyote")
modcd15 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + WeekFreqcoyote")
modcd16 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqdog")
modcd17 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqdog + WeekFreqcoyote")
modcd18 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqcoyote")
modcd19 <- coxModel(datacoydog, "Surv(Time_diff, Censor) ~ Tree_PA + Built_PA + WeekFreqcoyote + WeekFreqdog")

# model results
tab_model(list(modcd1$c.model, modcd2$c.model, modcd3$c.model, modcd4$c.model,
               modcd5$c.model, modcd6$c.model, modcd7$c.model, modcd8$c.model,
               modcd9$c.model, modcd10$c.model, modcd11$c.model, modcd12$c.model,
               modcd13$c.model, modcd14$c.model, modcd15$c.model, modcd16$c.model,
               modcd17$c.model,modcd18$c.model, modcd19$c.model, modcd20$c.model,
               modcd21$c.model,modcd22$c.model,modcd23$c.model,modcd24$c.model,
               modcd25$c.model, modcd26$c.model, modcd27$c.model, 
               modcd28$c.model, modcd29$c.model),
          show.aicc = TRUE,
          show.dev = TRUE, show.loglik = TRUE,
          file = "C:/Users/Jelan/OneDrive/Desktop/Courses, notes and assignments/Fourth Year Fall 2021/BIOD98-DirectedResearch in Biology/ImageDataExports/coyotedogmodels.html")

# model comparison table 
mod_tablecd <- model.sel(list(modcd1$c.model, modcd2$c.model, modcd3$c.model, modcd4$c.model,
                              modcd5$c.model, modcd6$c.model, modcd7$c.model, modcd8$c.model,
                              modcd9$c.model, modcd10$c.model, modcd11$c.model, modcd12$c.model,
                              modcd13$c.model, modcd14$c.model, modcd15$c.model, modcd16$c.model,
                              modcd17$c.model, modcd18$c.model, modcd19$c.model,
                              modcdnull$c.model))

#========== Test significance of overlap coefficient ===========================

# Create samples for the detection data 
T28fox <- dataset[((dataset$Species == "fox") & (dataset$RelativePath == "TUW28")),]$SunTime
T26fox <- dataset[((dataset$Species == "fox") & (dataset$RelativePath == "TUW26")),]$SunTime
T28coy <- dataset[((dataset$Species == "coyote") & (dataset$RelativePath == "TUW28")),]$SunTime
T26coy <- dataset[((dataset$Species == "coyote") & (dataset$RelativePath == "TUW26")),]$SunTime
T28dog <- dataset[((dataset$Species == "dog") & (dataset$RelativePath == "TUW28")),]$SunTime
T26dog <- dataset[((dataset$Species == "dog") & (dataset$RelativePath == "TUW26")),]$SunTime

T28foxdogoverlap <- overlapEst(T28fox, T28dog, type =  "Dhat4") 
T26foxdogoverlap <- overlapEst(T26fox, T26dog, type =  "Dhat4") 
T28foxdog <- bootstrap(T28fox, T28dog, 1000, type = "Dhat4")
T26foxdog <- bootstrap(T26fox, T26dog, 1000, type = "Dhat4")

T28foxdogCI <- bootCI(T28foxdogoverlap, T28foxdog, 0.84)
T26foxdogCI <- bootCI(T26foxdogoverlap, T28foxdog, 0.84)

T28coydogoverlap <- overlapEst(T28coy, T28dog, type =  "Dhat4") 
T26coydogoverlap <- overlapEst(T26coy, T26dog, type =  "Dhat4") 
T28coydog <- bootstrap(T28coy, T28dog, 1000, type = "Dhat4")
T26coydog <- bootstrap(T26coy, T26dog, 1000, type = "Dhat4")

T28coydogCI <- bootCI(T28coydogoverlap, T28coydog, 0.84)
T26coydogCI <- bootCI(T26coydogoverlap, T28coydog, 0.84)

#========================== hazard ratio Plots   ===============================
modsp <- coxModel(alltimes, "Surv(Time_diff, Censor) ~ Species")

# For both species 
fit <- survfit(Surv(Time_diff, Censor)~ Built_PA, data = alltimes)

ggsurvplot(fit,fun = "event", xlim = c(0, max(alltimes$Time_diff)), conf.int = TRUE) 

#by built area for foxes
fit2 <- survfit(Surv(Time_diff, Censor)~ Built_PA, data = datafoxdog)

ggsurvplot(fit2,fun = "event", xlim = c(0, max(datafoxdog$Time_diff)), 
           legend.labs = c("13%", "22%", "35%", "37%"),
           xlab = "Time after fox detection (hours)", ylab = "Probability of dog visit",
           legend.title = "Percent built area") 

#By built area for coyotes 
fit3 <- survfit(Surv(Time_diff, Censor)~ Built_PA, data = datacoydog)

ggsurvplot(fit3,fun = "event", xlim = c(0, max(datacoydog$Time_diff)), 
           legend.labs = c("13%", "22%", "35%", "37%"),
           xlab = "Time after coyote detection (hours)", ylab = "Probability of dog visit",
           legend.title = "Percent built area")


# Simulate and plot Hazard Ratios 
Sim1 <- coxsimLinear(modfd23$c.model, b = "Built_PA", 
                     qi = "Hazard Ratio", ci = 0.95,
                     Xj = seq(10, 40, by = 5), spin = TRUE)

simGG(Sim1, xlab = "", ylab = "", rug = FALSE,  pcolour = "black", lcolour = "black") +
  font_size(axis_title.x = 18, axis_title.y = 18, labels.x = 16, labels.y = 16, base.theme = theme_classic()) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank())



Sim2 <- coxsimLinear(modfd23$c.model, b = "WeekFreqdog", 
                     qi = "Hazard Ratio", ci = 0.95,
                     Xj = seq(0, 300, by = 20), spin = TRUE)

simGG(Sim2, xlab = "", ylab = "", rug = FALSE,  pcolour = "black", lcolour = "black") +
  font_size(axis_title.x = 18, axis_title.y = 18, labels.x = 16, labels.y = 16, base.theme = theme_classic()) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank())


Sim3 <- coxsimLinear(modcd9$c.model, b = "Built_PA", 
                     qi = "Hazard Ratio", ci = 0.95,
                     Xj = seq(10, 40, by = 5), spin = TRUE)

simGG(Sim3, xlab = "", ylab = "", pcolour = "black", lcolour = "black") + font_size(axis_title.x = 18, axis_title.y = 18, labels.x = 16, labels.y = 16, base.theme = theme_classic())



Sim4 <- coxsimLinear(modcd9$c.model, b = "WeekFreqdog", 
                     qi = "Hazard Rate", ci = 0.95,
                     Xj = seq(0,300, by = 20), spin = TRUE)

simGG(Sim4, xlab = "", ylab = "", pcolour = "black", lcolour = "black") + font_size(axis_title.x = 18, axis_title.y = 18, labels.x = 16, labels.y = 16, base.theme = theme_classic())


#====================== model diagnostics ======================================

#Plots of the continuous explanatory variable against martingale residuals of null cox proportional hazards model
ggcoxfunctional(modfd8$c.model, data = datafoxdog, ylim = c(-1,1)) 

#residuals by type
ggcoxdiagnostics(modfd5$c.model, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

#test for hazard proportionality 
ggcoxzph(cox.zph(modcd14$c.model))
plot(cox.zph(modcd14$c.model)[1])

model.avg(list(modcd6$c.model, modcd13$c.model, modcd8$c.model, modcd2$c.model))
