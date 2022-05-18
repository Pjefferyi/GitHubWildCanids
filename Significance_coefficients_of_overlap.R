library(overlap)


# This code was used to calculate the 84% confidence interval of the coefficients
# of overlap.These intervals were compared to infer whether differences in the 
# coefficients of overlap across sites were significant.


# load dataset 
dataset <- read.csv("C:/Users/Jelan/OneDrive/Desktop/Research/Professor_Molnar/Canid_temporal_segregation_manuscript/GitHubWildCanids/Datasets/TimelapseDataTUWall_WaitTimeDataset.csv")

#===========Test significance of changes in the overlap coefficient ============

# Create samples of sun times for the detection data 
T28fox <- dataset[((dataset$Species == "fox") & (dataset$RelativePath == "TUW28")),]$SunTime
T26fox <- dataset[((dataset$Species == "fox") & (dataset$RelativePath == "TUW26")),]$SunTime
T28coy <- dataset[((dataset$Species == "coyote") & (dataset$RelativePath == "TUW28")),]$SunTime
T26coy <- dataset[((dataset$Species == "coyote") & (dataset$RelativePath == "TUW26")),]$SunTime
T28dog <- dataset[((dataset$Species == "dog") & (dataset$RelativePath == "TUW28")),]$SunTime
T26dog <- dataset[((dataset$Species == "dog") & (dataset$RelativePath == "TUW26")),]$SunTime

# calculate the overlap for foxes and dogs in TUW28 (TU28 = low development) 
# and TUW26 (TUW26 = medium development)
T28foxdogoverlap <- overlapEst(T28fox, T28dog, type =  "Dhat4") 
T26foxdogoverlap <- overlapEst(T26fox, T26dog, type =  "Dhat4") 
T28foxdog <- bootstrap(T28fox, T28dog, 1000, type = "Dhat4")
T26foxdog <- bootstrap(T26fox, T26dog, 1000, type = "Dhat4")

#calculate the 84% confidence interval for the overlap coefficient at each site 
T28foxdogCI <- bootCI(T28foxdogoverlap, T28foxdog, 0.84)
T26foxdogCI <- bootCI(T26foxdogoverlap, T28foxdog, 0.84)

# calculate the overlap for coyotes and dogs in TUW28 (TU28 = low development) 
# and TUW26 (TUW26 = medium development)
T28coydogoverlap <- overlapEst(T28coy, T28dog, type =  "Dhat4") 
T26coydogoverlap <- overlapEst(T26coy, T26dog, type =  "Dhat4") 
T28coydog <- bootstrap(T28coy, T28dog, 1000, type = "Dhat4")
T26coydog <- bootstrap(T26coy, T26dog, 1000, type = "Dhat4")

#calculate the 84% confidence interval for the overlap coefficient at each site 
T28coydogCI <- bootCI(T28coydogoverlap, T28coydog, 0.84)
T26coydogCI <- bootCI(T26coydogoverlap, T28coydog, 0.84)