
# Library -----------------------------------------------------------------

#Libraries
library(lattice)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(Matrix)
library(foreach)
library(magrittr)
library(parallel)
library(sp)
library(INLA)
library(ggplot2)
#library(ggregplot) # not available for this version of r (4.3.1)
library(ggmap)
library(gstat)
library(ncdf4)
library(raster)
library(sf)
library(spdep)
library(tidyverse)
library(rworldmap)
library(maptools)
library(rworldxtra)
library(grDevices)
library(colorRamps)
library(RColorBrewer)
library(plyr)
library(fields)
library(corrplot)
library(dplyr)

# Loading in Data ---------------------------------------------------------

setwd("C:/Users/walkera2019/Dropbox/Recovery Gamma/Manuscript Review 1")


data = read.csv("Site_Coeff.10.24.23.csv")

data$Ocean_Name[which(data$Ocean_Name == "Red Sea")] <- "Indian"
data$Ocean_Name[which(data$Ocean_Name == "Arabian Gulf")] <- "Indian"
#hist(data$Historical_SST_Max, breaks = 10)

data$Years_Since_HW = as.numeric(as.Date(data$Start_Date, format = "%Y-%m-%d") - as.Date(data$Last_Heatwave_Date, format = "%Y-%m-%d"))/365

data$Longitude_Degrees[which(data$Longitude_Degrees < 0)] = data$Longitude_Degrees[which(data$Longitude_Degrees < 0)] + 360 # Shifting coordinates

data$Mean_Windspeed_DR[which(is.na(data$Mean_Windspeed_DR) == T)] <- 0


# Grabbing vars -----------------------------------------------------------

data_raw <- read.csv("C:/Users/walkera2019/Dropbox/Data Directory/All_Cover_All_Variables_10_2_23_density_turb.csv")
hist(data_raw$Density_Distance_100km)
hist(data_raw$Density_Distance_10km)
hist(data_raw$Density_Distance_1km)

data_raw$Reef_Density_100km[which(data_raw$Density_Distance_100km > 1)] <- NA
data_raw$Reef_Density_10km[which(data_raw$Density_Distance_10km > 1)] <- NA
data_raw$Reef_Density_1km[which(data_raw$Density_Distance_1km > 1)] <- NA


names(data_raw)
data$Population <- NA
data$Temperature_Maximum <- NA
data$Temperature_Kelvin_Standard_Deviation <- NA
data$Reef_Density_1km <- NA
data$Reef_Density_10km <- NA
data$Reef_Density_100km <- NA
data$Country_Name <- NA
data$Citation <- NA
blank <- data[0,]
for (i in levels(factor(data$Site_ID))) {
  x0 <- subset(data, Site_ID == i)
  y0 <- subset(data_raw, Site_ID == i)
  for (i1 in levels(factor(x0$New_ID))) {
    x <- subset(x0, New_ID == i1)
    y <- subset(y0, Date_Year <= x$Ending_Year)
    
    
    x$Population <- as.numeric(unique(y$Population))
    x$Reef_Density_1km <- as.numeric(unique(y$Reef_Density_1km))
    x$Reef_Density_10km <- as.numeric(unique(y$Reef_Density_10km))
    x$Reef_Density_100km <- as.numeric(unique(y$Reef_Density_100km))
    x$Country_Name <- unique(y$Country_Name)
    x$Citation <- unique(y$Citation)
    if (nlevels(factor(y$Temperature_Maximum)) > 0) {
      x$Temperature_Maximum <- max(as.numeric(na.omit(unique(y$Temperature_Maximum)))-273.15)
    }
    if (nlevels(factor(y$Temperature_Kelvin_Standard_Deviation)) > 0) {
      x$Temperature_Kelvin_Standard_Deviation <- mean(na.omit(unique(y$Temperature_Kelvin_Standard_Deviation)))
    }
    
    blank <- rbind(blank, x) 
  }
}

data <- blank


nlevels(factor(data$Country_Name))
# Corrplot ----------------------------------------------------------------


#
df = select_if(data, is.numeric)
names(df)#
df = subset(df, select = -c(r_coeff, Latitude_Degrees, Longitude_Degrees, Site_ID, New_ID, Ending_Year, N_samples, A_coeff, C_coeff, F_coeff, D_coeff, Velocity_distance, Skewness_distance, Kurtosis_distance, Hist_Max_SST_distance))
M = cor(df, use = "pairwise.complete.obs")

corrplot(M, method = "number")

png("Corplot_10.19.23.png", width = 2000, height = 2000)
corrplot(M, method = "number")
dev.off()

# Standardizing Variables across oceans -----------------------------------

standardize_function<-function(x){
  x.standardized=(x-mean(na.omit(x)))/sd(na.omit(x))
  return(x.standardized)
}


## Continuous Variables
data$Turbidity = standardize_function(data$Turbidity) #kd490
data$Depth = standardize_function(data$Depth) #depth in meters
data$Grav_tot = standardize_function(data$Grav_tot) # total gravity, human pressure metric
data$Velocity = standardize_function(data$Velocity) # SST Velocity, metric of SST change derived from Rob's prior work
data$Kurtosis = standardize_function(data$Kurtosis) # SST Kurtosis, derived from Rob's prior work
data$Skewness = standardize_function(data$Skewness) # SST Skewness, derived from Rob's prior work
data$Distance_to_shore = standardize_function(data$Distance_to_shore) # Distance to shore in meters
data$Population = standardize_function(data$Population) # Human population using nearest neighbor within buffer. Provided by Chelsey
data$Historical_Max_SST = standardize_function(data$Historical_Max_SST) # Max historical SST derived from Shannon's prior work
# data$DHW_1_Yr_Mean = standardize_function(data$DHW_1_Yr_Mean) # Mean DHW value at site within past year
# data$DHW_1_Yr_SD = standardize_function(data$DHW_1_Yr_SD) # SD of DHW at site within past year
# data$DHW_1_Yr_Max = standardize_function(data$DHW_1_Yr_Max) # Max DHW at site within past year
# data$DHW_5_Yr_Mean = standardize_function(data$DHW_5_Yr_Mean) # Mean DHW value at site within past five years
# data$DHW_5_Yr_SD = standardize_function(data$DHW_5_Yr_SD) # SD DHW value at site within past five years
# data$DHW_5_Yr_Max = standardize_function(data$DHW_5_Yr_Max) # Max DHW value at site within past five years

data$Cyclone_Freq = standardize_function(data$Cyclone_Freq) # Frequency of cyclones from 1970 to starting date
data$DHW_Freq_Prior = standardize_function(data$DHW_Freq_Prior) # Frequency of cyclones from 1970 to starting date

data$Initial_Cover = standardize_function(data$Initial_Cover) # Initial cover for recovery series
data$Initial_Algae = standardize_function(data$Initial_Algae) # Initial macroalgae cover for recovery series

data$Last_Heatwave_Intensity = standardize_function(data$Last_Heatwave_Intensity) # Intensity of most recent heatwave on reef prior to start date
data$Years_Since_HW = standardize_function(data$Years_Since_HW) # Number of years since last heatwave on reef prior to start date

data$Mean_DHW_DR <- standardize_function(data$Mean_DHW_DR)
data$Max_DHW_DR <- standardize_function(data$Max_DHW_DR)
data$SD_DHW_DR <- standardize_function(data$SD_DHW_DR)
data$Heatwaves_4_Number_DR <- standardize_function(data$Heatwaves_4_Number_DR)
data$Heatwaves_8_Number_DR <- standardize_function(data$Heatwaves_8_Number_DR)
data$Cyclone_Number_DR <- standardize_function(data$Cyclone_Number_DR)
data$Mean_Windspeed_DR <- standardize_function(data$Mean_Windspeed_DR)
data$Max_Windspeed_DR <- standardize_function(data$Max_Windspeed_DR)
data$Number_of_Disturbances_DR <- standardize_function(data$Number_of_Disturbances_DR)

data$Temperature_Maximum <- standardize_function(data$Temperature_Maximum)
data$Temperature_Kelvin_Standard_Deviation <- standardize_function(data$Temperature_Kelvin_Standard_Deviation)
data$Reef_Density_100km <- standardize_function(data$Reef_Density_100km)
data$Reef_Density_10km <- standardize_function(data$Reef_Density_10km)
data$Reef_Density_1km <- standardize_function(data$Reef_Density_1km)

data$High_DHW_Freq = standardize_function(data$High_DHW_Freq)
data$Low_DHW_Freq = standardize_function(data$Low_DHW_Freq) 

data$High_DHW_Number_DR = standardize_function(data$High_DHW_Number_DR)
data$Low_DHW_Number_DR = standardize_function(data$Low_DHW_Number_DR) 

data$Exposure[which(data$Exposure == "")] <- NA
data$Exposure <- factor(data$Exposure, levels = c("Sheltered", "Sometimes", "Exposed"))


# Testing artifacts that should be accounted for
data$N_samples = standardize_function(data$N_samples)
data$Analysis_Type <- factor(data$Analysis_Type)
# Subsetting by Ocean -----------------------------------------------------


Global <- data

# Establishing mesh for Oceans and Globe ----------------------------------

## Global
# Bind sites
Global_Loc = cbind(Global$Longitude_Degrees, Global$Latitude_Degrees)
#2 Setting up mesh
MeshPred = inla.mesh.2d(Global_Loc, max.edge = c(10, 24))
spde.pred= inla.spde2.matern(mesh=MeshPred, alpha=2)
s.index.p= inla.spde.make.index(name="sp.field.pred", n.spde=spde.pred$n.spde)
plot(MeshPred)
A2_global <- inla.spde.make.A(MeshPred, loc = Global_Loc)
dim(A2_global)  #
# 4 Create SPDE
clm.spde <- inla.spde2.pcmatern(mesh = MeshPred, #alpha = 2, #alpha is the default for 2d, not for time series
                                prior.range = c(50, 0.9), # P(range < 50) = 0.9
                                prior.sigma = c(1, 0.01) # P(sigma > 10) = 0.01       )
)
# Define the Matern correlation on the mesh
spde.globe   <- inla.spde2.matern(MeshPred, alpha=2)
#define spatial random field
w.index.globe <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde.globe$n.spde,
  n.group = 1,
  n.repl  = 1)



remove(MeshPred, spde.pred, clm.spde, Global_Loc)




# Creating X matrix for each group ----------------------------------------

## Global
# Creating dataframe for vars
X_Global = data.frame(Intercept = rep(1, nrow(Global) ), 
                      #Random Variables
                      Site_ID = Global$Site_ID,
                      Time = Global$Starting_Year,
                      
                      ## Continuous Variables
                      Initial_Cover = Global$Initial_Cover, # Initial cover for recovery series
                      Initial_Algae = Global$Initial_Algae, # Initial macroalgae cover for recovery series
                      #
                      Turbidity = Global$Turbidity, #kd490
                      Depth = Global$Depth, #depth in meters
                      #
                      Grav_tot = Global$Grav_tot, # total gravity, human pressure metric
                      Population = Global$Population, # Human population using nearest neighbor within buffer. Provided by Chelsey
                      Distance_to_shore = Global$Distance_to_shore, # Distance to shore in meters
                      #
                      Reef_Density_100km = data$Reef_Density_100km,
                      Reef_Density_10km = data$Reef_Density_10km,
                      Reef_Density_1km = data$Reef_Density_1km,
                      #
                      Velocity = Global$Velocity, # SST Velocity, metric of SST change derived from Rob's prior work
                      Kurtosis = Global$Kurtosis, # SST Kurtosis, derived from Rob's prior work
                      Skewness = Global$Skewness, # SST Skewness, derived from Rob's prior work
                      Hist_Max_SST = Global$Historical_Max_SST, # Max historical SST derived from Shannon's prior work
                      #
                      # DHW_1_Yr_Mean = Global$DHW_1_Yr_Mean, # Mean DHW value at site within past year
                      # DHW_1_Yr_SD = Global$DHW_1_Yr_SD, # SD of DHW at site within past year
                      # DHW_1_Yr_Max = Global$DHW_1_Yr_Max, # Max DHW at site within past year
                      # DHW_5_Yr_Mean = Global$DHW_5_Yr_Mean, # Mean DHW value at site within past five years
                      # DHW_5_Yr_SD = Global$DHW_5_Yr_SD, # SD DHW value at site within past five years
                      # DHW_5_Yr_Max = Global$DHW_5_Yr_Max, # Max DHW value at site within past five years
                      Cyclone_Freq = Global$Cyclone_Freq, # Frequency of cyclones from 1970 to starting date
                      Heatwave_4_Freq = Global$Heatwave_4_Freq, # Frequency of DHW >= 4 events from 1985 to starting year
                      Heatwave_8_Freq = Global$Heatwave_8_Freq, # Frequency of DHW >= 8 events from 1985 to starting year
                      Heatwave_4_8_Freq = Global$Heatwave_4_8_Freq, # Frequency of events between 4 and 8 DHW from 1985 to starting year
                      
                      Heatwaves_4_Number = Global$Heatwaves_4_Number, # Number of heatwaves with DHW >= 4
                      Heatwaves_8_Number = Global$Heatwaves_8_Number, # Number of heatwaves with DHW >= 8
                      Heatwaves_4_8_Number = Global$Heatwaves_4_8_Number,# Number of heatwaves between 4 and 8 DHW
                      
                      Low_DHW_Freq = Global$Low_DHW_Freq,
                      High_DHW_Freq = Global$High_DHW_Freq,
                      
                      High_DHW_Number_DR = Global$High_DHW_Number_DR,
                      Low_DHW_Number_DR = Global$Low_DHW_Number_DR,
                      
                      
                      Low_DHW_Freq_Prior = Global$Low_DHW_Freq_Prior,
                      High_DHW_Freq_Prior = Global$High_DHW_Freq_Prior,
                      DHW_Freq_Prior = Global$DHW_Freq_Prior,
                      
                      
                      Last_Heatwave_Intensity = Global$Last_Heatwave_Intensity, # Intensity of most recent heatwave on reef prior to start date
                      Years_Since_HW = Global$Years_Since_HW, # Number of years since last heatwave on reef prior to start date
                      
                      Temperature_Max = Global$Temperature_Maximum,
                      Temperature_SD = Global$Temperature_Kelvin_Standard_Deviation,
                      
                      Mean_DHW_DR = Global$Mean_DHW_DR,
                      Max_DHW_DR = Global$Max_DHW_DR,
                      SD_DHW_DR = Global$SD_DHW_DR,
                      
                      Heatwaves_4_Number_DR = Global$Heatwaves_4_Number_DR,
                      Heatwaves_8_Number_DR = Global$Heatwaves_8_Number_DR,
                      Cyclone_Number_DR = Global$Cyclone_Number_DR,
                      Mean_Windspeed_DR = Global$Mean_Windspeed_DR,
                      Max_Windspeed_DR = Global$Max_Windspeed_DR,

                      ## Categorical
                      Ocean_Name = factor(Global$Ocean_Name),
                      Ecoregion_Name = factor(Global$Ecoregion_Name),
                      Habitat_Type = factor(Global$Habitat_Type),
                      Exposure = Global$Exposure,
                      # Artifacts
                      Sample_Size = Global$N_samples,
                      Analysis = Global$Analysis_Type
)
names(X_Global)
Global_Stk <- inla.stack(
  tag  = "Est",
  data = list(y = Global$r),  
  A    = list(A2_global, 1),                      
  effects = list(                 
    w = w.index.globe,            #Spatial field  
    X = as.data.frame(X_Global)))  #Covariates


# Var.plot function -------------------------------------------------------

Var.plot = function(dataframe) {
  
  Var = NA
  mu = NA
  ci = NA
  df = data.frame(Var, mu, ci)
  df = df[0,]
  
  
  
  for (i in names(dataframe$marginals.fixed)) {
    Var = as.character(i)
    mu = inla.emarginal(function(x) x, dataframe$marginals.fixed[[as.character(i)]])             #posterior mean for fixed
    
    # 95% credibility interval
    cred_in95 = inla.hpdmarginal(0.95, dataframe$marginals.fixed[[as.character(i)]])   
    cred_in95[2]
    ci95 = abs(cred_in95[1] - mu)
    
    # 90% credibility interval
    cred_in90 = inla.hpdmarginal(0.90, dataframe$marginals.fixed[[as.character(i)]])   
    cred_in90[2]
    ci90 = abs(cred_in90[1] - mu)
    
    # 80% credibility
    cred_in80 = inla.hpdmarginal(0.80, dataframe$marginals.fixed[[as.character(i)]])   
    cred_in80[2]
    ci80 = abs(cred_in80[1] - mu)
    
    # credible intervarls for fixed
    df1 = data.frame(Var, mu, ci95, ci90, ci80)
    df = rbind(df, df1)
  }
  
  
  #Marking if a row is significant
  df$sig <- NA
  for (i1 in 1:nrow(df)) {
    low95 <- sign(df$mu[i1] - df$ci95[i1]) # Checks if it is positive or negative
    low90 <- sign(df$mu[i1] - df$ci90[i1]) # Checks if it is positive or negative
    low80 <- sign(df$mu[i1] - df$ci80[i1]) # Checks if it is positive or negative
    
    high95 <- sign(df$mu[i1] + df$ci95[i1]) # Checks if it is positive or negative
    high90 <- sign(df$mu[i1] + df$ci90[i1]) # Checks if it is positive or negative
    high80 <- sign(df$mu[i1] + df$ci80[i1]) # Checks if it is positive or negative
    
    df$sig[i1] <- sum(low95, low90, low80, high95, high90, high80)
    
    #-6 = 95% credible negative
    #-4 = 90% credible negative
    #-2 = 80% credible negative
    #0 = no credible
    #2 = 80 credible positive
    #4 = 90% credible positive
    #6 = 95% credible positive
    
    
  }
  cols <- c("-6" = "#EF1414", "-4" = "#EF1414", "-2" = "#9A5555", "0" = "#6F7070", "2" = "#457DB2", "4" = "#1485EF", "6" = "#1485EF")
  y <- ggplot(df, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
    geom_hline(yintercept = 0) + 
    #geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
    geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
    geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
    scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
    coord_flip() +
    xlab("") +
    ylab("Coefficient") +
    theme(legend.position = "none") +
    ggtitle("Coefficient Plot ")
  
  return(y)
}


# INLA.function only on global data -------------------------------

INLA.function = function(Variables) {
  prec.prior = list(prec = list(param = c(0.001, 0.001)))
  ## Global
  f0 = as.formula(paste0("y ~ -1 + ", paste0(Variables, collapse=" + ")))
  I_global = inla(f0,
                  family = "log normal",
                  data = inla.stack.data(Global_Stk),
                  control.compute = list(waic=TRUE),
                  control.predictor = list(A = inla.stack.A(Global_Stk), compute=TRUE))
  summary(I_global)
  print(summary(I_global)$waic$waic)
  
  Var.plot(dataframe = I_global)
  #return(I_global)
}


# Running Models ----------------------------------------------------------


names(X_Global)

All_Variables = c("f(Site_ID, model='iid')", # Random effect of site
                  "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                  "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                  "f(w, model=spde.globe)",  # Random effect of location
                  "Initial_Cover", "Initial_Algae", 
                  "Turbidity", 
                  "Depth", 
                  "Grav_tot", "Population", 
                  "Distance_to_shore", 
                  "Reef_Density_100km", "Reef_Density_10km", "Reef_Density_1km", 
                  "Velocity", "Kurtosis", "Skewness", "Hist_Max_SST", 
                  "Cyclone_Freq", "Heatwave_4_Freq", "Heatwave_8_Freq", "Heatwave_4_8_Freq",
                  "Heatwaves_4_Number", "Heatwaves_8_Number", "Heatwaves_4_8_Number",  
                  "Last_Heatwave_Intensity", "Years_Since_HW", 
                  "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR", "SD_DHW_DR",
                  "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR", "Number_of_Disturbances_DR",
                  "Ocean_Name", "Ecoregion_Name", "Habitat_Type", "Exposure")





INLA.function(Variables = c("f(Site_ID, model='iid')", # Random effect of site
                            "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                            "f(Time, model='rw1', hyper=prec.prior)"  # Random effect of time
))
# Running model only using random effects
INLA.function(Variables = c("f(Site_ID, model='iid')", # Random effect of site
                            "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                            "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                            "f(w, model=spde.globe)"))
# WAIC = -1290.126


# Adding in initial cover and initial macroalagae cover
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",
                "Initial_Cover", "Initial_Algae")) 
# WAIC = -1809.653
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "Initial_Cover", "Initial_Algae")) 
# WAIC = -1702.512





# Adding in disturbance events during recovery
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",
                "Initial_Cover", "Initial_Algae",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR"
)) 



# 
# > INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
#                   +                 "f(Ocean_Name, model='iid')",  # Random effect of Ocean
#                   +                 "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
#                   +                 "f(w, model=spde.globe)",  # Random effect of location
#                   +                 "Initial_Cover", "Initial_Algae", 
#                   +                 "Turbidity", 
#                   +                 "Depth", 
#                   +                 "Grav_tot", "Population", 
#                   +                 "Distance_to_shore", 
#                   +                 "Reef_Density_100km", "Reef_Density_10km", "Reef_Density_1km", 
#                   +                 "Velocity", "Kurtosis", "Skewness", "Hist_Max_SST", 
#                   +                 "Cyclone_Freq", "Heatwave_4_Freq", "Heatwave_8_Freq", "Heatwave_4_8_Freq",
#                   +                 "Heatwaves_4_Number", "Heatwaves_8_Number", "Heatwaves_4_8_Number",  
#                   +                 "Last_Heatwave_Intensity", "Years_Since_HW", 
#                   +                 "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR", "SD_DHW_DR",
#                   +                 "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR", "Number_of_Disturbances_DR"))
# [1] -2008.811






INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover", "Initial_Algae", 
                "Turbidity", 
                "Depth", 
                "Grav_tot", "Population", 
                "Distance_to_shore", 
                "Reef_Density_100km", "Reef_Density_10km", "Reef_Density_1km", 
                "Velocity", "Kurtosis", "Skewness", "Hist_Max_SST", 
                "Cyclone_Freq", 
                
                "Last_Heatwave_Intensity", "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))


# -1948.911




INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover", "Initial_Algae", 
                "Turbidity", 
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
#WAIC -1953.39



# Adding exposure as a random effect. Leads to reduced WAIC
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(Exposure, model='iid')",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
#-1949.47


# Adding exposure as a fixed effect
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))

# WAIC =-2037.139 but effects are all very negative. Does make other effects clearer though. Should probably not include in main plot and add to supps


# Adding habitat as fixed effect
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# -2041.635 but all categories are negative. Doesn't really enhance much


# Adding habitat as random effect
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "f(Habitat_Type, model='iid')",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
#WAIC is -2030.327


INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_4_Number_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
#-2043.439

# Swapping Number_DHW_4_DR for max_dhw_dr
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Years_Since_HW", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# -2092.824



# Removing time since last heatwave
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))

# WAIC = -2113.072


# Removing Temperature max and number of heatwaves
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Cyclone_Number_DR", "Mean_Windspeed_DR"))

# WAIC = -2108.531


# Removing algae coral cover interaction
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                
                "Initial_Cover","Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC = -2090.135. Good to keep interaction


# adding heatwave_4_8_freq
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Heatwave_4_8_Freq",
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Cyclone_Number_DR", "Mean_Windspeed_DR"))
#WAIC = -2093.503


# adding population, removing DHW 4 8 freq
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Population",
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Cyclone_Number_DR", "Mean_Windspeed_DR"))

# WAIC is -2095.879




# Going back to model with # WAIC = -2113.072 and removing Heatwaves_8_number_dr 
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC = -2107.08

# Adding Heatwaves_8_Number_DR back in and adding interaction between turbidity and depth
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                
                "Initial_Cover*Initial_Algae", 
                
                "Turbidity*Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC = -2097.364


# Removing interaction between turbidity and depth, adding in interaction between turbidity and max dhw 
# prior model was WAIC = -2113.072
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                
                "Initial_Cover*Initial_Algae", 
                
                "Turbidity*Max_DHW_DR",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC = -2093.396


# Removing interaction between turbidity and max dhw, adding interaction between algae and turbidity
# prior model was WAIC = -2113.072
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                
                "Initial_Cover*Initial_Algae", 
                
                "Initial_Algae*Turbidity", 
                "Max_DHW_DR",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC = -2093.3




INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))


# 10.4.23 -----------------------------------------------------------------

INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Max_DHW_DR*Mean_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))

# WAIC = -2113.072



INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Low_DHW_Freq", "Max_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC = -1994.583




INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Low_DHW_Freq", "High_DHW_Freq", "Max_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))



# Models without spatial effect -------------------------------------------

INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Max_DHW_DR*Mean_DHW_DR",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))
# WAIC is -1957

# Adding in DHW freq
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Max_DHW_DR", "Low_DHW_Freq", "High_DHW_Freq",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"
))

#-1897.787

# adding analysis type
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(Analysis, model='iid')",  # Random effect of Analysis type
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Max_DHW_DR", "Low_DHW_Freq", "High_DHW_Freq",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"
))
#-2005.402

# Adding sample size and DHW Freq Prior
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(Analysis, model='iid')",  # Random effect of Analysis type
                "f(Sample_Size, model='iid')",  # Random effect of sample size
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_Max", "Temperature_SD", "Max_DHW_DR", "DHW_Freq_Prior",
                "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"
))
#-2156.829


# Adding 
INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(w, model=spde.globe)",  # Random effect of location
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(Analysis, model='iid')",  # Random effect of Analysis type
                "f(Sample_Size, model='iid')",  # Random effect of sample size
                #
                #"factor(Exposure)",
                "factor(Habitat_Type)",
                #
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", 
                "Kurtosis", 
                "Skewness", 
                "Cyclone_Freq", 
                "Temperature_SD", 
                "Max_DHW_DR", 
                "Mean_Windspeed_DR"
))
#-2281.087

INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "factor(Exposure)",
                "factor(Habitat_Type)",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Temperature_SD", "DHW_Freq_Prior", "Max_DHW_DR",
                "Mean_Windspeed_DR"))

## WAIC = -2006.783



INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                #"f(Exposure, model='iid')",
                "factor(Habitat_Type)",
                "f(Analysis, model='iid')",
                "f(Sample_Size, model='iid')",
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(w, model=spde.globe)",  # Random effect of location
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Temperature_SD", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", "Kurtosis", "Skewness", 
                "Cyclone_Freq", 
                "Heatwaves_8_Number", "Max_DHW_DR",
                "Mean_Windspeed_DR"))

# WAIC -2291.889

# Models 10.19.23 ---------------------------------------------------------

INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                "f(w, model=spde.globe)",  # Random effect of location
                "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                "f(Analysis, model='iid')",  # Random effect of Analysis type
                "f(Sample_Size, model='iid')",  # Random effect of sample size
                #
                #"factor(Exposure)",
                "factor(Habitat_Type)",
                #
                "Initial_Cover*Initial_Algae", 
                "Turbidity",
                "Depth", 
                "Grav_tot",
                "Distance_to_shore", 
                "Reef_Density_10km",
                "Velocity", 
                "Kurtosis", 
                "Skewness", 
                "Cyclone_Freq", 
                "DHW_Freq_Prior",
                #"Temperature_SD", 
                "Max_DHW_DR", 
                "Max_Windspeed_DR"
))

# WAIC = -2300.299

# Best Model --------------------------------------------------------------




INLA.function = function(Variables) {
  prec.prior = list(prec = list(param = c(0.001, 0.001)))
  ## Global
  f0 = as.formula(paste0("y ~ -1 + ", paste0(Variables, collapse=" + ")))
  I_global = inla(f0,
                  family = "log normal",
                  data = inla.stack.data(Global_Stk),
                  control.compute = list(waic=TRUE),
                  control.predictor = list(A = inla.stack.A(Global_Stk), compute=TRUE))
  summary(I_global)
  print(summary(I_global)$waic$waic)
  
  #Var.plot(dataframe = I_global)
  return(I_global)
}

#WAIC = -2113.072
I_best <- INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                          "factor(Exposure)",
                          "factor(Habitat_Type)",
                          "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                          "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                          "f(w, model=spde.globe)",  # Random effect of location
                          "Initial_Cover*Initial_Algae", 
                          "Turbidity",
                          "Depth", 
                          "Grav_tot",
                          "Distance_to_shore", 
                          "Reef_Density_10km",
                          "Velocity", "Kurtosis", "Skewness", 
                          "Cyclone_Freq", 
                          "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR",
                          "Heatwaves_8_Number_DR", "Cyclone_Number_DR", "Mean_Windspeed_DR"))








Var = NA
mu = NA
ci = NA
df = data.frame(Var, mu, ci)
df = df[0,]

for (i in names(I_best$marginals.fixed)) {
  Var = as.character(i)
  mu = inla.emarginal(function(x) x, I_best$marginals.fixed[[as.character(i)]])             #posterior mean for fixed
  
  # 95% credibility interval
  cred_in95 = inla.hpdmarginal(0.95, I_best$marginals.fixed[[as.character(i)]])   
  cred_in95[2]
  ci95 = abs(cred_in95[1] - mu)
  
  # 90% credibility interval
  cred_in90 = inla.hpdmarginal(0.90, I_best$marginals.fixed[[as.character(i)]])   
  cred_in90[2]
  ci90 = abs(cred_in90[1] - mu)
  
  # 80% credibility
  cred_in80 = inla.hpdmarginal(0.80, I_best$marginals.fixed[[as.character(i)]])   
  cred_in80[2]
  ci80 = abs(cred_in80[1] - mu)
  
  # credible intervarls for fixed
  df1 = data.frame(Var, mu, ci95, ci90, ci80)
  df = rbind(df, df1)
}
#Marking if a row is significant
df$sig <- NA
for (i1 in 1:nrow(df)) {
  low95 <- sign(df$mu[i1] - df$ci95[i1]) # Checks if it is positive or negative
  low90 <- sign(df$mu[i1] - df$ci90[i1]) # Checks if it is positive or negative
  low80 <- sign(df$mu[i1] - df$ci80[i1]) # Checks if it is positive or negative
  
  high95 <- sign(df$mu[i1] + df$ci95[i1]) # Checks if it is positive or negative
  high90 <- sign(df$mu[i1] + df$ci90[i1]) # Checks if it is positive or negative
  high80 <- sign(df$mu[i1] + df$ci80[i1]) # Checks if it is positive or negative
  
  df$sig[i1] <- sum(low95, low90, low80, high95, high90, high80)
  
  #-6 = 95% credible negative
  #-4 = 90% credible negative
  #-2 = 80% credible negative
  #0 = no credible
  #2 = 80 credible positive
  #4 = 90% credible positive
  #6 = 95% credible positive
  
  
}

# Establishing Plot
cols <- c("-6" = "#EF1414", "-4" = "#EF1414", "-2" = "#9A5555", "0" = "#6F7070", "2" = "#457DB2", "4" = "#1485EF", "6" = "#1485EF") # Colors


# Subsetting Habitat Type and Exposure
unique(df$Var)
df_factors <- df[which(df$Var %in% c("factor(Exposure)", "factor(Exposure)Exposed", "factor(Exposure)Sheltered",
                                     "factor(Exposure)Sometimes", "factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]
df_cont <- df[which(df$Var %in% c("Initial_Cover", "Initial_Algae", "Turbidity", "Depth", "Grav_tot", "Distance_to_shore", "Reef_Density_10km", "Velocity", "Kurtosis", "Skewness",
                                  "Cyclone_Freq", "Temperature_Max", "Temperature_SD", "Mean_DHW_DR", "Max_DHW_DR", "Heatwaves_8_Number_DR", "Cyclone_Number_DR",
                                  "Mean_Windspeed_DR", "Initial_Cover:Initial_Algae")),]


df_cont$Var
df_cont$Var <- c("Initial Coral Cover", "Initial Macroalgae", "Turbidity", "Depth", "Total Gravity", "Distance to Shore", 
                 "Reef Density (10km)", "SST Velocity", "SST Kurtosis", "SST Skewness", "Prior Cyclone Frequency",
                 "Max SST", "SST SD", "Mean DHW During Recovery", "Max DHW During Recovery", "N Heatwaves During Recovery", "N Cyclones During Recovery", "Mean Windspeed", "Initial Coral Cover:Initial Macroalgae")
ggplot(df_cont, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0) + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab("Coefficient") +
  theme(legend.position = "none") +
  ggtitle("Coefficient Plot")



ggplot(df_factors, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0) + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab("Coefficient") +
  theme(legend.position = "none") +
  ggtitle("Coefficient Plot")




# Best Model 10.6.23 --------------------------------------------------------------




INLA.function = function(Variables) {
  prec.prior = list(prec = list(param = c(0.001, 0.001)))
  ## Global
  f0 = as.formula(paste0("y ~ -1 + ", paste0(Variables, collapse=" + ")))
  I_global = inla(f0,
                  family = "log normal",
                  data = inla.stack.data(Global_Stk),
                  control.compute = list(waic=TRUE),
                  control.predictor = list(A = inla.stack.A(Global_Stk), compute=TRUE))
  summary(I_global)
  print(summary(I_global)$waic$waic)
  
  #Var.plot(dataframe = I_global)
  return(I_global)
}

#WAIC = -2113.072
I_best <- INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                          "f(w, model=spde.globe)",  # Random effect of location
                          "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                          "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                          "f(Analysis, model='iid')",  # Random effect of Analysis type
                          "f(Sample_Size, model='iid')",  # Random effect of sample size
                          #
                          #"factor(Exposure)",
                          "factor(Habitat_Type)",
                          #
                          "Initial_Cover*Initial_Algae", 
                          "Turbidity",
                          "Depth", 
                          "Grav_tot",
                          "Distance_to_shore", 
                          "Reef_Density_10km",
                          "Velocity", 
                          "Kurtosis", 
                          "Skewness", 
                          "Cyclone_Freq", 
                          "Temperature_SD", 
                          "Max_DHW_DR", 
                          "Mean_Windspeed_DR"
))
#-2281.087







Var = NA
mu = NA
ci = NA
df = data.frame(Var, mu, ci)
df = df[0,]

for (i in names(I_best$marginals.fixed)) {
  Var = as.character(i)
  mu = inla.emarginal(function(x) x, I_best$marginals.fixed[[as.character(i)]])             #posterior mean for fixed
  
  # 95% credibility interval
  cred_in95 = inla.hpdmarginal(0.95, I_best$marginals.fixed[[as.character(i)]])   
  cred_in95[2]
  ci95 = abs(cred_in95[1] - mu)
  
  # 90% credibility interval
  cred_in90 = inla.hpdmarginal(0.90, I_best$marginals.fixed[[as.character(i)]])   
  cred_in90[2]
  ci90 = abs(cred_in90[1] - mu)
  
  # 80% credibility
  cred_in80 = inla.hpdmarginal(0.80, I_best$marginals.fixed[[as.character(i)]])   
  cred_in80[2]
  ci80 = abs(cred_in80[1] - mu)
  
  # credible intervarls for fixed
  df1 = data.frame(Var, mu, ci95, ci90, ci80)
  df = rbind(df, df1)
}
#Marking if a row is significant
df$sig <- NA
for (i1 in 1:nrow(df)) {
  low95 <- sign(df$mu[i1] - df$ci95[i1]) # Checks if it is positive or negative
  low90 <- sign(df$mu[i1] - df$ci90[i1]) # Checks if it is positive or negative
  low80 <- sign(df$mu[i1] - df$ci80[i1]) # Checks if it is positive or negative
  
  high95 <- sign(df$mu[i1] + df$ci95[i1]) # Checks if it is positive or negative
  high90 <- sign(df$mu[i1] + df$ci90[i1]) # Checks if it is positive or negative
  high80 <- sign(df$mu[i1] + df$ci80[i1]) # Checks if it is positive or negative
  
  df$sig[i1] <- sum(low95, low90, low80, high95, high90, high80)
  
  #-6 = 95% credible negative
  #-4 = 90% credible negative
  #-2 = 80% credible negative
  #0 = no credible
  #2 = 80 credible positive
  #4 = 90% credible positive
  #6 = 95% credible positive
  
  
}

# Establishing Plot
cols <- c("-6" = "#EF1414", "-4" = "#EF1414", "-2" = "#9A5555", "0" = "#6F7070", "2" = "#457DB2", "4" = "#1485EF", "6" = "#1485EF") # Colors


`%notin%` <- Negate(`%in%`)

# Subsetting Habitat Type and Exposure
unique(df$Var)
df_factors <- df[which(df$Var %in% c("factor(Habitat_Type)", "factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]
df_cont <- df[which(df$Var %notin% c("factor(Habitat_Type)","factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]


df_cont$Var
df_cont$Var <- c("Initial Coral Cover", "Initial Macroalgae", "Turbidity", "Depth", "Total Gravity", "Distance to Shore", 
                 "Reef Density (10km)", "SST Velocity", "SST Kurtosis", "SST Skewness", "Prior Cyclone Frequency",
                 "SST SD", "Max DHW During Recovery",  "Mean Cyclone Windspeed", "Initial Coral Cover:Initial Macroalgae")
ggplot(df_cont, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0) + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab(substitute(paste(italic('β'), " coefficient"))) +
  theme(legend.position = "none") +
  ggtitle("")



ggplot(df_factors, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0) + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab(substitute(paste(italic('β'), " coefficient"))) +
  theme(legend.position = "none") +
  ggtitle("")



# Best Model 10.19.23 -----------------------------------------------------




INLA.function = function(Variables) {
  prec.prior = list(prec = list(param = c(0.001, 0.001)))
  ## Global
  f0 = as.formula(paste0("y ~ -1 + ", paste0(Variables, collapse=" + ")))
  I_global = inla(f0,
                  family = "log normal",
                  data = inla.stack.data(Global_Stk),
                  control.compute = list(waic=TRUE),
                  control.predictor = list(A = inla.stack.A(Global_Stk), compute=TRUE))
  summary(I_global)
  print(summary(I_global)$waic$waic)
  
  #Var.plot(dataframe = I_global)
  return(I_global)
}

#WAIC = -2113.072
I_best <- INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                          "f(w, model=spde.globe)",  # Random effect of location
                          "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                          "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                          "f(Analysis, model='iid')",  # Random effect of Analysis type
                          "f(Sample_Size, model='iid')",  # Random effect of sample size
                          #
                          #"factor(Exposure)",
                          "factor(Habitat_Type)",
                          #
                          "Initial_Cover*Initial_Algae", 
                          "Turbidity",
                          "Depth", 
                          "Grav_tot",
                          "Distance_to_shore", 
                          "Reef_Density_10km",
                          "Velocity", 
                          "Kurtosis", 
                          "Skewness", 
                          "Cyclone_Freq", 
                          "DHW_Freq_Prior",
                          #"Temperature_SD", 
                          "Max_DHW_DR", 
                          "Max_Windspeed_DR"
))
#-2298.596







Var = NA
mu = NA
ci = NA
df = data.frame(Var, mu, ci)
df = df[0,]

for (i in names(I_best$marginals.fixed)) {
  Var = as.character(i)
  mu = inla.emarginal(function(x) x, I_best$marginals.fixed[[as.character(i)]])             #posterior mean for fixed
  
  # 95% credibility interval
  cred_in95 = inla.hpdmarginal(0.95, I_best$marginals.fixed[[as.character(i)]])   
  cred_in95[2]
  ci95 = abs(cred_in95[1] - mu)
  
  # 90% credibility interval
  cred_in90 = inla.hpdmarginal(0.90, I_best$marginals.fixed[[as.character(i)]])   
  cred_in90[2]
  ci90 = abs(cred_in90[1] - mu)
  
  # 80% credibility
  cred_in80 = inla.hpdmarginal(0.80, I_best$marginals.fixed[[as.character(i)]])   
  cred_in80[2]
  ci80 = abs(cred_in80[1] - mu)
  
  # credible intervarls for fixed
  df1 = data.frame(Var, mu, ci95, ci90, ci80)
  df = rbind(df, df1)
}
#Marking if a row is significant
df$sig <- NA
for (i1 in 1:nrow(df)) {
  low95 <- sign(df$mu[i1] - df$ci95[i1]) # Checks if it is positive or negative
  low90 <- sign(df$mu[i1] - df$ci90[i1]) # Checks if it is positive or negative
  low80 <- sign(df$mu[i1] - df$ci80[i1]) # Checks if it is positive or negative
  
  high95 <- sign(df$mu[i1] + df$ci95[i1]) # Checks if it is positive or negative
  high90 <- sign(df$mu[i1] + df$ci90[i1]) # Checks if it is positive or negative
  high80 <- sign(df$mu[i1] + df$ci80[i1]) # Checks if it is positive or negative
  
  df$sig[i1] <- sum(low95, low90, low80, high95, high90, high80)
  
  #-6 = 95% credible negative
  #-4 = 90% credible negative
  #-2 = 80% credible negative
  #0 = no credible
  #2 = 80 credible positive
  #4 = 90% credible positive
  #6 = 95% credible positive
  
  
}

# Establishing Plot
cols <- c("-6" = "#EF1414", "-4" = "#EF1414", "-2" = "#9A5555", "0" = "#6F7070", "2" = "#457DB2", "4" = "#1485EF", "6" = "#1485EF") # Colors


`%notin%` <- Negate(`%in%`)

# Subsetting Habitat Type and Exposure
unique(df$Var)
df_factors <- df[which(df$Var %in% c("factor(Habitat_Type)", "factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]

df_cont <- df[which(df$Var %notin% c("factor(Habitat_Type)","factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]


df_cont$Var
df_cont$Var <- c("Initial Coral Cover", "Initial Macroalgae", "Turbidity", "Depth", "Total Gravity", "Distance to Shore", 
                 "Reef Density (10km)", "SST Velocity", "SST Kurtosis", "SST Skewness", "Prior Cyclone Frequency", "Prior Heatwave Frequency",
                 "Max DHW During Recovery",  "Max Cyclone Windspeed", "Initial Coral Cover:Initial Macroalgae")

df_factors$Var <- c("Habitat NA", "Back Reef Slope", "Deep Lagoon", "Inner Reef Flat", 
                    "Outer Reef Flat", "Plateau", "Reef Crest", "Reef Slope", "Shallow Lagoon",
                    "Sheltered Reef Slope", "Terrestrial Reef Flat")


coef_cont <- ggplot(df_cont, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0, size = 1.5, linetype = "dotdash", color = "#404040") + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2.5) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 3.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab(substitute(paste(italic('β'), " coefficient"))) +
  theme(legend.position = "none") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.4, "cm"),

    panel.border = element_rect(colour = "black", fill=NA, size=1),
  )



ggsave("C:/Users/walkera2019/Dropbox/Recovery Gamma/Manuscript Review 1/INLA Plot/INLA Continuous.png", coef_cont, 
       width = 8, height = 8)



coef_hab <- ggplot(df_factors, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0, size = 1.5, linetype = "dotdash", color = "#404040") + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab(substitute(paste(italic('β'), " coefficient"))) +
  theme(legend.position = "none") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.4, "cm"),
    
    panel.border = element_rect(colour = "black", fill=NA, size=1),
  )

ggsave("C:/Users/walkera2019/Dropbox/Recovery Gamma/Manuscript Review 1/INLA Plot/INLA Habitat.png", coef_hab, 
       width = 8, height = 8)



# Best Model 10.24.23 -----------------------------------------------------




INLA.function = function(Variables) {
  prec.prior = list(prec = list(param = c(0.001, 0.001)))
  ## Global
  f0 = as.formula(paste0("y ~ -1 + ", paste0(Variables, collapse=" + ")))
  I_global = inla(f0,
                  family = "log normal",
                  data = inla.stack.data(Global_Stk),
                  control.compute = list(waic=TRUE),
                  control.predictor = list(A = inla.stack.A(Global_Stk), compute=TRUE))
  summary(I_global)
  print(summary(I_global)$waic$waic)
  
  #Var.plot(dataframe = I_global)
  return(I_global)
}

#WAIC = -2113.072
I_best <- INLA.function(c("f(Site_ID, model='iid')", # Random effect of site
                          "f(w, model=spde.globe)",  # Random effect of location
                          "f(Ocean_Name, model='iid')",  # Random effect of Ocean
                          "f(Time, model='rw1', hyper=prec.prior)",  # Random effect of time
                          "f(Analysis, model='iid')",  # Random effect of Analysis type
                          "f(Sample_Size, model='iid')",  # Random effect of sample size
                          #
                          "factor(Habitat_Type)",
                          #
                          "Initial_Cover*Initial_Algae", 
                          "Turbidity",
                          "Depth", 
                          "Grav_tot",
                          "Distance_to_shore", 
                          "Reef_Density_100km",
                          "Velocity", 
                          "Kurtosis", 
                          "Skewness", 
                          "Cyclone_Freq", 
                          "DHW_Freq_Prior",
                          "Max_DHW_DR", 
                          "Max_Windspeed_DR"
))
#-2298.596







Var = NA
mu = NA
ci = NA
df = data.frame(Var, mu, ci)
df = df[0,]

for (i in names(I_best$marginals.fixed)) {
  Var = as.character(i)
  mu = inla.emarginal(function(x) x, I_best$marginals.fixed[[as.character(i)]])             #posterior mean for fixed
  
  # 95% credibility interval
  cred_in95 = inla.hpdmarginal(0.95, I_best$marginals.fixed[[as.character(i)]])   
  cred_in95[2]
  ci95 = abs(cred_in95[1] - mu)
  
  # 90% credibility interval
  cred_in90 = inla.hpdmarginal(0.90, I_best$marginals.fixed[[as.character(i)]])   
  cred_in90[2]
  ci90 = abs(cred_in90[1] - mu)
  
  # 80% credibility
  cred_in80 = inla.hpdmarginal(0.80, I_best$marginals.fixed[[as.character(i)]])   
  cred_in80[2]
  ci80 = abs(cred_in80[1] - mu)
  
  # credible intervarls for fixed
  df1 = data.frame(Var, mu, ci95, ci90, ci80)
  df = rbind(df, df1)
}
#Marking if a row is significant
df$sig <- NA
for (i1 in 1:nrow(df)) {
  low95 <- sign(df$mu[i1] - df$ci95[i1]) # Checks if it is positive or negative
  low90 <- sign(df$mu[i1] - df$ci90[i1]) # Checks if it is positive or negative
  low80 <- sign(df$mu[i1] - df$ci80[i1]) # Checks if it is positive or negative
  
  high95 <- sign(df$mu[i1] + df$ci95[i1]) # Checks if it is positive or negative
  high90 <- sign(df$mu[i1] + df$ci90[i1]) # Checks if it is positive or negative
  high80 <- sign(df$mu[i1] + df$ci80[i1]) # Checks if it is positive or negative
  
  df$sig[i1] <- sum(low95, low90, low80, high95, high90, high80)
  
  #-6 = 95% credible negative
  #-4 = 90% credible negative
  #-2 = 80% credible negative
  #0 = no credible
  #2 = 80 credible positive
  #4 = 90% credible positive
  #6 = 95% credible positive
  
  
}

#write.csv(x = df, file = "INLA Model Output 12.1.23.csv")


# Establishing Plot
cols <- c("-6" = "#EF1414", "-4" = "#EF1414", "-2" = "#9A5555", "0" = "#6F7070", "2" = "#457DB2", "4" = "#1485EF", "6" = "#1485EF") # Colors


`%notin%` <- Negate(`%in%`)

# Subsetting Habitat Type and Exposure
unique(df$Var)
df_factors <- df[which(df$Var %in% c("factor(Habitat_Type)", "factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]

df_cont <- df[which(df$Var %notin% c("factor(Habitat_Type)","factor(Habitat_Type)Back Reef Slope", "factor(Habitat_Type)Deep Lagoon",
                                     "factor(Habitat_Type)Inner Reef Flat", "factor(Habitat_Type)Outer Reef Flat", "factor(Habitat_Type)Plateau",
                                     "factor(Habitat_Type)Reef Crest", "factor(Habitat_Type)Reef Slope", "factor(Habitat_Type)Shallow Lagoon", 
                                     "factor(Habitat_Type)Sheltered Reef Slope", "factor(Habitat_Type)Terrestrial Reef Flat")),]


df_cont$Var
df_cont$Var <- c("Initial Coral Cover", "Initial Macroalgae", "Turbidity", "Depth", "Total Gravity", "Distance to Shore", 
                 "Reef Density (100km)", "SST Velocity", "SST Kurtosis", "SST Skewness", "Prior Cyclone Frequency", "Prior Heatwave Frequency",
                 "Max DHW During Recovery",  "Max Windspeed During Recovery", "Initial Coral Cover:Initial Macroalgae")

df_factors$Var <- c("Habitat NA", "Back Reef Slope", "Deep Lagoon", "Inner Reef Flat", 
                    "Outer Reef Flat", "Plateau", "Reef Crest", "Reef Slope", "Shallow Lagoon",
                    "Sheltered Reef Slope", "Terrestrial Reef Flat")


coef_cont <- ggplot(df_cont, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0, size = 1.5, linetype = "dotdash", color = "#404040") + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2.5) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 3.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab(substitute(paste(italic('β'), " coefficient"))) +
  theme(legend.position = "none") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.4, "cm"),
    
    panel.border = element_rect(colour = "black", fill=NA, size=1),
  )



ggsave("C:/Users/walkera2019/Dropbox/Recovery Gamma/Manuscript Review 1/INLA Plot/INLA Continuous 12.1.23.pdf", coef_cont, 
       width = 8, height = 8, dpi = 300)



coef_hab <- ggplot(df_factors, aes(x = reorder(Var, mu), y = mu, color = factor(sig))) +
  geom_hline(yintercept = 0, size = 1.5, linetype = "dotdash", color = "#404040") + 
  #geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = mu - ci95, ymax = mu + ci95), width = 0, size = 1) +
  geom_errorbar(aes(ymin = mu - ci90, ymax = mu + ci90), width = 0, size = 2) +
  geom_errorbar(aes(ymin = mu - ci80, ymax = mu + ci80), width = 0, size = 2.5) +
  scale_colour_manual(values = cols, breaks = c("-6", "-4", "-2", "0", "2", "4", "6")) +
  coord_flip() +
  xlab("") +
  ylab(substitute(paste(italic('β'), " coefficient"))) +
  theme(legend.position = "none") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.4, "cm"),
    
    panel.border = element_rect(colour = "black", fill=NA, size=1),
  )

ggsave("C:/Users/walkera2019/Dropbox/Recovery Gamma/Manuscript Review 1/INLA Plot/INLA Habitat 10.24.23.png", coef_hab, 
       width = 8, height = 8)


