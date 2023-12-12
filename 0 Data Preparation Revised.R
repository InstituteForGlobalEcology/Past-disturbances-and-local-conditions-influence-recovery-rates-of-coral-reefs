
### Library -----------------------------------------------------------------
library(minpack.lm)
library(dplyr)
library(rjags)
library(coda)
library(drc) # used for gompertz
### Loading data in ---------------------------------------------------------

#What data do I need?
##Coral time series data
##Diff.Year function, but consider day and month as well to get more exact time


Data_Directory <- "C:/Users/Andy Walker/Dropbox/Data Directory"
setwd(Data_Directory)

data_sites <- read.csv("All_Cover_All_Variables_3_28_23.csv")
data <- read.csv("All_Cover_All_Variables_9_29_23_macroalgae.csv")
data <- data[which(data$Site_ID %in% data_sites$Site_ID),]

remove(data_sites)
##  Data Cleaning
#Removing NA values for X.Hard_Coral_Cover and Year
data <- subset(data, X.Hard_Coral_Cover != "NA") # 6 rows removed
data <- subset(data, Date_Year != "NA") # No rows removed

data <- subset(data, select = c("Site_ID", "Sample_ID", "Latitude_Degrees", "Longitude_Degrees", "Ocean_Name", "Habitat_Type", 
                                "Habitat_Distance", "Date_Day", "Date_Month", "Date_Year", "Exposure", "Turbidity", "Grav_tot", "Depth", "Distance_to_shore", "X.Hard_Coral_Cover", "X.Macroalgae_Cover", "Ecoregion_Name", "State_Island_Province_Name"))


Sites <- levels(factor(data$Site_ID))
min(data$Date_Year)
max(data$Date_Year)

# Full Date Tool ----------------------------------------------------------

data$Missing_Month <- F
data$Missing_Day <- F

for (i in 1:nrow(data)) {
  #specifying year month and day from row
  year <- data$Date_Year[i]
  month <- data$Date_Month[i]
  day <- data$Date_Day[i]
  if (is.na(day) == TRUE) { #If there is no day this block assigns a day Not an ideal solution but for now it works
    day <- 27
    data$Date_Day[i] <- 27 #setting arbitrary day for NA
    data$Missing_Day[i] <- T
  }
  if (is.na(month) == TRUE) {#If there is no month this block assigns a month Not an ideal solution but for now it works
    if (data$Latitude_Degrees[i] < 0 ) { # southern hemisphere
      data$Date_Month[i] <- 2 #setting arbitrary month for NA
      month <- 2
      data$Missing_Month[i] <- T
    }
    if (data$Latitude_Degrees[i] > 0) { # Northern hemisphere
      data$Date_Month[i] <- 8 #setting arbitrary month for NA
      month <- 8
      data$Missing_Month[i] <- T
    }
  }
  #creating date format
  date <- paste(year, "-", month, "-", day, sep = "") #takes the three date columns and makes one date out of them
  survey_date <- as.Date(date, format = "%Y-%m-%d") #converts date to date format for date maths
  data$Date_Full[i] <- date #creates column in data for full date for simpler future use
}

remove(day, month, year, i)

data$Date_Full <- as.Date(data$Date_Full, format = "%Y-%m-%d")

nrow(subset(data, Missing_Month == T)) #1246 rows missing month data
nrow(subset(data, Missing_Day == T)) #1386 rows missing month data

nlevels(factor(subset(data, Missing_Month == T)$Site_ID)) #300 sites
nlevels(factor(subset(data, Missing_Day == T)$Site_ID)) #354 sites

### Removing Sites that have no change in coral cover --------------------

blank <- data[0,]
for (x in Sites) {
  x1 <- subset(data, Site_ID == x)
  if (nlevels(factor(x1$X.Hard_Coral_Cover)) != 1) { # If there is only one value in whole site it is excluded
    blank <- rbind(blank, x1)
  }
}
data <- blank
remove(x1, blank)


data$Date_Full <- as.Date(data$Date_Full, format = "%Y-%m-%d")
data.class(data$Date_Full)
data$Date_Full <- format(data$Date_Full, format = "%Y-%m-%d")
### Collapsing Data into One Row per Date per Site -----------------------------------

site <- levels(factor(data$Site_ID))
blank <- data[0,]
for (i in site) {
  x <- subset(data, Site_ID == i)
  dates <- levels(factor(x$Date_Full))
  for (o in dates) {
    d <- as.Date(o, format = "%Y-%m-%d")
    if (identical(which(dates == d-1), integer(0)) == TRUE) { # If there wasn't a date the day prior
      y <- subset(x, (Date_Full == d) | (Date_Full == d + 1))
      y1 <- y[1,]
      y1$X.Hard_Coral_Cover <- mean(y$X.Hard_Coral_Cover)
      y1$Depth <- mean(y$Depth)
      blank <- rbind(blank, y1)
    }
  }
} #creates 18362 rows 
data <- blank
remove(blank, x, dates, site, y, y1, i, o)


## Filtering for sites with 2 or more points -----------------------------

# Filtering for sites with 3 or more points
blank <- data[0,]
for (i in levels(factor(data$Site_ID))) {
  x <- subset(data, Site_ID == i)
  if (nrow(x) >= 2) {
    blank <- rbind(blank, x)
  }
}
data <- blank
remove(blank)



### Difference in Time Between Surveys -----------------------------------

site <- levels(factor(data$Site_ID))
blank <- data[0,]
for (i in site) {
  x <- subset(data, Site_ID == i)
  start_date <- min(as.Date(x$Date_Full, format = "%Y-%m-%d"))
  for (o in 1:nrow(x)) {
    x$Date_Diff[o] <- round((difftime(as.Date(x$Date_Full[o], format = "%Y-%m-%d"), start_date, unit = "days"))/365, digits = 3)
  }
  x <- x[order(x$Date_Diff),] #ordering by date_diff
  blank <- rbind(blank, x)
}
data <- blank
remove(blank, i, o, start_date, x, site)




### Within Site Vector ------------------------------------------------------

data$index <- ave(data$Date_Diff, # Create numbering variable that counts up by
                  data$Site_ID,   #one integer each time a site is surveyed
                  FUN = seq_along)



# Loading in disturbances -------------------------------------------------

## Cyclone data has been extracted for sites using Cyclone Extractor WIP.R
#cyclones <- read.csv("C:/Users/Andy Walker/Dropbox/Disturbance History/Full_Date_Severe_Cyclones_100km.8.14.23.csv")
cyclones <- read.csv("C:/Users/Andy Walker/Dropbox/Disturbance History/Full_Date_Severe_Cyclones_200km.8.24.23.csv")

cyclones$Full_Date <- as.Date(cyclones$Full_Date, format = "%m/%d/%Y")
data.class(cyclones$Full_Date)
range(cyclones$m_per_s) # Making sure minimum is above 33
### NOTE
## Must filter cyclone data based on windspeed and distance from point
##


## DHW data
# Here I am loading in all DHW data to the same df so that in the loop later on I can subset for a site, and if there isn't data on it I can skip it without an error occurring
dhw_list <- list.files(path = "C:/Users/Andy Walker/Dropbox/Disturbance History/DHW 9.21.23/DHW by Site/", pattern="*.csv") # Grabbing all CSVs in the folder. Each file should only contain values of 4 or more
DHW_Sites <- NA
for (i in 1:length(dhw_list)) {
  x <- strsplit(dhw_list[i], "_")
  x <- x[[1]]
  x <- x[2]
  DHW_Sites <- c(DHW_Sites, x)
}
DHW_Sites <- na.omit(DHW_Sites)
remove(x, i, dhw_list)



# Calculating Frequency of Disturbances -----------------------------------

data$Cyclone_Freq <- NA
data$Heatwave_Freq <- NA

data$Dist_Type <- NA# Type of prior disturbance
data$Dist_Date <- NA # last disturance date
data$Dist_Intensity <- NA # For cyclone or heatwave that was prior
data$Last_Heatwave_Date <- NA  #Date of last heatwave
data$Last_Heatwave_Intensity <- NA # Intensity of prior heatwave
data$Heatwave_4_Freq <- NA # Frequency of heatwaves 4 DHW or greater
data$Heatwave_8_Freq <- NA # Frequency of heatwaves 8 DHW or greater
data$Heatwave_4_8_Freq <- NA # Frequency of heatwaves between 4 and 8 DHW

data$Heatwaves_4_Number <- NA # Number of heatwaves 4 DHW or greater
data$Heatwaves_8_Number <- NA # Number of heatwaves 8 DHW or greater
data$Heatwaves_4_8_Number <- NA # Number of heatwaves between 4 and 8 DHW


DHW_Cutoff <- 4 # DHW minimum value for loop

blank <- data[0,]
Max_Dist <- 1 # In arc degrees. Max distance for cyclones to occur at site
for (i in levels(factor(data$Site_ID))) {
  cover0 <- subset(data, Site_ID == i)
  cover0$Date_Full <- as.Date(cover0$Date_Full, format = "%Y-%m-%d")
  ##
  ## Loading in cyclones
  site_cyclone <- subset(cyclones, Site_ID == i)
  site_cyclone <- subset(site_cyclone, Distance <= Max_Dist)
  if (nrow(site_cyclone) >= 1) {
    site_cyclone$Date<- site_cyclone$Full_Date
    site_cyclone$event[1] <- 1
    if (nrow(site_cyclone) > 1) {
      for (j1 in 2:nrow(site_cyclone)) {
        if (site_cyclone$Date[j1-1]+1 == site_cyclone$Date[j1]) {
          site_cyclone$event[j1] <- site_cyclone$event[j1-1]
        }
        if (site_cyclone$Date[j1-1]+1 != site_cyclone$Date[j1]) {
          site_cyclone$event[j1] <- site_cyclone$event[j1-1]+1
        }
      } 
    }
    m_per_s <- NA
    Date <- NA
    Distance <- NA
    df0 <- data.frame(Date, m_per_s, Distance)
    df0 <- df0[0,]
    for (j2 in unique(site_cyclone$event)) {
      event1 <- subset(site_cyclone, event == j2)
      m_per_s <- max(event1$m_per_s)
      Distance <- min(event1$Distance)
      Date <- as.character(max(event1$Date[which(event1$m_per_s == max(event1$m_per_s))]))
      cyc_df <- data.frame(Date, m_per_s, Distance)
      df0 <- rbind(df0, cyc_df)
    }
    site_cyclone <- df0
    site_cyclone$Type <- "Cyclone"
    site_cyclone$Date <- as.Date(site_cyclone$Date, format = "%Y-%m-%d")
  }
  if (nrow(site_cyclone) == 0) {
    Date <- NA
    Type <- NA
    site_cyclone <- data.frame(Date, Type)
    site_cyclone <- site_cyclone[0,]
  }
  ##
  ## Loading in DHW events
  if (length(which(DHW_Sites == i)) == 1) { # If there is a file for DHW events, load it in
    site_DHW <- read.csv(paste("C:/Users/Andy Walker/Dropbox/Disturbance History/DHW 9.21.23/DHW by Site/Site_", i, "_DHW.9.21.23.csv", sep = ""))
    site_DHW$event <- NA
    site_DHW$Date <- as.Date(site_DHW$Date, format = "%Y-%m-%d")
    site_DHW <- subset(site_DHW, DHW >= DHW_Cutoff)
    if (nrow(site_DHW) >= 1) {
      site_DHW$event[1] <- 1
      if (nrow(site_DHW) > 1) {
        for (j1 in 2:nrow(site_DHW)) {
          if (site_DHW$Date[j1-1]+1 == site_DHW$Date[j1]) {
            site_DHW$event[j1] <- site_DHW$event[j1-1]
          }
          if (site_DHW$Date[j1-1]+1 != site_DHW$Date[j1]) {
            site_DHW$event[j1] <- site_DHW$event[j1-1]+1
          }
        } 
      }
      DHW <- NA
      Date <- NA
      df0 <- data.frame(Date, DHW)
      df0 <- df0[0,]
      for (j2 in unique(site_DHW$event)) {
        event1 <- subset(site_DHW, event == j2)
        DHW <- max(event1$DHW)
        Date <- as.character(max(event1$Date[which(event1$DHW == max(event1$DHW))]))
        dhw_df <- data.frame(Date, DHW)
        df0 <- rbind(df0, dhw_df)
      }
      site_DHW <- df0
      site_DHW$Type <- "Heatwave"
      site_DHW$Date <- as.Date(site_DHW$Date, format = "%Y-%m-%d")
    }
    if (nrow(site_DHW) == 0) {
      Date <- NA
      Type <- NA
      site_DHW <- data.frame(Date, Type)
      site_DHW <- site_DHW[0,]
    }
    
  }
  if (length(which(DHW_Sites == i)) == 0) { # If there is no DHW events for the site, add blank dataframe
    Date <- NA
    Type <- NA
    site_DHW <- data.frame(Date, Type)
    site_DHW <- site_DHW[0,]
  }
  ##
  ## Combining disturbances
  if ((nrow(site_cyclone) >= 1) & (nrow(site_DHW) >= 1)) {
    Disturbances <- rbind(site_cyclone[,c("Date", "Type")], site_DHW[,c("Date", "Type")])
  }
  if ((nrow(site_cyclone) == 0) & (nrow(site_DHW) >= 1)) {
    Disturbances <- site_DHW[,c("Date", "Type")]
    
  }
  if ((nrow(site_cyclone) >= 1) & (nrow(site_DHW) == 0)) {
    Disturbances <- site_cyclone[,c("Date", "Type")]
  }
  if ((nrow(site_cyclone) == 0) & (nrow(site_DHW) == 0)) {
    Date <- NA
    Type <- NA
    Disturbances <- data.frame(Date, Type)
    Disturbances <- Disturbances[0,]
  }
  Disturbances <- Disturbances[order(Disturbances$Date),]
  if (nrow(Disturbances) >= 1) {
    #Disturbances1 <- subset(Disturbances, (Date >= min(cover0$Date_Full)) & (Date <= max(cover0$Date_Full)))
    #if (which(Disturbances$Date == min(Disturbances1$Date)) > 1) {
    #  Disturbances1 <- rbind(Disturbances[which(Disturbances$Date == min(Disturbances1$Date))-1,], Disturbances1)
    #}
    
    #Disturbances <- Disturbances1
    
    #Disturbances <- Disturbances[((min(which(Disturbances$Date >= min(cover0$Date_Full))))-1):max(which(Disturbances$Date <= max(cover0$Date_Full))),] # Subsetting for disturbances that occur during survey period at site and one priod
    Disturbances$Event <- 1:nrow(Disturbances)
  }
  
  
  
  for (i1 in 1:nrow(cover0)) {
    Disturbances_1 <- subset(Disturbances, Date <= cover0$Date_Full[i1])
    cover0$Cyclone_Freq[i1] <- nrow(subset(Disturbances_1, Type == "Cyclone"))/(as.numeric(format(cover0$Date_Full[i1], format = "%Y")) - 1970)
 
    if (nrow(Disturbances_1) >= 1) { # if condition for if there are disturbances
      cover0$Dist_Type[i1] <- Disturbances_1$Type[which(Disturbances_1$Date == max(Disturbances_1$Date))]# Type of prior disturbance
      cover0$Dist_Date[i1] <- as.character(max(Disturbances_1$Date)) # last disturance date
      
      
      if (Disturbances_1$Type[which(Disturbances_1$Date == max(Disturbances_1$Date))] == "Cyclone") {
        cover0$Dist_Intensity <- max(site_cyclone$m_per_s[which(site_cyclone$Date == max(Disturbances_1$Date))])
      }
      if (Disturbances_1$Type[which(Disturbances_1$Date == max(Disturbances_1$Date))] == "Heatwave") {
        cover0$Dist_Intensity <- site_DHW$DHW[which(site_DHW$Date == max(Disturbances_1$Date))]
      } 
    }
   
    
    Disturbances_2 <- subset(site_DHW[which(site_DHW$Date %in% Disturbances_1$Date),], DHW >= 4)
    Disturbances_2 <- na.omit(Disturbances_2)
    if (nrow(Disturbances_2) >= 1) {# If condition for last heatwave
      cover0$Last_Heatwave_Date[i1] <- as.character(max(Disturbances_2$Date))  #Date of last heatwave
      cover0$Last_Heatwave_Intensity[i1] <- site_DHW$DHW[which(site_DHW$Date == max(Disturbances_2$Date))] # Intensity
      
      cover0$Heatwave_4_Freq[i1] <- nrow(Disturbances_2)/(as.numeric(format(cover0$Date_Full[i1], format = "%Y")) - 1985)

      Disturbances_3 <- subset(site_DHW[which(site_DHW$Date %in% Disturbances_2$Date),], DHW >= 8)
      Disturbances_4 <- subset(site_DHW[which(site_DHW$Date %in% Disturbances_2$Date),], DHW < 8)
      
      cover0$Heatwave_8_Freq[i1] <- nrow(Disturbances_3)/(as.numeric(format(cover0$Date_Full[i1], format = "%Y")) - 1985)
      cover0$Heatwave_4_8_Freq[i1] <- nrow(Disturbances_4)/(as.numeric(format(cover0$Date_Full[i1], format = "%Y")) - 1985)
      
      cover0$Heatwaves_4_Number[i1] <- nrow(Disturbances_2)
      cover0$Heatwaves_8_Number[i1] <- nrow(Disturbances_3)
      cover0$Heatwaves_4_8_Number[i1] <- nrow(Disturbances_4)
      
    }
    if (nrow(Disturbances_2) == 0) {
      cover0$Last_Heatwave_Date[i1] <- NA
      cover0$Last_Heatwave_Intensity[i1] <- NA
      cover0$Heatwave_4_Freq[i1] <- NA
      cover0$Heatwave_8_Freq[i1] <- NA
      cover0$Heatwave_4_8_Freq[i1] <- NA
      cover0$Heatwaves_4_Number[i1] <- NA
      cover0$Heatwaves_8_Number[i1] <- NA
      cover0$Heatwaves_4_8_Number[i1] <- NA
    }
  }
  
  blank <- rbind(blank, cover0)
}

plot(blank$Heatwave_4_Freq, blank$Heatwave_8_Freq)
plot(blank$Heatwave_4_Freq, blank$Heatwaves_4_Number)
plot(blank$Heatwave_8_Freq, blank$Heatwaves_8_Number)
plot(blank$Heatwave_4_8_Freq, blank$Heatwaves_4_8_Number)


data <- blank





# Splitting dataframe -----------------------------------------------------

# Filtering for sites with 3 or more points
Two_Points <- data[0,]
blank <- data[0,]
for (i in levels(factor(data$Site_ID))) {
  x <- subset(data, Site_ID == i)
  if (nrow(x) >= 3) {
    blank <- rbind(blank, x)
  }
  if (nrow(x) == 2) {
    if ((x$X.Hard_Coral_Cover[2]-x$X.Hard_Coral_Cover[1]) > 0) {
      Two_Points <- rbind(Two_Points, x)
    }
  }
}
data <- blank
remove(blank)

### Recovery Detector V5 (NEW)---------------------

Sites <- levels(factor(data$Site_ID))
data$recovery_number <- NA
data$disturbance_number <- NA
Norm_Coral_Cover <- data$X.Hard_Coral_Cover
blank <- data[0,]

## Recovery Series Detector
#This set of code recognizes increases in time series, and will mark increases as recovery series

for (x in Sites) {
  
  
  ### Subsetting to the Site_ID
  series <- subset(data, Site_ID == x)
  series$recovery_number <- NA # Marking the recovery number as NA just in case
  
  
  
  ### Creating relative coral cover
  # These lines of code stretch the coral cover so that the upper value in the
  # series is equal to 100%, and the lower value in the series is equal to 0%.
  # This allows me to add a relative 10% error to the data when analyzing.
  # The error has to be relative, because some sites such as the Atlantic time series
  # have much lower upper coral cover values than sites in the Pacific, so adding
  # an absolute 10% value would cause bias across sites.
  series$X.Hard_Coral_Cover <- series$X.Hard_Coral_Cover - min(series$X.Hard_Coral_Cover) # Setting minimum cover to 0 by subtracting the minimum coral cover from all values in the vector
  series$X.Hard_Coral_Cover <- (series$X.Hard_Coral_Cover/(max(series$X.Hard_Coral_Cover)))*100 # Setting upper cover to 100% by dividing by the max coral cover in series and then multiplying by 100
  
  
  
  ### Finding vertices
  # These lines detect local minima and maxima in the time series.
  # This is a nice little shortcut to find the sets of points with the greatest 
  # magnitude of change, since the greatest magnitude of change will always occur between
  # minima and maxima. This just allows us to calculate the values faster.
  ## Finding the maxima
  point <- which(diff(c(TRUE,diff(series$X.Hard_Coral_Cover)>=0,FALSE))<0) # This line detects which rows contain "peaks"
  type <- rep("peak", times = length(point)) # Creating labels for these points as peaks
  peaks <- data.frame(point, type) # Creating a dataframe for the peaks
  ## Finding the minima
  point <- which(diff(c(FALSE,diff(series$X.Hard_Coral_Cover)>0,TRUE))>0) # This line detects which rows contain "valleys"
  type <- rep("valley", times = length(point)) # Creating labels for points as valleys
  valleys <- data.frame(point, type) # Creating dataframe for valleys
  ## Combining the dataframes
  vertices <- rbind(peaks, valleys) # Combining peaks and valleys
  vertices <- vertices[order(vertices$point),] # Sorting by the original row number (this sorts the points chronologically
  
  
  
  ### Expanded grid
  # Expand.grid is a function that creates a dataframe with all possible 
  # combinations of the vectors you provide it.
  # In this case, I am giving it a vector that contains all the points that are
  # marked as minima and maxima. Using this vector, it will create a dataframe with
  # all possible combinations of points. Then I just filter out the sets that 
  # don't make sense (e.g. the first point takes place chronologically after the second point,
  # or if the starting and ending point are the same).
  points <- vertices$point # Creating a vector of all points from the prior dataframe
  point_grid <- expand.grid(points, points) # Creating a grid of all possible combinations of points
  point_grid$Duplicate <- "N" # Marking all points as "N" for duplication (This is just a default that will be altered in the next lines)
  
  
  
  ### Filtering the expanded grid
  # Here I sort the expanded grid to remove all sets of points where the starting
  # point has a higher value than the ending point (takes place chronologically after the ending point),
  # or where the starting and ending point are the same.
  ## Marking Bad Rows of Point Grid
  for (z in 1:nrow(point_grid)) { # Start of duplicate removal loop
    x1 <- point_grid[z, 1] #Specifying starting column
    x2 <- point_grid[z, 2] #Specifying ending column
    x3 <- x1 >= x2  # If the value in the starting column is greater than or equal to the value in the ending column..                
    if (x3 == TRUE) {
      point_grid$Duplicate[z] <- "Y" #If the row is bad we mark it with a "Y" so we can remove it later
    }
  } # End of duplicate removal loop
  point_grid <- subset(point_grid, Duplicate != "Y") # Subsetting to remove duplicates marked by prior loop
  point_grid <- point_grid[order(point_grid$Var1, decreasing = FALSE),] #O rdering by the value of the first column (value of starting row of Point Grid)
  point_grid$sigma <- NA # Creating the sigma column and marking all values as NA for now. This will become the row where we calculate the change between points.
  
  
  
  ### Calculating Sigma
  # Here I calculate sigma between the points in the point grid.
  # I do this by taking a row in the point grid, and determining its starting and
  # ending value. Then I made two vectors from it. The first goes from the starting
  # point to the point before the end point. The second goes from one past the starting
  # point to the end point. From here, I get the coral cover values for all these points.
  # The change in coral cover between each point it just the first vector minus the second vector.
  # Once the change in coral cover has been determined, I take the sum of all differences and divide
  # by the number of years that this series took place. I divide by the N of years 
  # between the starting and ending point so that the sigma is not biased towards longer time series.
  # We also add an extra 10% to the sum of coral cover to account for sampling error
  for (i in 1:nrow(point_grid)) {
    start_var <- point_grid$Var1[i] # Starting point
    end_var <- point_grid$Var2[i] # Ending point
    start_cc <- series$X.Hard_Coral_Cover[start_var:(end_var-1)] # Starting coral cover
    end_cc <- series$X.Hard_Coral_Cover[(start_var+1):end_var] # Ending coral cover
    point_grid$sigma[i] <- (sum(end_cc - start_cc)+10)/(series$Date_Diff[end_var] - series$Date_Diff[start_var]) # Calculating sigma for set of points
  }
  
  
  
  ### Determining Recovery Series
  # In this loop we are determining which sets of points are recovery series.
  # We do this by starting at the first point in the series, finding out at what
  # point it reaches its max magnitude, and whether the magnitude is negative or positive.
  # If the magnitude is negative it is skipped and the ending point of that series
  # is used as the starting point for the next series. If the magnitude is positive
  # then the set of points is marked as a recovery series, and the end point is marked as the
  # starting point for the next iteration of the loop. 
  # Also, if the first point in the time series is part of a recovery series, we check to see if
  # the change between the first point and the second point is positive or negative.
  # If the change is negative we mark the start of the time series as the second point.
  # This is done for the first point because we have a bias to include the first point
  # as a part of a recovery series since it is the only point in the time series that we consistently
  # start with. All other points are determined based on magnitude and where a prior series ended, 
  # so we don't need to worry about those. 
  end_point <- min(point_grid$Var1) # Marking the end_point object as the first point in the series
  for (y in as.numeric(levels(factor(point_grid$Var1)))) {
    if (y == 1 | y == end_point) { # If y is equal to the starting point or the end_point object we run this part. This if condition prevents all lines of the point grid from being run
      
      #Subsetting to the desired row of point_grid
      y1 <- subset(point_grid, Var1 == y)
      
      
      #Making sure the first time step is positive
      change_1 <- series$X.Hard_Coral_Cover[y+1] - series$X.Hard_Coral_Cover[y] # Change in coral cover in first time step
      if (change_1 > 0) { # If it is a positive change it is recovery
        
        # Finding the maximum magnitude of y1, and subsetting to that row (y2)
        y1mag <- max(abs(y1$sigma)) # Taking absolute value to find maximum magnitude and ignore any negative symbols
        y2 <- subset(y1, sigma == y1mag) # Subsetting to the row with the maximum magnitude
        if (nrow(y2) == 0) { # If the maximum magnitude is actually negative, then the nrow will equal 0, in which case this condition is run
          y2 <- subset(y1, sigma == -y1mag) # Subsetting for a negative magnitude
        }
        if (nrow(y2) > 1) { # If there happen to be two rows with identical magnitude
          y2 <- y2[nrow(y2),] # Subset to the last row chronologically
        }
        
        
        
        for (m in 1:nrow(y1)) { # This loop makes sure the series doesn't dip below the original starting cover
          sigma <- y1$sigma[m] # Coral cover change in row m
          
          if (sigma < 0) { # If change goes below starting cover
            break_point <- y1$Var2[m-1] # Marks the break point as the Var2 parameter (will become end_point)
            z1 <- which(y1$Var2 == break_point) # The row in which y1 sigma passes 0 (going negative)
            z2 <- max(which(y1$sigma == y2$sigma)) # The row in which y1 has max magnitude
            
            if (z2 > z1) { # If the max magnitude row is past the row where the sum dipped negative..
              end_point <- break_point # Marks end_point as break point
              y2$sigma <- -9999 # Marks y2 sigma as -9999 for an if condition later
            }
            remove(break_point) # Break_point is removed to prevent future iterations of the loop from using it
            break # Loop is broken
          }
        }
        
        # NOTE! I am setting recovery to be > 1% relative change in coral cover. This is an arbetrary choice, and can be tweaked. The main reason is for sites like 764 where the trend is declining, but it sees a step thats > 0 and says yeah thats recovery
        if (y2$sigma > 1) { # If this condition is true the change is positive, now we need to mark as a recovery series
          rec_num <- max(as.integer(levels(factor(series$recovery_number)))) # Checks what the prior recovery number was
          if (is.infinite(rec_num)) { # If there wasn't a recovery number, it'll come up as infinite
            series$recovery_number[y2$Var1:y2$Var2] <- 1 # Here we assign the first recovery number
          }
          else { # If there were prior recovery numbers
            rec_num <- rec_num + 1 # We assign a new recovery number that is one above the last one
            series$recovery_number[y2$Var1:y2$Var2] <- rec_num # Here we set the recovery number in the proper rows
          }
          end_point <- y2$Var2 # And we mark the end_point as the last point in the series so that we can continue the loop
        }
        
        if (y2$sigma == -9999) {
          # This part is the continuation of what we do if the recovery series dipped below the starting coral cover.
          # We mark it as a recovery series up to that break_point
          rec_num <- max(as.integer(levels(factor(series$recovery_number)))) # Checks what the prior recovery number was
          if (is.infinite(rec_num)) { # If there wasn't a recovery number, it'll come up as infinite
            series$recovery_number[y:end_point] <- 1 # Here we assign the first recovery number
          }
          else { # If there were prior recoveries we do the same steps as seen above
            rec_num <- rec_num + 1
            series$recovery_number[y:end_point] <- rec_num
          }
        }
        #
        #
      }
      if (change_1 < 0) { # If the change was less than 0 it was a disturbance
        end_point <- y1$Var2[1] # We simply mark the end_point and move on
      }
    }
  }
  
  
  
  ### Recovery Combiner
  # This loop checks if the row prior to a recovery series was also a recovery series
  # If it is, then the loop checks the difference between the end of the prior recovery series and the current recovery series
  # If the difference is greater than -5% (still a relative cover term, not absolute), 
  # then the loop asks if the starting value of current recovery series is greater than the starting value of the prior recovery series
  # If the starting value of the current recovery series is greater than the starting value of the prior series, then it combines the series into one
  series$recovery_number[which(is.na(series$recovery_number) == TRUE)] <- 0
  rec_numbers <- rev(as.numeric(levels(factor(series$recovery_number))))
  rec_numbers <- rec_numbers[1:length(rec_numbers)-1]
  if (length(rec_numbers) > 1) {
    for (v in rec_numbers) {
      if (v > 1) { #first recovery will not have a prior recovery series
        line <- min(which(series$recovery_number == v)) #says which rows in series have recovery number v
        rec_check <- series$recovery_number[line - 1]
        if (rec_check != 0) { #if the prior line is a recovery number
          v1 <- series$X.Hard_Coral_Cover[line-1] #end of prior rec series
          v2 <- series$X.Hard_Coral_Cover[line] #start of current rec series
          
          v_diff <- v2 - v1
          if (v_diff > -10) { #if greater than -10% error (more positive)
            v3 <- series$X.Hard_Coral_Cover[min(which(series$recovery_number == v)-1)] #starting value for coral cover of prior rec series #### THIS LINE ISNT WORKING BECAUSE WE OVERWROTE A PRIOR RECOVERY NUMBER
            if (v2 > v3) {
              v4 <- which(series$recovery_number == v)
              series$recovery_number[v4] <- v-1
            }
          }
        }
      }
    } #end of recovery combiner
  }
  
  
  
  ### Recovery Renamer
  # Because of the possibility of overwriting recovery numbers in the prior step,
  # this loop renames all recovery numbers to ensure that the numbers are linear
  # and in order.
  s0 <- as.numeric(levels(factor(series$recovery_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      series$recovery_number[which(series$recovery_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  
  blank <- rbind(blank, series)
}

### Disturbance detector V2 (NEW) -------------------------------------------------

data <- blank
data$disturbance_number <- NA
blank <- data[0,]

for (x in Sites) {
  series <- subset(data, Site_ID == x)
  
  
  ### Creating relative coral cover
  # These lines of code stretch the coral cover so that the upper value in the
  # series is equal to 100%, and the lower value in the series is equal to 0%.
  # This allows me to subtract a relative 10% error to the data when analyzing.
  # The error has to be relative, because some sites such as the Atlantic time series
  # have much lower upper coral cover values than sites in the Pacific, so subtracting
  # an absolute 10% value would cause bias across sites.
  series$X.Hard_Coral_Cover <- series$X.Hard_Coral_Cover - min(series$X.Hard_Coral_Cover) #setting minimum cover to 0
  series$X.Hard_Coral_Cover <- (series$X.Hard_Coral_Cover/(max(series$X.Hard_Coral_Cover)))*100 #creating relative term where 100 is the max cover (minimum cover is 0 from prior line)
  
  
  
  ### Finding vertices
  # These lines detect local minima and maxima in the time series.
  # This is a nice little shortcut to find the sets of points with the greatest 
  # magnitude of change, since the greatest magnitude of change will always occur between
  # minima and maxima. This just allows us to calculate the values faster.
  ## Finding the maxima
  point <- which(diff(c(TRUE,diff(series$X.Hard_Coral_Cover)>=0,FALSE))<0) #this is for finding peaks
  type <- rep("peak", times = length(point))
  peaks <- data.frame(point, type)
  ## Finding minima
  point <- which(diff(c(FALSE,diff(series$X.Hard_Coral_Cover)>0,TRUE))>0) #this is for finding valleys
  type <- rep("valley", times = length(point))
  valleys <- data.frame(point, type)
  ## Making dataframe
  vertices <- rbind(peaks, valleys) #binding dataframe
  vertices <- vertices[order(vertices$point),] #ordering by point
  
  
  
  
  ### Expanded grid
  # Expand.grid is a function that creates a dataframe with all possible 
  # combinations of the vectors you provide it.
  # In this case, I am giving it a vector that contains all the points that are
  # marked as minima and maxima. Using this vector, it will create a dataframe with
  # all possible combinations of points. Then I just filter out the sets that 
  # don't make sense (e.g. the first point takes place chronologically after the second point,
  # or if the starting and ending point are the same).
  points <- vertices$point
  point_grid <- expand.grid(points, points)
  point_grid$Duplicate <- "N"
  
  
  
  ### Filtering the expanded grid
  # Here I sort the expanded grid to remove all sets of points where the starting
  # point has a higher value than the ending point (takes place chronologically after the ending point),
  # or where the starting and ending point are the same.
  ## Marking Bad Rows of Point Grid
  for (z in 1:nrow(point_grid)) { #start of duplicate removal loop
    x1 <- point_grid[z, 1] #Specifying starting column
    x2 <- point_grid[z, 2] #Specifying ending column
    x3 <- x1 >= x2                  
    if (x3 == TRUE) { #
      point_grid$Duplicate[z] <- "Y" #If the row is bad we mark it with a Y
    }
  } #end of duplicate removal loop
  point_grid <- subset(point_grid, Duplicate != "Y") #subsetting to remove duplicates marked by prior loop
  point_grid <- point_grid[order(point_grid$Var1, decreasing = FALSE),] #Ordering by the value of the first column (value of starting row of Change Grid)
  point_grid$sigma <- NA
  
  
  
  ### Calculating Sigma
  # Here I calculate sigma between the points in the point grid.
  # I do this by taking a row in the point grid, and determining its starting and
  # ending value. Then I made two vectors from it. The first goes from the starting
  # point to the point before the end point. The second goes from one past the starting
  # point to the end point. From here, I get the coral cover values for all these points.
  # The change in coral cover between each point it just the first vector minus the second vector.
  # Once the change in coral cover has been determined, I take the sum of all differences and divide
  # by the number of years that this series took place. I divide by the N of years 
  # between the starting and ending point so that the sigma is not biased towards longer time series.
  # We also subtract 10% from the sum of coral cover to account for sampling error. This
  # is opposite of how we did it for the recovery detector, since we are trying to find disturbances here.
  for (i in 1:nrow(point_grid)) {
    start_var <- point_grid$Var1[i]
    end_var <- point_grid$Var2[i]
    start_cc <- series$X.Hard_Coral_Cover[start_var:(end_var-1)]
    end_cc <- series$X.Hard_Coral_Cover[(start_var+1):end_var]
    point_grid$sigma[i] <- (sum(end_cc - start_cc)-10)/(series$Date_Diff[end_var] - series$Date_Diff[start_var]) #adding 10% negative error in coral cover here
  }
  
  
  
  ### Determining Disturbance Series
  # In this loop we are determining which sets of points are disturbance series.
  # We do this by starting at the first point in the series, finding out at what
  # point it reaches its max magnitude, and whether the magnitude is negative or positive.
  # If the magnitude is positive it is skipped and the ending point of that series
  # is used as the starting point for the next series. If the magnitude is negative,
  # then the set of points is marked as a disturbance series, and the end point is marked as the
  # starting point for the next iteration of the loop. 
  # Also, if the first point in the time series is part of a disturbance series, we check to see if
  # the change between the first point and the second point is positive or negative.
  # If the change is positive we mark the start of the time series as the second point.
  # This is done for the first point because we have a bias to include the first point
  # as a part of a disturbance series since it is the only point in the time series that we consistently
  # start with. All other points are determined based on magnitude and where a prior series ended, 
  # so we don't need to worry about those. 
  end_point <- min(point_grid$Var1)
  for (y in as.numeric(levels(factor(point_grid$Var1)))) {
    if (y == 1 | y == end_point) {
      #Subsetting to the desired row
      y1 <- subset(point_grid, Var1 == y)
      
      change_1 <- series$X.Hard_Coral_Cover[y+1] - series$X.Hard_Coral_Cover[y] #change in first time step
      if (change_1 < 0) { #disturbance ###
        
        #finding max mag of y1, and also subsetting to that row (y2)
        y1mag <- max(abs(y1$sigma))
        y2 <- subset(y1, sigma == y1mag)
        if (nrow(y2) == 0) { #if sigma is negative ... 
          y2 <- subset(y1, sigma == -y1mag)
        }
        if (nrow(y2) > 1) { #subsetting if there are more than 1 row that are identical
          y2 <- y2[nrow(y2),]
        }
        
        
        #
        #
        for (m in 1:nrow(y1)) { #this loop mkes sure the series doesn't dip below the original starting cover
          sigma <- y1$sigma[m]
          if (sigma > 0) { #if change goes above starting cover...
            if (m == 1) {
              break_point <- y1$Var2[m] #end Var2 parameter (should become end_point)
            }
            if (m != 1) {
              break_point <- y1$Var2[m-1] #end Var2 parameter (should become end_point)
            }
            z1 <- which(y1$Var2 == break_point) #the row in which y1 sigma passes 0 (going negative)
            z2 <- max(which(y1$sigma == y2$sigma)) #the row in which y1 has max magnitude
            if (z2 > z1) { #if the max mag row is past the row where the sum dipped positive (make end point the break point and restart series),
              #also mark as recovery series up to breakpoint
              end_point <- break_point
              y2$sigma <- 9999
            }
            remove(break_point)
            break
          }
        }
        #
        #
        
        #
        # NOTE! I am setting recovery to be > 1% relative change in coral cover. This is an arbetrary choice, and can be tweaked. The main reason is for sites like 764 where the trend is declining, but it sees a step thats > 0 and says yeah thats recovery
        if (y2$sigma < 1) { #If true it is negative, now we need to mark as a d disturbance
          rec_num <- max(as.integer(levels(factor(series$disturbance_number))))
          if (is.infinite(rec_num)) {
            series$disturbance_number[y2$Var1:y2$Var2] <- 1
          }
          else {
            rec_num <- rec_num + 1
            series$disturbance_number[y2$Var1:y2$Var2] <- rec_num
          }
          #
          #
          end_point <- y2$Var2
        }
        
        if (y2$sigma == 9999) {
          rec_num <- max(as.integer(levels(factor(series$disturbance_number))))
          if (is.infinite(rec_num)) {
            series$disturbance_number[y:end_point] <- 1
          }
          else {
            rec_num <- rec_num + 1
            series$disturbance_number[y:end_point] <- rec_num
          }
        }
        #
        #
      }
      if (change_1 > 0) {#just specify the next endpoint as Var2 ###
        end_point <- y1$Var2[1]
      }
    } #nothing should be outside this bracket
  }
  
  
  
  ### Disturbance Combiner
  # This loop checks if the row prior to a disturbance series was also a disturbance series
  # If it is, then the loop checks the difference between the end of the prior disturbance series and the current disturbance series
  # If the difference is greater than +5% (still a relative cover term, not absolute), 
  # then the loop asks if the starting value of current disturbance series is greater than the starting value of the prior disturbance series
  # If the starting value of the current disturbance series is less than the starting value of the prior series, then it combines the series into one
  series$disturbance_number[which(is.na(series$disturbance_number) == TRUE)] <- 0
  dist_numbers <- rev(as.numeric(levels(factor(series$disturbance_number))))
  dist_numbers <- dist_numbers[1:length(dist_numbers)-1]
  if (length(dist_numbers) > 1) {
    for (p in dist_numbers) {
      if (p > 1) { #first recovery will not have a prior recovery series
        line <- min(which(series$disturbance_number == p)) #says which rows in series have disturbance number p
        dist_check <- series$disturbance_number[line - 1]
        if (dist_check != 0) { #if the prior line is a recovery number
          v1 <- series$X.Hard_Coral_Cover[line-1] #end of prior rec series
          v2 <- series$X.Hard_Coral_Cover[line] #start of current rec series
          
          v_diff <- v2 - v1
          if (v_diff < 10) { #if the difference is less than +10%
            v3 <- series$X.Hard_Coral_Cover[min(which(series$disturbance_number == p)-1)] #ending value for coral cover of prior rec series
            if (v2 < v3) { #if v2 is less than v3 because we are looking at disturbances, so it is inverse of recovery
              v4 <- which(series$disturbance_number == p)
              series$disturbance_number[v4] <- p-1
            }
          }
        }
      }
    } #end of disturbance combiner
  }
  
  
  
  ### Disturbance Renamer
  # Because of the possibility of overwriting disturbance numbers in the prior step,
  # this loop renames all disturbance numbers to ensure that the numbers are linear
  # and in order.
  s0 <- as.numeric(levels(factor(series$disturbance_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      series$disturbance_number[which(series$disturbance_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  blank <- rbind(blank, series)
}
data <- blank
data$disturbance_number[which(is.na(data$disturbance_number) == TRUE)] <- 0
data$recovery_number[which(is.na(data$recovery_number) == TRUE)] <- 0


data_backup <- data

### Site cleaner (NEW) ------------------------------------------------------------

### Site Cleaner
# This block checks the heads and tails of recovery and disturbance series to
# see if it can extend the recovery or disturbance series at all.
# The reason it has to do this is because the prior detectors find the series
# with the greatest magnitude, so it often leaves out the plateaus

blank <- data[0,]
Sites <- levels(factor(data$Site_ID))

for (i in Sites) {
  
  # Subsetting to the Site_ID and creating relative coral cover
  # See prior annotations for more context
  x <- subset(data, Site_ID == i)
  x$X.Hard_Coral_Cover <- x$X.Hard_Coral_Cover - min(x$X.Hard_Coral_Cover) 
  x$X.Hard_Coral_Cover <- (x$X.Hard_Coral_Cover/(max(x$X.Hard_Coral_Cover)))*100
  
  
  ## Recovery End Checker
  # This block of code checks if the point directly following the recovery series falls within +-10% relative cover of the series maximum cover.
  # If it does, it adds the point to the series and tries the next point.
  # This loop will continue until the if condition is false, at which point it breaks.
  
  recovery_numbers <- as.numeric(levels(factor(x$recovery_number))) # Creating a vector of recovery numbers
  recovery_numbers <- recovery_numbers[2:length(recovery_numbers)] # Removing the 0 from the list
  recovery_numbers <- as.numeric(na.omit(recovery_numbers)) # Removing any NA values and making sure the list is numeric. Shouldn't ever have NA's, but just in case.
  if (length(recovery_numbers) >= 1 & sum(recovery_numbers) != 0) { # This condition makes sure that there is a recovery series. If the sum is 0 then the only number in the recovery column is 0, in which case we skip
    for (z in recovery_numbers) {
      x1 <- subset(x, recovery_number == z) # Subsetting the recovery series
      max_index <- max(x$index[which(x$recovery_number == z)]) # Looking at the last point in series
      max_cover <- x1$X.Hard_Coral_Cover[nrow(x1)] # Looking at the end cover in the series
      if (max_index != max(x$index)) {
        if (length(which(x$recovery_number == z)) != 0) { # If the length of the vector containing row numbers for the recovery number (z) is not equal to 0..
          sequence <- seq(max(which(x$recovery_number == z)) + 1, max(x$index), 1) # This sequence created contains every row within the recovery series
          
          for (z1 in min(sequence):max(sequence)) { # For every row within the sequence
            max_index <- max(x$index[which(x$recovery_number == z)]) # Looking at the last point in series
            x2 <- x[z1,]
            
            if (max_index == max(x$index)) { # If the max point is the same as the last point in the full series we break the loop
              break
            }
            
            if (max_index != max(x$index)) { # If the max point in the recovery series is not equal to the max point in the full series.. 
              if (between(x2$X.Hard_Coral_Cover, max_cover - 10, max_cover + 10) == TRUE) { # If the next point falls within the +- 10% relative coral cover for the max of series
                x$recovery_number[max_index+1] <- z # Marking the point as a part of the recovery series
              }
              
              if ((between(x2$X.Hard_Coral_Cover, max_cover - 10, max_cover + 10) == FALSE)) { # If the next point doesn't fall within the +- 10% relative coral cover for the max of series
                recovery_numbers <- levels(factor(x$recovery_number))
                recovery_numbers <- as.integer(recovery_numbers)
                recovery_numbers <- recovery_numbers[2:length(recovery_numbers)]
                break
              }  
            }
          }
        }
      }
    }
  }
  
  
  ## Disturbance End Checker
  # This block of code checks if the point directly following the disturbance series falls within +-10% relative cover of the series minimum cover.
  # If it does, it adds the point to the series and tries the next point.
  # This loop will continue until the if condition is false, at which point it breaks.  
  #This loop is the same as the one above, but for disturbances. 
  #There are some minor differences, such as this one using minimum cover instead of maximum, but overall the function is the same
  
  disturbance_numbers <- as.numeric(levels(factor(x$disturbance_number))) # Getting levels of disturbance numbers
  if (disturbance_numbers[1] == 0) {
    disturbance_numbers <- disturbance_numbers[2:length(disturbance_numbers)] # Removing the 0 from the list
  }
  disturbance_numbers <- as.numeric(na.omit(disturbance_numbers)) # Omitting NA values
  if (length(disturbance_numbers) >= 1 & sum(disturbance_numbers) != 0) { # If the length of the disturbance number vectors is greater than 0 and the sum does not equal 0..
    for (g in disturbance_numbers) {
      x1 <- subset(x, disturbance_number == g) # Subsetting to the disturbance series
      max_index <- max(x$index[which(x$disturbance_number == g)]) # Looking at the last point in series
      min_cover <- x1$X.Hard_Coral_Cover[nrow(x1)] # Looking at the end cover in the series
      
      if (length(which(x$disturbance_number == g)) != 0) { # If the length of the disturbance number vector does not equal 0..
        if ((max(which(x$disturbance_number == g)) + 1 <= max(x$index))) { # If the max row with the disturbance number g + 1 is less than the max index
          
          sequence <- seq(max(which(x$disturbance_number == g)) + 1, max(x$index), 1) # Sequence for row values within the disturbance number
          for (g1 in min(sequence):max(sequence)) {
            max_index <- max(x$index[which(x$disturbance_number == g)]) # Looking at the last point in series
            x2 <- x[g1,]
            
            if (max_index == max(x$index)) { # If the max point is the same as the last point in the series we break the loop
              break
            }
            
            if (max_index != max(x$index)) {
              if (between(x2$X.Hard_Coral_Cover, min_cover - 10, min_cover + 10) == TRUE) { # If the next point falls within the + or - 5% relative coral cover for the max of series
                x$disturbance_number[max_index+1] <- g # Marking the point as a part of the disturbance series
              }
              
              if ((between(x2$X.Hard_Coral_Cover, min_cover - 10, min_cover + 10) == FALSE)) {
                disturbance_numbers <- levels(factor(x$disturbance_number))
                disturbance_numbers <- as.integer(disturbance_numbers)
                disturbance_numbers <- disturbance_numbers[2:length(disturbance_numbers)]
                break
              }  
            }
          }
        }
      }
    }
  }
  
  ## Recovery Beginning Checker
  # This loop looks at the beginning of a recovery series to see if there is a
  # plateau before the series that can be incorporated into the series itself.
  # It does this by considering each point prior to the start of the series and
  # seeing if it falls within the +-10% bounds of the minimum cover of the series.
  # If it does, it adds the point to the series and looks at the next point that
  # is one time step prior to the point that was just observed.
  # It is worth noting that this checks the values in reverse chronological order,
  # so it will never include points that are disconnected from the main series 
  # chronologically.
  recovery_numbers <- as.numeric(levels(factor(x$recovery_number))) # Creating a vector of recovery numbers
  if (recovery_numbers[1] == 0) {
    recovery_numbers <- recovery_numbers[2:length(recovery_numbers)] # Removing the 0 from the list
  }
  recovery_numbers <- as.numeric(na.omit(recovery_numbers)) # Removing any NA values and making sure the list is numeric
  if (length(recovery_numbers) >= 1 & sum(recovery_numbers) != 0) { # This condition makes sure that there is a recovery series. If the sum is 0 then the only number in the recovery column is 0, in which case we skip
    for (z in recovery_numbers) {
      x1 <- subset(x, recovery_number == z) # Subsetting to the recovery series
      min_index <- min(x$index[which(x$recovery_number == z)]) # Looking at the first point in series
      min_cover <- x1$X.Hard_Coral_Cover[1] # Minimum cover is used for recovery because we are looking at the starting point
      
      if (min_index != min(x$index)) {
        if (length(which(x$recovery_number == z)) != 0) {
          sequence <- rev(seq(min(x$index), min(which(x$recovery_number == z)) - 1, 1))
          for (z1 in max(sequence):min(sequence)) {
            min_index <- min(x$index[which(x$recovery_number == z)]) # Looking at the first point in series
            x2 <- x[z1,]
            
            if (min_index == min(x$index)) { # If the max point is the same as the last point in the series we break the loop
              break
            }
            
            if (min_index != min(x$index)) {
              if (between(x2$X.Hard_Coral_Cover, min_cover - 10, min_cover + 10) == TRUE) { # If the next point falls within the + or - 5% relative coral cover for the max of series
                x$recovery_number[min_index-1] <- z # Marking the point as a part of the recovery series
              }
              if ((between(x2$X.Hard_Coral_Cover, min_cover - 10, min_cover + 10) == FALSE)) {
                recovery_numbers <- levels(factor(x$recovery_number))
                recovery_numbers <- as.integer(recovery_numbers)
                recovery_numbers <- recovery_numbers[2:length(recovery_numbers)]
                break
              }  
            }
          }
        }
      }
    }
  }
  
  
  ## Disturbance Beginning Checker
  # This loop looks at the beginning of a disturbance series to see if there is a
  # plateau before the series that can be incorporated into the series itself.
  # It does this by considering each point prior to the start of the series and
  # seeing if it falls within the +-10% bounds of the maximum cover of the series.
  # If it does, it adds the point to the series and looks at the next point that
  # is one time step prior to the point that was just observed.
  # It is worth noting that this checks the values in reverse chronological order,
  # so it will never include points that are disconnected from the main series 
  # chronologically.
  disturbance_numbers <- as.numeric(levels(factor(x$disturbance_number))) # Creating a vector of recovery numbers
  if (disturbance_numbers[1] == 0) {
    disturbance_numbers <- disturbance_numbers[2:length(disturbance_numbers)] # Removing the 0 from the list
  }
  disturbance_numbers <- as.numeric(na.omit(disturbance_numbers)) # Removing any NA values and making sure the list is numeric
  if (length(disturbance_numbers) >= 1 & sum(disturbance_numbers) != 0) { # This condition makes sure that there is a recovery series. If the sum is 0 then the only number in the recovery column is 0, in which case we skip
    for (z in disturbance_numbers) {
      x1 <- subset(x, disturbance_number == z) # Subsetting to the recovery series
      min_index <- min(x$index[which(x$disturbance_number == z)]) # Looking at the first point in series
      max_cover <- x1$X.Hard_Coral_Cover[1]
      
      if (min_index != min(x$index)) {
        if (length(which(x$disturbance_number == z)) != 0) {
          sequence <- rev(seq(min(x$index), min(which(x$disturbance_number == z)) - 1, 1))
          for (z1 in max(sequence):min(sequence)) {
            min_index <- min(x$index[which(x$disturbance_number == z)]) # Looking at the first point in series
            x2 <- x[z1,]
            
            if (min_index == min(x$index)) { # If the max point is the same as the last point in the series we break the loop
              break
            }
            
            if (min_index != min(x$index)) {
              if (between(x2$X.Hard_Coral_Cover, max_cover - 10, max_cover + 10) == TRUE) { # If the next point falls within the + or - 5% relative coral cover for the max of series
                x$disturbance_number[min_index-1] <- z # Marking the point as a part of the recovery series
              }
              if ((between(x2$X.Hard_Coral_Cover, max_cover - 10, max_cover + 10) == FALSE)) {
                disturbance_numbers <- levels(factor(x$disturbance_number))
                disturbance_numbers <- as.integer(disturbance_numbers)
                disturbance_numbers <- disturbance_numbers[2:length(disturbance_numbers)]
                break
              }  
            }
          }
        }
      }
    }
  }
  
  
  
  series <- x # Renaming x to series
  
  ## Recovery Renamer
  # See prior annotations
  s0 <- as.numeric(levels(factor(series$recovery_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      series$recovery_number[which(series$recovery_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  ## Disturbance Renamer
  # See prior annotations
  s0 <- as.numeric(levels(factor(series$disturbance_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      series$disturbance_number[which(series$disturbance_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  
  ## Recovery/Disturbance Combiners
  # This block of code checks if the recovery or disturbance series can be
  # combined with each other into longer time series. 
  # It does this by checking the row prior to the start of a recovery or disturbance
  # series. If the prior line is also a recovery or disturbance series, it checks 
  # if the coral cover has a difference of -10% for recovery series, and 10% for disturbance
  # series. If it falls within these bounds, the series are combined, and renamed.
  # These combiners rerun over and over to ensure that all possible series are combined.
  
  ## Recovery Combiner
  # This loop checks if the prior row is part of a recovery series. If it is, it checks if they are similar enough to be combined
  for (int in 1:(max(series$recovery_number)+2)) {
    rec_numbers <- rev(as.numeric(levels(factor(series$recovery_number))))
    if (rec_numbers[length(rec_numbers)] == 0) {
      rec_numbers <- rec_numbers[1:length(rec_numbers)-1]
    }
    for (v in rec_numbers) {
      if (v > 1) { # First recovery will not have a prior recovery series
        line <- min(which(series$recovery_number == v)) # Saying which rows in series have recovery number v
        rec_check <- series$recovery_number[line - 1]
        if (rec_check != 0) { # If the prior line is a recovery number
          v1 <- series$X.Hard_Coral_Cover[line-1] # End of prior rec series
          v2 <- series$X.Hard_Coral_Cover[line] # Start of current rec series
          
          v_diff <- v2 - v1
          if (v_diff > -10) { # If greater than -10% error (more positive)
            v3 <- series$X.Hard_Coral_Cover[min(which(series$recovery_number == v)-1)] # Starting value for coral cover of prior rec series #### THIS LINE ISNT WORKING BECAUSE WE OVERWROTE A PRIOR RECOVERY NUMBER
            if (v2 > v3 - 10) {
              v4 <- which(series$recovery_number == v)
              series$recovery_number[v4] <- v-1
            }
          }
        }
      }
    } # End of recovery combiner
  }
  
  ## Disturbance Combiner
  # See annotations above
  for (int in 1:(max(series$disturbance_number)+2)) {
    # This loop checks if the prior row is part of a disturbance series. If it is, it checks if they are similar enough to be combined  
    dist_numbers <- rev(as.numeric(levels(factor(series$disturbance_number))))
    if (dist_numbers[length(dist_numbers)] == 0) {
      dist_numbers <- dist_numbers[1:length(dist_numbers)-1]
    }
    if (length(dist_numbers) > 1) {
      for (p in dist_numbers) {
        if (p > 1) { # First disturbance will not have a prior disturbance series
          line <- min(which(series$disturbance_number == p)) # Saying which rows in series have disturbance number p
          dist_check <- series$disturbance_number[line - 1]
          if (dist_check != 0) { # If the prior line is a recovery number
            v1 <- series$X.Hard_Coral_Cover[line-1] # End of prior rec series
            v2 <- series$X.Hard_Coral_Cover[line] # Start of current rec series
            
            v_diff <- v2 - v1
            if (v_diff < 10) { # If the difference is less than +10%
              v3 <- series$X.Hard_Coral_Cover[min(which(series$disturbance_number == p)-1)] # Ending value for coral cover of prior rec series
              if (v2 < v3 + 10) { # If v2 is less than v3 because we are looking at disturbances, so it is inverse of recovery
                v4 <- which(series$disturbance_number == p)
                series$disturbance_number[v4] <- p-1
              }
            }
          }
        }
      } # End of disturbance combiner
    }
  }
  
  
  
  ## Recovery Renamer
  # See prior annotations
  s0 <- as.numeric(levels(factor(series$recovery_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      series$recovery_number[which(series$recovery_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  ## Disturbance Renamer
  # See prior annotations
  s0 <- as.numeric(levels(factor(series$disturbance_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      series$disturbance_number[which(series$disturbance_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  blank <- rbind(blank, series)
}

data <- blank

data_backup_2 <- data
#removing disturbances inside recovery series and vice versa
#loop would start here



### Internal series remover (NEW) -------------------------------------------------

#This loop removes disturbance and recovery series that are within another series (e.g. a recovery series occurs within a disturbance series)
#This removes miniseries that occur within overarching series of the opposite type if the change in coral cover within the miniseries is under a absolute magnitude of 10%

## Internal Series Remover
# This block removes any disturbance or recovery series that are nested inside
# a series of the opposite type if the change in relative coral cover is less than
# +10% for recovery, and -10% for distrubances.
# The logic behind this is that we are considering a +-10% human error in surveying
# So if the change is less than this error within a larger trend, we should remove
# the smaller trend to reduce noise in our data.
blank <- data[0,]
for (i in Sites) {
  
  # Subsetting and making relative coral cover
  x <- subset(data, Site_ID == i)
  x$X.Hard_Coral_Cover <- x$X.Hard_Coral_Cover - min(x$X.Hard_Coral_Cover)
  x$X.Hard_Coral_Cover <- (x$X.Hard_Coral_Cover/(max(x$X.Hard_Coral_Cover)))*100 
  
  ## Recovery Remover
  recovery_numbers <- as.numeric(levels(factor(x$recovery_number)))
  if (sum(recovery_numbers) != 0) {
    if (recovery_numbers[1] == 0) { # Removing 0 from the mix
      recovery_numbers <- recovery_numbers[2:length(recovery_numbers)]
    }
    for (i1 in recovery_numbers) {
      # What are the starting and ending points for recovery series
      start_rec_point <- min(which(x$recovery_number == i1))
      end_rec_point <- max(which(x$recovery_number == i1))
      # Do they have a disturbance number
      start_rec_point_dist_num <- x$disturbance_number[start_rec_point] # Pulling disturbance number from start point
      end_rec_point_dist_num <- x$disturbance_number[end_rec_point] # Pulling from end point
      # Is the disturbance number the same
      if (start_rec_point_dist_num == end_rec_point_dist_num & end_rec_point_dist_num != 0 & start_rec_point_dist_num != 0) {
        
        if (start_rec_point != 1) {
          prior_point_dist <- x$disturbance_number[start_rec_point - 1]
        }
        if (start_rec_point == 1) {
          prior_point_dist <- x$disturbance_number[start_rec_point]
        }
        
        if (end_rec_point == max(x$index)) {
          next_point_dist <- x$disturbance_number[end_rec_point]
        }
        if (end_rec_point != max(x$index)) {
          next_point_dist <- x$disturbance_number[end_rec_point + 1]
        }
        
        
        if (is.na(next_point_dist) == TRUE) {
          next_point_dist <- prior_point_dist
        }
        if (is.na(prior_point_dist) == TRUE) {
          prior_point_dist <- next_point_dist
        }
        if (prior_point_dist == next_point_dist & next_point_dist != 0 & prior_point_dist != 0) { #if the points have the same value and are not equal to 0
          CC_change <- x$X.Hard_Coral_Cover[end_rec_point] - x$X.Hard_Coral_Cover[start_rec_point] 
          if (CC_change < 10) {
            x$recovery_number[which(x$recovery_number == i1)] <- 0
          }
        }
      }
    }
    
    ## Recovery Renamer
    # See prior annotations
    s0 <- as.numeric(levels(factor(x$recovery_number)))
    if (s0[1] == 0) {
      s0 <- na.omit(s0[2:length(s0)])
    }
    if (length(s0) != 0 & sum(s0) != 0) {
      s2 <- 1:length(s0)
      s3 <- cbind(s0, s2)
      for (s in 1:nrow(s3)) {
        x$recovery_number[which(x$recovery_number == s3[s,1])] <- s3[s,2]
      }
    }
  }
  
  ## Disturbance Remover
  disturbance_numbers <- as.numeric(levels(factor(x$disturbance_number)))
  if (sum(disturbance_numbers) != 0) {
    if (disturbance_numbers[1] == 0) { # Removing 0 from the mix
      disturbance_numbers <- disturbance_numbers[2:length(disturbance_numbers)]
    }
    for (i2 in disturbance_numbers) {
      # What are the starting and ending points for recovery series
      start_dist_point <- min(which(x$disturbance_number == i2))
      end_dist_point <- max(which(x$disturbance_number == i2))
      
      # Do they have a recovery number
      start_dist_point_rec_num <- x$recovery_number[start_dist_point] #pulls disturbance number from start point
      end_dist_point_rec_num <- x$recovery_number[end_dist_point] #pulls from end point
      
      if (start_dist_point_rec_num == end_dist_point_rec_num & end_dist_point_rec_num != 0 & start_dist_point_rec_num != 0) {
        if (start_dist_point != 1) {
          prior_point_rec <- x$recovery_number[start_dist_point - 1]
        }
        if (start_dist_point == 1) {
          prior_point_rec <- x$recovery_number[start_dist_point]
        }
        
        if (end_dist_point == max(x$index)) {
          next_point_rec <- x$recovery_number[end_dist_point]
        }
        if (end_dist_point != max(x$index)) {
          next_point_rec <- x$recovery_number[end_dist_point + 1]
        }
        
        if (is.na(next_point_rec) == TRUE) {
          next_point_rec <- prior_point_rec
        }
        if (is.na(prior_point_rec) == TRUE) {
          prior_point_rec <- next_point_rec
        }
        if (prior_point_rec == next_point_rec & next_point_rec != 0 & prior_point_rec != 0) { #if the points have the same value and are not equal to 0
          CC_change <- x$X.Hard_Coral_Cover[end_dist_point] - x$X.Hard_Coral_Cover[start_dist_point] 
          if (CC_change > -10) {
            x$disturbance_number[which(x$disturbance_number == i2)] <- 0
          }
        }
      }
    }
    
    ## Disturbance Renamer
    s0 <- as.numeric(levels(factor(x$disturbance_number)))
    if (s0[1] == 0) {
      s0 <- na.omit(s0[2:length(s0)])
    }
    if (length(s0) != 0 & sum(s0) != 0) {
      s2 <- 1:length(s0)
      s3 <- cbind(s0, s2)
      for (s in 1:nrow(s3)) {
        x$disturbance_number[which(x$disturbance_number == s3[s,1])] <- s3[s,2]
      }
    }
  }
  
  blank <- rbind(blank, x)
}

data <- blank

## Adding back in the standard coral cover since we replaced it with relative cover previously
data$X.Hard_Coral_Cover <- Norm_Coral_Cover


### Start/End Error Checker (NEW) -----------------------------------------------------------
#added in loop that checks normal coral cover for recovery and disturbance series
#if the end coral cover is lower than the starting cover of a recovery series, remove it
#if the end coral cover is higher than the starting coral cover of a disturbance series, remove it

## Start/End Error Checker
# This block checks the recovery and disturbance series for any logical errors
# that may have occured while running the prior blocks. 
# Specifically, if the end coral cover is lower than the starting coral cover
# for a recovery series, the block goes back and removes that end point and 
# removes the last point in the series up to the minimum point in the recovery series
# For disturbance series, if the coral cover at the end is higher than the beginning coral cover,
# then it removed the last point in the series up to the maximum coral cover point within the disturbance series

blank <- data[0,]
for (i in Sites) {
  
  x <- subset(data, Site_ID == i)
  
  recovery_numbers <- as.numeric(levels(factor(x$recovery_number)))
  if (sum(recovery_numbers) != 0) {
    if (recovery_numbers[1] == 0) { #removing 0 from the mix
      recovery_numbers <- recovery_numbers[2:length(recovery_numbers)]
    }
    
    for (i1 in recovery_numbers) {
      CC_1 <- x$X.Hard_Coral_Cover[min(which(x$recovery_number == i1))]
      CC_2 <- x$X.Hard_Coral_Cover[max(which(x$recovery_number == i1))]
      if (CC_1 > CC_2) {
        x1 <- subset(x, recovery_number == i1)
        x$recovery_number[min(which(x$recovery_number == i1)):which(x$X.Hard_Coral_Cover == min(x1$X.Hard_Coral_Cover)) - 1] <- 0
      }
    }
  }
  
  ## Recovery Renamer
  # See prior annotations
  s0 <- as.numeric(levels(factor(x$recovery_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      x$recovery_number[which(x$recovery_number == s3[s,1])] <- s3[s,2]
    }
  }
  
  
  disturbance_numbers <- as.numeric(levels(factor(x$disturbance_number)))
  if (sum(disturbance_numbers) != 0) {
    if (disturbance_numbers[1] == 0) { #removing 0 from the mix
      disturbance_numbers <- disturbance_numbers[2:length(disturbance_numbers)]
    }
    
    for (i2 in disturbance_numbers) {
      CC_1 <- x$X.Hard_Coral_Cover[min(which(x$disturbance_number == i2))]
      CC_2 <- x$X.Hard_Coral_Cover[max(which(x$disturbance_number == i2))]
      if (CC_1 < CC_2) {
        x2 <- subset(x, disturbance_number == i2)
        x$disturbance_number[min(which(x$disturbance_number == i2)):which(x$X.Hard_Coral_Cover == max(x2$X.Hard_Coral_Cover)) - 1] <- 0     
      }
    }
  }
  
  ## Disturbance Renamer
  # See prior annotations
  s0 <- as.numeric(levels(factor(x$disturbance_number)))
  if (s0[1] == 0) {
    s0 <- na.omit(s0[2:length(s0)])
  }
  if (length(s0) != 0 & sum(s0) != 0) {
    s2 <- 1:length(s0)
    s3 <- cbind(s0, s2)
    for (s in 1:nrow(s3)) {
      x$disturbance_number[which(x$disturbance_number == s3[s,1])] <- s3[s,2]
    }
  }
  blank <- rbind(blank, x)
}

data <- blank



###### Site.Plotter Function ---------------------------------------------------

#Site.Plotter function
site.plotter <- function(site) {
  x <- subset(data, Site_ID == site)
  plot(x$Date_Diff, x$X.Hard_Coral_Cover, main = paste("Graph of Site", site), type = "l")
  
  y <- subset(x, is.na(recovery_number) == FALSE)
  y1 <- levels(factor(y$recovery_number))
  for (i in y1) {
    if (i != 0) {
      y2 <- subset(y, recovery_number == i)
      lines(y2$Date_Diff, y2$X.Hard_Coral_Cover, col = "blue", type = "l", pch = 19, lwd = 3)
    }
  }
  
  y3 <- subset(x, is.na(disturbance_number) == FALSE)
  y4 <- levels(factor(y3$disturbance_number))
  for (z in y4) {
    if (z != 0) {
      y5 <- subset(y3, disturbance_number == z)
      lines(y5$Date_Diff, y5$X.Hard_Coral_Cover, col = "red", type = "l", pch = 19, lty = 2, lwd = 3)
    }
  }
  lines(x$Date_Diff, x$X.Hard_Coral_Cover, type = "p", pch = 16, cex = 1.5)
}

#Examples
site.plotter(1015)
site.plotter(1016)

site.plotter(1039)

site.plotter(1033)
site.plotter(1034)
site.plotter(16621)
site.plotter(1255)
site.plotter(1256)
site.plotter(1257)


# Adding back in two point data -------------------------------------------


## Combining Two_Point and data
# calculating date difference for Two_Point
site <- levels(factor(Two_Points$Site_ID))
blank <- Two_Points[0,]
for (i in site) {
  x <- subset(Two_Points, Site_ID == i)
  start_date <- min(as.Date(x$Date_Full, format = "%Y-%m-%d"))
  for (o in 1:nrow(x)) {
    x$Date_Diff[o] <- round((difftime(as.Date(x$Date_Full[o], format = "%Y-%m-%d"), start_date, unit = "days"))/365, digits = 3)
  }
  x <- x[order(x$Date_Diff),] #ordering by date_diff
  x$index <- c(1,2)
  blank <- rbind(blank, x)
}
Two_Points <- blank
Two_Points$recovery_number <- 1 # All sites are recoveries, and all sites only have two points, meaning only one recovery series
Two_Points$disturbance_number <- 0
data <- rbind(data, Two_Points)
data$New_ID <- NA


###### Determining impact for each recovery series (Site_coeff_Improved)-----------------------------

## Creating Site_Coeff df scaffolding
# Making data.frame
## Creating Site_Coeff df scaffolding
# Making data.frame
Site_ID <- NA
New_ID <- NA
Latitude_Degrees <- NA
Longitude_Degrees <- NA
Ocean_Name <- NA
Ecoregion_Name <- NA
Habitat_Type <- NA
Starting_Year <- NA
Ending_Year <- NA
Start_Date <- NA
End_Date <- NA
N_samples <- NA
Depth <- NA
Distance_to_shore <- NA
Exposure <- NA
Turbidity <- NA
Grav_tot <- NA
Prior_Cover <- NA
Initial_Cover <- NA
Initial_Algae <- NA
Impact <- NA
Relative_Impact <- NA
Cyclone_Freq <- NA
Last_Heatwave_Date <- NA
Last_Heatwave_Intensity <- NA

Heatwave_4_8_Freq <- NA
Heatwave_8_Freq <- NA
Heatwave_4_Freq <- NA
Heatwaves_4_Number <- NA
Heatwaves_8_Number <- NA
Heatwaves_4_8_Number <- NA


Site_Coeff <- data.frame(Site_ID, New_ID, Latitude_Degrees, Longitude_Degrees, Ocean_Name, Ecoregion_Name, Habitat_Type, 
                         Starting_Year, Ending_Year, Start_Date, End_Date, N_samples, Depth, Distance_to_shore, Exposure, 
                         Turbidity, Grav_tot, Initial_Cover, Initial_Algae, Prior_Cover, Impact, Relative_Impact,
                         Cyclone_Freq, Last_Heatwave_Date, Last_Heatwave_Intensity, Heatwave_4_8_Freq, Heatwave_4_Freq, Heatwave_8_Freq,
                         Heatwaves_4_8_Number, Heatwaves_4_Number, Heatwaves_8_Number)

Site_Coeff <- Site_Coeff[0,]


data$New_ID <- NA
blank <- data[0,]

for (i in levels(factor(data$Site_ID))) {
  s <- subset(data, Site_ID == i)
  rec_nums <- as.numeric(levels(factor(s$recovery_number)))
  if ((rec_nums[1] == 0) & (length(rec_nums) > 1)) {
    rec_nums <- rec_nums[2:length(rec_nums)]
  }
  for (i1 in rec_nums) {
    r <- subset(s, recovery_number == i1)
    
    ### Creating New_ID for blank df
    r$New_ID <- paste(r$Site_ID, "000", r$recovery_number, sep = "")
    
    ## Combining r back into blank df to become rec_data
    blank <- rbind(blank, r)
    
    #Creating row for Site_Coeff
    Site_ID <- unique(r$Site_ID)
    New_ID <- unique(r$New_ID)
    Latitude_Degrees <- r$Latitude_Degrees[1]
    Longitude_Degrees <- r$Longitude_Degrees[1]
    Ocean_Name <- r$Ocean_Name[1]
    Ecoregion_Name <- r$Ecoregion_Name[1]
    Habitat_Type <- r$Habitat_Type[1]
    Starting_Year <- min(r$Date_Year)
    Ending_Year <- max(r$Date_Year)
    Start_Date <- min(r$Date_Full)
    N_samples <- nrow(r)
    End_Date <- max(r$Date_Full)
    Depth <- mean(r$Depth, na.rm = T)
    Distance_to_shore <- unique(r$Distance_to_shore)
    Exposure <- unique(r$Exposure)
    Turbidity <- as.numeric(unique(r$Turbidity))
    Grav_tot <- as.numeric(unique(r$Grav_tot))
    Initial_Cover <- min(r$X.Hard_Coral_Cover)
    Initial_Algae <- r$X.Macroalgae_Cover[1]
    
    Cyclone_Freq <- r$Cyclone_Freq[1]
    Last_Heatwave_Date <- r$Last_Heatwave_Date[1]
    Last_Heatwave_Intensity <- r$Last_Heatwave_Intensity[1]
    
    Heatwave_4_8_Freq <- r$Heatwave_4_8_Freq[1]
    Heatwave_8_Freq <- r$Heatwave_8_Freq[1]
    Heatwave_4_Freq <- r$Heatwave_4_Freq[1]
    Heatwaves_4_Number <- r$Heatwaves_4_Number[1]
    Heatwaves_8_Number <- r$Heatwaves_8_Number[1]
    Heatwaves_4_8_Number <- r$Heatwaves_4_8_Number[1]
    
    
    #Prior Cov, Impact, Relative Impact
     if (min(r$index) > 1) {
       p_cov <- subset(s, index <= min(r$index))
       
       ### Finding vertices
       # These lines detect local minima and maxima in the time series.
       # This is a nice little shortcut to find the sets of points with the greatest 
       # magnitude of change, since the greatest magnitude of change will always occur between
       # minima and maxima. This just allows us to calculate the values faster.
       ## Finding the maxima
       point <- which(diff(c(TRUE,diff(p_cov$X.Hard_Coral_Cover)>=0,FALSE))<0) # This line detects which rows contain "peaks"
       type <- rep("peak", times = length(point)) # Creating labels for these points as peaks
       peaks <- data.frame(point, type) # Creating a dataframe for the peaks
       
       Prior_Cover <- p_cov$X.Hard_Coral_Cover[which(p_cov$index == max(peaks$point))] # absolute CC value
       
    
       # Pulling disturbance
       # getting lowest first point in recovery
       point <- which(diff(c(FALSE,diff(r$X.Hard_Coral_Cover)>0,TRUE))>0) # This line detects which rows contain "valleys"
       type <- rep("valley", times = length(point)) # Creating labels for points as valleys
       valleys <- data.frame(point, type) # Creating dataframe for valleys
       #
       r$index[min(valleys$point)] # Initial lowest point for R
       #p_cov <- s[s$index[max(peaks$point)]:r$index[min(valleys$point)],] # Range for impact and disturbance
     
       ## Pickup here
       # Determining disturbance
       # Was it the cause of impact?
      # Was there a disturbance in the middle of recovery?
       # Intensity?
      # Etc
     
       #Impact <- Prior_Cover - Initial_Cover
       #Relative_Impact <- ((Prior_Cover/Prior_Cover) - (Initial_Cover/Prior_Cover))*100
       
       
     } else { # If index is 1, indicating first point in series
       Prior_Cover <- NA
       Impact <- NA
       Relative_Impact <- NA
    }
    
    
    Site_Coeff_1 <- data.frame(Site_ID, New_ID, Latitude_Degrees, Longitude_Degrees, Ocean_Name, Ecoregion_Name, Habitat_Type, 
                               Starting_Year, Ending_Year, Start_Date, End_Date, N_samples, Depth, Distance_to_shore, Exposure, 
                               Turbidity, Grav_tot, Initial_Cover, Initial_Algae, Prior_Cover, Impact, Relative_Impact,
                               Cyclone_Freq, Last_Heatwave_Date, Last_Heatwave_Intensity, Heatwave_4_8_Freq, Heatwave_4_Freq, Heatwave_8_Freq,
                               Heatwaves_4_8_Number, Heatwaves_4_Number, Heatwaves_8_Number)
    Site_Coeff <- rbind(Site_Coeff, Site_Coeff_1)
  }
  
  
}

rec_data <- subset(blank, recovery_number != 0)
Site_Coeff <- subset(Site_Coeff, endsWith(as.character(Site_Coeff$New_ID), "0") == F)
Site_Coeff$r_coeff <- NA
Site_Coeff$A_coeff <- NA
Site_Coeff$C_coeff <- NA
Site_Coeff$F_coeff <- NA
Site_Coeff$D_coeff <- NA



## Remaking index for within new recovery series ---------------------------
blank <- rec_data[0,]

for (i in levels(factor(rec_data$New_ID))) {
  x <- subset(rec_data, New_ID == i)
  for (i1 in 1:nrow(x)) {
    x1 <- x$index[i1]
    x$index[i1] <- which(x$index == x1)
  }
  blank <- rbind(blank, x)
}
rec_data <- blank

remove(i, i1, x1, x, blank)


###Setting Each New Site ID Start Date to 0 (need to relabel) --------------------------------

#Creating Vector of them for a loop
New_ID <- levels(factor(rec_data$New_ID))

blank <- rec_data[0,]
#Setting each recovery series start date to 0
for (i in New_ID) {
  x <- subset(rec_data, New_ID == i)
  
  start_date <- min(x$Date_Diff)
  x$Date_Diff <- x$Date_Diff - start_date
  blank <- rbind(blank, x)
}

rec_data <- blank #removes all rows that are not a part of a recovery series

remove(x, start_date, i, blank)




### Sorting for Gompertz and Exponential NLS --------------------------------------------

# Three categories
# 1. Sites that have 2 points so we can only run a base exponential growth equation on
# 2. Sites that are fit for a gompertz curve
# 3. Sites that have 3 or more points and only exponential growth

Site_Coeff$Analysis_Type <- NA


Site_Coeff_backup <- Site_Coeff
Exp_Equation_Only <- subset(Site_Coeff, N_samples == 2)
Site_Coeff <- subset(Site_Coeff, N_samples >= 3)
New_IDs <- levels(factor(Site_Coeff$New_ID))
Gomp_df <- Site_Coeff[0,]
for (i in New_IDs) {
  x <- subset(rec_data, New_ID == i)
  
  # This condition checks to see if the last slope has the greatest magnitude. if it does then it is likely exponential
  slopes <- (x$X.Hard_Coral_Cover[2:nrow(x)] - x$X.Hard_Coral_Cover[1:(nrow(x)-1)])/(x$Date_Diff[2:nrow(x)] - x$Date_Diff[1:(nrow(x)-1)])
  if (slopes[length(slopes)] != max(slopes)) { # If the last slope does not have the biggest magnitude
    
    x1 <- subset(Site_Coeff, New_ID == i)
    Gomp_df <- rbind(Gomp_df, x1)
  }
}
nrow(Gomp_df)
#949 sites

# Sorting out exponential equations

Exp_df <- Gomp_df[0,]
for (i in New_IDs) {
  if (nrow(subset(Gomp_df, New_ID == i)) == 0) {
    x <- subset(Site_Coeff, New_ID == i)
    Exp_df <- rbind(Exp_df, x)
  }
}
nrow(Exp_df) #675

# Removing 3 sample sites from gomp_df and moving them to Exp_Equation df
# This is because Gomp cannot use 3 samples, but Exponential NLS can. So we sort out the 3 sample sites for the Exp NLS
# and the 3 sample sites that do not have clean exponential growth are places in the Exp_Equation df
x <- subset(Gomp_df, N_samples == 3)
Gomp_df <- subset(Gomp_df, N_samples > 3)
nrow(x) + nrow(Gomp_df)
nrow(Exp_Equation_Only)
Exp_Equation_Only <- rbind(Exp_Equation_Only, x)
nrow(Exp_Equation_Only)



nrow(Exp_Equation_Only) + nrow(Gomp_df) + nrow(Exp_df) #3262
nrow(Site_Coeff_backup) # 3262


### Running Exponential Equation for Sites with 2 or 3 Samples -------------------
Exp_Equation_Only$Analysis_Type <- "Exp_Equation"
for (i in 1:nrow(Exp_Equation_Only)) {
  d <- subset(rec_data, New_ID == Exp_Equation_Only$New_ID[i])
  if (nrow(d) == 2) {
    d$X.Hard_Coral_Cover[1] # t coral cover
    d$X.Hard_Coral_Cover[2] # t + 1 coral cover
    d$Date_Diff[2] - d$Date_Diff[1]
    
    Exp_Equation_Only$r_coeff[i] <- log(((d$X.Hard_Coral_Cover[2] + 1)/(d$X.Hard_Coral_Cover[1] + 1)))/(d$Date_Diff[2] - d$Date_Diff[1])
  }
  if (nrow(d) == 3) {
    d$X.Hard_Coral_Cover[1:nrow(d)-1] # t coral cover
    d$X.Hard_Coral_Cover[2:nrow(d)] # t + 1 coral cover
    d$Date_Diff[2:nrow(d)] - d$Date_Diff[1:nrow(d)-1]
    
    Exp_Equation_Only$r_coeff[i] <- max(log(((d$X.Hard_Coral_Cover[2:nrow(d)] + 1)/(d$X.Hard_Coral_Cover[1:nrow(d)-1] + 1)))/(d$Date_Diff[2:nrow(d)] - d$Date_Diff[1:nrow(d)-1]))
  }
  
  
}

plot(Exp_Equation_Only$Starting_Year, Exp_Equation_Only$r_coeff, ylim = c(0, 4))

remove(d)



### Running Gompertz Equation for Gomp_df -----------------------------------

Gomp_df$Analysis_Type <- "Gompertz_NLS"

for (i in 1:nrow(Gomp_df)) {
  
  df <- subset(rec_data, New_ID == Gomp_df$New_ID[i])
  
  x <- df$Date_Diff - df$Date_Diff[1]
  y <- df$X.Hard_Coral_Cover
  
  ## min cover coefficient
  min_cov <- min(df$X.Hard_Coral_Cover)
  
  Gomp <- drm(y ~ x, fct = G.4(fixed = c(NA, min_cov, NA, NA), names = c("b", "c", "d", "e")),
              lowerl = c(NA, max(df$X.Hard_Coral_Cover) - 0.05*(max(df$X.Hard_Coral_Cover) - min(df$X.Hard_Coral_Cover)), -10),
              upperl = c(NA, 100, NA))
  Gomp_df$C_coeff[i] <- -as.numeric(Gomp$coefficients[1]) #b coefficient (rate). Need to remove negative sign
  Gomp_df$A_coeff[i] <- as.numeric(Gomp$coefficients[2]) # d coefficient (max cover)
  Gomp_df$F_coeff[i] <- as.numeric(Gomp$coefficients[3]) #e coefficient (f, inflection point)
  Gomp_df$D_coeff[i] <- min_cov
  
}
Gomp_df <- subset(Gomp_df, C_coeff > 0) # removing series with negative rates
#summary(Gomp)

plot(Gomp_df$Starting_Year, Gomp_df$C_coeff, ylim = c(0, 50))


### Converting Gompertz Equation to Exponential Equation --------------------
Gomp_df <- subset(Gomp_df, F_coeff >= 0)
Gomp_df$Lag <- NA
for (i in 1:nrow(Gomp_df)) {
  df <- Gomp_df[i,]
  a <- df$A_coeff
  c <- df$C_coeff
  f <- df$F_coeff
  d <- df$D_coeff
  t_step <- 0.001
  
  t <- seq(0, f, t_step)
  
  y <- (a-d)*exp(-exp(-c*(t-f)))+d
  
  t0 <- t[max(which(round(y, digits = 2) == round(min(y), digits = 2)))]
  y0 <- y[max(which(round(y, digits = 2) == round(min(y), digits = 2)))]
  
  
  df_1 <- subset(rec_data, New_ID == df$New_ID)
  t1 <- df$F_coeff
  yt <- (a-d)*exp(-exp(-c*(t1-f)))+d #t1
  
  r <- log(yt/y0)/(t1-t0)
  
  Gomp_df$r_coeff[i] <- r
  Gomp_df$Lag[i] <- t0
  
}


Gomp_df <- subset(Gomp_df, r_coeff >= 0) # removing rows that have a negative r coeff. These rows are weird time series that my recovery detector didnt do a good job on
### Running Exponential Equation for Exp_df ---------------------------------
Exp_df$Analysis_Type <- "Exp_NLS"

for (i in 1:nrow(Exp_df)) {
  site <- Exp_df$New_ID[i]
  df <- subset(rec_data, New_ID == site)
  y <- df$X.Hard_Coral_Cover
  
  yt_1 <- df$X.Hard_Coral_Cover[2:nrow(df)]
  yt <-  df$X.Hard_Coral_Cover[1:nrow(df)-1]
  t_t <- df$Date_Diff[2:nrow(df)] - df$Date_Diff[1:nrow(df)-1]
  delta_y <- yt_1 - yt
  df1 <- cbind(t_t, yt, yt_1, delta_y)
  df1 <- subset(df1, delta_y > 0)
  
  if (nrow(df1) == 1) { # if there is only one positive row in df1
    Exp_df$Analysis_Type[i] <- "Exp_Equation"
    Exp_df$r_coeff[i] <- max(log(((df$X.Hard_Coral_Cover[2:nrow(df)] + 1)/(df$X.Hard_Coral_Cover[1:nrow(df)-1] + 1)))/(df$Date_Diff[2:nrow(df)] - df$Date_Diff[1:nrow(df)-1]))
  }
  
  if (nrow(df1) >= 2) { # if there's two or more positive rowsin df1
    Exp <-  nlsLM(yt_1 ~ yt*exp(r*t_t),
                  start = list(r = 0.01),
                  trace = F,
                  control=nls.control(maxiter=1024)
    )
    Exp_df$r_coeff[i] <- as.numeric(coefficients(Exp)[1])
  }
}

Exp_df <- subset(Exp_df, r_coeff >= 0)


### Combining all data ------------------------------------------------------

Exp_Equation_Only$Lag <- NA
Exp_df$Lag <- NA
Filtered_Site_Coeff<- rbind(Exp_Equation_Only, Gomp_df, Exp_df)



Site_Coeff <- subset(Filtered_Site_Coeff, r_coeff <= 3)

nlevels(factor(Site_Coeff$Site_ID)) #1925
nrow(Site_Coeff)
nrow(subset(Site_Coeff, is.na(Initial_Algae) == F)) # 2994 rows with algae
#Comparing to old CSV
x <- read.csv("C:/Users/Andy Walker/Dropbox/Recovery Gamma/Manuscript Review 1/Filtered_Rate_Coeff_DHW_9.20.2023.csv")
nrow(subset(x, is.na(Initial_Algae) == F)) #460 rows with algae
remove(x)

#write.csv(Site_Coeff, "Filtered_Rate_Coeff_DHW_9.20.2023.csv")



# Detecting disturbances during recovery ----------------------------------

# Using start and end date



## Cyclone data has been extracted for sites using Cyclone Extractor WIP.R
#cyclones <- read.csv("C:/Users/Andy Walker/Dropbox/Disturbance History/Full_Date_Severe_Cyclones_100km.8.14.23.csv")
cyclones <- read.csv("C:/Users/Andy Walker/Dropbox/Disturbance History/Full_Date_Severe_Cyclones_200km.8.24.23.csv")

cyclones$Full_Date <- as.Date(cyclones$Full_Date, format = "%m/%d/%Y")
data.class(cyclones$Full_Date)
range(cyclones$m_per_s) # Making sure minimum is above 33
### NOTE
## Must filter cyclone data based on windspeed and distance from point
##


## DHW data
# Here I am loading in all DHW data to the same df so that in the loop later on I can subset for a site, and if there isn't data on it I can skip it without an error occurring
dhw_list <- list.files(path = "C:/Users/Andy Walker/Dropbox/Disturbance History/DHW 9.21.23/DHW by Site/", pattern="*.csv") # Grabbing all CSVs in the folder. Each file should only contain values of 4 or more
DHW_Sites <- NA
for (i in 1:length(dhw_list)) {
  x <- strsplit(dhw_list[i], "_")
  x <- x[[1]]
  x <- x[2]
  DHW_Sites <- c(DHW_Sites, x)
}
DHW_Sites <- na.omit(DHW_Sites)
remove(x, i, dhw_list)







Site_Coeff$Start_Date <- as.Date(Site_Coeff$Start_Date, format = "%Y-%m-%d")
Site_Coeff$End_Date <- as.Date(Site_Coeff$End_Date, format = "%Y-%m-%d")
data.class(Site_Coeff$Start_Date)

Site_Coeff$Mean_DHW_DR <- NA# Mean of DHW during recovery (DR)
Site_Coeff$Max_DHW_DR <- NA# Max of DHW during recovery
Site_Coeff$SD_DHW_DR <- NA # SD of DHW during recovery
Site_Coeff$Heatwaves_4_Number_DR <- NA # Number of heatwaves 4 DHW or greater during recovery
Site_Coeff$Heatwaves_8_Number_DR <- NA # Number of heatwaves 8 DHW or greater during recovery
Site_Coeff$Cyclone_Number_DR <- NA # Number of cyclones during recovery
Site_Coeff$Mean_Windspeed_DR <- NA # Mean windspeed in m per s of cyclones during recovery
Site_Coeff$Max_Windspeed_DR <- NA # Mean windspeed in m per s of cyclones during recovery
Site_Coeff$Number_of_Disturbances_DR <- NA # total number of disturbances during recovery 


DHW_Cutoff <- 0 # DHW minimum value for loop
Max_Dist <- 1 # In arc degrees. Max distance for cyclones to occur at site
blank <- Site_Coeff[0,]

for (i in levels(factor(Site_Coeff$Site_ID))) {
  
  Site <- subset(Site_Coeff, Site_ID == i)
  
  
  ## Loading in cyclones
  site_cyclone <- subset(cyclones, Site_ID == i)
  site_cyclone <- subset(site_cyclone, Distance <= Max_Dist)
  if (nrow(site_cyclone) >= 1) {
    site_cyclone$Date<- site_cyclone$Full_Date
    site_cyclone$event[1] <- 1
    if (nrow(site_cyclone) > 1) { # Marking events
      for (j1 in 2:nrow(site_cyclone)) {
        if (site_cyclone$Date[j1-1]+1 == site_cyclone$Date[j1]) {
          site_cyclone$event[j1] <- site_cyclone$event[j1-1]
        }
        if (site_cyclone$Date[j1-1]+1 != site_cyclone$Date[j1]) {
          site_cyclone$event[j1] <- site_cyclone$event[j1-1]+1
        }
      } 
    }
  }
  


  ##
  ## Loading in DHW events
  if (length(which(DHW_Sites == i)) == 1) { # If there is a file for DHW events, load it in
    site_DHW <- read.csv(paste("C:/Users/Andy Walker/Dropbox/Disturbance History/DHW 9.21.23/DHW by Site/Site_", i, "_DHW.9.21.23.csv", sep = ""))
    site_DHW$event <- NA
    site_DHW$Date <- as.Date(site_DHW$Date, format = "%Y-%m-%d")
    

    #
    #
    
    site_DHW_4 <- subset(site_DHW, DHW >= DHW_Cutoff)
    if (nrow(site_DHW_4) >= 1) {
      site_DHW_4$event[1] <- 1
      if (nrow(site_DHW_4) > 1) {
        for (j1 in 2:nrow(site_DHW_4)) {
          if (site_DHW_4$Date[j1-1]+1 == site_DHW_4$Date[j1]) {
            site_DHW_4$event[j1] <- site_DHW_4$event[j1-1]
          }
          if (site_DHW_4$Date[j1-1]+1 != site_DHW_4$Date[j1]) {
            site_DHW_4$event[j1] <- site_DHW_4$event[j1-1]+1
          }
        } 
      }
    }
  }
  
  
  for (i1 in levels(factor(Site$New_ID))) {
    x <- subset(Site, New_ID == i1)
    time_diff <- as.numeric(x$End_Date - x$Start_Date)
    
    DHW_4_new <- subset(site_DHW_4, Date >= x$Start_Date & Date < x$End_Date)
    DHW_new <- subset(site_DHW, Date >= x$Start_Date & Date < x$End_Date)
    cyclone <- subset(site_cyclone, Full_Date >= x$Start_Date & Full_Date < x$End_Date)
    

    #Extract mean, max, sd here
    x$Mean_DHW_DR <- mean(c(DHW_new$DHW, rep(0, time_diff - nrow(DHW_new))), na.rm = T)
    x$Max_DHW_DR <- max(c(DHW_new$DHW, rep(0, time_diff - nrow(DHW_new))), na.rm = T)
    x$SD_DHW_DR<- sd(c(DHW_new$DHW, rep(0, time_diff - nrow(DHW_new))), na.rm = T)
  

    
    x$Heatwaves_4_Number_DR <- nlevels(factor(subset(DHW_4_new, DHW >= 4)$event)) # Number of events with DHW 4 or more
    x$Heatwaves_8_Number_DR <- nlevels(factor(subset(DHW_4_new, DHW >= 8)$event)) # Number of events with DHW 8 or more
    x$Cyclone_Number_DR <- nlevels(factor(cyclone$event))
    x$Mean_Windspeed_DR <- mean(cyclone$m_per_s)
    x$Max_Windspeed_DR <- max(cyclone$m_per_s)
    x$Number_of_Disturbances_DR <- sum(nlevels(factor(cyclone$event)) + nlevels(factor(subset(DHW_4_new, DHW >= 4)$event)))
    
    blank <- rbind(blank, x)
  }
}


# Convert date back to character

Site_Coeff <- blank
Site_Coeff$Start_Date <- as.character(Site_Coeff$Start_Date)
Site_Coeff$End_Date <- as.character(Site_Coeff$End_Date)


#write.csv(Site_Coeff, file = "C:/Users/Andy Walker/Dropbox/Data Directory/Site_Coeff_w_DR_Vars.10.2.23.csv")




# Frequency of low intensity dhw ------------------------------------------

Site_Coeff = read.csv("C:/Users/Andy Walker/Dropbox/Recovery Gamma/Manuscript Review 1/Site_Coeff_All_Poly_vars.10.2.23.csv")

Site_Coeff$Start_Date <- as.Date(Site_Coeff$Start_Date, format = "%Y-%m-%d")
Site_Coeff$End_Date <- as.Date(Site_Coeff$End_Date, format = "%Y-%m-%d")

Site_Coeff$Low_DHW_Freq <- NA
Site_Coeff$High_DHW_Freq <-  NA

Site_Coeff$Low_DHW_Freq_Prior <- NA
Site_Coeff$High_DHW_Freq_Prior <-  NA
Site_Coeff$DHW_Freq_Prior <-  NA

Site_Coeff$High_DHW_Number_DR <- NA
Site_Coeff$Low_DHW_Number_DR <- NA

blank <- Site_Coeff[0,]
DHW_Cutoff <- 4
for (i in levels(factor(Site_Coeff$Site_ID))) {
  
  Site <- subset(Site_Coeff, Site_ID == i)

  ##
  # ## Loading in DHW events
  # site_DHW <- read.csv(paste("C:/Users/Andy Walker/Dropbox/Disturbance History/DHW 9.21.23/DHW by Site/Site_", i, "_DHW.9.21.23.csv", sep = ""))
  # site_DHW$Date <- as.Date(site_DHW$Date, format = "%Y-%m-%d")
  # 
  if (length(which(DHW_Sites == i)) == 1) { # If there is a file for DHW events, load it in
    site_DHW <- read.csv(paste("C:/Users/Andy Walker/Dropbox/Disturbance History/DHW 9.21.23/DHW by Site/Site_", i, "_DHW.9.21.23.csv", sep = ""))
    site_DHW$event <- NA
    site_DHW$Date <- as.Date(site_DHW$Date, format = "%Y-%m-%d")
    
    
    #
    #
    
    site_DHW_4 <- subset(site_DHW, DHW >= DHW_Cutoff)
    if (nrow(site_DHW_4) >= 1) {
      site_DHW_4$event[1] <- 1
      if (nrow(site_DHW_4) > 1) {
        for (j1 in 2:nrow(site_DHW_4)) {
          if (site_DHW_4$Date[j1-1]+1 == site_DHW_4$Date[j1]) {
            site_DHW_4$event[j1] <- site_DHW_4$event[j1-1]
          }
          if (site_DHW_4$Date[j1-1]+1 != site_DHW_4$Date[j1]) {
            site_DHW_4$event[j1] <- site_DHW_4$event[j1-1]+1
          }
        } 
      }
    }
  }
  
  
  for (i1 in levels(factor(Site$New_ID))) {
    x <- subset(Site, New_ID == i1)
    time_diff <- as.numeric(x$End_Date - x$Start_Date)/365
    
    ## Subsetting for prior events
    time_diff_prior <- as.numeric(x$Start_Date - as.Date("25/3/1985", format = "%d/%m/%Y"))/365
    
    DHW_old <- subset(site_DHW_4, Date < x$Start_Date)
    
    #Low_Intensity_Prior <- subset(DHW_old, DHW < 4)
    #x$Low_DHW_Freq_Prior <- nrow(Low_Intensity_Prior)/time_diff_prior
    
    #High_Intensity_Prior <- subset(DHW_old, DHW >= 8)
    #x$High_DHW_Freq_Prior <-  nrow(High_Intensity_Prior)/time_diff_prior
    
    x$DHW_Freq_Prior <-  nlevels(factor(DHW_old$event))/time_diff_prior
    
    
    
    
    
    
    ## Subsetting for events during recovery
    # 
    # DHW_new <- subset(site_DHW, Date >= x$Start_Date & Date < x$End_Date)
    # 
    # Low_Intensity <- subset(DHW_new, DHW < 4)
    # 
    # x$Low_DHW_Freq <- nrow(Low_Intensity)/time_diff
    # 
    # High_Intensity <- subset(DHW_new, DHW > 8)
    # 
    # x$High_DHW_Freq <-  nrow(High_Intensity)/time_diff
    # 
    # x$High_DHW_Number_DR <- nrow(High_Intensity)
    # x$Low_DHW_Number_DR <- nrow(Low_Intensity)
    #   
    blank <- rbind(blank, x)
  }
}

blank$Max_Windspeed_DR[which(is.infinite(blank$Max_Windspeed_DR) == T)] <- 0
#hist(blank$Low_DHW_Number_DR, breaks = 50)
hist(blank$DHW_Freq_Prior, breaks = 10)

#hist(blank$Low_DHW_Freq, breaks = 50)
#hist(blank$High_DHW_Freq, breaks = 100)

#hist(blank$Low_DHW_Freq_Prior, breaks = 50)
#hist(blank$High_DHW_Freq_Prior, breaks = 100)
#hist(blank$DHW_Freq_Prior, breaks = 50)

#blank$High_DHW_Freq[which(blank$High_DHW_Freq > 0.2)] <- NA

hist(blank$Mean_DHW_DR, breaks = 100)
hist(blank$Max_DHW_DR, breaks = 100)


#write.csv(x = blank, file = "Site_Coeff_All_Poly_vars.10.19.23.csv")

