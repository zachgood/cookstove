#' Load Accelerometer Data
#'
#' This function loads in accelerometer data.  It cleans the raw CSV files and
#' combines multiple files if needed. It assumes that the data is in
#' 1 or 2 .CSV files in the specific format that comes from the accelerometer.
#' Also assumes that fileName1 is data recoreded before fileName2. If there is
#' only one file, fileName2 and printHrs can be omited.
#'
#' @param fileName1 Path to first input file or filename if in current directory
#' @param fileName2 Path to second input file or filename if in current directory. DEFAULT NULL
#' @param printHrs boolean varable. True will print time recorded to concle. DEFAULT FALSE
#' @return A dataframe with 4 columns: 1. time stamps and 2., 3., 4. raw accelerometer data on 3 axis
#' @export


loadAccData <- function(fileName1, fileName2 = NULL, printHrs = FALSE){
  acceldata <- utils::read.table(file = fileName1, sep = ",", skip = 8,
                                 nrows = length(readLines(fileName1)) - 9)
  #lms <- acceldata[nrow(acceldata),1] #last measured second
  #hours <- lms/60/60 #Hours measured

  if(!is.null(fileName2)){ # If there are two files passed, combine them
    acceldata2 <- utils::read.table(file = fileName2, sep = ",", skip = 8,
                                    nrows = length(readLines(fileName2)) - 9)
    lms2 <- acceldata[nrow(acceldata2),1] #last measured second
    hours2 <- lms2/60/60 #Hours measured

    # Add seconds already observed to seconds in second file b/c it starts at 0
    acceldata2[,1] <- acceldata2[,1] + acceldata[nrow(acceldata),1]

    # Merge data together
    #lms <- lms + lms2
    #hours <- hours + hours2
    acceldata <- rbind(acceldata[,0:4], acceldata2[,0:4])
  }
  lms <- acceldata[nrow(acceldata),1] #last measured second
  hours <- lms/60/60 #Hours measured

  if(printHrs){
    cat('Seconds: ', lms)
    cat('\nHours: ', hours)
  }
  #acceldata[-nrow(acceldata), 2:4]
  colnames(acceldata) <- c("time","Ax","Ay","Az")
  acceldata[-nrow(acceldata), 1:4]
}


#' Convert Accelerometer Data to Vector Magnitude
#'
#' This function converts acceldata to 1-minute vecter magnatude data.
#' This makes the data much lighter and easier to work with or plot.
#' The function takes 3 steps: 1. Vector magnitude, calculatd as square root
#' of sum of squares, 2. Vector magnitude, after subtracting out mode and
#' taking abs value, & 3. convert to 1-minute block average.
#' Each data point of the data returned represents 1 minute of raw acc data.
#'
#' @param acceldata Dataframe object returned by loadAccData function
#' @return 'vm3' - A vector 1 minute data points.
#' @export

vecMag <- function(acceldata){
  # Vector magnitude, calculatd as square root of sum of squares
  # vm <- sqrt(acceldata[, 1]^2 + acceldata[, 2]^2 + acceldata[, 3]^2)
  vm <- sqrt(acceldata[, 2]^2 + acceldata[, 3]^2 + acceldata[, 4]^2)

  # Vector magnitude, after subtracting out mode and then making it all positive
  vm2 <- vm - as.numeric(names(sort(table(vm), decreasing = TRUE))[1])  #get mode
  vm2 <- abs(vm2)

  # Vector magnitude, after converting to 1-minute data
  # vm3 <- blockaves(x = vm2, window = window) * 12 * 60

  # points recoreded / window = num blocks --> want num blocks == minutes rec
  # window = pts rec / mins rec
  window_len <- nrow(acceldata) / (acceldata[nrow(acceldata),1]/60)
  vm3 <- accelerometry::blockaves(x = vm2, window = window_len) * 12 * 60
}


#' Convert Vector Magnitude to Movement Data
#'
#' This function converts vector magnitude data to movement data.
#' Vector magnitude data is the data returned by the vecMag function.
#' The functions applys a rolling standard deviation to the data input
#' to smooth the data and more closely resemble continuous movement.
#'
#' @param vm3 A vector of data points such as that returned by vecMag
#' @param rollLength Length of data to apply on rolling std dev. DEFAULT 5
#' @return 'movement' - A vector of smooth movement data points
#' @export

makeMovement <- function(vm3, rollLength = 5){
  `%>%` <- magrittr::`%>%`
  movement <- vm3 %>% as.matrix %>% roll::roll_sd(rollLength)
  movement <- movement[rollLength:nrow(movement)] #First 4 rows are NA b/c rolling of 5
}


#' Convert Movement Data to Compliance data
#'
#' This function converts movement data to compliance binary data.
#' 1 means the subject was wearing the bag at that minute, 0 means not.
#' Thus, each point is defined as compliant - 1, or noncompliant - 0.
#' The threshold variable set the level of septeration between which
#' points are compliant and which are noncompliant.
#'
#' NOTE: Various threshold testing was performed to select the optimal level of 400
#'
#' @param movement A vector of smooth movement points such as that returned by makeMovement
#' @param threshold Septeration point btw compliant & noncompliant. DEFAULT 400
#' @return A vector of 1s and 0s to indicate compliance
#' @export

makeCompliance <- function(movement, threshold = 400){
  #Set compliance/movement threshold at 300
  movement[movement <= threshold] <- 0
  movement[movement > 0] <- 1
  movement
}


#Function to change sequences of non-compliant points to compliant points
#mv = list of accel data indicating compliant (1) or non-compliant (2)
#seqLen = minimum number of acceptable sequeciencal points
#nonToComp = boolean indicating change of non-comp to comp (TRUE) or comp to non-comp (FALSE)


#' Change Sequences of Compliance Data to Opposite Classification
#'
#' This function will convert non-compliant points to compliant points
#' or vice versa.  It is used by the switchFake function to further increase
#' the accuracy of the data on each subject.  Thus, the data is further smoothed.
#'
#' Motivation: a subject might set the bag down for a few minutes and cause
#' the data to be non-compliant.  If the bag is only set down for a small period
#' of time, we dont want that period to be misclassified as non-compliant.
#' A similar situation exists for misclassifying compliant periods.
#'
#' @param mv List of accel data indicating compliant or non-compliant
#' @param seqLen Minimum number of acceptable sequeciencal points. DEFAULT 10
#' @param nonToComp Boolean indicating change of non-comp to comp - TRUE or comp to non-comp - FALSE
#' @return A vector of 1s and 0s without sequences of 1s or 0s smaller than 'seqLen' argument
#' @export

switchCompliance <- function(mv, seqLen = 10, nonToComp = TRUE){
  `%>%` <- magrittr::`%>%`
  #Get list of repeated observations
  r <- rle(mv)
  #Select the observations that are comp/non-comp and sequentially less than seqLen
  if(nonToComp) w <- which(!r$values & r$lengths < seqLen)
  else w <- which(r$values & r$lengths < seqLen)
  #Make a list of indexs that should be switched
  switch <- lapply(w, function(x){
    before <- sum(r$lengths[1:(x-1)])
    (before+1):(before+ r$lengths[x])
  }) %>% unlist
  #Switch those indexes
  mv[switch] = (nonToComp*1)
  mv
}


#Function to remove "fake" compliance/non-compliance observations
#mv: list of accel data indicating compliant (1) or non-compliant (2)
#seqLen: minimum number of acceptable sequeciencal points
#nonCompFirst: boolean indicating if "fake" non-compliance should be switched first (TRUE) or second (FALSE)


#' Remove Misclassified Compliance Data Points
#'
#' This function will convert misclassified non-compliant points to compliant
#' points and misclassified compliant points to non-compliant.
#' It uses the switchCompliance function twice: once with nonToComp = TRUE &
#' once with nonToComp = FALSE.  nonToComp = TRUE is preformed first to error
#' on the side of compliance but this can be controled by the 'nonCompFirst' boolean
#'
#' Motivation: a subject might set the bag down for a few minutes and cause
#' the data to be non-compliant.  If the bag is only set down for a small period
#' of time, we dont want that period to be misclassified as non-compliant.
#' A similar situation exists for misclassifying compliant periods.
#'
#' @param mv List of accel data indicating compliant or non-compliant
#' @param seqLen Minimum number of acceptable sequeciencal points. DEFAULT 10
#' @param nonCompFirst Boolean indicating order of switching compliance type
#' @return A vector of 1s and 0s without sequences of 1s AND 0s smaller than 'seqLen' argument
#' @export

switchFake <- function(mv, seqLen = 10, nonCompFirst = TRUE){
  `%<>%` <- magrittr::`%<>%`
  mv %<>% switchCompliance(seqLen, nonCompFirst)
  mv %<>% switchCompliance(seqLen, !nonCompFirst)
  mv
}


# This function converts raw accelerometer data to compliance data

#' Convert Raw Accelerometer Data to Compliance Data
#'
#' This function will convert Raw Accelerometer Data to Compliance Data
#' using all of the accelerometer functions using their default arguements
#'
#' @param fileName1 Path to first input file or filename if in current directory
#' @param fileName2 Path/name of second input file if present. DEFAULT NULL
#' @param cut24 boolean varable. True will cut the data at 24 hours before returning. DEFAULT FALSE
#' @param attachHours boolean varable. True will include number of raw hours recorded in returned value. DEFALUT FALSE
#' @return A vector of 1s and 0s to indicate compliance and non-compliance. if attachHours = TRUE; list with 2 elements: 1. raw hours recorede, 2. comp data
#' @export

accToComp <- function(fileName1, fileName2 = NULL, cut24 = FALSE, attachHours = FALSE){
  acceldata = loadAccData(fileName1, fileName2)
  lms <- acceldata[nrow(acceldata),1] #last measured second
  hours <- lms/60/60 #Hours measured
  if(cut24) { # Cut at 24 hours
    if(acceldata[nrow(acceldata),1] > 86700) { # At least 24 hours recorded
      acceldata <- acceldata[1:min(which(acceldata[,1] > 86700)),]
    } else {
      cat("Unable to cut acc data at 24 hours: less than 24 hours recorded")
    }
  }
  vm3 = vecMag(acceldata)
  movement = makeMovement(vm3)
  mv = makeCompliance(movement)
  compData = switchFake(mv)
  if(attachHours){
    return(list(hours,compData))
  } else {
    return(compData)
  }
}
