#List of functions to pull data


#Code from Trinh Nguyen's OMR Data Pull 8/4/2022; fixed incorrect end_date url
calc_OMR <- function(dateStart, dateEnd = NULL, timing, extrap = F, proof = F, showQAQC = F) {
  # Do you have the two required packages?
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("zoo package is required", call. = FALSE)
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr package is required", call. = FALSE)
  }
  
  # The dateStart date has to dateStart 14 days before the specified dateStart date.
  # This is to account for NA padding during the rolling mean calculations.
  # The 13 is correct here.
  dateStartPad <- as.Date(dateStart) - lubridate::days(13)
  
  # If no dateEnd date is given, use today
  if (is.null(dateEnd)) {
    dateEnd <- lubridate::today()
  }
  
  # For hourly data, timing = "uv"; for daily it's "dv"
  if (!timing %in% c("daily", "hourly")) {
    stop("Timing variable should either be daily or hourly.")
  }
  
  if (timing %in% "daily") {
    timing <- "dv"
    columnNames <- c("agency", "site", "date", "flow", "qaqc")
  } else {
    if (timing %in% "hourly") {
      timing <- "uv"
      columnNames <- c("agency", "site", "date", "time","TZ", "flow", "qaqc")
    }
  }
  
  # Old R
  oldURL <- paste0("https://nwis.waterdata.usgs.gov/nwis/",
                   timing,
                   "?cb_72137=on&format=rdb&site_no=11313405&referred_module=sw&period=&begin_date=",
                   dateStartPad,
                   "&end_date=", 
                   dateEnd)
  
  oldR <- readr::read_table2(oldURL, 
                             col_names = columnNames,
                             skip = 31) %>%
    rename(oldR = site, oldFlow = flow)
  cat("Old River URL: \n", oldURL, "\n")
  
  # Middle R
  middleURL <- paste0("https://nwis.waterdata.usgs.gov/nwis/",
                      timing,
                      "?cb_72137=on&format=rdb&site_no=11312676&referred_module=sw&period=&begin_date=",
                      dateStartPad,
                      "&end_date=", 
                      dateEnd)
  
  middleR <- readr::read_table2(middleURL, 
                                col_names = columnNames,
                                skip = 31) %>%
    rename(middleR = site, middleFlow = flow)
  cat("Middle R URL: \n", middleURL, "\n")
  
  # Is there data for the time frame that was specified?
  if (nrow(middleR) == 0 | nrow(oldR) == 0) {
    cat("Middle R has", nrow(middleR), "rows. \n" )
    cat("Old R has", nrow(oldR), "rows. \n")
    stop("One of the gage had no data for this time frame")
  }
  
  # QAQC range
  if (isTRUE(showQAQC)) {
    
    cat("For Old River \n")
    
    oldR %>% 
      group_by(qaqc) %>% 
      slice(1, n()) %>% 
      print()
    
    cat("For Middle River \n")
    warning("Function will stop here. Disable QAQC to run further.", call. = F)
    return(middleR %>% 
             group_by(qaqc) %>% 
             slice(1, n()) %>% 
             print())
  }
  
  # Creating OMR table from the tidally filtered means
  if (timing %in% "uv") {
    omrDF <- full_join(oldR, middleR, by = c("agency", "date", "time", "TZ"))
  } else {
    if (timing %in% "dv") {
      omrDF <- full_join(oldR, middleR, by = c("agency", "date"))
    }
  }
  
  omrDF <- omrDF %>%
    transmute(date,
              week = lubridate::week(date),
              month = lubridate::month(date),
              year = lubridate::year(date),
              waterYear = year + (month > 9),
              oldFlow, middleFlow,
              # You want this to not have na.rm
              omrFlow = oldFlow + middleFlow) %>% 
    # If timing = daily, then this gives the same results
    group_by(date, week, month, year, waterYear) %>%
    summarise(oldFlow = mean(oldFlow),
              middleFlow = mean(middleFlow),
              omrFlow = mean(omrFlow),
              .groups = "drop")
  # mutate(omrFlow5 = zoo::rollmean(omrFlow, 5, na.pad = T, align = "right"),
  #        omrFlow14 = zoo::rollmean(omrFlow, 14, na.pad = T, align = "right")) %>%
  # # Filtering out the padded dates
  # filter(date >= dateStart)
  
  # Did you want to extrapolate missing values from the other?
  if (extrap) {
    if (sum(is.na(omrDF$oldFlow), is.na(omrDF$middleFlow)) == 0) {
      message("Extrapolation wasn't needed. Full data for both rivers.")
    }
    # Extrapolation occuring on DAILY levels
    oldLM <- lm(oldFlow ~ middleFlow, data = omrDF)
    middleLM <- lm(middleFlow ~ oldFlow, data = omrDF)
    # Need to predict at only the NA values.
    omrDF <- omrDF %>%
      mutate(oldRExtrap = case_when(is.na(oldFlow) ~ predict(oldLM, omrDF),
                                    TRUE ~ oldFlow),
             middleRExtrap = case_when(is.na(middleFlow) ~ predict(middleLM, omrDF),
                                       TRUE ~ middleFlow)) %>%
      mutate(omrFlowExtrap = oldRExtrap + middleRExtrap) %>%
      ungroup()
    # mutate(omrFlow5Extrap = zoo::rollmean(omrFlowExtrap, 5, na.pad = T, align = "right"),
    #        omrFlow14Extrap = zoo::rollmean(omrFlowExtrap, 14, na.pad = T, align = "right"))
    
    if (proof) {
      proofPlot <- ggplot2::ggplot(omrDF, ggplot2::aes(oldFlow, middleFlow)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", formula = "y ~ x") +
        ggplot2::labs(title = paste0("Regression of daily Middle ~ Old flow", ", adjR = ", 
                                     round(summary(middleLM)$adj.r.squared, 2)),
                      subtitle = paste0("From ", min(omrDF$date), " to ", max(omrDF$date)),
                      x = "Old River Flow (cfs)", y = "Middle River Flow (cfs)") +
        ggplot2::theme_classic(base_size = 18)
      
      print(proofPlot)
    }
  }
  
  return(omrDF)
}


#Code to pull Dayflow data from Pascale Goertler https://github.com/goertler/inundation/blob/main/R/get_dayflow.R
get_dayflow <- function(){
  # get metadata
  m <- jsonlite::fromJSON("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536.jsonld")
  
  file_table <- m$`@graph`
  
  file_table <- subset(file_table, `dct:format` == "CSV")
  
  urls <- grep("results", file_table$`dcat:accessURL`$`@id`, value = TRUE)
  
  
  # read in the data
  col_types <- readr::cols(.default = readr::col_character())
  dat <- lapply(urls, readr::read_csv, col_types=col_types, show_col_types = FALSE, progress = FALSE)
  suppressWarnings(dat <- lapply(dat, function(x){
    if (is.null(x$SJR)){
      x$SJR <- NA
    }
    if (is.null(x$SAC)){
      x$SAC <- NA
    }
    if (is.null(x$EXPORTS)){
      x$EXPORTS <- NA
    }
    return(x[, c("Date", "SJR","SAC","EXPORTS")])
  }))
  
  
  # bind data
  dayflow <- do.call(rbind, dat)
  
  # rename columns
  dayflow$Date <- lubridate::parse_date_time(dayflow$Date, orders = c("mdy", "ymd", "dmy"))
  dayflow$SJR <- as.numeric(dayflow$SJR)
  dayflow$SAC <- as.numeric(dayflow$SAC)
  dayflow$EXPORTS <- as.numeric(dayflow$EXPORTS)
  # clean names
  dayflow <- janitor::clean_names(dayflow)
  
  # remove duplicates
  i <- which(!duplicated(dayflow))
  dayflow <- dayflow[i, ]
  
  # start when data starts
  dayflow_clean <- na.omit(dayflow)
  
  return(dayflow_clean)
  
}


####function for predicting with models
add_predictions_with_ci <- function(model, newdata) {
  preds <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$pred <- preds$fit
  newdata$se <- preds$se.fit
  newdata$lower_ci <- preds$fit - 1.96 * preds$se.fit
  newdata$upper_ci <- preds$fit + 1.96 * preds$se.fit
  return(newdata)
}

add_predictions_with_ci_80 <- function(model, newdata) {
  preds <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$pred <- preds$fit
  newdata$se <- preds$se.fit
  newdata$lower_ci <- preds$fit - 1.28 * preds$se.fit
  newdata$upper_ci <- preds$fit + 1.28 * preds$se.fit
  return(newdata)
}
