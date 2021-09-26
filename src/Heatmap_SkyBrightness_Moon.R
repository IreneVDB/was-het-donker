library(dplyr)
library(purrr)
library(lubridate)
library(stringr)
library(tidyr)

locdata <- "/home/pi/data"
locinput <- "/home/pi/data"
locout <- "/home/pi/output"

# 1) Get the data when connected to raspberry pi based on station name, year, month:
get.data.from.file <- function(station, jaar, maand){
  filename <- paste0("SQM_", station, "_", jaar, "-", str_pad(maand, 2, "left", "0"))
  filepath <- paste0(file.path(locdata, station, filename), ".dat")
  
  df <- read.table(file = filepath, sep = ";", stringsAsFactors = FALSE) %>%
    'colnames<-'(c("UTCDateTime", "LocalDateTime","Temperature",
                   "Counts", "Frequency","MSAS")) %>%
    mutate(UTCDateTime = as_datetime(UTCDateTime, tz = "UTC"),
           LocalDateTime = as_datetime(LocalDateTime),
           StartDate = case_when(substr(
             LocalDateTime, 12, 13) > 12 ~ as_date(LocalDateTime),
             TRUE ~ as_date(LocalDateTime) - days(1))) 
  
  return(df)
  }

# 2) Make table with daily MSAS and moon phase: ----
make.day.table <- function(df, value = "MSAS"){
  require(lunar)
  # get all dates to fill up missing days:
  first_date <- paste0(substr(df$StartDate[1], 1, 8), "01")
  all_days <- data.frame(StartDate = seq(as_date(first_date), 
                                         as_date(paste0(substr(first_date,1,8), days_in_month(first_date))), by = "days")) 
  df_day <- df %>%
    group_by(StartDate) %>%
    summarise(mean = mean(!!as.name(value), na.rm = TRUE),
              min = min(!!as.name(value), na.rm = TRUE),
              max = max(!!as.name(value), na.rm = TRUE),
              .groups = "drop") %>%
    left_join(all_days, ., by = "StartDate") %>%
    mutate(moonphase_noon = lunar.illumination(StartDate, shift = -1),
           moonphase_midnight = lunar.illumination(StartDate, shift = 11)) %>%
    arrange(desc(StartDate))
  }

# 3) Make table with mean MSAS per time intreval (row) and day (column): ----
make.time.table <- function(df, value = "MSAS", timewindow_min = 10){
  Timebin <- data.frame(Time = paste0(rep(c(12:23, str_pad(0:11, 2, "left", "0")), each = 60 / timewindow_min), ":",
                                      rep(str_pad(seq(0, 60 - timewindow_min, by = timewindow_min), 
                                                  2, "left", "0"), 24)))
  
  first_date <- paste0(substr(df$StartDate[1], 1, 8), "01")
  all_days <- data.frame(StartDate = seq(as_date(first_date), 
                                         as_date(paste0(substr(first_date, 1, 8), 
                                                        days_in_month(first_date))), by = "days")) 
  
  df_time <- df %>%
    mutate(Time = as.character(round_date(LocalDateTime, unit = paste(
      timewindow_min, "mins")))) %>%
    group_by(Time) %>% 
    summarise(StartDate = first(StartDate),
              value = mean(!!as.name(value)), .groups = "drop") %>%
    mutate(Time = substr(Time, 12, 16)) %>%
    left_join(all_days, ., by = "StartDate") %>%
    pivot_wider(., names_from = StartDate, values_from = value) %>%
    left_join(Timebin, ., by = "Time") %>%
    'rownames<-'(.$Time) %>%
    select(-1)
  
}

# 4) Make the plot:
  # a) title and subtitle
make.header <- function(station){
  par(mai = c(0, 0, 0, 0), xpd = TRUE, col = "grey40", family = "sans", xpd = NA)
  plot(c(0,1), c(0,1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes = FALSE, asp = NA)
  text(x = 1, y = 0.8, labels = paste(station), font = 2, adj = c(1, 0.5), cex = 3)
  }

  # b) heatmap of MSAS:
make.heatmap <- function(time_df, value = "MSAS", range = c(5, 25),
                         low_col = "white", high_col = "black",
                         station){
  
  # normalized data:
  mat <- (t(as.matrix(rev(time_df))) - range[1]) / (range[2] - range[1])
  mat[is.na(mat) | mat < 0] <- 0
  mat[mat > 1] <- 1
  # colors:
  col.vector <- colorRamp(c(low_col, high_col))(as.vector(t(mat)))
  fill <- rgb(col.vector, maxColorValue = 255)
  # labels:
  maandlabel <- c("Januari", "Februari", "Maart", "April", "Mei", "Juni",
                  "Juli", "Augustus", "September", "Oktober", "November", 
                  "December")[month(as_date(colnames(time_df)[1]))]
  xlab <- "Tijd (UTC +1)"
  ylab <- paste(maandlabel, year(as_date(colnames(time_df)[1])))
  
  # plot:
  par(mai = c(1.27, 1.27, 0.6, 0.2) / 2.54, xpd = TRUE, fg = "grey40", col = "grey40", family = "sans", lwd = 0.5)
  plot(c(0, ncol(mat)), c(0, nrow(mat)), type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", axes = FALSE, asp = NA)
  symbols(x = rep(1:ncol(mat) - 0.5, nrow(mat)),
          y = rep(1:nrow(mat) - 0.5, each = ncol(mat)),
          rectangles = matrix(rep(c(1, 0.95), ncol(mat) * nrow(mat)), byrow = TRUE, ncol = 2),
          bg = fill, fg = NULL, add = TRUE, inches = FALSE)
  
  # y-axis = Day of month
  axis(side = 2, line = 0.2, at = 1:nrow(mat) - 0.5, lwd = 0.5, tcl = -0.2, #gap.axis = 0,
       labels = rev(1:nrow(mat)), las = 2, hadj = 0.5, cex.axis = 0.75, col.axis = "grey40")
  mtext(text = ylab, side = 2, line = 3, at = 0.5 * nrow(mat), font = 2, padj = 1)
  # x-axis = Time
  xlab.pos <- which(substr(rownames(time_df), 3, 5) == ":00") - 0.5
  axis(side = 1, line = 0.2, at = 1:ncol(mat) - 0.5, lwd = 0.5, tcl = -0.2,
       labels = rep("", ncol(mat)), col.axis = "grey40")
  axis(side = 1, line = 0.2, at = xlab.pos, lwd = 0, lwd.ticks = 0.8, tcl = -0.4, #gap.axis = 0, 
       padj = -1.2, labels = c(12:23, 0:11), col.axis = "grey40", cex.axis = 0.8)
  mtext(text = xlab, side = 1, line = 2.5, at = 0.5 * ncol(mat), font = 2, padj = 0)
  
  # get lat and lon values for station - if available
  loc <- read.table(file = paste0(locinput, "/loc.txt"), sep = " ", stringsAsFactors = FALSE) %>%
    'colnames<-'(c("station", "lat","lon", "color"))
 
  if(station %in% loc$station) {
    require(suncalc)
    lat <- loc[which(loc$station == station), "lat"]
    lon <- loc[which(loc$station == station), "lon"]
    
    # scaling factor (default = 24 hours over 144 units)
    hours_x <- as.numeric(difftime(as_datetime(paste0(colnames(time_df)[2], " ", rev(rownames(time_df))[1], ":00")),
                                   as_datetime(paste0(colnames(time_df)[1], " ", rownames(time_df)[1], ":00")),
                                   units = "hours"))
    x_per_h <- hours_x / (ncol(mat) - 1)
    
    # Times for morning and evening civil twilight (correct to UTC + 1):
    Sunset <- map_dfr(as_date(colnames(time_df)), suncalc::getSunlightTimes,
                      lat = lat , lon = lon, keep = "sunset") %>% 
      mutate(x_sunset = as.numeric(difftime(sunset + hours(1), as_datetime(paste0(
        date, " 11:55:00")), units = "hours")) / x_per_h)
    # Daw: use time of next day: 
    Dawn <- map_dfr(as_date(colnames(time_df)) + days(1), suncalc::getSunlightTimes,
                    lat = lat , lon = lon, keep = "dawn") %>%
      mutate(x_dawn = as.numeric(difftime(dawn + hours(1), as_datetime(paste0(
        date - days(1), " 11:55:00")), units = "hours")) / x_per_h)
    
    # add lines: 
    points(x = Sunset$x_sunset, y = ncol(time_df):1 - 0.5, type = "l", col = low_col, lty = 3, lwd = 1.35)
    points(x = Dawn$x_dawn, y = ncol(time_df):1 - 0.5, type = "l", col = low_col, lty = 3, lwd  = 1.35)
    # annotate lines:
    mtext(text = c("Zonsopkomst\n(BCMS)", "Zonsondergang\n(BCAS)"),
          at = c(Dawn$x_dawn[1], Sunset$x_sunset[1]), side = 3, line = 0.5, adj = 0.5, cex = 0.8)
    segments(x0 = c(Dawn$x_dawn[1], Sunset$x_sunset[1]), y0 = nrow(mat) + 0.1,
             x1 = c(Dawn$x_dawn[1], Sunset$x_sunset[1]), y1 = nrow(mat) + 0.5)
  }
}

  # c) moonphase plot:
make.moonphase.plot <- function(day_df, low_col = "white", high_col = "black"){
  
  # y-axis coordinates for moonphases (4 quarts):
  mooncycle <- 29.530587981
  ncycles <- ceiling(nrow(day_df) / mooncycle)
  
  moonphase <- day_df$moonphase_midnight
  ypos <- 1:nrow(day_df) - 0.5
  fit <- nls(moonphase ~ 0.5 * (1 - cos(2 * pi * (ypos - shift) / mooncycle)), 
             start = list(shift = 14))
  shift <- as.numeric(coef(fit))
 
  NM <- seq(shift - mooncycle, shift + ncycles * mooncycle, by = mooncycle)
  EK <- NM - 0.25 * mooncycle
  VM <- NM - 0.5 * mooncycle
  LK <- NM - 0.75 * mooncycle
  NM <- NM[which(NM > 0 & NM <= nrow(day_df))]
  EK <- EK[which(EK > 0 & EK <= nrow(day_df))]
  VM <- VM[which(VM > 0 & VM <= nrow(day_df))]
  LK <- LK[which(LK > 0 & LK <= nrow(day_df))]
  
  # layout for plot:
  ma_cm = c(1.27, 0.2, 0.6, 0)
  par(mai = ma_cm / 2.54, xpd = TRUE, family = "sans", lwd = 0.5, col = "grey40")
  xlim <- c(-0.2, 1.2)
  ylim <- c(0, nrow(day_df))
  px_x <- xlim[2] - xlim[1]
  px_y <- ylim[2] - ylim[1]
  
  plot(xlim, ylim, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", axes = FALSE, asp = NA)
  symbols(x = 0.5, y = 0.5 * nrow(day_df),
          rectangles = matrix(c(px_x, px_y), ncol = 2),
          bg = high_col, fg = NULL, add = TRUE, inches = FALSE)
  points(x = moonphase, y = ypos, type = "l", col = low_col, lty = 1, lwd = 2)
  
  # add circle on every quart:
  symbols(x = c(rep(0, length(NM)), rep(0.5, length(c(EK, LK))), rep(1, length(VM))),
          y = c(NM, EK, LK, VM),
          circles = rep(0.1, length(c(NM, EK, VM, LK))), 
          bg = c(rep(high_col, length(c(NM, EK, LK))), rep(low_col, length(VM))), 
          fg = low_col, lwd = 1, add = TRUE, inches = FALSE)
  # add label on every quart:
  text(x = rep(0, length(NM)), y = NM, labels = "NM", pos = 4, offset = 1.2, col = low_col, font = 2, cex = 1.2)
  text(x = rep(0.5, length(EK)), y = EK, labels = "EK", pos = 4, offset = 1.2, col = low_col, font = 2, cex = 1.2)
  text(x = rep(0.5, length(LK)), y = LK, labels = "LK", pos = 4, offset = 1.2, col = low_col, font = 2, cex = 1.2)
  text(x = rep(1, length(VM)), y = VM, labels = "VM", pos = 2, offset = 1.2, col = low_col, font = 2, cex = 1.2)
  
  # y unit size:
  in_x <- par("fin")[1] - (ma_cm[2] + ma_cm[4])/2.54
  in_y <- par("fin")[2] - (ma_cm[1] + ma_cm[3])/2.54
  in_px_x <- in_x / px_x
  in_px_y <- in_y / px_y
  
  # polygon for 1st and 3rd quart:
  polygon(x = rep(0.5 + 0.1 * cos(seq(1.5, 2.5, by = 0.1) * pi), length(EK)),
          y = rep(EK, each = 11) + 0.1 * in_px_x/in_px_y * sin(seq(1.5, 2.5, by = 0.1) * pi),
          border = NA, col = low_col) 
  polygon(x = rep(0.5 + 0.1 * cos(seq(0.5, 1.5, by = 0.1) * pi), length(LK)),
          y = rep(LK, each = 11) + 0.1 * in_px_x/in_px_y * sin(seq(0.5, 1.5, by = 0.1) * pi), 
          border = NA, col = low_col) 
  
  # x-axis for moon phase:
  axis(side = 1, line = 0.2, at = seq(0,1, by = 0.1), lwd = 0.5, tcl = -0.2,
       labels = rep("", 11), col.axis = "grey40")
  axis(side = 1, line = 0.2, at = seq(0,1, by = 0.5), lwd = 0.8, tcl = -0.4, #gap.axis = 0, 
       padj = -1.2, labels = c(0, 0.5, 1), col.axis = "grey40", cex.axis = 0.8)
  mtext(text = "Maan", side = 1, line = 2.5, at = 0.5, font = 2, padj = 0)
}

  # d) legend heatmap:
make.lgd.hm <- function(low_col = "white", high_col = "black", 
                        value = "MSAS", unit = "mag/arcsec^2", range = c(5, 25)){
  colfunc <- colorRampPalette(c(low_col, high_col))
  legend_image <- as.raster(matrix(colfunc(50), nrow = 1))
  
  par(mai = c(0.5, 2.5, 0.5, 2) / 2.54, xpd = TRUE,col = "grey40", family = "sans")
  plot(range, c(0, 1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes = FALSE, asp = NA)
  rasterImage(legend_image, range[1], 0, range[2], 1)
  axis(side = 1, line = 0.25, at = seq(range[1], range[2], by = 2.5), lwd = 0.8, tcl = -0.2, #gap.axis = 0, 
       padj = -2, labels = seq(range[1], range[2], by = 2.5), col.axis = "grey50", cex.axis = 0.8)
  text(x = range[1], y = 0.5, labels = paste0(value, "\n(", unit, ")"), 
       font = 2, cex = 1.1, pos = 2, offset = 0.5)
  
}

# e) references:
add.source <- function(){
  par(mai = c(0, 0, 0, 0), xpd = TRUE, col = "grey50", family = "sans")
  plot(c(0,1), c(0,1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes = FALSE, asp = NA)
  text(x = 1, y = 0, xpd = NA, font = 3, adj = c(1, 0),
       labels = "Bron: washetdonker.nl\nDatavisualisatie: jebentwatjemeet.nl")
}

#5 ) print plot(s):
print.plot <- function(station, jaar, maand, print,
                       timewindow_min = 10, 
                       value = "MSAS", range = c(5, 25), unit = "mag/arcsec^2",
                       low_col = "white", high_col = "black"){
  
  data <- get.data.from.file(station = station, jaar = jaar, maand = maand)
  
  day_df <- make.day.table(df = data, value = value)
  time_df <- make.time.table(df = data, value = value, timewindow_min = timewindow_min)
  
  if(is.na(print) == FALSE){
    path <- file.path(locout, station, jaar)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    file <- paste0(station, "_", jaar, "_", maand, ".", print)
    print_name <- file.path(path, file)
    # print on A4 landscape (11.64 x 8.27 inch)
    print_width <- 11.64
    print_height <- 8.27 
    if(print == "png"){
      png(print_name, width = print_width, height = print_height, res = 400, units = "in")
    } else if(print == "jpeg" | print == "jpg"){
      jpeg(print_name, width = print_width, height = print_height, res = 400, units = "in")
    } else if(print == "pdf"){
      pdf(print_name, width = print_width, height = print_height) 
    }
  }
  
  par(omi = rep(0.5, 4))
  layout(mat = matrix(c(1, 1, 2, 3, 4, 5), byrow = TRUE, ncol = 2), 
         widths = c(6, 1), heights = c(3, 28, 5))
  make.header(station = station)
  make.heatmap(time_df = time_df, station = station, value = value, range = range, 
               low_col = low_col, high_col = high_col)
  make.moonphase.plot(day_df = day_df, low_col = low_col, high_col = high_col)
  make.lgd.hm(low_col = low_col, high_col = high_col, value = value, unit = unit, range = range)
  add.source()
  if(is.na(print) == FALSE){
    dev.off()
  }
}
