################################################################################
##
## INFERRING ACTIVITY PATTERNS FROM TIMES-OF-ARRIVAL (TOA) OF ACOUSTIC SIGNALS
##
## Identify behavioural patterns in multiple species from acoustic detections
## using the TOA-based method
##
################################################################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set global options to display time with millisecond precision
options(digits.secs = 3)

# # Install "activityTOA" package from source
# install.packages("activityTOA_1.0.tar.gz", type = "source")
# 
# # The most recent version of the package can be installed from GitHub
# devtools::install_github("aspillaga/activityTOA", build_vignettes = TRUE)

# Load libraries
library(activityTOA)
library(data.table)
library(lubridate)
library(suncalc)
library(momentuHMM)
library(parallel)

# Create results directory
result_path <- paste0("./results/")
dir.create(result_path, showWarnings = FALSE)


# 1. LOAD AND PREPARE TRACKING DATA ============================================

# Load detection and sensor data
detect <- fread("./data/fish_detection_data.csv.gz")
sensor <- fread("./data/fish_sensor_data.csv.gz")

# Load transmitter and deployment metadata
tag_metadata <- fread("./data/fish_tag_metadata.csv")
deploy <- fread("./data/fish_deployment_metadata.csv")

# Check detection counts per fish and deployment
tab <- table(detect$fish_id, detect$deployment)
tab

# Convert fish_id to factor
detect[, fish_id := factor(fish_id, levels = tag_metadata$fish_id)]

# Assign temperature to each detection by matching 10-minute intervals
indx_t <- match(paste0(floor_date(detect$date_time, "10min"),
                       detect$deployment_id),
                paste0(sensor$date_time, sensor$deployment_id))
detect[, temperature_receiver := sensor$temperature[indx_t]]

# Assign salinity and depth of each deployment
indx_d <- match(detect$deployment_id, deploy$deployment_id)
detect[, salinity := deploy$salinity_psu[indx_d]]
detect[, depth := deploy$depth[indx_d]]

# Estimate sound speed using temperature, salinity and depth
detect[, sound_speed := soundSpeed(temperature_receiver, salinity, depth)]

# Number of detections per individual
table(detect$fish_id)


# 2. ESTIMATE TOAs AND DISTANCES ===============================================

# Extract fractional seconds from detection timestamps
detect[, fs := fractSec(date_time)]

# Add signal granularity to each detection to later unwrap the fs drift
detect[, gran := tag_metadata$granularity[match(fish_id, tag_metadata$fish_id)]]

# Unwrap time drift per individual and deployment using known granularity
detect[, drift := unwrapDrift(fs, date_time, granularity = gran[1]),
       by = list(fish_id, deployment_id)]

# Helper function to extract TOAs by fitting a cubic smoothing spline to the
# drift
extractTOA <- function(drift, time.stamps) {
  df <- ceiling(length(unique(as.Date(time.stamps)))*0.5)
  spl <- smooth.spline(time.stamps, drift, df = df)
  resid <- drift - predict(spl, x = as.numeric(time.stamps))$y
  return(resid)
}

# Extract TOAs
detect[, toa := extractTOA(drift, date_time), by = list(fish_id, deployment_id)]

# Calculate relative distances by multiplying TOA by the sound speed
detect[, dist := toa * sound_speed]

# Set the minimum distance to 0 m as a reference
detect[, dist := dist - min(dist), by = list(fish_id, deployment_id)]


# 2.1. PLOT TOA ESTIMATION EXAMPLE ---------------------------------------------

ind <- "CORHIP-892" # Example individual to illustrate the method
data <- detect[fish_id == ind]
data

# Define date range for X axis
date_range <- c(floor_date(min(data$date_time), "days") - 24*3600,
                ceiling_date(max(data$date_time), "days") + 24*3600)
dates <- seq(date_range[1], date_range[2], "days")

# Set point and spline colours
col_pt <- "#0C3096"
col_spl <- "#FC3508"

# Open plot device
png(paste0(result_path, "fig1_toa_estimation.png"), width = 2400, height = 2000, 
    res = 300, pointsize = 11)
par(mfrow = c(2, 2), mar = c(3.1, 4.8, 3.1, 2), xpd = FALSE)

# A. Raw fractional seconds
plot(fs ~ date_time, data = data, col = col_pt, pch = 16, cex = 0.4,
     ann = FALSE, axes = FALSE)
axis(1, at = dates, lwd = 0, lwd.ticks = 1, labels = format(dates, "%b-%d"))
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
title(ylab = "Seconds", line = 3.2)
title(main = "Raw drift (0 – 1 constrained)")
box()

# B. Unwrapped drift + spline
df <- ceiling(length(unique(as.Date(data$date_time)))*0.5)
spl <- smooth.spline(data$date_time, data$drift, df = df)
plot(drift ~ date_time, data = data, col = col_pt, pch = 16, cex = 0.4,
     ann = FALSE, axes = FALSE)
axis(1, at = dates, lwd = 0, lwd.ticks = 1, labels = format(dates, "%b-%d"))
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
lines(spl, col = col_spl, lwd = 2)
title(ylab = "Seconds", line = 3.2)
title(main = "Unwrapped drift")
legend("topright", legend = paste0("Spline (df=", df, ")"), lwd = 2,
       col = col_spl, bty = "n")
box()

# C. TOA values (spline residuals)
plot(toa ~ date_time, data = data, col = col_pt, pch = 16,  cex = 0.4,
     ann = FALSE, axes = FALSE)
axis(1, at = dates, lwd = 0, lwd.ticks = 1, labels = format(dates, "%b-%d"))
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
title(ylab = "TOA (seconds)", line = 3.2)
usr <- par("usr")
xpos <- usr[2] + diff(usr[1:2])*0.04
ydiff <- diff(usr[3:4])
par(xpd = TRUE)
text(rep(xpos), usr[3:4] + c(0.1, -0.1) * ydiff, labels = c("–", "+"),
     pos = c(1, 3), cex = 1.6, font = 2)
arrows(x0 = xpos, x1 = xpos, y0 = usr[3] + ydiff*0.07, y1 = usr[4] - ydiff*0.1,
       code = 2, length = 0.08)
text(xpos, mean(usr[3:4]), "Distance to receiver", srt = 270, xpd = TRUE,
     adj = c(0.5, -1))
title(main = "Time-of-arrival")
box()

dev.off()


# 3. HMM FOR BEHAVIOURAL STATE ESTIMATION ======================================

# Aggregate data into 5-minute intervals and calculate mean relative distances
# for each fish and deployment
detect[, time_int := floor_date(date_time, "5min")]
data_t <- detect[, .(n = .N, dist = mean(dist, na.rm = TRUE)),
                 by = list(date_time = time_int, deployment_id, fish_id)]

# Calculate step lengths
data_t[, step := c(NA, abs(diff(dist))), by = list(fish_id, deployment_id)]

# Extract sunrise/sunset times for each deployment duration
sun <- lapply(deploy$deployment_id, function(depl) {
  # Date sequence
  dates <- seq(min(as.Date(data_t$date_time[data_t$deployment_id == depl])-1),
               max(as.Date(data_t$date_time[data_t$deployment_id == depl])+1),
               "days")

  # Deployment coordinates
  coord <- deploy[deployment_id == depl, c("latitude", "longitude")]

  # Get sunrise and sunset times
  sun_t <- getSunlightTimes(dates, lat = coord$latitude, lon = coord$longitude,
                            keep = c("sunrise", "sunset"))
  sun_t$deployment_id <- depl

  return(sun_t)
})
sun <- rbindlist(sun)

# Assign day/night phase to each observation
indx_dn <-  match(paste(data_t$deployment_id, as.Date(data_t$date_time)),
                  paste(sun$deployment_id, sun$date))
data_t[, daynight := ifelse(date_time >= sun$sunrise[indx_dn] &
                              data_t$date_time <= sun$sunset[indx_dn],
                            "day", "night")]
data_t[, daynight := factor(daynight, levels = c("day", "night"))]

# Multiple HMMs will be fitted for each individual, initialising each model
# with randomised starting parameters to find the global likelihood maximum

# Number of cores to use for parallel computation
ncores <- detectCores() - 1

# Set global seed for reproducibility
set.seed(42)

# Loop over each fish ID
hmm_list <- lapply(tag_metadata$fish_id, function(ind) {

  cat("Processing individual ", ind, " (", which(tag_metadata$fish_id==ind),
      "/", nrow(tag_metadata), ")\n", sep = "")

  # Subset the data
  data_hmm <- data_t[fish_id == ind]

  # In the case of the "RAJUND-516" individual, add an "ID" variable with the
  # deployment ID so that each deployment is considered a random factor in the
  # same HMM fit
  if (ind == "RAJUND-516") {
    data_t[, ID := deployment_id]
    data_t <- data_t[order(ID, date_time)]
  }

  # Prepare data of current individual for HMM ("step" is set as covariate)
  data_hmm <- prepData(data_hmm, coordNames = NULL,covNames = "step")

  # Define the names of behavioural states
  states <- c("Resting", "Active")

  # Define the distribution of step lengths
  distrib <- list(step = "gamma")

  # Number of HMM fits to perform with different random initialisations
  n_runs <- 20

  # Use step length quantiles to guide initialisation ranges to generate
  # plausible starting parameters
  val <- quantile(data_t$step, c(0.05, 0.5, 0.95), na.rm = TRUE)

  # Pre-generate random seeds for each parallel run to ensure reproducibility
  run_seeds <- sample.int(10^6, n_runs)

  # Fit the models in parallel using mclapply
  hmm_runs <- mclapply(1:n_runs, mc.cores = ncores, function(i) {

    # Set a unique seed for this run
    set.seed(run_seeds[i])

    # Randomise starting values (mean and sd) for each state
    step_mean <- runif(2, min = c(val[1], val[2]), max = c(val[2], val[3]))
    step_sd <- runif(2, min = c(val[1], val[2]), max = c(val[2], val[3]))
    par0 <- list(step = c(step_mean, step_sd))

    # Fit HMM with a day-night covariate for state transition probabilities
    hmm <- fitHMM(data = data_hmm,
                  nbStates = length(states),
                  stateNames = states,
                  dist = distrib,
                  Par0 = par0,
                  formula = ~daynight)
    if (class(hmm)[1] == "try-error") hmm <- NULL

    # Return the fitted model with the initialisation used
    return(list(hmm, par0))

  })

  return(hmm_runs)

})

# Add individual names to the result list
names(hmm_list) <- tag_metadata$fish_id

# Extract negative log-likelihoods to select the best-fitting model (lowest
# nlog_lk value)
nlog_lk <- lapply(hmm_list, function(x) {
  unlist(lapply(x, function(y) y[[1]]$mod$minimum))
})
nlog_lk # Consistent values across runs
best_fits <- lapply(nlog_lk, which.min)

# Extract mean step length estimates for active and resting states
mean_steps <- lapply(hmm_list, function(x) {
  t(sapply(x, function(y) y[[1]]$mle$step[1, ]))
})
mean_steps # Consistent values across runs

# Select best models for each individual (smaller nlog-lk)
hmms <- lapply(seq_along(hmm_list), function(i) {
  hmm_list[[i]][[best_fits[[i]]]][[1]]
})

# Apply Viterbi decoding to assign the most likely behavioural state to each
# detection based on the best model
data_t[, state := viterbi(hmms[[unique(fish_id)]]), by = fish_id]


# 3.1. PLOT HMM RESULTS --------------------------------------------------------

# Function to shade nighttime periods in plots using sunrise/sunset data
plotSun <- function(sun, night_col = "#D3DAE2") {
  usr <- par("usr")
  rect(xleft = sun$sunset[-nrow(sun)], xright = sun$sunrise[-1],
       ybottom = usr[3], ytop = usr[4],
       col = night_col, border = "transparent")
}

# Define colours used for plotting points by behavioural state and lines
col_pt <- c("#04AA96", "#EA9D03")
col_line <- "#01415B"


# PLOT RELATIVE DISTANCES

# Open plotting device
png(paste0(result_path, "fig3_distance_individuals.png"), width = 2200,
    height = 2800, res = 300)

# Plot one row per individual
par(mfrow = c(5, 1), mar = c(2.8, 4.8, 2.3, 1.1))
plot_ind <- tag_metadata$fish_id[!grepl("RAJUND", tag_metadata$fish_id)]

for (ind in plot_ind) {

  # Subset data for the individual
  data_sub <- data_t[fish_id == ind]
  sun_sub <- sun[deployment_id == unique(data_sub$deployment_id)]
  sun_sub$date <- ymd(sun_sub$date, tz = "UTC")

  # Set limits for the Y axis
  ylim <- c(0, max(data_sub$dist, na.rm = TRUE))

  # For the "Control Tag", centre the distance in 0 (stationary tag)
  if (ind == "CONTAG-106") {
    data_sub$dist <- data_sub$dist - mean(data_sub$dist, na.rm = TRUE)
    ylim <- c(-60, 60)
  }

  # Create an empty plot and add night shading
  plot(dist ~ date_time, data = data_sub, ylim = ylim, type = "n",
       axes = FALSE, ann = FALSE)
  plotSun(sun_sub)

  # Add distance lines and points
  lines(dist ~ date_time, data = data_sub, col = col_line, lwd = 1)
  points(dist ~ date_time, data = data_sub, cex = 0.5, pch = 16,
         col = col_pt[state])

  # Plot axes
  axis(1, at = sun_sub$date, labels = format(sun_sub$date, "%b-%d"), lwd = 0,
       lwd.ticks = 1)
  axis(2, lwd = 0, lwd.ticks = 1, las = 1)
  title(ylab = "Distance to receiver (m)", line = 3.1)
  box()

  # Add legend and individual and panel labels
  legend("topright", legend = c("Resting", "Active"), col = col_pt, pch = 16,
         inset = c(0, -0.16), cex = 0.8, bg = "white", horiz = T, xpd = TRUE,
         bty = "n")
  mtext(tag_metadata$common_name[tag_metadata$fish_id == ind],
        side = 3, at = par("usr")[1], font = 2,
        adj = 0, line = 0.2, cex = 0.7)
  mtext(LETTERS[which(tag_metadata$fish_id == ind)], side = 3,
        at = par("usr")[1] - diff(par("usr")[1:2]/20), font = 2, adj = 0,
        line = 0.8, cex = 1)
}

dev.off()


# PLOT STEP LENGTHS

# Open plotting device
png(paste0(result_path, "steps_individuals.png"), width = 2200, height = 2800, 
     res = 300)

# Plot one row per individual
par(mfrow = c(5, 1), mar = c(2.8, 4.8, 2.3, 1.1))

for (ind in plot_ind) {

  # Subset data for the individual
  data_sub <- data_t[fish_id == ind]
  sun_sub <- sun[deployment_id == unique(data_sub$deployment_id)]
  sun_sub$date <- ymd(sun_sub$date, tz = "UTC")

  # Set limits for the Y axis
  ylim <- c(0, max(data_sub$step, na.rm = TRUE))

  # For the "Control Tag", increase the Y range to match the other individuals
  if (ind == "CONTAG-106") ylim <- c(-10, 20)

  # Create an empty plot and add night shading
  plot(step ~ date_time, data = data_sub, ylim = ylim, type = "n",
       axes = FALSE, ann = FALSE)
  plotSun(sun_sub)

  # Add step length lines and points
  lines(step ~ date_time, data = data_sub, col = col_line, lwd = 1)
  points(step ~ date_time, data = data_sub, cex = 0.5, pch = 16,
         col = col_pt[state])

  # Plot axes
  axis(1, at = sun_sub$date, labels = format(sun_sub$date, "%b-%d"), lwd = 0,
       lwd.ticks = 1)
  axis(2, lwd = 0, lwd.ticks = 1, las = 1)
  title(ylab = "Step length (m)", line = 3.1)
  box()

  # Add legend and individual and panel labels
  legend("topright", legend = c("Resting", "Active"), col = col_pt, pch = 16,
         inset = c(0, -0.16), cex = 0.8, bg = "white", horiz = T, xpd = TRUE,
         bty = "n")
  mtext(tag_metadata$common_name[tag_metadata$fish_id == ind],
        side = 3, at = par("usr")[1], font = 2,
        adj = 0, line = 0.2, cex = 0.7)
  mtext(LETTERS[which(tag_metadata$fish_id == ind)], side = 3,
        at = par("usr")[1] - diff(par("usr")[1:2]/20), font = 2, adj = 0,
        line = 0.8, cex = 1)
}

dev.off()


# 4. COMPARE MULTIPLE DEPLOYMENTS AND POSITIONED TRACK =========================

# 4.1. Adjust HMM to positioned trajectory -------------------------------------

track <- fread("./data/fish_rajund_trajectory.csv")

# Aggregate data in 5-minute intervals and calculate mean X and Y coordinates
track[, time_int := floor_date(date_time, "5min")]
track_t <- track[, .(n = length(tag_id),
                     x = mean(x, na.rm = TRUE),
                     y = mean(y, na.rm = TRUE)),
                 by = list(date_time = time_int, fish_id)]

# Add day/night phase to each position (use sun data from deployment_id 14)
sun_sub <- sun[deployment_id == 14]
indx_dn <- match(as.Date(track_t$date_time), sun_sub$date)
track_t$daynight <- ifelse(track_t$date_time >= sun_sub$sunrise[indx_dn] &
                             track_t$date_time <= sun_sub$sunset[indx_dn],
                           "day", "night")
track_t$daynight <- factor(track_t$daynight, levels = c("day", "night"))

# Plot trajectory
plot(y ~ x, data = track_t, asp = 1, type = "b")

# Prepare data for HMM
track_hmm <- prepData(track_t, type = "UTM", coordNames = c("x", "y"))

# State labels for HMM
states <- c("Resting", "Active")

# Parameter distributions
distrib <- list(step = "gamma", angle = "vm")

# Number of HMM fits to perform with different random initialisations
n_runs <- 20

# Use step length quantiles to guide initialisation ranges to generate
# plausible starting parameters
step_val <- quantile(track_hmm$step, c(0.05, 0.5, 0.95), na.rm = TRUE)
angle_val <- c(0, 0.1, 15)

# Pre-generate random seeds for each parallel run to ensure reproducibility
set.seed(27) # Global seed
run_seeds <- sample.int(10^6, n_runs)

hmm_tr_list <- mclapply(1:n_runs, mc.cores = ncores, function(i) {

  # Set a unique seed for this run
  set.seed(run_seeds[i])

  # Randomise starting values (mean and sd) for each state
  step_mean <- runif(2, min = c(step_val[1], step_val[2]),
                     max = c(step_val[2], step_val[3]))
  step_sd <- runif(2, min = c(step_val[1], step_val[2]),
                   max = c(step_val[2], step_val[3]))
  angle_con <- runif(2, min = c(angle_val[1], angle_val[2]),
                     max = c(angle_val[2], angle_val[3]))
  par0 <- list(step = c(step_mean, step_sd), angle = angle_con)

  # Fit HMM with a day/night covariate
  hmm <- fitHMM(data = track_hmm, nbStates = length(states),
                stateNames = states, dist = distrib,
                Par0 = par0, formula = ~daynight)
  if (class(hmm)[1] == "try-error") hmm <- NULL

  # Return the fitted model with the initialisation used
  return(list(hmm, par0))
})

# Extract negative log-likelihoods to select the best-fitting model (lowest
# nlog_lk value)
nlog_lk_tr <- unlist(lapply(hmm_tr_list, function(y) y[[1]]$mod$minimum))
nlog_lk_tr # Consistent values across runs
best_fit <- which.min(nlog_lk_tr)

# Extract mean step length estimates for active and resting states
mean_steps <- t(sapply(hmm_tr_list, function(y) y[[1]]$mle$step[1, ]))
mean_steps # Consistent values across runs

# Select best models for each individual (smaller nlog-lk)
hmm_tr <- hmm_tr_list[[best_fit]][[1]]

# Apply Viterbi decoding to assign the most likely behavioural state to each
# detection based on the best model
track_t$state <- viterbi(hmm_tr)
track_t$step <- track_hmm$step


# 4.2. Plot results ------------------------------------------------------------

png(paste0(result_path, "fig4_deployments_track.png"), width = 2200, 
    height = 2900, res = 300)
layout(mat = matrix(c(1, 2), ncol =  1), heights = c(0.8, 0.2))
par(mar = c(0.5, 4.5, 2.5, 1.1))

# PLOT RELATIVE DISTANCES

data_ru <- data_t[fish_id == "RAJUND-516"]
sun_sub <- sun[deployment_id == data_ru$deployment_id[1]]

# Maximum observed distance and number of deployments
max_dist <- round(max(data_ru$dist) * 1.1, -1)
n_depl <- length(unique(data_ru$deployment_id))

# Define Y and X limits
ymax <- max_dist * (n_depl)
xlim <- range(data_ru$date_time)

# Empty plot
plot(dist ~ date_time, data = data_ru, type = "n", axes = FALSE, ann = FALSE,
     ylim = c(-50, ymax), yaxs = "i", xlim = xlim)
plotSun(sun_sub)

# Positions for Y labels
lab <- c(0, 150, 300)

# Loop for each deployment
n <- 0
at <- numeric(0)
depl_list <- sort(unique(data_ru$deployment_id))
st_list <- deploy$station[match(depl_list, deploy$deployment_id)]

for (d in rev(depl_list)) {
  abline(h = n[n!=0], lty = 2, col = "gray30")
  sub <- data_t[deployment_id == d]
  lines(sub$dist + n ~ sub$date_time, col = col_line, lwd = 1)
  points(sub$dist + n ~ sub$date_time, cex = 0.5, pch = 16,
         col = col_pt[sub$state])
  at <- c(at, lab + n)
  n <- n + max_dist
}

# Add Y axis
axis(2, at = at, labels = rep(lab, n_depl), las = 1, lwd = 0, lwd.ticks = 1)
title(ylab = "Distance to receiver (m)", line = 3.1)

# Add legend
legend("topright", legend = c("Resting", "Active"), pch = 16, col = col_pt,
       inset = c(0.01, -0.005), cex = 0.8, horiz = TRUE, bty = "n")
text(rep(par("usr")[1], n_depl), y = seq(max_dist, ymax+max_dist, max_dist),
     labels = paste("", rev(st_list)), adj = c(-0.15, 2), cex = 0.9, font = 2)
mtext("A", side = 3,
      at = par("usr")[1] - diff(par("usr")[1:2]/10), font = 2, adj = 0,
      line = -0.8, cex = 1.2)
box()


# PLOT TRACK RESULTS

par(mar = c(3.2, 4.5, 0.5, 1.1))

track_t$y_pos <- track_t$y - min(track_t$y)
plot(y_pos ~ date_time, data = track_t, type = "n", axes = FALSE, ann = FALSE,
     xlim = xlim)
plotSun(sun_sub)

lines(y_pos ~ date_time, data = track_t, col = col_line, lwd = 1)
points(y_pos ~ date_time, data = track_t, cex = 0.5, pch = 16,
       col = col_pt[track_t$state])
axis(2, las = 1, lwd = 0, at = c(0, 150, 300, 450), lwd.ticks = 1)
axis(1, at = as.POSIXct(sun_sub$date), labels = format(sun_sub$date, "%b-%d"),
     lwd = 0, lwd.ticks = 1)
title(ylab = "Y coord. (m)", line = 3.1)
text(par("usr")[1], y = par("usr")[4],
     labels = "Trajectory",
     adj = c(-0.15, 2), cex = 0.9, font = 2)
mtext("B", side = 3,
      at = par("usr")[1] - diff(par("usr")[1:2]/10), font = 2, adj = 0,
      line = 0, cex = 1.2)
box()

dev.off()


# 4.3. Compare results based on TOA-distances and the positioned track ---------

# Add state estimated from trajectory to detection data (5-minute bins)
data_ru$state_track <- track_t$state[match(data_ru$date_time, 
                                           track_t$date_time)]

# Calculate number of observations per deployment that agree on behavioural
# classification
tr_summ <- data_ru[, .(n = .N,
                       agree_tr = sum(state == state_track, na.rm = TRUE),
                       total_tr = sum(!is.na(state) & !is.na(state_track))),
                   by = deployment_id]

# Mean number of observations per deployment
round(mean(tr_summ$n), 1)
round(sd(tr_summ$n), 1)

# Agreement percentage mean and sd
round(100*mean(tr_summ$agree_tr / tr_summ$total_tr), 1)
round(100*sd(tr_summ$agree_tr / tr_summ$total_tr), 1)
nrow(track_t)

# Extract model parameters from the HMM fitted to the positioned track
param_tr <- CIreal(hmm_tr, alpha = 0.95, parms = c("step"))
step_mean_tr <- paste0(round(param_tr$step$est["mean", ], 1), "±",
                       round(param_tr$step$se["mean", ], 1))
step_sd_tr <- paste0(round(param_tr$step$est["sd", ], 1), "±",
                     round(param_tr$step$se["sd", ], 1))

# Percentage of activity at day and night
act_day_tr <- round(100*(sum(track_t$state == 2 & track_t$daynight == "day") /
                           sum(track_t$daynight == "day")), 1)
act_night_tr <- round(100*(sum(track_t$state == 2 &
                                 track_t$daynight == "night") /
                             sum(track_t$daynight == "night")), 1)


# 5. RESULTS TABLE =============================================================

# Extract general movement metrics
mov_summary <- data_t[, .(
  n = .N,
  dist_mean = round(mean(dist, na.rm = TRUE), 1),
  dist_sd =  round(sd(dist, na.rm = TRUE), 1),
  dist_max = round(max(dist,  na.rm = TRUE), 1),
  dist_min_95 = round(quantile(dist, probs = 0.025, na.rm = TRUE), 1),
  dist_max_95 = round(quantile(dist, probs = 0.975, na.rm = TRUE), 1),
  dist_range_95 = NA,
  step_mean = round(mean(abs(step), na.rm = TRUE), 1),
  step_sd = round(sd(abs(step), na.rm = TRUE), 1),
  step_max = round(max(abs(step),   na.rm = TRUE), 1),
  step_max_95 = round(quantile(abs(step), probs = 0.975, na.rm = TRUE), 1),
  active_perc_day = round(100*sum(state == 2 & daynight == "day") /
                            sum(daynight == "day"), 1),
  active_perc_night = round(100*sum(state == 2 & daynight == "night") /
                              sum(daynight == "night"), 1),
  active_perc_total = round(100*sum(state == 2) / length(state), 1)
  ), by = list(fish_id)][order(fish_id)]

# Calculate 95% distance range
mov_summary[, dist_range_95 := dist_max_95 - dist_min_95]


# Extract step parameter estimates from HMMs (mean ± sd for each state)
hmm_results <- rbindlist(lapply(seq_along(hmms), function(i) {
  param <- CIreal(hmms[[i]], alpha = 0.95, parms = c("step"))
  data.table(
    id = names(hmms)[i],
    rest_step_mean = paste0(round(param$step$est["mean", "Resting"], 1), "±",
                            round(param$step$se["mean", "Resting"], 1)),
    rest_step_sd = paste0(round(param$step$est["sd", "Resting"], 1), "±",
                          round(param$step$se["sd", "Resting"], 1)),
    act_step_mean = paste0(round(param$step$est["mean", "Active"], 1), "±",
                            round(param$step$se["mean", "Active"], 1)),
    act_step_sd = paste0(round(param$step$est["sd", "Active"], 1), "±",
                          round(param$step$se["sd", "Active"], 1))
  )
}))

# Merge HMM results into summary table
mov_summary <- cbind(mov_summary, hmm_results)

# Print summary table
mov_summary

# Final version of "Table 2" with common names and formatted metrics
results <- data.table(
  fish_id = tag_metadata$common_name,
  mean_dist = paste0(mov_summary$dist_mean, " [", mov_summary$dist_min_95, ", ",
                     mov_summary$dist_max_95, "]"),
  max_dist = mov_summary$dist_max,
  rest_step_mean = mov_summary$rest_step_mean,
  rest_step_sd = mov_summary$rest_step_sd,
  act_step_mean = mov_summary$act_step_mean,
  act_step_sd = mov_summary$act_step_sd,
  active_perc_day = mov_summary$active_perc_day,
  active_perc_night = mov_summary$active_perc_night
)

results

# Save results to CSV
fwrite(results, file = "./results/table3_species_results.csv")


# 6. EFFECT OF BODY SIZE AND TRANSMITTER POWER ON MOVEMENT PARAMETERS ==========

# Test the effect of body size (weight) and transmitter power output on
# TOA-based step lengths. Body weight is estimated from body length using 
# specific length–weight relationships available in FishBase.


# Select fish individuals and data for the analysis (exclude control tag)
fish_ref <- tag_metadata[species != "" & !is.na(species)]
data_f <- data_t[fish_id %in% fish_ref$fish_id]

# Rename measurement names to match FishBase codes (TL, FL, WD)
types <- c(TL = "Total length", FL = "Furcal length", WD = "Disk width")
fish_ref[, type_fb := names(types)[match(length_type, types)]]


# 6.1. Extract length–weight coefficients from FishBase ------------------------

# library(rfishbase)
# 
# # Extract "a" and "b" parameter estimates for each species
# lw <- data.table(length_weight(fish_ref$species))
# 
# # Keep only the rows matching the length type used for each species
# keep <- c()
# for (i in 1:nrow(fish_ref)) {
#   keep <- c(keep, which(lw$Species == fish_ref$species[i] &
#                           lw$Type == fish_ref$type_fb[i]))
# }
# lw <- lw[keep, ]
# table(lw$Species, lw$Type)
# 
# # Estimate mean "a" and "b" values for each species
# param <- lw[, .(a = mean(a, na.rm = TRUE), b = mean(b, na.rm = TRUE)),
#             by = Species][match(fish_ref$species, Species)]
# 
# param


# 6.2. Estimate fish weights and adjust GLM ------------------------------------

# Parameters from FishBase (clean version kept for reproducibility)
param <- data.table(Species = c("Epinephelus marginatus", "Coryphaena hippurus",
                                "Gymnura altavela", "Dasyatis pastinaca", 
                                "Raja undulata"),
                    a = c(0.01102500, 0.03208417, 0.00772000, 0.02424857, 
                          0.00432800),
                    b = c(3.097725, 2.834250, 3.094933, 3.230986, 3.145000))


# Estimate body weight (kg) using length–weight formula
fish_ref$weight <- param$a * fish_ref$length_cm^param$b / 1000

# Add transmitter power and fish weight to detection data
indx <- match(data_f$fish_id, fish_ref$fish_id)
data_f[, power := fish_ref$power_output_db[indx]]
data_f[, weight := fish_ref$weight[indx]]
data_f[, species := fish_ref$species[indx]]

# Plot step lengths by species and state
boxplot(step ~ species, data = data_f[state == 1, ], outline = FALSE)
boxplot(step ~ species, data = data_f[state == 2, ], outline = FALSE)


# Fit Gamma GLMs (log link) separately for Resting (1) and Active (2) states
glm_rest <- glm(step ~ power + weight, family = Gamma(link = "log"),
                data = data_f[state == 1])
glm_activ <- glm(step ~ power + weight, family = Gamma(link = "log"), 
                 data = data_f[state == 2])

# Diagnostic plots of model fit
par(mfrow = c(2, 2))
plot(glm_rest, ask = FALSE)
plot(glm_activ, ask = FALSE)

# Summaries of model results
summary(glm_rest)
summary(glm_activ)

