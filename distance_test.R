################################################################################
##
## INFERRING ACTIVITY PATTERNS FROM TIMES-OF-ARRIVAL (TOA) OF ACOUSTIC SIGNALS
##
## Distance test for method validation
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

# Create results directory
result_path <- paste0("./results/")
dir.create(result_path, showWarnings = FALSE)


# 1. LOAD AND PREPARE TRACKING DATA ============================================

# Load detection and sensor data
detect <- fread("./data/range_test_detection_data.csv")
sensor <- fread("./data/range_test_sensor_data.csv")

# Load transmitter and deployment metadata
tag_metadata <- fread("./data/range_test_tag_metadata.csv")
deploy <- fread("./data/range_test_deployment_metadata.csv")

# Keep only detections assigned to a transmitter deployment
detect <- detect[deployment_id != ""]

# Assign temperature to each detection by matching 10-minute intervals
indx_t <- match(floor_date(detect$date_time, "10min"), sensor$date_time)
detect$temperature_receiver <- sensor$temperature[indx_t]

# Estimate sound speed using temperature and a salinity of 37 PSU
detect[, sound_speed := soundSpeed(temperature_receiver, sal = 37, depth = 5)]

# Extract fractional seconds from detection timestamps
detect[, fs := fractSec(date_time)]

# Calculate test time as minutes since the start of deployment D1
start_time <- deploy$time_in[deploy$deployment_id == "D1"]
detect[, test_time := difftime(date_time, start_time, units = "mins")]

# Define deployment colours for plotting
depl_list <- unique(detect$deployment_id)
col <- c("#4D4D4D","#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
         "#E6AB02", "#A6761D", "#364B99")
names(col) <- depl_list
detect[, color := col[deployment_id]]


# 2. APPLICATION OF THE METHOD TO A SINGLE TRANSMITTER =========================

# Show number of detections per transmitter and deployment
tab <- table(detect$tag_id, detect$deployment_id)
tab

# Select data from a single transmitter
id_sel <- "OPi-2004"
data <- detect[tag_id == id_sel]

# Plot raw fractional seconds
plot(fs ~ test_time, data = data, pch = 16, cex = 0.8, col = color,
     ann = FALSE, axes = FALSE)
abline(v = 0, col = "gray60", lty = 2)
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
title(xlab = "Test time (min)", line = 2.6)
title(ylab = "Seconds", line = 3)
title(main = "Raw drift (0 – 1 constrained)")
box()

# Unwrap drift using tag granularity
gran <- tag_metadata$granularity[tag_metadata$tag_id == id_sel]
data[, drift := unwrapDrift(fs, date_time, gran)]

# Fit linear model to calibration data only
lm <- lm(drift ~ test_time, data = data[deployment_id == "CAL"])

# Plot unwrapped drift and linear model
plot(drift ~ test_time, data = data, pch = 16, cex = 0.8, col = color,
     ann = FALSE, axes = FALSE)
abline(v = 0, col = "gray60", lty = 2)
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
abline(lm, col = "red")
title(xlab = "Test time (min)", line = 2.6)
title(ylab = "Seconds", line = 3)
title(main = "Unwrapped time drift")
legend("topright", legend = "Linear model", lwd = 1, col = "red", bty = "n")
box()

# Calculate TOAs as the difference between the observed fractional second values
# and the values predicted by the linear model
data[, toa := drift - predict(lm, data)]

# Convert TOA to distance using sound speed
data[, dist := toa * sound_speed]

# Offset distances relative to the mean of deployment D1 (dist = 0 m)
data[, dist := dist - mean(dist[deployment_id == "D1"])]

# Plot TOA
plot(toa ~ test_time, data = data, pch = 16, cex = 0.8, col = color,
     ann = FALSE, axes = FALSE)
abline(v = 0, col = "gray60", lty = 2)
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
title(xlab = "Test time (min)", line = 2.6)
title(ylab = "TOA (seconds)", line = 3)
title(main = "Time-of-arrival")
box()

# Plot distances
plot(dist ~ test_time, data = data, pch = 16, cex = 0.8, col = color,
     ann = FALSE, axes = FALSE)
abline(h = unique(deploy$dist), col = "gray60", lty = 2)
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
title(xlab = "Test time (min)", line = 2.6)
title(ylab = "Distance (m)", line = 3)
title(main = "Distance estimations")
box()

# Plot signal-to-noise ratio
plot(snr ~ dist, data = data[deployment_id != "CAL", ], pch = 16, cex = 0.8,
     col = color, ann = FALSE, axes = FALSE)
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
title(xlab = "Distance (m)", line = 2.6)
title(ylab = "SNR", line = 3)
title(main = "SNR vs Distance")
box()


# 3. APPLY THE METHOD TO ALL TRANSMITTERS ======================================

# Initialise empty variables for drift, TOA and distances
detect[, c("drift", "toa", "dist") := .(NA, NA, NA)]

# Loop through all transmitters
for (tag in tag_metadata$tag_id) {

  # Subset data
  indx_tag <- which(detect$tag_id == tag)
  data_sub <- detect[indx_tag, ]

  # Unwrap drift based on tag granularity
  gran <- tag_metadata$gran[tag_metadata$tag_id == tag]
  data_sub[, drift := unwrapDrift(fs, date_time, gran)]

  # Fit linear model to calibration deployment and extract residuals
  lm_t <- lm(drift ~ date_time, data = data_sub[deployment_id == "CAL"])
  data_sub[, toa := drift - predict(lm_t, data_sub)]

  # Calculate distance and apply a baseline correction based on the mean of
  # deployment D1 (dist = 0 m)
  data_sub[, dist := toa * sound_speed]
  data_sub[, dist := dist - mean(dist[deployment_id == "D1"])]

  # Add values to the original dataset
  detect$drift[indx_tag] <- data_sub$drift
  detect$toa[indx_tag] <- data_sub$toa
  detect$dist[indx_tag] <- data_sub$dist
}


# Plot results for the manuscript
png(paste0(result_path, "fig2_distance_test_results.png"),
    width = 2700, height = 2200, res = 300, pointsize = 14)
par(mfrow = c(2, 2), mar = c(4.1, 4.8, 1, 1))

# Legend labels
legend_labels <-  c("Calibration", names(col)[-1])

# A: Unwrapped drift of one transmitter
plot(drift ~ test_time, data = data, type = "n", ann = FALSE, axes = FALSE)
abline(v = 0, col = "gray60", lty = 2)
points(drift ~ test_time, data = data, pch = 16, cex = 0.8, col = color)
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, lwd = 0, lwd.ticks = 1, las = 1)
abline(lm, col = "red")
title(xlab = "Test time (min)", line = 2.6)
title(ylab = "Seconds", line = 3)
legend("topright", legend = "Linear model", lwd = 1, col = "red", bty = "n")

legend("topright", legend = legend_labels[c(1:5, NA, 6:9)],
       pch = 16, col = col[c(1:5, NA, 6:9)], inset = c(0, 0.1),
       ncol = 2, bty = "n", text.width = 40)
legend("bottomleft", legend = c("No. of tags = 1", paste("Tag ID = ", id_sel)),
       bty = "n", cex = 0.8, inset = c(-0.03, 0.02))
mtext("A", 2, las = 1, at = par("usr")[4], font = 2, cex = 1.4, adj = 4.1)

box()

# B: Distance vs test time
data_test <- detect[deployment_id != "CAL"] # Remove calibration period

plot(dist ~ test_time, data = data_test, type = "n", axes = FALSE, ann = FALSE)
abline(h = deploy$dist, col = "gray50", lty = 2)
points(dist ~ test_time, data = data_test, pch = 16,
       col = adjustcolor(data_test$color, 0.5))
axis(2, las = 1, lwd = 0, lwd.ticks = 1)
axis(1, las = 1, lwd = 0, lwd.ticks = 1)
title(xlab = "Test time (mins)", line = 2.6)
title(ylab = "Distance to receiver (m)", line = 3)
legend("bottomright", legend = legend_labels[-1], pch = 16, col = col[-1],
       ncol = 2, bg = "white", box.lty  = 0)
legend("topleft", legend = c("No. of tags = 8"),
       bty = "n", cex = 0.8, inset = c(-0.03, 0.07))
mtext("B", 2, las = 1, at = par("usr")[4], font = 2, cex = 1.4, adj = 4.1)
box()

# C: SNR vs Distance
plot(snr ~ dist, data = data_test, type = "n", axes = FALSE, ann = FALSE)
points(snr ~ dist, data = data_test, pch = 16,
       col = adjustcolor(data_test$color, 0.5))
axis(2, las = 1, lwd = 0, lwd.ticks = 1)
axis(1, las = 1, lwd = 0, lwd.ticks = 1)
title(xlab = "Distance to receiver (m)", line = 2.6)
title(ylab = "SNR", line = 3)
legend("topright", legend = legend_labels[-1], pch = 16, col = col[-1],
       ncol = 2, bty = "n")
legend("topright", legend = c("No. of tags = 8"),
       bty = "n", cex = 0.8, inset = c(0.02, 0.35))
mtext("C", 2, las = 1, at = par("usr")[4], font = 2, cex = 1.4, adj = 4.1)
box()

dev.off()


# 4. EVALUATE DISTANCE ESTIMATES ===============================================

# Accuracy: difference between estimated and true distance
indx_d <- match(data_test$deployment_id, deploy$deployment_id)
real_dist <- deploy$distance_m[indx_d]
data_test[, accuracy := dist - real_dist]
cat("Mean accuracy:", round(mean(data_test$accuracy, na.rm = TRUE), 2),
    "±", round(sd(data_test$accuracy, na.rm = TRUE), 2), "m")

# Precision: deviation from tag-specific mean distance
data_test$precision <- NA
for (test in unique(data_test$deployment_id)) {
  for (tag in unique(data_test$tag_id)) {
    indx <- data_test$deployment_id == test & data_test$tag_id == tag
    mean <- mean(data_test$dist[indx])
    data_test$precision[indx] <- data_test$dist[indx] - mean
  }
}
cat("Mean precision:", round(mean(abs(data_test$precision), na.rm = TRUE), 2),
    "±", round(sd(abs(data_test$precision), na.rm = TRUE), 2), "m")


# 5. EXPORT RESULTS TABLE ======================================================

# Mean distance and sd per tag and deployment
mean_dist <- tapply(data_test$dist,
                    list(data_test$tag_id, data_test$deployment_id),
                    mean, na.rm = TRUE)
sd_dist <- tapply(data_test$dist,
                  list(data_test$tag_id, data_test$deployment_id),
                  sd, na.rm = TRUE)
results <- paste0(round(mean_dist, 1), "±", round(sd_dist, 1))
results <- gsub("NA±NA|±NA", "", results)

# Format results as a data frame
results <- data.frame(matrix(results, ncol = ncol(mean_dist),
                             nrow = nrow(mean_dist)))
col_names <- colnames(mean_dist)
c_indx <- match(col_names, deploy$deployment_id)
col_names <- paste0(col_names, " (", deploy$distance_m[c_indx], " m)")
colnames(results) <- col_names

# Add transmitter metadata
indx_m <- match(rownames(mean_dist), tag_metadata$tag_id)
results <- cbind(manufacturer = tag_metadata$tag_manufacturer[indx_m],
                 tag_id = rownames(mean_dist),
                 tag_model = tag_metadata$tag_model[indx_m], results)

results

# Export results
fwrite(results, paste0(result_path, "table2_distance_test_results.csv"))

