## ----echo=FALSE,results='hide',message=FALSE----------------------------------
######################################################################
### Author: Michael Höhle <http://www.math.su.se/~hoehle>
### Course: The Mathematics and Statistics of Infectious Disease Outbreaks
###         (MT3002 - Summer 2020) at the Department of Mathematics,
###         Stockholm University.
###         (https://kurser.math.su.se/course/view.php?id=911)
###         given by Tom Britton and Michael Höhle
###    
### License (for the slides):
### This work is licensed under a <a rel="license"
### href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
### Attribution-ShareAlike 4.0 International License</a>.
### License (for the code):
### GNU General Public License v3 - https://www.gnu.org/licenses/gpl-3.0.en.html
###
### Description:
###  Rnw File for lecture 10
###
### History:
###  -- 09 Aug 2020 file created 
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(RColorBrewer)
library(lubridate)


## -----------------------------------------------------------------------------
options(width=90,prompt="R> ")


## ----results="hide"-----------------------------------------------------------
######################################################################
## R code generating graphs and numbers for the multivariate and scan
## statistics sections of the book chapter "Prospective Detection of
## Outbreaks" by B. Allévius and M. Höhle in the Handbook of
## Infectious Disease Data Analysis.
##
## Author: Benjamin Allévius <http://www.su.se/english/profiles/bekj9674-1.194276>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 2017-11-10
######################################################################

library(tidyverse)
library(lubridate)
library(magrittr)
library(scales)
library(surveillance)


# Munge data ===================

# Menigococcal data from surveillance
library(surveillance)
data(imdepi)

# Aggregate across event types, age and gender
meningo_cases <- imdepi$events@data %>% as_tibble %>% mutate(count = 1L)

# Grab coordinates
meningo_coords <- coordinates(imdepi$events)
meningo_cases %<>%
  mutate(x_coord = meningo_coords[, 1], y_coord = meningo_coords[, 2])

# Add dates and remove unneeded columns
meningo_cases %<>%
  mutate(day = ceiling(time),
         date = as.Date("2001-12-31") + lubridate::days(day), # Exact start date is unknown
         year = as.integer(lubridate::year(date)),
         month = as.integer(month(date))) %>%
  select(-eps.t, -eps.s, -.obsInfLength, -.sources, -.bdist, -.influenceRegion)

# Munge district-level data
district_data <- imdepi$stgrid %>% as_tibble %>%
  mutate(population = as.integer(popdensity * area))

tile_location <- tibble(tile = unique(district_data$tile),
                        location = 1:length(unique(district_data$tile)))
district_data <- left_join(district_data, tile_location, by = "tile")

# Get the total population (population is constant across time)
total_pop <- sum((district_data %>% filter(BLOCK == 1))$population)

# Monthly counts with covariates
meningo_monthly <- meningo_cases %>% select(tile, BLOCK) %>%
  mutate(count = 1L) %>%
  group_by(tile, BLOCK) %>%
  summarize(count = n()) %>%
  ungroup %>%
  right_join(district_data %>% select(tile, location, BLOCK, area, popdensity, population),
             by = c("tile", "BLOCK")) %>%
  mutate(count = ifelse(is.na(count), 0L, count)) %>%
  mutate(total_pop = total_pop) %>%
  rename(time = BLOCK) %>%
  arrange(location, time) %>%
  mutate(year = 2002 + floor((time - 1) / 12),
         month = ifelse(time %% 12 == 0, 12, time %% 12),
         date = as.Date(paste(year, month, "01", sep = "-")))



## ----warning=FALSE------------------------------------------------------------
# Scan statistics ==============================================================
suppressPackageStartupMessages(library(scanstatistics))

# Shapefile for the districts of Germany
load(system.file("shapes", "districtsD.RData", package = "surveillance"))

## ----compute_satscan, cache=TRUE----------------------------------------------
# Extract dates
dates <- (meningo_monthly %>% filter(tile == "01001") %>% select(date))$date

# Parameters for the scan statistic
zones <- coordinates(districtsD) %>% coords_to_knn(k = 15) %>% knn_zones
scan_length <- 6 # Scanning window covers last 6 months

# Define surveillance period
t2_start <- 12 * 2 + 1
t2_end <- 12 * 4
t2_length <- t2_end - t2_start + 1

# Parameters for the surveillance period
scan_start <- t2_start
scan_end <- t2_end
scan_mc <- 999
scan_alpha <- 1 / (12 * 5)

# Store replicate scan statistics
replicates <- rep(NA, (scan_end - scan_start + 1) * scan_mc)

# Kulldorff's scan statistic
scan_df <- tibble(date = dates[scan_start:scan_end],
                  score = NA,
                  crit = NA,
                  pval = NA,
                  zone = NA,
                  duration = NA,
                  relrisk_in = NA,
                  relrisk_out = NA)

relrisk_support <- seq(1, 15, by = 0.1)
prev_relrisk_prob <- rep(1, length(relrisk_support))

# Run the scan statistics
idx <- 1
for (i in scan_start:scan_end) {
  time_window <- seq(max(1, i - scan_length + 1), i, by = 1)
  obs_counts <- meningo_monthly %>% filter(time %in% time_window)

  # Kulldorff's scan statistic
  scan <- scan_pb_poisson(obs_counts, zones, n_mcsim = scan_mc)
  repl_idx <- ((idx-1) * scan_mc + 1):(idx * scan_mc)
  replicates[repl_idx] <- scan$replicates$score

  scan_df$score[idx] <- scan$MLC$score
  scan_df$crit[idx] <- quantile(replicates[1:tail(repl_idx, 1)],
                                1 - scan_alpha,
                                type = 8)
  scan_df$pval[idx] <- scanstatistics:::mc_pvalue(scan$MLC$score, replicates[1:tail(repl_idx, 1)])
  scan_df$zone[idx] <- scan$MLC$zone_number
  scan_df$duration[idx] <- scan$MLC$duration
  scan_df$relrisk_in[idx] <- scan$MLC$relrisk_in
  scan_df$relrisk_out[idx] <- scan$MLC$relrisk_out

  print(paste0("idx = ", idx))
  idx <- idx + 1
}


## ----scan_score_plot, fig.align="center", warning=FALSE-----------------------
# Plot the score of the MLC over time
scan_score_plot <- ggplot(gather(scan_df,
                                 key = "type", value = "value",
                                 score, crit)) +
  geom_line(aes(x = date, y = value, color = type, linetype = type)) +
  scale_x_date(date_breaks = "6 month",
               date_minor_breaks = "1 month",
               labels = date_format("%b-%Y")) +
  scale_color_manual(name  = "",
                     breaks=c("score", "crit"),
                     labels=c("Score", "Critical value"),
                     values = c("gray47", "black")) +
  scale_linetype_manual(name  = "",
                        breaks=c("score", "crit"),
                        labels=c("Score", "Critical value"),
                        values = c("dashed", "solid")) +
  xlab("Date") + ylab(expression(lambda[W])) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank())
scan_score_plot


## ----scan_mlc_map, fig.align="center", warning=FALSE--------------------------
# Calculate the overlap between clusters
zone_olap <- function(z1, z2) {
  length(base::intersect(z1, z2)) / length(base::union(z1, z2))
}

vec_zone_olap <- Vectorize(zone_olap, c("z1", "z2"))

zone_overlap <- rep(NA, nrow(scan_df) - 1)
for (i in 2:nrow(scan_df)) {
  zone_overlap[i-1] <- zone_olap(zones[[scan_df$zone[i - 1]]],
                                 zones[[scan_df$zone[i]]])
}

# Extract the MLC
MLC_zone <- zones[[scan_df$zone[which.max(scan_df$score)]]]

# Incidences per 100,000 people
meningo_incidence <- meningo_monthly %>%
  group_by(tile, location) %>%
  summarise(count = sum(count),
            population = population[1]) %>%
  ungroup %>%
  mutate(incidence = count * 100000 / population)

# Label the cluster
scan_clust <- districtsD@data %>%
  as_tibble %>%
  mutate(tile = as.factor(KEY),
         id = KEY,
         popdensity = POPULATION / AREA) %>%
  left_join(tile_location, by = "tile") %>%
  mutate(MLC = ifelse(location %in% MLC_zone, "Yes", "No")) %>%
  left_join(meningo_incidence, by = "tile")


# Make map data plotable
district_map <- fortify(districtsD) %>%
  as.tbl %>%
  left_join(scan_clust, by = "id") %>%
  mutate(state = substr(KEY, 1, 2))

cluster_state <- filter(district_map, MLC == "Yes" & order == 1)$state

district_map %<>%
  mutate(MLC_in_state = (state == cluster_state),
         MLC = factor(MLC, levels = c("Yes", "No")))


# Time series plot of observed scan statistic
scan_mlc_map <- ggplot(district_map %>% filter(MLC_in_state)) +
  theme_minimal() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = MLC),
               color = "black") +
  labs(x = "", y = "") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("gray37", "white")) +
  theme(legend.position = "none") +
  # theme(legend.position=c(0.9, 0.35)) +
  coord_equal()

library(ggspatial)
scan_mlc_map + 
    annotation_north_arrow(which_north = "grid", location="br") 

