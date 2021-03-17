##----------------------
## Analysis of three pairing methodologies
## SM 2020-03-17
##----------------------

## Load required packages
library(RODBC)
library(ggplot2)
library(reshape)
library(reshape2)
library(mgcv)
library(grid)
library(viridis)

## Connect to database
con <- odbcConnectAccess(
  database_path
  )

## ---Foyle Fisher T90 90mm and diamond 80mm 120mm SMP-----------------------------------

dat_lengths <- sqlQuery(con, "SELECT *
                               FROM 201905_T90_90mm_V_T0_80mm_120mmSMP_lengths;")

dat_weights <- sqlQuery(con, "SELECT *
                               FROM 201905_T90_90mm_V_T0_80mm_120mmSMP_weights;")

dat_haul <- sqlQuery(con, "SELECT *
                            FROM 201905_T90_90mm_V_T0_80mm_120mmSMP_haul;")

## Work on valid hauls

dat_lengths <- merge(dat_lengths,
                     dat_haul[, c("tow_id", "tow_valid")])
dat_lengths <- subset(dat_lengths, tow_valid == 1)
dat_weights <- merge(dat_weights,
                     dat_haul[, c("tow_id", "tow_valid")])
dat_weights <- subset(dat_weights, tow_valid == 1)


# ## Work on the following species 
spp <- c("Cod", "Whiting", "Haddock", "Megrim")
# ## Assign length of vector of species to be worked on
n_spp <- length(spp)

## Subset data to species of interest named above
dat_lengths <- subset(dat_lengths, species %in% spp)
dat_weights <- subset(dat_weights, species %in% spp)

## Get subsratios
dat_weights$subsratio <- with(dat_weights, subsample / weight)

## If no subsratio is given, assume subsample to have subsratio of 1
dat_weights$subsratio[is.na(dat_weights$subsratio) |
                        dat_weights$subsratio == 0] <- 1

## Drop unused levels
dat_lengths <- droplevels(dat_lengths)
dat_weights <- droplevels(dat_weights)

## Merge subsratios from weights dataframe
dat <- merge(dat_lengths,
             dat_weights[, c(
               "tow_id", "species", "compartment", "subsratio")])
dat$tow_valid <- NULL

## Calculate raised counts
dat$raising_factor <- with(dat, 1/subsratio)
dat$count_raised <- with(dat, count_n * raising_factor)

## Reshape the data
## Measured counts
vars2keep <- c(
  "species", "compartment", "length_id", "tow_id", "count_n"
)

## Long data by haul
dat_long <- melt(data = dat[, vars2keep], 
                 id.vars = c("species", "compartment", "length_id", "tow_id"))

## Wide format of counts at length in compartment by haul
dat_wide_haul <- cast(
  dat_long, species + length_id + tow_id ~ compartment + variable, sum
)

## Factor species
dat_wide_haul$species <- factor(as.character(dat_wide_haul$species), levels = spp)

## Prepare subsratio dataframe
vars2keep_subs <- c("species", "compartment", "tow_id", "subsratio")

## Reshape the data to get pivot table views of subsratio
## for species per compartment and haul. 
dat_long_subs <- unique(melt(data = dat[, vars2keep_subs], 
                             id.vars = c("species", "compartment", "tow_id")))

dat_wide_subs <- cast(dat_long_subs, species + tow_id ~ compartment)

## No whiting found in Experimental tows 1 & 2, replace NA's with 1
dat_wide_subs$Experimental[is.na(dat_wide_subs$Experimental == 0)] <- 1

## Rename to standardize
colnames(dat_wide_subs)[colnames(dat_wide_subs) == "Control"] <- "Control_subsratio"
colnames(dat_wide_subs)[colnames(dat_wide_subs) == "Experimental"] <- "Experimental_subsratio"


## Bootstrap GAM with within haul replacement 
run_boot_fn <- function(sppi, n_boot){
  ## function runs the bootstrap for a given species with
  ## a given number of bootstrap replicates defined by n_boot
  ## predicted lengths
  length_range <- range(sub_dat$length_id)
  length_sequence <- seq(length_range[1], length_range[2])
  ## Prepare dfs to hold predictions
  boot_pred <- matrix(NA, nrow = n_boot, ncol = length(length_sequence))
  pred_df <- data.frame(
    length_id = length_sequence)
  ## run the bootstrap
  for(j in 1:n_boot){
    tows <- as.integer(unique(dat$tow_id))
    J <- length(tows)  
    ## experimental tows - NB replace = TRUE here
    exp_hauls <- as.integer(sample(tows, length(tows), replace = TRUE))
    ## control tows
    con_hauls <- as.integer(sample(tows, length(tows), replace = TRUE))
    ## bootstrap dataframe
    boot_dat <- NULL
    ## Go fish with above sampled tow pairs
    for(k in 1:J){
      ## Experimental
      ## Subset to sampled tow 
      exp_tow_boot <- subset(sub_dat, tow_id == exp_hauls[k])
      ## Remove control counts
      exp_tow_boot$Control_count_n <- NULL
      ## Remove lengths not found in given tow experimental compartment
      ## This gives us the lengths available to sample from and the number
      ## of samples to take from them
      ## resample lengths within the haul
      ## Available lengths to sample from
      exp_lengths_boot <- with(exp_tow_boot, rep(length_id, times = Experimental_count_n))  
      ## Sample lengths the same number that were sampled in that tow
      resamp_exp_lengths_boot <- sample(exp_lengths_boot, 
                                        replace = TRUE, 
                                        size = sum(exp_tow_boot$Experimental_count_n))
      if(length(resamp_exp_lengths_boot) == 0){
        resamp_exp_lengths_boot <- NA
      }
      ## Standardize data names
      exp_resampled <- as.data.frame(table(resamp_exp_lengths_boot, useNA = "always"))
      exp_resampled$length_id <- as.numeric(as.character(exp_resampled$resamp_exp_lengths_boot))
      exp_resampled$Experimental_count_n <- exp_resampled$Freq
      ## Assign tow id
      exp_resampled$tow_id <- exp_hauls[k]
      ## Assign subsratio
      exp_resampled <- merge(exp_resampled, 
                             subset(dat_wide_subs, species == sppi)[, c("tow_id", "Experimental_subsratio")])
      ## Raise counts
      exp_resampled$Experimental_raising_factor <- 
        with(exp_resampled, 1/Experimental_subsratio)
      exp_resampled$Experimental_count_raised <- 
        round(with(exp_resampled, Experimental_count_n * Experimental_raising_factor))
      ## Repeat process for control tow
      ## Subset to sampled tow 
      con_tow_boot <- subset(sub_dat, tow_id == con_hauls[k])
      ## Remove experimental counts
      con_tow_boot$Experimental_count_n <- NULL
      ## Remove lengths not found in given tow control compartment
      ## This gives us the lengths available to sample from and the number
      ## of samples to take from them
      ## Available lengths to sample from
      con_lengths_boot <- with(con_tow_boot, rep(length_id, times = Control_count_n))
      ## Sample lengths the same number that were sampled in that tow
      resamp_con_lengths_boot <- sample(con_lengths_boot, 
                                        replace = TRUE, 
                                        size = sum(con_tow_boot$Control_count_n))
      if(length(resamp_con_lengths_boot) == 0){
        resamp_con_lengths_boot <- NA
      }        
      ## Standardize data names
      con_resampled <- as.data.frame(table(resamp_con_lengths_boot, useNA = "always"))
      con_resampled$length_id <- as.numeric(as.character(con_resampled$resamp_con_lengths_boot))
      con_resampled$Control_count_n <- con_resampled$Freq
      ## Assign tow id
      con_resampled$tow_id <- con_hauls[k]
      ## Assign subsratio
      con_resampled <- merge(con_resampled, 
                             subset(dat_wide_subs, species == sppi)[, c("tow_id", "Control_subsratio")])                             
      ## Raise counts
      con_resampled$Control_raising_factor <- 
        with(con_resampled, 1/Control_subsratio)
      con_resampled$Control_count_raised <- 
        round(with(con_resampled, Control_count_n * Control_raising_factor))
      ## Merge dataframes
      merge_resampled_df <- merge(exp_resampled[, c("length_id", "Experimental_count_raised")], 
                                  con_resampled[, c("length_id", "Control_count_raised")], all = TRUE)
      merge_resampled_df <- subset(merge_resampled_df, !is.na(length_id))
      ## NA's found for lengths not found in a compartment, replace with count 0
      merge_resampled_df[is.na(merge_resampled_df)] <- 0
      ## Populate boot dataframe
      if(nrow(merge_resampled_df) > 0){
        boot_dat <- rbind(boot_dat, merge_resampled_df)
      }
    }
    # Make matrix of counts
    spp_count_mat <- as.matrix(boot_dat[, c(
      "Experimental_count_raised", "Control_count_raised")])
    ## fit GAM
    if(length(unique(boot_dat$length_id)) >= 5){
      gam_fit <- gam(spp_count_mat ~ s(length_id, k = 5),
                     data = boot_dat, family = binomial,
                     method = "ML")
    }else{
      gam_fit <- gam(spp_count_mat ~
                       s(length_id,
                         k = length(unique(boot_dat$length_id))),
                     data = boot_dat, family = binomial,
                     method = "ML")
    }
    ## get predictions overall
    pred <- predict(gam_fit, newdata = pred_df)
    boot_pred[j,] <- plogis(pred) 
    
  }
  
  pred_df <- data.frame(species = sppi,
                        length_id = length_sequence,
                        proportion = apply(boot_pred, 2, mean, na.rm = TRUE),
                        upper = apply(boot_pred, 2, quantile, probs = 0.975, na.rm = TRUE),
                        lower = apply(boot_pred, 2, quantile, probs = 0.025, na.rm = TRUE),
                        se_fit = apply(qlogis(boot_pred), 2, sd, na.rm = TRUE),
                        method = "Random")
  
  
  return(pred_df)
}

## Prepare dataframe to hold proportional data
boot_all_pred_df <- NULL

## Run 10000 times
for(i in 1:n_spp){
  sub_dat <- subset(dat_wide_haul, species == spp[i])
  print(paste("Bootstrapping", spp[i]))
  boot_all_pred_df <- rbind(boot_all_pred_df,
                            run_boot_fn(sppi = spp[i], n_boot = 10000))
}

## Overall analysis for proportion plot
vars2keep_overall <- c(
  "species", "compartment", "length_id", "count_raised", "tow_id", "count_n", "subsratio"
)

dat_long_overall <- melt(data = dat[, vars2keep_overall], 
                         id.vars = c("species", "compartment", "length_id", "tow_id"))

## Cast dataframe for overall analysis
dat_wide_overall <- cast(dat_long_overall, species + length_id ~ compartment + variable, sum)

## assign total raised counts to dataframes
dat_wide_overall$count_raised_total <- with(
  dat_wide_overall, Control_count_raised + Experimental_count_raised
)

## assign proportion in alternate compartment to dataframes
dat_wide_overall$prop_experimental <- with(
  dat_wide_overall, Experimental_count_raised / count_raised_total
)

## Plot
# png("double_bootstrap.png", width = 190, height = 140, units='mm', res = 300)
ggplot(dat_wide_overall, aes(x = length_id, y = prop_experimental)) +
  # geom_point(aes(size = log(n)), pch = 1) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  ylim(c(0, 1)) +
  geom_ribbon(data = boot_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "#48494B") +
  geom_line(data = boot_all_pred_df, aes(y = proportion),
            colour = "black") +
  geom_hline(yintercept = 0.5, linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Length (cm)") + ylab("Proportion retained in test gear") +
  facet_wrap(~ species, ncol = 2, scales = "free_x") 
# dev.off()




##-----------------

## Twin rig paired 

## Select variables to keep
vars2keep <- c(
  "species", "compartment", "length_id", "count_raised", "tow_id", "count_n", "subsratio"
)

## Foyle Fisher T90 90mm and diamond 80mm 120mm SMP

## Work on valid hauls

twin_dat_lengths <- merge(dat_lengths, 
                          dat_haul[, c("tow_id", "tow_valid")])
twin_dat_lengths <- subset(dat_lengths, tow_valid == 1)
twin_dat_weights <- merge(dat_weights, 
                          dat_haul[, c("tow_id", "tow_valid")])
twin_dat_weights <- subset(dat_weights, tow_valid == 1)

## Get subsratios

twin_dat_weights$subsratio <- with(twin_dat_weights, subsample / weight)

## If no subsratio is given, assume subsample to have subsratio of 1

twin_dat_weights$subsratio[is.na(twin_dat_weights$subsratio) | 
                             twin_dat_weights$subsratio == 0] <- 1

## Drop unused levels

twin_dat_lengths <- droplevels(twin_dat_lengths)
twin_dat_weights <- droplevels(twin_dat_weights)

## Merge subsratios from weights dataframe

twin_dat <- merge(twin_dat_lengths, 
                  twin_dat_weights[, c(
                    "tow_id", "species", "compartment", "subsratio")])
twin_dat$tow_valid <- NULL

## Subset data to species of interest named above

twin_dat <- subset(twin_dat, species %in% spp)

## use own calculated raising factors for totals

## twin
twin_dat$raising_factor <- with(twin_dat, 1/subsratio)
twin_dat$count_raised <- with(twin_dat, count_n * raising_factor)

## Reshape the data to get pivot table views of raised count
## and proportions per length class. 

twin_dat_long <- melt(data = twin_dat[, vars2keep], 
                      id.vars = c("species", "compartment", "length_id", "tow_id"))

## check the numbers for duplicates. 
## If there are no duplicates, <0 rows> will be displayed

twin_check <- cast(twin_dat_long, species + length_id + tow_id ~ 
                     compartment + variable, length)
twin_check[apply(twin_check, 1, FUN = function(z){any(z[-(1:3)] > 1)}),]

## Cast dataframe for overall analysis

twin_dat_wide <- cast(twin_dat_long, species + length_id ~ compartment + variable, sum)

## assign total raised counts to dataframes

twin_dat_wide$count_raised_total <- with(
  twin_dat_wide, Control_count_raised + Experimental_count_raised
)

## assign proportion in alternate compartment to dataframes

twin_dat_wide$prop_experimental <- with(
  twin_dat_wide, Experimental_count_raised / count_raised_total
)

## Factor species and tow id

twin_dat_wide$species <- factor(as.character(twin_dat_wide$species), levels = spp)

## Paired
twin_pair_dat_wide <- cast(twin_dat_long, species + length_id + tow_id ~ 
                             compartment + variable, sum)

## CM CHECK
## replace 0 subsratios with 1
twin_pair_dat_wide$Control_subsratio[twin_pair_dat_wide$Control_subsratio == 0] <- 1
twin_pair_dat_wide$Experimental_subsratio[
  twin_pair_dat_wide$Experimental_subsratio == 0] <- 1


## assign total raised counts

twin_pair_dat_wide$count_raised_total <- with(
  twin_pair_dat_wide, Control_count_raised + Experimental_count_raised
)

## assign proportion in experimental

twin_pair_dat_wide$prop_experimental <- with(
  twin_pair_dat_wide, Experimental_count_raised / count_raised_total
)

## Factor species and tow id

twin_pair_dat_wide$species <- factor(as.character(twin_pair_dat_wide$species), 
                                     levels = spp)
twin_pair_dat_wide$tow_id <- factor(twin_pair_dat_wide$tow_id)

## create container for predictions
twin_pair_all_pred_df <- NULL
## Loop through species fitting GAM and
## populating predictions dataframe
for(i in 1:length(spp)){
  #print(i)
  sppi <- spp[i]
  dat0 <- subset(twin_pair_dat_wide, species == sppi)
  spp_count_mat <- as.matrix(dat0[, c("Experimental_count_n", "Control_count_n")])
  offset_gam <-
    log(dat0[, "Experimental_subsratio"] / dat0[, "Control_subsratio"])
  # GAM
  gam_fit <- gam(spp_count_mat ~ s(length_id, k = 5) +
                   offset(offset_gam) + s(tow_id, bs="re"),
                 data = dat0, family = binomial,
                 method = "ML")
  ## get predictions overall
  length_range <- range(dat0$length_id)
  pred_df <- data.frame(
    species = sppi,
    length_id = length_range[1]:length_range[2],
    offset_gam = 0, tow_id = dat0$tow_id[1])
  pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE, exclude = "s(tow_id)")
  pred_df$proportion <- plogis(pred$fit)
  pred_df$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  pred_df$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  pred_df$se_fit <- pred$se.fit
  ## get predictions by haul
  twin_pair_all_pred_df <- rbind(twin_pair_all_pred_df, pred_df)
  rm(list = c("pred_df", "gam_fit", "pred"))
}

plot_twin_pair_prop_gam <- ggplot(twin_dat_wide,
                                  aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = twin_pair_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "darkblue") +
  geom_line(data = twin_pair_all_pred_df, aes(y = proportion),
            colour = "darkblue")+ ggtitle("GAM \n predicted proportions") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Length (cm)") + ylab("Proportion in T90 90mm") +
  ggtitle("Foyle Fisher T90 90mm and diamond 80mm 120mm SMP \n twinrig paired") +
  theme(plot.title = element_text(hjust = 0.5))
print(plot_twin_pair_prop_gam)







## -------------------------------

## Alt paired

## Prepare data for alternate pairs

twin_alt_dat_long <- twin_dat_long
twin_alt_dat_long$tow_id_alt <- twin_alt_dat_long$tow_id

## C1 - E2
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 2 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 1
## C2 - E1
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 1 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 2
## C3 - E5
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 5 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 3
## C5 - E3
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 3 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 5
## C6 - E7
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 7 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 6
## C7 - E6
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 6 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 7
## C8 - E9
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 9 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 8
## C9 - E8
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 8 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 9
## C10 - E11
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 11 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 10
## C11 - E10
twin_alt_dat_long$tow_id_alt[twin_alt_dat_long$tow_id == 10 & 
                               twin_alt_dat_long$compartment == "Experimental"] <- 11

## alt
twin_alt_dat_wide <- cast(twin_alt_dat_long, species + length_id + tow_id_alt ~ 
                            compartment + variable, sum)

## replace 0 subsratios with 1
twin_alt_dat_wide$Control_subsratio[twin_alt_dat_wide$Control_subsratio == 0] <- 1
twin_alt_dat_wide$Experimental_subsratio[
  twin_alt_dat_wide$Experimental_subsratio == 0] <- 1


## assign total raised counts to dataframes

twin_alt_dat_wide$count_raised_total <- with(
  twin_alt_dat_wide, Control_count_raised + Experimental_count_raised
)

## assign proportion in alternate compartment to dataframes

twin_alt_dat_wide$prop_experimental <- with(
  twin_alt_dat_wide, Experimental_count_raised / count_raised_total
)

## Factor species and tow id

twin_alt_dat_wide$species <- factor(as.character(twin_alt_dat_wide$species), levels = spp)
twin_alt_dat_wide$tow_id_alt <- factor(twin_alt_dat_wide$tow_id_alt)

## create container for predictions
twin_alt_all_pred_df <- NULL
## Loop through species fitting GAM and
## populating predictions dataframe
for(i in 1:length(spp)){
  #print(i)
  sppi <- spp[i]
  dat0 <- subset(twin_alt_dat_wide, species == sppi)
  spp_count_mat <- as.matrix(dat0[, c("Experimental_count_n", "Control_count_n")])
  offset_gam <-
    log(dat0[, "Experimental_subsratio"] / dat0[, "Control_subsratio"])
  # GAM
  gam_fit <- gam(spp_count_mat ~ s(length_id, k = 5) +
                   offset(offset_gam) + s(tow_id_alt, bs="re"),
                 data = dat0, family = binomial,
                 method = "ML")
  ## get predictions overall
  length_range <- range(dat0$length_id)
  pred_df <- data.frame(
    species = sppi,
    length_id = length_range[1]:length_range[2],
    offset_gam = 0, tow_id_alt = dat0$tow_id_alt[1])
  pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE, exclude = "s(tow_id_alt)")
  pred_df$proportion <- plogis(pred$fit)
  pred_df$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  pred_df$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  pred_df$se_fit <- pred$se.fit
  ## get predictions by haul
  twin_alt_all_pred_df <- rbind(twin_alt_all_pred_df, pred_df)
  rm(list = c("pred_df", "gam_fit", "pred"))
}

plot_twin_alt_prop_gam <- ggplot(twin_dat_wide,
                                 aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  geom_ribbon(data = twin_alt_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "darkblue") +
  geom_line(data = twin_alt_all_pred_df, aes(y = proportion),
            colour = "darkblue")+ ggtitle("GAM \n predicted proportions") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Length (cm)") + ylab("Proportion in T90 90mm") +
  ggtitle("Foyle Fisher T90 90mm and diamond 80mm 120mm SMP \n alternate tow paired") +
  theme(plot.title = element_text(hjust = 0.5))
print(plot_twin_alt_prop_gam)




##-------------

## Plots!

## Set plot layout and margins
par(mfrow = c(4, 3), mar = c(2, 2.2, 1, 1), oma = c(2, 3.5, 1, .5))
library(ggplot2); theme_set(theme_bw())
## Twinrig paired
plot_twin_pair_prop_gam_fig <- ggplot(twin_dat_wide,
                                      aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = twin_pair_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              # alpha = 0.2, fill = "#8DC73F60") +
              alpha = 0.2, fill = "#48494B") +
  geom_line(data = twin_pair_all_pred_df, aes(y = proportion),
            colour = "black") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none") +
  # ggtitle("Twin-rig paired") +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0.5))

## Alternate paired
plot_twin_alt_prop_gam_fig <- ggplot(twin_dat_wide,
                                     aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = twin_alt_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              # alpha = 0.2, fill = "#8DC73F60") +
              alpha = 0.2, fill = "#48494B") +
  geom_line(data = twin_alt_all_pred_df, aes(y = proportion),
            colour = "black") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none") +
  # ggtitle("Alternate paired") +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0.5))

## Bootstrap paired
plot_twin_boot_prop_gam_fig <- ggplot(twin_dat_wide, 
                                      aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = boot_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              # alpha = 0.2, fill = "#8DC73F60") +
              alpha = 0.2, fill = "#48494B") +
  geom_line(data = boot_all_pred_df, aes(y = proportion),
            colour = "black") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none") +
  # ggtitle("Bootstrap random paired") +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0.5))

# png("Fig4_twin_vs_alt_vs_bootstrap_GAM_SM_10k.png", width = 190, height = 210, units='mm', res = 300)

grid.arrange(plot_twin_pair_prop_gam_fig, plot_twin_alt_prop_gam_fig, plot_twin_boot_prop_gam_fig, 
             ncol=3, nrow = 1, 
             bottom = "Length (cm)", left = "Proportion retained in test gear")
# top = "Overall proportion retained in T90 90mm")

# dev.off()




## Standard error plots

twin_pair_all_pred_df$pairing_method <- "Simultaneous" 
twin_alt_all_pred_df$pairing_method <- "Alternate"
boot_all_pred_df$pairing_method <- "Random bootstrap"

twin_pair_all_pred_df <- subset(twin_pair_all_pred_df, select = -c(offset_gam, tow_id))
twin_alt_all_pred_df <- subset(twin_alt_all_pred_df, select = -c(offset_gam, tow_id_alt))
boot_all_pred_df <- subset(boot_all_pred_df, select = -c(method))
se_fit_df <- rbind(twin_pair_all_pred_df, twin_alt_all_pred_df, boot_all_pred_df)



# png("se_draft_SM_10k.png", width = 190, height = 140, units='mm', res = 300)
ggplot(se_fit_df, aes(x = length_id, y = se_fit)) + 
  geom_line(aes(color = pairing_method)) + 
  facet_wrap(~ species, ncol = 2, scales = "free") +
  theme(legend.position = "none") +
  xlab("Length (cm)") + ylab("Standard error") +
  scale_color_manual(values = viridis(3, end = 0.9), name = "Pairing method") 
# dev.off()

# png("se_draft_legend_SM_10k.png", width = 190, height = 150, units='mm', res = 300)
se_fit_df$pairing_method <- factor(se_fit_df$pairing_method, 
                                   levels = c("Simultaneous", "Alternate", "Random bootstrap"))
ggplot(se_fit_df, aes(x = length_id, y = se_fit)) + 
  geom_line(aes(color = pairing_method)) + 
  facet_wrap(~ species, ncol = 2, scales = "free") +
  theme(legend.position = "none") +
  xlab("Length (cm)") + ylab("Standard error") +
  scale_color_manual(values = viridis(3, end = 0.9), name = "Pairing method") + 
  theme(legend.position = "bottom") 
# dev.off()



## Figure 1
## Get mcrs data
dat_mcrs <- sqlQuery(con, ('SELECT * FROM irish_sea_celtic_seas_weight_at_length'))
dat_mcrs <- subset(unique(dat_mcrs[, c("species_code","species", "mcrs")]))
dat_mcrs <- subset(dat_mcrs, species %in% spp)

twin_dat_wide <- merge(twin_dat_wide,
                       dat_mcrs[, c("species", "mcrs")])


inches2mm <- 25.4

pdf("Figure_3_no_embed.pdf",
    height = 120/inches2mm,
    width = 170/inches2mm)

library(ggplot2); theme_set(theme_bw(base_size = 12))
plot_lf <- ggplot(twin_dat_wide, aes(x = length_id, y = Control_count_raised)) +
  geom_line() +
  geom_line(aes(y = Experimental_count_raised), lty = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  facet_wrap( ~ species, scales = "free", ncol = 2) +
  ## Plot mcrs length 
  geom_vline(aes(xintercept=mcrs), col = "red", data=dat_mcrs, show.legend=FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Length (cm)") +
  ylab("Raised count")
print(plot_lf)


dev.off()

embed_fonts("Figure_3_no_embed.pdf", outfile = "Figure_3_embed.pdf")


## 4 x 1

pdf("Figure_3_4x1_no_embed.pdf",
    height = 220/inches2mm,
    width = 85/inches2mm)

library(ggplot2); theme_set(theme_bw(base_size = 12))
plot_lf <- ggplot(twin_dat_wide, aes(x = length_id, y = Control_count_raised)) +
  geom_line() +
  geom_line(aes(y = Experimental_count_raised), lty = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  facet_wrap( ~ species, scales = "free", ncol = 1) +
  ## Plot mcrs length 
  geom_vline(aes(xintercept=mcrs), col = "red", data=dat_mcrs, show.legend=FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Length (cm)") +
  ylab("Raised count")
print(plot_lf)


dev.off()

embed_fonts("Figure_3_4x1_no_embed.pdf", outfile = "Figure_3_4x1_embed.pdf")





## Figure 4

library(ggplot2); theme_set(theme_bw(base_size = 12))
plot_prop_bubble_fig <- ggplot(twin_dat_wide_haul, 
                               aes(x = length_id, y = prop_Experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(alpha = 0.35, aes(colour = factor(tow_id), size = count_raised_total)) +
  # facet_grid(~ species, ncol = 4, scales = "free_x", space = "free") +
  facet_grid(~ species, scales = "free_x") +
  ylab("A") + xlab(NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = "none", axis.title.x=element_blank(),
        axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5))
# print(plot_prop_bubble_fig)

library(ggplot2); theme_set(theme_bw(base_size = 12))
plot_alt_prop_bubble_fig <- ggplot(twin_alt_dat_wide_haul, 
                                   aes(x = length_id, y = prop_Experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(alpha = 0.35, aes(colour = factor(tow_id_alt), size = count_raised_total)) +
  # facet_grid(~ species, scales = "free_x", labeller = NULL) +
  facet_grid(~ species, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank()) +
  xlab(NULL) + ylab("B")
# print(plot_alt_prop_bubble_fig)

inches2mm <- 25.4
pdf("Figure_4_no_embed.pdf",
    height = 80/inches2mm,
    width = 170/inches2mm)

grid.arrange(plot_prop_bubble_fig, plot_alt_prop_bubble_fig, 
             ncol=1, nrow = 2, left = "Proportion retained in T90 90mm", bottom = "Length (cm)")
dev.off()


embed_fonts("Figure_4_no_embed.pdf", outfile = "Figure_4_embed.pdf")



## figure 5

## Set plot layout and margins
par(mfrow = c(4, 3), mar = c(2, 2.2, 1, 1), oma = c(2, 3.5, 1, .5))
library(ggplot2); theme_set(theme_bw(base_size = 12))
## Twinrig paired
plot_twin_pair_prop_gam_fig <- ggplot(twin_dat_wide,
                                      aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = twin_pair_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              # alpha = 0.2, fill = "#8DC73F60") +
              alpha = 0.35, fill = "#48494B") +
  geom_line(data = twin_pair_all_pred_df, aes(y = proportion),
            colour = "black") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        legend.position = "none") +
  # ggtitle("Twin-rig paired") +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0.5))

## Alternate paired
plot_twin_alt_prop_gam_fig <- ggplot(twin_dat_wide,
                                     aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = twin_alt_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              # alpha = 0.2, fill = "#8DC73F60") +
              alpha = 0.35, fill = "#48494B") +
  geom_line(data = twin_alt_all_pred_df, aes(y = proportion),
            colour = "black") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none") +
  # ggtitle("Alternate paired") +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0.5))

## Bootstrap paired
plot_twin_boot_prop_gam_fig <- ggplot(twin_dat_wide, 
                                      aes(x = length_id, y = prop_experimental)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(aes(size = log(count_raised_total)), pch = 1) +
  # geom_point(aes(size = count_raised_total), pch = 1, colour = "slategrey") +
  geom_ribbon(data = boot_all_pred_df,
              aes(y = proportion, ymin = lower, ymax = upper),
              # alpha = 0.2, fill = "#8DC73F60") +
              alpha = 0.35, fill = "#48494B") +
  geom_line(data = boot_all_pred_df, aes(y = proportion),
            colour = "black") +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none") +
  # ggtitle("Bootstrap random paired") +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0.5))

pdf("Figure_5_no_embed.pdf",
    height = 180/inches2mm,
    width = 170/inches2mm)

grid.arrange(plot_twin_pair_prop_gam_fig, plot_twin_alt_prop_gam_fig, plot_twin_boot_prop_gam_fig,
             ncol=3, nrow = 1,
             bottom = "Length (cm)", left = "Proportion retained in T90 90mm")

dev.off()


embed_fonts("Figure_5_no_embed.pdf", outfile = "Figure_5_embed.pdf")

# pdf("Figure_5_no_embed_v2.pdf",
#     height = 180/inches2mm,
#     width = 170/inches2mm)
# 
# grid.arrange(plot_twin_pair_prop_gam_fig, plot_twin_alt_prop_gam_fig, plot_twin_boot_prop_gam_fig, 
#              ncol=3, nrow = 1, 
#              bottom = "Length (cm)", left = "Proportion retained in T90 90mm")
# 
# dev.off()
# 
# 
# embed_fonts("Figure_5_no_embed_v2.pdf", outfile = "Figure_5_embed_v2.pdf")
















## figure 6

pdf("Figure_6_no_embed.pdf",
    height = 125/inches2mm,
    width = 170/inches2mm)

library(ggplot2); theme_set(theme_bw(base_size = 12))
# png("se_draft_SM_10k.png", width = 190, height = 140, units='mm', res = 300)
ggplot(se_fit_df, aes(x = length_id, y = se_fit, linetype = pairing_method)) + 
  geom_line(aes(color = pairing_method)) +
  facet_wrap(~ species, ncol = 2, scales = "free") +
  theme(legend.position = "none") +
  xlab("Length (cm)") + ylab("Standard error of estimated log-odds") +
  scale_color_manual(values = viridis(3, end = 0.9), name = "Pairing method") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dev.off()


embed_fonts("Figure_6_no_embed.pdf", outfile = "Figure_6_embed.pdf")


pdf("Figure_6_4x1_no_embed.pdf",
    height = 225/inches2mm,
    width = 85/inches2mm)

library(ggplot2); theme_set(theme_bw(base_size = 12))
# png("se_draft_SM_10k.png", width = 190, height = 140, units='mm', res = 300)
ggplot(se_fit_df, aes(x = length_id, y = se_fit, linetype = pairing_method)) + 
  geom_line(aes(color = pairing_method)) +
  facet_wrap(~ species, ncol = 1, scales = "free") +
  theme(legend.position = "none") +
  xlab("Length (cm)") + ylab("Standard error of estimated log-odds") +
  scale_color_manual(values = viridis(3, end = 0.9), name = "Pairing method") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dev.off()


embed_fonts("Figure_6_4x1_no_embed.pdf", outfile = "Figure_6_4x1_embed.pdf")




## figure 7 random effects
# ## Load data
# load("20200514_random_effects_workspace.RData")

pdf("Figure_7_no_embed.pdf",
    height = 120/inches2mm,
    width = 170/inches2mm)

## Plot random effect by tow, facet wrapped by species
library(ggplot2); theme_set(theme_bw(base_size = 12))
ggplot(data = df_random_effects, aes(tow, random_effect)) +
  geom_point(size = 3, shape = 1) + 
  geom_point(shape = 17, data = df_alt_random_effects, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Log-odds random effect", x = "Haul") + 
  theme( # remove the horizontal grid lines
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    # explicitly set the vertical lines
    panel.grid.major.x = element_line( size=.01, color="#D3D3D3"),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.1, fill = NA)
  ) +
  facet_wrap(~ species, drop = TRUE)

dev.off()


embed_fonts("Figure_7_no_embed.pdf", outfile = "Figure_7_embed.pdf")

pdf("Figure_7_4x1_no_embed.pdf",
    height = 225/inches2mm,
    width = 85/inches2mm)

## Plot random effect by tow, facet wrapped by species
library(ggplot2); theme_set(theme_bw(base_size = 12))
ggplot(data = df_random_effects, aes(tow, random_effect)) +
  geom_point(size = 3, shape = 1) + 
  geom_point(shape = 17, data = df_alt_random_effects, size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Log-odds random effect", x = "Haul") + 
  theme( # remove the horizontal grid lines
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    # explicitly set the vertical lines
    panel.grid.major.x = element_line( size=.01, color="#D3D3D3"),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.1, fill = NA)
  ) +
  facet_wrap(~ species, drop = TRUE, ncol = 1)

dev.off()


embed_fonts("Figure_7_4x1_no_embed.pdf", outfile = "Figure_7_4x1_embed.pdf")



