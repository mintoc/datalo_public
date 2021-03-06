##-----------------------------
## Simple pairing simulation
## CM: 1/7/20
## based on earlier sim for Nephrops 
##-----------------------------
library(ggplot2); theme_set(theme_bw())
library(reshape)
library(mgcv)

## number of paired hauls
J <- 10

scenario <- c("Exchangeable", "Unbalanced", "Balanced", "Sequential", "Random")

##all_scen_pred_df <- NULL
##scen_dat <- NULL

low <- 100
high <- 1000

######################
## parallel component
######################
library(doSNOW)
## number of simulations
nsim <- 1000
## register 10 nodes
cl <- makeSOCKcluster(10)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## parallel loop
all_sim_res <- foreach(sim = 1:nsim, .combine = rbind, .options.snow = opts) %dopar% {
    library(mgcv)
    library(reshape)
    ## dataframe container for predictions
    all_pred_df <- NULL
    for(scen in scenario){
        print(scen)
        ## abundances - what we change for each scenario
        if(scen == "Exchangeable"){
            nj <- rep(high, J)
        }
        if(scen == "Unbalanced"){
            nj <- c(rep(low, 3), rep(high, J-3))
        }
        if(scen == "Balanced"){
            nj <- c(rep(low, 4), rep(high, J-4))
        }
        if(scen == "Sequential"){
            nj <- round(seq(low, high, length = J))
        }
        if(scen == "Random"){
            nj <- sample(c(low, high), size = J, replace = TRUE)
        }
        ## to keep all the same except abundance
        set.seed(sim)
        dat <- NULL
        ## standard deviation of lognormal
        sd <- 0.3
        for(j in 1:J){
            ## number of fish at length contacting gear
            n_table <- table(round(rlnorm(nj[j], meanlog = log(30) - sd^2/2, sdlog = sd)))
            n <- as.numeric(n_table)
            l <- as.numeric(names(n_table))
            ## assume equal fishing power between the gears
            pC <- 0.5
            pT <- 1 - pC
            ## contact selection
            ## control
            rlC <- plogis(0.2 * (l - 30))
            rlT <- plogis(0.5 * (l - 35))
            ## counts
            nT <- rpois(length(n), lambda = n * pT * rlT)
            nC <- rpois(length(n), lambda = n * pC * rlC)
            dat <- rbind(dat, data.frame(l = l, haul = j, nT = nT, nC = nC))
            dat <- subset(dat, nT > 0 | nC > 0)
        }
        dat$scenario <- scen
        ##scen_dat <- rbind(scen_dat, dat)
        ##
        dat_long <- melt(dat, id.vars = c("haul", "l"))
        dat_long$variable <- factor(as.character(dat_long$variable), levels = c("nC", "nT"))
        ## ggplot(dat_long, aes(x = l, y = value)) +
        ##     geom_bar(stat = "identity") +
        ##     facet_grid(haul ~ variable)
        dat$fhaul <- factor(paste0("h", dat$haul), levels = paste0("h", 1:J))
        ## now apply the three methods to see what uncertainty they give
        ## predicted length range approximately 95% of density of lognormal
        ## qlnorm(c(0.025, 0.975), meanlog = log(30) - sd^2/2, sdlog = sd)
        l_pred <- 16:52
        ##---------------
        ## SIMULTANEOUS
        ##---------------
        Y <- with(dat, cbind(nT, nC))
        gam_fit <- gam(Y ~ s(l, k = 10) + s(fhaul, bs="re"),
                       data = dat, family = binomial,
                       method = "ML")
        ## get predictions overall
        ##length_range <- range(dat$l)
        pred_df <- data.frame(
            l = l_pred,
            fhaul = dat$fhaul[1])
        ##
        pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE, exclude = "s(fhaul)")
        pred_df$fhaul <- NULL
        pred_df$proportion <- plogis(pred$fit)
        pred_df$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
        pred_df$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
        pred_df$method <- "Simultaneous pairs"
        pred_df$scenario <- scen
        pred_df$sim <- sim
        ##
        all_pred_df <- rbind(all_pred_df, pred_df)
        ##-----------
        ## ALTERNATE 
        ##-----------
        ## make alternate data
        pairs <- matrix(1:J, ncol = 2, byrow = TRUE)
        ##
        alt_dat <- NULL
        for(i in 1:nrow(pairs)){
            haul_i <- subset(dat, haul == pairs[i, 1])
            haul_j <- subset(dat, haul == pairs[i, 2])
            ## merge 1st test with second control
            tmp_i <- merge(haul_i[, c("l", "nT")], haul_j[, c("l", "nC")], all = TRUE)
            tmp_i[is.na(tmp_i)] <- 0
            tmp_i$haul <- pairs[i, 1]
            alt_dat <- rbind(alt_dat, tmp_i)
            ## merge second test with 1st control 
            tmp_j <- merge(haul_j[, c("l", "nT")], haul_i[, c("l", "nC")], all = TRUE)
            tmp_j[is.na(tmp_j)] <- 0
            tmp_j$haul <- pairs[i, 2]
            alt_dat <- rbind(alt_dat, tmp_j)
        }
        alt_dat_long <- melt(alt_dat, id.vars = c("haul", "l"))
        alt_dat_long$variable <- factor(as.character(alt_dat_long$variable), levels = c("nC", "nT"))
        ##
        ##ggplot(alt_dat_long, aes(x = l, y = value)) +
        ##    geom_bar(stat = "identity") +
        ##    facet_grid(haul ~ variable)
        alt_dat$fhaul <- factor(paste0("h", alt_dat$haul), levels = paste0("h", 1:J))
        Y <- with(alt_dat, cbind(nT, nC))
        gam_fit <- gam(Y ~ s(l, k = 10) + s(fhaul, bs="re"),
                       data = alt_dat, family = binomial,
                       method = "ML")
        ## get predictions overall
        length_range <- range(dat$l)
        pred_df <- data.frame(
            l = l_pred,
            fhaul = dat$fhaul[1])
        pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE, exclude = "s(fhaul)")
        pred_df$fhaul <- NULL
        pred_df$proportion <- plogis(pred$fit)
        pred_df$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
        pred_df$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
        pred_df$method <- "Alternate pairs"
        pred_df$scenario <- scen
        pred_df$sim <- sim
        ##
        all_pred_df <- rbind(all_pred_df, pred_df)
        ##--------
        ## RANDOM 
        ##--------
        l <- l_pred
        n_boot <- 1e3
        ##n_boot <- 10
        boot_pred <- matrix(NA, nrow = n_boot, ncol = length(l))
        pred_df <- data.frame(l = l)
        for(j in 1:n_boot){
            ## Assign hauls from subsetted data to avoid hauls with no data
            ## experimental tows - NB replace = TRUE here
            T_hauls <- sample(1:J, replace = TRUE)
            ## control tows
            C_hauls <- sample(1:J, replace = TRUE)
            boot_dat <- NULL
            for(i in 1:J){
                haul_i <- subset(dat, haul == T_hauls[i])
                haul_j <- subset(dat, haul == C_hauls[i])
                ## resample lengths within the haul - "s" for sample
                ## test
                l_T <- with(haul_i, rep(l, times = nT))
                ##sl_T <- sample(l_T, replace = TRUE)
                sl_T <- sample(l_T)
                shaul_i <- as.data.frame(table(sl_T))
                shaul_i$l <- as.numeric(as.character(shaul_i$sl_T))
                shaul_i$nT <- shaul_i$Freq
                ## control
                l_C <- with(haul_j, rep(l, times = nC))
                ##sl_C <- sample(l_C, replace = TRUE)
                sl_C <- sample(l_C)
                shaul_j <- as.data.frame(table(sl_C))
                shaul_j$l <- as.numeric(as.character(shaul_j$sl_C))
                shaul_j$nC <- shaul_j$Freq
                ##
                tmp_i <- merge(shaul_i[, c("l", "nT")], shaul_j[, c("l", "nC")], all = TRUE)
                tmp_i[is.na(tmp_i)] <- 0
                tmp_i$haul <- i
                ##if(nrow(tmp_i) > 0){
                boot_dat <- rbind(boot_dat, tmp_i)
                ##}
            }
            boot_dat$fhaul <- factor(paste0("h", boot_dat$haul), levels = paste0("h", 1:J))
            Y <- with(boot_dat, cbind(nT, nC))
            gam_fit <- gam(Y ~ s(l, k = 10),
                           data = boot_dat, family = binomial,
                           method = "ML")
            pred <- predict(gam_fit, newdata = pred_df)
            boot_pred[j,] <- plogis(pred)
        }
        pred_df <- data.frame(l = l,
                              proportion = apply(boot_pred, 2, mean, na.rm = TRUE),
                              upper = apply(boot_pred, 2, quantile, probs = 0.975, na.rm = TRUE),
                              lower = apply(boot_pred, 2, quantile, probs = 0.025, na.rm = TRUE),
                              method = "Random pairs",
                              scenario = scen,
                              sim = sim)
        all_pred_df <- rbind(all_pred_df, pred_df)
    }
    all_pred_df
}

close(pb)
stopCluster(cl)

##save(all_sim_res, file = paste0("all_res.RData"))

##all_sim_res <- unique(all_sim_res)

load("all_res.RData")

## plot
all_sim_res$method[all_sim_res$method == "Simultaneous pairs"] <- "Simultaneous match"
all_sim_res$method[all_sim_res$method == "Alternate pairs"] <- "Consecutive match"
all_sim_res$method[all_sim_res$method == "Random pairs"] <- "Random match"

all_sim_res$method <- factor(as.character(all_sim_res$method), levels = c("Simultaneous match", "Consecutive match", "Random match"))

all_sim_res$scenario <- factor(as.character(all_sim_res$scenario), levels = c("Exchangeable", "Balanced", "Unbalanced", "Sequential", "Random"))

## aggregate
mean_df <- aggregate(proportion ~ l + method + scenario,
                     FUN = mean,
                     data = all_sim_res)
##

lwr_df <- aggregate(lower ~ l + method + scenario,
                    FUN = mean,
                    data = all_sim_res)
##names(lwr_df)[names(lwr_df) == "proportion"] <- "lwr"
##
upr_df <- aggregate(upper ~ l + method + scenario,
                    FUN = mean,
                    data = all_sim_res)
##names(upr_df)[names(upr_df) == "proportion"] <- "upr"

tmp <- merge(mean_df, lwr_df)

plot_df <- merge(tmp, upr_df)

true_df <- data.frame(l = 16:52)

true_df$rlC <- plogis(0.2 * (true_df$l - 30))
true_df$rlT <- plogis(0.5 * (true_df$l - 35))

true_df$proportion <- with(true_df, rlT / (rlC + rlT))

label_df <- unique(plot_df[, c("method", "scenario")])
label_df <- subset(label_df, method == "Simultaneous match")
label_df$l <- 17
label_df$proportion <- 0.7
label_df$lab <- NA
label_df$lab[label_df$scenario == "Exchangeable"] <- "A"
label_df$lab[label_df$scenario == "Balanced"] <- "B"
label_df$lab[label_df$scenario == "Unbalanced"] <- "C"
label_df$lab[label_df$scenario == "Sequential"] <- "D"
label_df$lab[label_df$scenario == "Random"] <- "E"

## label_df$lab[label_df$scenario == "Exchangeable"] <- "(a)"
## label_df$lab[label_df$scenario == "Balanced"] <- "(b)"
## label_df$lab[label_df$scenario == "Unbalanced"] <- "(c)"
## label_df$lab[label_df$scenario == "Sequential"] <- "(d)"
## label_df$lab[label_df$scenario == "Random"] <- "(e)"

inches2mm <- 25.4
library(extrafont)

theme_set(theme_bw(base_size = 12))

pdf("Figure_2_08_03_21_no_embed.pdf",
    height = 170/inches2mm,
    width = 170/inches2mm)
ggplot(plot_df, aes(x = l, y = proportion)) +
    geom_ribbon(aes(y = proportion, ymin = lower, ymax = upper),
                alpha = 0.2, fill = "darkgrey", col = grey(0.1), lwd = 0.2) +
    geom_line(aes(y = proportion)) +
    geom_hline(yintercept = 0.5, linetype = 3) +
    geom_point(data = true_df, size = 0.5, pch = 1) + 
    facet_grid(scenario ~ method) +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Length (cm)") +
    ylab("Proportion retained in test gear") +
    geom_text(data = label_df, aes(label = lab))
dev.off()

embed_fonts("Figure_2_08_03_21_no_embed.pdf", outfile = "Figure_2_08_03_21_embed.pdf")


## all means
png("simulation_all_means.png", height = 8, width = 8, units = "in", res = 400)
ggplot(all_sim_res, aes(x = l, y = proportion, group = sim)) +
    geom_line(alpha = 0.05, col = "forestgreen", lwd = 0.2) +
    facet_grid(scenario ~ method) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size=14)) +
    xlab("Length (cm)") +
    ylab("Proportion retained in test gear")
dev.off()

png("simulation_all_CIs.png", height = 8, width = 8, units = "in", res = 400)
ggplot(all_sim_res, aes(x = l, y = upper, group = sim)) +
    geom_line(alpha = 0.05, col = "turquoise", lwd = 0.2) +
    geom_line(aes(y = lower), alpha = 0.05, col = "turquoise", lwd = 0.2) +
    facet_grid(scenario ~ method) +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Length (cm)") +
    ylab("Proportion retained in test gear")
dev.off()


## Brier score at length
true_df$p_true <- true_df$proportion

all_sim_res <- merge(all_sim_res, true_df[, c("l", "p_true")])

all_sim_res$e2 <- with(all_sim_res, (proportion - p_true)^2)

brier_df <- aggregate(e2 ~ l + method + scenario,
                      FUN = mean,
                      ##FUN = median,
                      data = all_sim_res)

names(brier_df)[names(brier_df) == "e2"] <- "brier"

pdf("brier.pdf", height = 9, width = 5)
ggplot(brier_df, aes(x = l, y = brier, colour = method)) +
    geom_line() +
    facet_wrap(~ scenario, ncol = 1) +
    scale_colour_manual(values = c("cadetblue", "orange2", "purple")) +
    geom_vline(xintercept = c(20, 35, 50), linetype = 2) +
    theme(legend.position = "bottom")
dev.off()

## average confidence interval widths
all_sim_res$width <- with(all_sim_res, (upper - lower))

width_df <- aggregate(width ~ l + method + scenario,
                      FUN = mean,
                      data = all_sim_res)

pdf("width.pdf", height = 9, width = 5)
ggplot(width_df, aes(x = l, y = width, colour = method)) +
    geom_line() +
    facet_wrap(~ scenario, ncol = 1) +
    scale_colour_manual(values = c("cadetblue", "orange2", "purple")) +
    geom_vline(xintercept = c(20, 35, 50), linetype = 2) +
    theme(legend.position = "bottom")
dev.off()


lengths <- c(20, 35, 50)

## output precision table
library(reshape)
wide <- subset(width_df, l %in% lengths)

longd <- melt(wide, id.vars = c("l", "method", "scenario"))
names(longd)[names(longd) == "method"] <- "methodd"

out <- cast(data = longd, scenario + l ~ methodd)

write.csv(out, file = "precision_table.csv", row.names = FALSE)
