# _____________________________________________________________________________
# Minimum required imports
# _____________________________________________________________________________
library(rddensity)
library(rdrop2)
library(lpdensity)
library(ggplot2)
library(rdd)
library(car)
library(dplyr)
library(lubridate)
library(dummies)
library(extrafont)
library(gridExtra)

#For generation of Times New Roman fonts (optional for replication)
#font_import()

# _____________________________________________________________________________
# Requisite Data
# _____________________________________________________________________________

station.dat.4 <-read.csv("station_data_dataverse.csv")
external.sample <-read.csv("external_sample.csv")


#Data cleaning
#Eliminating transactions with 0 kWh 
sub<-which(station.dat.4$kwhTotal==0.000)
station.dat.4<-station.dat.4[-sub,]

#Subset free transactions
free <- station.dat.4[station.dat.4$dollars == 0,]


# _____________________________________________________________________________
# DESCRIPTIVE STATISTICS
# _____________________________________________________________________________

#Placing paid transactions in a separate dataset (376 observations)
sub<-which(station.dat.4$dollars==0)
station.dat.4.2<-station.dat.4[-sub,]

#Generate count of transactons by user ID for histogram
station.dat.agg2<-tally(group_by(station.dat.4, userId))

#Eliminating observations with no reported zipcode for distance calculations
station.dat.4.3<-station.dat.4[station.dat.4$reportedZip == 1,]

#Calculating total distance per user
station.dat.agg<-aggregate(station.dat.4.3$distance ~ station.dat.4.3$userId, data = station.dat.4.3, sum)
station.dat.agg$distance<-station.dat.agg$`station.dat.4.3$distance`

#Create matrix to display descriptive statistics
table_2 <- matrix(nrow = 6,ncol = 5, dimnames = list(c("Charge time (hours)",
                                                       "Total consumption (kWh)",
                                                       "Repeat transactions per user (count)",
                                                       "Session revenue ($)",
                                                       "Estimated daily commute distance - one way (miles)",
                                                       "Electric vechicle miles traveled per user (miles)"),
                                                     c("M", "SD", "Min", "Max", "Active sessions")))


#3 datasets used in calculating descriptive stats:
#all observations, observations with cost > 0, and observations reporting zip code

#Calculating mean of each variable
table_2[,"M"] <- c(round(mean(station.dat.4$chargeTimeHrs),2),
                   round(mean(station.dat.4$kwhTotal),2),
                   round(mean(station.dat.agg2$n),2),
                   round(mean(station.dat.4.2$dollars),2),
                   round(mean(station.dat.4$distance[!is.na(station.dat.4$distance)]),2),
                   round(mean(station.dat.agg$distance),2))

#Calculating standard deviation of each variable
table_2[,"SD"] <- c(round(sd(station.dat.4$chargeTimeHrs),2),
                    round(sd(station.dat.4$kwhTotal),2),
                    round(sd(station.dat.agg2$n),2),
                    round(sd(station.dat.4.2$dollars),2),
                    round(sd(station.dat.4$distance[!is.na(station.dat.4$distance)]),2),
                    round(sd(station.dat.agg$distance),2))

#Calculating minimum of each variable
table_2[,"Min"] <- c(round(min(station.dat.4$chargeTimeHrs),2),
                     round(min(station.dat.4$kwhTotal),2),
                     round(min(station.dat.agg2$n),2),
                     round(min(station.dat.4.2$dollars),2),
                     round(min(station.dat.4$distance[!is.na(station.dat.4$distance)]),2),
                     round(min(station.dat.agg$distance),2))

#Calculating maximum of each variable
table_2[,"Max"] <- c(round(max(station.dat.4$chargeTimeHrs),2),
                     round(max(station.dat.4$kwhTotal),2),
                     round(max(station.dat.agg2$n),2),
                     round(max(station.dat.4.2$dollars),2),
                     round(max(station.dat.4$distance[!is.na(station.dat.4$distance)]),2),
                     round(max(station.dat.agg$distance),2))

#Calculating total sessions used in analysis for each variable
table_2[,"Active sessions"] <- c(nrow(station.dat.4),
                                 nrow(station.dat.4),
                                 nrow(station.dat.4),
                                 nrow(station.dat.4.2),
                                 nrow(station.dat.4.3),
                                 nrow(station.dat.4.3))

table_2


#Calculating statistic on number of users who are treated at 4 hours
#station.dat.test <- station.dat.4[station.dat.4$chargeTimeHrs >= 4,]
#a <- c(station.dat.test$userId)
#station.dat.4 <- station.dat.4[station.dat.4$userId %in% a,]
#length(unique(station.dat.4$userId))

# _____________________________________________________________________________
# Creation of function for creating McCrary Figures
# _____________________________________________________________________________

#adapts output of DCdensity function for more flexible aesthetic changes to figures
#function adapted from GitHub user mikedecr to calculate at 95% CI
#https://gist.github.com/mikedecr/6ae9c63b6d28c43b068ddc0d85e8897b
mccrary <-
  function (runvar, cutpoint, bin = NULL, bw = NULL, verbose = FALSE,
            plot = TRUE, ext.out = FALSE, htest = FALSE)
  {
    library(rdd)
    runvar <- runvar[complete.cases(runvar)]
    rn <- length(runvar)
    rsd <- sd(runvar)
    rmin <- min(runvar)
    rmax <- max(runvar)
    if (missing(cutpoint)) {
      if (verbose)
        cat("Assuming cutpoint of zero.\n")
      cutpoint <- 0
    }
    if (cutpoint <= rmin | cutpoint >= rmax) {
      stop("Cutpoint must lie within range of runvar")
    }
    if (is.null(bin)) {
      bin <- 2 * rsd * rn^(-1/2)
      if (verbose)
        cat("Using calculated bin size: ", sprintf("%.3f",
                                                   bin), "\n")
    }
    l <- floor((rmin - cutpoint)/bin) * bin + bin/2 + cutpoint
    r <- floor((rmax - cutpoint)/bin) * bin + bin/2 + cutpoint
    lc <- cutpoint - (bin/2)
    rc <- cutpoint + (bin/2)
    j <- floor((rmax - rmin)/bin) + 2
    binnum <- round((((floor((runvar - cutpoint)/bin) * bin +
                         bin/2 + cutpoint) - l)/bin) + 1)
    cellval <- rep(0, j)
    for (i in seq(1, rn)) {
      cnum <- binnum[i]
      cellval[cnum] <- cellval[cnum] + 1
    }
    cellval <- (cellval/rn)/bin
    cellmp <- seq(from = 1, to = j, by = 1)
    cellmp <- floor(((l + (cellmp - 1) * bin) - cutpoint)/bin) *
      bin + bin/2 + cutpoint
    if (is.null(bw)) {
      leftofc <- round((((floor((lc - cutpoint)/bin) * bin +
                            bin/2 + cutpoint) - l)/bin) + 1)
      rightofc <- round((((floor((rc - cutpoint)/bin) * bin +
                             bin/2 + cutpoint) - l)/bin) + 1)
      if (rightofc - leftofc != 1) {
        stop("Error occurred in bandwidth calculation")
      }
      cellmpleft <- cellmp[1:leftofc]
      cellmpright <- cellmp[rightofc:j]
      P.lm <- lm(cellval ~ poly(cellmp, degree = 4, raw = T),
                 subset = cellmp < cutpoint)
      mse4 <- summary(P.lm)$sigma^2
      lcoef <- coef(P.lm)
      fppleft <- 2 * lcoef[3] + 6 * lcoef[4] * cellmpleft +
        12 * lcoef[5] * cellmpleft * cellmpleft
      hleft <- 3.348 * (mse4 * (cutpoint - l)/sum(fppleft *
                                                    fppleft))^(1/5)
      P.lm <- lm(cellval ~ poly(cellmp, degree = 4, raw = T),
                 subset = cellmp >= cutpoint)
      mse4 <- summary(P.lm)$sigma^2
      rcoef <- coef(P.lm)
      fppright <- 2 * rcoef[3] + 6 * rcoef[4] * cellmpright +
        12 * rcoef[5] * cellmpright * cellmpright
      hright <- 3.348 * (mse4 * (r - cutpoint)/sum(fppright *
                                                     fppright))^(1/5)
      bw = 0.5 * (hleft + hright)
      if (verbose)
        cat("Using calculated bandwidth: ", sprintf("%.3f",
                                                    bw), "\n")
    }
    if (sum(runvar > cutpoint - bw & runvar < cutpoint) == 0 |
        sum(runvar < cutpoint + bw & runvar >= cutpoint) == 0)
      stop("Insufficient data within the bandwidth.")
    if (plot) {
      d.l <- data.frame(cellmp = cellmp[cellmp < cutpoint],
                        cellval = cellval[cellmp < cutpoint], dist = NA,
                        est = NA, lwr = NA, upr = NA)
      pmin <- cutpoint - 2 * rsd
      pmax <- cutpoint + 2 * rsd
      for (i in 1:nrow(d.l)) {
        d.l$dist <- d.l$cellmp - d.l[i, "cellmp"]
        w <- kernelwts(d.l$dist, 0, bw, kernel = "triangular")
        newd <- data.frame(dist = 0)
        pred <- predict(lm(cellval ~ dist, weights = w, data = d.l),
                        interval = "confidence", level = 0.95, newdata = newd)
        d.l$est[i] <- pred[1]
        d.l$lwr[i] <- pred[2]
        d.l$upr[i] <- pred[3]
      }
      d.r <- data.frame(cellmp = cellmp[cellmp >= cutpoint],
                        cellval = cellval[cellmp >= cutpoint], dist = NA,
                        est = NA, lwr = NA, upr = NA)
      for (i in 1:nrow(d.r)) {
        d.r$dist <- d.r$cellmp - d.r[i, "cellmp"]
        w <- kernelwts(d.r$dist, 0, bw, kernel = "triangular")
        newd <- data.frame(dist = 0)
        pred <- predict(lm(cellval ~ dist, weights = w, data = d.r),
                        interval = "confidence", newdata = newd)
        d.r$est[i] <- pred[1]
        d.r$lwr[i] <- pred[2]
        d.r$upr[i] <- pred[3]
      }
      plot(d.l$cellmp, d.l$est, lty = 1, lwd = 2, col = "black",
           type = "l", xlim = c(pmin, pmax), ylim = c(min(cellval[cellmp <=
                                                                    pmax & cellmp >= pmin]), max(cellval[cellmp <=
                                                                                                           pmax & cellmp >= pmin])), xlab = NA, ylab = NA,
           main = NA)
      lines(d.l$cellmp, d.l$lwr, lty = 2, lwd = 1, col = "black",
            type = "l")
      lines(d.l$cellmp, d.l$upr, lty = 2, lwd = 1, col = "black",
            type = "l")
      lines(d.r$cellmp, d.r$est, lty = 1, lwd = 2, col = "black",
            type = "l")
      lines(d.r$cellmp, d.r$lwr, lty = 2, lwd = 1, col = "black",
            type = "l")
      lines(d.r$cellmp, d.r$upr, lty = 2, lwd = 1, col = "black",
            type = "l")
      points(cellmp, cellval, type = "p", pch = 20)
    }
    cmp <- cellmp
    cval <- cellval
    padzeros <- ceiling(bw/bin)
    jp <- j + 2 * padzeros
    if (padzeros >= 1) {
      cval <- c(rep(0, padzeros), cellval, rep(0, padzeros))
      cmp <- c(seq(l - padzeros * bin, l - bin, bin), cellmp,
               seq(r + bin, r + padzeros * bin, bin))
    }
    dist <- cmp - cutpoint
    w <- 1 - abs(dist/bw)
    w <- ifelse(w > 0, w * (cmp < cutpoint), 0)
    w <- (w/sum(w)) * jp
    fhatl <- predict(lm(cval ~ dist, weights = w), newdata = data.frame(dist = 0))[[1]]
    w <- 1 - abs(dist/bw)
    w <- ifelse(w > 0, w * (cmp >= cutpoint), 0)
    w <- (w/sum(w)) * jp
    fhatr <- predict(lm(cval ~ dist, weights = w), newdata = data.frame(dist = 0))[[1]]
    thetahat <- log(fhatr) - log(fhatl)
    sethetahat <- sqrt((1/(rn * bw)) * (24/5) * ((1/fhatr) +
                                                   (1/fhatl)))
    z <- thetahat/sethetahat
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    if (verbose) {
      cat("Log difference in heights is ", sprintf("%.3f",
                                                   thetahat), " with SE ", sprintf("%.3f", sethetahat),
          "\n")
      cat("  this gives a z-stat of ", sprintf("%.3f", z),
          "\n")
      cat("  and a p value of ", sprintf("%.3f", p), "\n")
    }
    
    estimates <- data.frame(dhat = c(d.l$est, d.r$est), dlower = c(d.l$lwr, d.r$lwr), dupper = c(d.l$upr, d.r$upr),
                            force = c(rep(0, length(d.l$est)), rep(1, length(d.r$est)))
    )
    
    if (ext.out)
      return(
        list(
          theta = thetahat, se = sethetahat, z = z,
          p = p, binsize = bin, bw = bw, cutpoint = cutpoint,
          data = data.frame(cellmp, cellval, force = c(rep(0, length(d.l$est)), rep(1, length(d.r$est)))),
          estimates = estimates
        )
      )
    else if (htest) {
      structure(list(statistic = c(z = z), p.value = p, method = "McCrary (2008) sorting test",
                     parameter = c(binwidth = bin, bandwidth = bw, cutpoint = cutpoint),
                     alternative = "no apparent sorting"), class = "htest")
    }
    else return(p)
  }
# _____________________________________________________________________________
# Generating McCrary Figures
# _____________________________________________________________________________

#Transform variable to date for later calculations
station.dat.4$created<-as.Date(station.dat.4$created)

#Create variable for week of program
station.dat.4$weeknos <- (interval(min(station.dat.4$created), station.dat.4$created) %/% weeks(1)) + 1

#Generate total sessions by userID
station.dat.4 <-
  station.dat.4 %>%
  dplyr::group_by(userId) %>%
  dplyr::mutate(total = length(unique(sessionId)))

#Generate subsets needed for high and low-volume users
high_volume <- station.dat.4[station.dat.4$total >= 20,]
all_volume <- station.dat.4

# -----------------------------------------------------------------------------
# (a) 4-hour cutoff
# -----------------------------------------------------------------------------

#see comments from generation of figure (a); follows parallel proess

mccrary_b <- mccrary(station.dat.4$chargeTimeHrs, 4, bin = 4.15*sd(station.dat.4$chargeTimeHrs)*length(station.dat.4$chargeTimeHrs)^(-.5),
                     bw = 0.275, verbose = TRUE,
                     plot = TRUE, ext.out = TRUE, htest = FALSE)

sub_mb <- mccrary_b$estimates
sub_mb_pre <- sub_mb
sub_mb_post <- sub_mb
sub_mb_pre[sub_mb_pre$force == 1,] <- NA
sub_mb_post[sub_mb_post$force == 0,] <- NA 
sub_mb1 <- mccrary_b$data

figure_s2a <- ggplot(sub_mb1, aes(y = cellval, x = cellmp)) +
  geom_line(aes(y = sub_mb_post$dlower, col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5)+
  geom_line(aes(y = sub_mb_post$dhat, col = "est"), colour="grey30",size = 0.5) + 
  geom_line(aes(y = sub_mb_post$dupper), linetype= "dashed",colour="grey30",size = 0.5)+
  geom_line(aes(y = sub_mb_pre$dlower, col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5)+
  geom_line(aes(y = sub_mb_pre$dhat, col = "est"), colour="grey30",size = 0.5) + 
  geom_line(aes(y = sub_mb_pre$dupper, col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5)+
  labs(x="Duration of Charge", y = "Density of Transcations") +
  coord_cartesian(xlim=c(0,7), ylim = c(0,.5)) +
  geom_vline(xintercept=4, size = 1.5) +
  geom_ribbon(aes(ymin=sub_mb_post$dlower, ymax=sub_mb_post$dupper), linetype=2, alpha=0.2, fill = "skyblue2")+
  geom_ribbon(aes(ymin=sub_mb_pre$dlower, ymax=sub_mb_pre$dupper), linetype=2, alpha=0.2, fill = "skyblue2")+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=11, face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

figure_s2a

# -----------------------------------------------------------------------------
# (b) 2-hour cutoff
# -----------------------------------------------------------------------------

#Get output of estimates and CIs
mccrary_a <- mccrary(station.dat.4$chargeTimeHrs, 2, bin = 5*sd(station.dat.4$chargeTimeHrs)*length(station.dat.4$chargeTimeHrs)^(-.5),
                     bw = 0.6, verbose = TRUE,
                     plot = TRUE, ext.out = TRUE, htest = FALSE)

#isolate estimates
sub_ma <- mccrary_a$estimates

#create replicas to serve as lines for pre and post-cutoff
sub_ma_pre <- sub_ma
sub_ma_post <- sub_ma

#replace estimates before or after cutoff with NA to maintain series length as required by ggplot
sub_ma_pre[sub_ma_pre$force == 1,] <- NA
sub_ma_post[sub_ma_post$force == 0,] <- NA 

#isolate axis data
sub_ma1 <- mccrary_a$data

#note: ggplot will throw warning for removed NA values. This is by design
#and is necessary for accurate output
figure_s2b <- ggplot(sub_ma1, aes(y = cellval, x = cellmp)) +
  geom_line(aes(y = sub_ma_post$dlower, col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5)+
  geom_line(aes(y = sub_ma_post$dhat, col = "est"), colour="grey30",size = 0.5) + 
  geom_line(aes(y = sub_ma_post$dupper), linetype= "dashed",colour="grey30",size = 0.5)+
  geom_line(aes(y = sub_ma_pre$dlower, col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5)+
  geom_line(aes(y = sub_ma_pre$dhat, col = "est"), colour="grey30",size = 0.5) + 
  geom_line(aes(y = sub_ma_pre$dupper, col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5)+
  labs(x="Duration of Charge", y = "Density of Transcations") +
  coord_cartesian(xlim=c(0,7), ylim = c(0,.5)) +
  geom_vline(xintercept=2, size = 1.5) +
  geom_ribbon(aes(ymin=sub_ma_post$dlower, ymax=sub_ma_post$dupper), linetype=2, alpha=0.2, fill = "salmon")+
  geom_ribbon(aes(ymin=sub_ma_pre$dlower, ymax=sub_ma_pre$dupper), linetype=2, alpha=0.2, fill = "salmon")+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=11, face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

figure_s2b

figure_s2 <- grid.arrange(figure_s2a, figure_s2b, ncol=2)

# _____________________________________________________________________________
# Figure on additions of new users over time
# _____________________________________________________________________________

#Take one row per user for use in generating cumulative frequency
high_volume <- station.dat.4[match(unique(high_volume$userId), high_volume$userId),]
high_volume$dummy <- 1
high_volume <- high_volume[order(high_volume$created),]

#Running sum of high-volume users over time
high_volume$cumsum <- cumsum(high_volume$dummy)

#Find max running total by week for plotting
high_volume <- 
  high_volume %>%
  dplyr::group_by(weeknos) %>%
  dplyr::mutate(finalcumf = max(cumsum))

#Take only first observation per week for plotting
high_volume <- high_volume[match(unique(high_volume$weeknos), high_volume$weeknos),]

#Take one row per user for use in generating cumulative frequency
all_volume <- station.dat.4[match(unique(all_volume$userId), all_volume$userId),]
all_volume$dummy <- 1
all_volume <- all_volume[order(all_volume$created),]

#Running sum of high-volume users over time
all_volume$cumsum <- cumsum(all_volume$dummy)

#Find max running total by week for plotting
all_volume <- 
  all_volume %>%
  dplyr::group_by(weeknos) %>%
  dplyr::mutate(finalcumf = max(cumsum))

#Take only first observation per week for plotting
all_volume <- all_volume[match(unique(all_volume$weeknos), all_volume$weeknos),]

#Extract only relevant columns for plotting
all_vector <- all_volume[,c('weeknos','finalcumf')]
high_vector <- high_volume[,c('weeknos','finalcumf')]

#Pad weeks that do not have a new user added with the cumulative number of users from
#previous week (46 total weeks in program)
for (i in 1:46) {
  if (!(i %in% all_vector$weeknos)) {
    all_vector[nrow(all_vector) + 1,] <- c(i, all_vector[i-1,][2])
    all_vector <- all_vector[order(all_vector$weeknos),]
  }
  if (!(i %in% high_vector$weeknos)) {
    high_vector[nrow(high_vector) + 1,] <- c(i, high_vector[i-1,][2])
    high_vector <- high_vector[order(high_vector$weeknos),]
  }
}

#flags to allow plot to distinguish between low and high-volume users
all_vector$flag <- c("all")
high_vector$flag <- c("high")
final_frame <- rbind(all_vector,high_vector)

figure_s5a <- ggplot(final_frame, aes(x=weeknos, y=finalcumf, col=flag)) +
  geom_line(size = 0.5) +
  labs(title = "", x="Week of program", y=  expression("Total number of users"), color = "", fill = "white") +
  scale_color_manual(labels = c("All users", "High-repeat users"), values = c("grey5", "orange2")) +
  scale_fill_identity(name='',guide = 'legend',labels = c('Plug-out Times','Plug-in Times'))+
  theme(legend.position = c(.1, 1),legend.justification = c(0, 1),legend.text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans")) +
  scale_x_continuous(breaks = seq(0, 45, by = 5)) +
  theme(legend.key=element_blank())

figure_s5a

# _____________________________________________________________________________
# Figure transactions per station by week
# _____________________________________________________________________________

#Find total transactions per week
station.dat.s4 <- 
  station.dat.4 %>%
  dplyr::group_by(weeknos) %>%
  dplyr::mutate(transactions = n())

#Take one row per user for use in generating cumulative frequency
cumulative_stats <- station.dat.s4[match(unique(station.dat.s4$stationId), station.dat.s4$stationId),]
cumulative_stats$dummy <- 1
cumulative_stats <- cumulative_stats[order(cumulative_stats$created),]

#Running sum of stations
cumulative_stats$cumsum <- cumsum(cumulative_stats$dummy)

#Find max running total by week for plotting
cumulative_stats <- 
  cumulative_stats %>%
  dplyr::group_by(weeknos) %>%
  dplyr::mutate(finalcumf = max(cumsum))

#Take only first observation per week for plotting
cumulative_stats <- cumulative_stats[match(unique(cumulative_stats$weeknos), cumulative_stats$weeknos),]

#Extract only relevant columns for plotting
final <- cumulative_stats[,c('weeknos','finalcumf')]

#Pad weeks that do not have a new user added with the cumulative number of users from
#previous week (46 total weeks in program)
for (i in 1:46) {
  if (!(i %in% final$weeknos)) {
    final[nrow(final) + 1,] <- c(i, final[i-1,][2])
    final <- final[order(final$weeknos),]
  }
}

#Replace single week that did not have station added (week 6)
cumulative_stats <- station.dat.s4[match(unique(station.dat.s4$weeknos), station.dat.s4$weeknos),]
cumulative_stats_last <- cumulative_stats[c('weeknos','transactions')]
row <- c(6,0)
cumulative_stats_last[nrow(cumulative_stats_last) + 1,] <- list(6,0)

#Merge for finak plotting
last <- merge(cumulative_stats_last, final, by = "weeknos")

#Generate transactions per week for plot
last$plot = last$transactions / last$finalcumf

figure_s5b <- ggplot(last, aes(x=weeknos, y=plot, col = )) +
  geom_line(size = 0.5) +
  labs(title = "", x="Week of program", y=  expression("Transactions per station"), fill = "white") +
  geom_line(aes(x = weeknos, y = plot), color = 'orange2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans")) +
  scale_x_continuous(breaks = seq(0, 45, by = 5)) +
  scale_y_continuous(breaks = seq(0, 2, by = 1)) +
  theme(legend.key=element_blank())

figure_s5b

figure_s5 <- grid.arrange(figure_s5a, figure_s5b, ncol=2)

# _____________________________________________________________________________
# Data Processing & Variable Creation
# _____________________________________________________________________________

#Create variable for lag in kWh per session by user
station.dat.4<-
  station.dat.4 %>%
  dplyr::group_by(userId) %>% 
  dplyr::mutate(Lag1 = dplyr::lag(kwhTotal,n = 1, default = NA))

#Natural log of the difference between natural log of the current and previous transaction
station.dat.4$delta.kwh.lag.ln<- log(station.dat.4$kwhTotal)-log(station.dat.4$Lag1)

#Remove all observartions of NA
sub<-which(is.na(station.dat.4$delta.kwh.lag.ln))
station.dat.4<-station.dat.4[-sub,]

#station.dat.4 <- station.dat.4[!duplicated(station.dat.4$userId),]
station.dat.4 <- ungroup(station.dat.4)

#Square and cubic terms
station.dat.4$charge3<-station.dat.4$chargeTimeHrs^3
station.dat.4$charge2<-station.dat.4$chargeTimeHrs^2

#Pull month out of datetime column
station.dat.4$month <- month(station.dat.4$created, label=TRUE)

#Dummies of months
month_new <- dummy(station.dat.4$month)
new <- as.data.frame(month_new)
new$sessionId <- station.dat.4$sessionId

#Merge dummies back to master
station.dat.4 <- merge(station.dat.4, new, by = "sessionId")




# _____________________________________________________________________________
# Robustness Checks
# _____________________________________________________________________________

#Calculate optimal bandwidth at cutpoints of 2 and 4 hours
obw.2<-IKbandwidth(X=station.dat.4$chargeTimeHrs, Y=station.dat.4$delta.kwh.lag.ln, 
                   cutpoint = 2,verbose =TRUE, kernel = "triangular")
obw.4<-IKbandwidth(X=station.dat.4$chargeTimeHrs, Y=station.dat.4$delta.kwh.lag.ln, 
                   cutpoint = 4,verbose =TRUE, kernel = "triangular")

#----------------------------------
# NEXT: Sharp RD, 4hrs
#----------------------------------

#No clusertering, no covariates (sharp, 4hrs)
RDest4_1 <-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = station.dat.4, cutpoint = 4,
                      verbose = TRUE, se.type='HC0')

#Clustering at facility with day of week covariates (sharp, 4hrs)
RDest4_3<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri , data = station.dat.4, cutpoint = 4, 
                     verbose = TRUE, cluster=station.dat.4$facilityType, se.type='HC0')

#No clustering with day of week covariates/month covariates (sharp, 4hrs)
RDest4_6<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 4, 
                     verbose = TRUE, se.type='HC0')

#No clustering with day of week covariates/month covariates, cubic term (sharp, 4hrs)
RDest4_7<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3 + Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 4, 
                     verbose = TRUE, se.type='HC0')

#Clustering at facility type, day of week covariates/month covariates (sharp, 4hrs)
RDest4_8<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 4, 
                     verbose = TRUE,cluster=station.dat.4$facilityType, se.type='HC0')

#Clustering at location ID, day of week covariates/month covariates (sharp, 4hrs)
RDest4_9<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 4, 
                     verbose = TRUE,cluster=station.dat.4$locationId, se.type='HC0')

#Clustering at station ID, day of week covariates/month covariates (sharp, 4hrs)
RDest4_10<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 4, 
                     verbose = TRUE,cluster=station.dat.4$stationId, se.type='HC0')

#----------------------------------
# NEXT: Sharp RD, 2hrs            
#----------------------------------

#No clusertering, no covariates (sharp, 2hrs)
RDest2_1 <-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = station.dat.4, cutpoint = 2,
                      verbose = TRUE, se.type='HC0')

#Clustering at facility with day of week covariates (sharp, 2hrs)
RDest2_3<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri , data = station.dat.4, cutpoint = 2, 
                     verbose = TRUE, cluster=station.dat.4$facilityType, se.type='HC0')

#No clustering, day of week covariates/month covariates (sharp, 2hrs)
RDest2_6<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 2, 
                     verbose = TRUE, se.type='HC0')

#No clustering, day of week covariates/month covariates, cubic term (sharp, 2hrs)
RDest2_7<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3 + Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 2, 
                     verbose = TRUE, se.type='HC0')

#Clustering at facility type, day of week covariates/month covariates (sharp, 2hrs)
RDest2_8<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 2, 
                     verbose = TRUE, cluster=station.dat.4$facilityType, se.type='HC0')

#Clustering at location ID, day of week covariates/month covariates (sharp, 2hrs)
RDest2_9<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 2, 
                     verbose = TRUE, cluster=station.dat.4$locationId, se.type='HC0')

#Clustering at station ID, day of week covariates/month covariates (sharp, 2hrs)
RDest2_10<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 2, 
                     verbose = TRUE, cluster=station.dat.4$stationId, se.type='HC0')



#----------------------------------
# NEXT: Sharp RD (Managers), 2hrs
#----------------------------------

#Subset by users who are charging with car given to managers
station.dat.5<-subset(station.dat.4, station.dat.4$managerVehicle == 1)

#Calculate optimal bandwidth for these users
obw.2_Managers<-IKbandwidth(X=station.dat.5$chargeTimeHrs, Y=station.dat.5$delta.kwh.lag.ln, 
                            cutpoint = 2,verbose =TRUE, kernel = "triangular")

#No clusertering, no covariates (sharp, 2hrs/managers only)
RDest_Man_1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = station.dat.5, cutpoint = 2,
                        verbose = TRUE, se.type='HC0')

#Clustering at facility with day of week covariates (sharp, 2hrs/managers only)
RDest_Man_3<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri , data = station.dat.5, cutpoint = 2, 
                        verbose = TRUE, cluster=station.dat.5$facilityType, se.type='HC0')

#No clustering with day of week covariates/monthly dummies (sharp, 2hrs/managers only)
RDest_Man_6<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.5, cutpoint = 2, 
                        verbose = TRUE, se.type='HC0')

#No clustering with day of week covariates/monthly dummies, cubic term (sharp, 2hrs/managers only)
RDest_Man_7<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3 + Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.5, cutpoint = 2, 
                        verbose = TRUE, se.type='HC0')

#Clustering at facility type, day of week covariates/monthly dummies (sharp, 2hrs/managers only)
RDest_Man_8<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.5, cutpoint = 2, 
                        verbose = TRUE, cluster=station.dat.5$facilityType, se.type='HC0')

#Clustering at location ID, day of week covariates/monthly dummies (sharp, 2hrs/managers only)
RDest_Man_9<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.5, cutpoint = 2, 
                        verbose = TRUE, cluster=station.dat.5$locationId, se.type='HC0')

#Clustering at station ID, day of week covariates/monthly dummies (sharp, 2hrs/managers only)
RDest_Man_10<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.5, cutpoint = 2, 
                        verbose = TRUE, cluster=station.dat.5$stationId, se.type='HC0')

#----------------------------------
# NEXT: Sharp RD (Managers), 4hrs
#----------------------------------

#Subset by users who are charging with car given to managers
station.dat.5<-subset(station.dat.4, station.dat.4$managerVehicle == 1)

#Calculate optimal bandwidth for these users
obw.4_Managers<-IKbandwidth(X=station.dat.5$chargeTimeHrs, Y=station.dat.5$delta.kwh.lag.ln, 
                            cutpoint = 4,verbose =TRUE, kernel = "triangular")

#No clusertering, no covariates (sharp, 2hrs/managers only)
RDest_Man4_1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = station.dat.5, cutpoint = 4,
                        verbose = TRUE, se.type='HC0')


#Clustering at facility with day of week covariates (sharp, 2hrs/managers only)
RDest_Man4_3<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri , data = station.dat.5, cutpoint = 4, 
                        verbose = TRUE, cluster=station.dat.5$facilityType, se.type='HC0')


#No clustering with day of week covariates/monthly covariates (sharp, 2hrs/managers only)
RDest_Man4_6<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov, data = station.dat.5, cutpoint = 4, 
                         verbose = TRUE, se.type='HC0')

#No clustering with day of week covariates/monthly covariates, cubic term (sharp, 2hrs/managers only)
RDest_Man4_7<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3 + Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov, data = station.dat.5, cutpoint = 4, 
                         verbose = TRUE, se.type='HC0')

#Clustering at facility type, day of week covariates/monthly covariates (sharp, 2hrs/managers only)
RDest_Man4_8<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov, data = station.dat.5, cutpoint = 4, 
                         verbose = TRUE, cluster=station.dat.5$facilityType, se.type='HC0')

#Clustering at location ID, day of week covariates/monthly covariates (sharp, 2hrs/managers only)
RDest_Man4_9<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov, data = station.dat.5, cutpoint = 4, 
                         verbose = TRUE, cluster=station.dat.5$locationId, se.type='HC0')

#Clustering at station ID, day of week covariates/monthly covariates (sharp, 2hrs/managers only)
RDest_Man4_10<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov, data = station.dat.5, cutpoint = 4, 
                         verbose = TRUE, cluster=station.dat.5$stationId, se.type='HC0')

#----------------------------------
# NEXT: Placebo (Sharp)
#----------------------------------

#Calculate optimal bandwidth at cutpoints of 3 hours
obw.3<-IKbandwidth(X=station.dat.4$chargeTimeHrs, Y=station.dat.4$delta.kwh.lag.ln, 
                   cutpoint = 3,verbose =TRUE, kernel = "triangular")

#No clusertering, no covariates (sharp, 3hrs/placebo)
RDest.p1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = station.dat.4, cutpoint = 3,
                     verbose = TRUE, se.type='HC0')

#Clustering at facility with day of week covariates (sharp, 3hrs/placebo)
RDest.p3<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri , data = station.dat.4, cutpoint = 3, 
                     verbose = TRUE, cluster=station.dat.4$facilityType, se.type='HC0')

#No clustering with day of week covariates/month covariates (sharp, 3hrs/placebo)
RDest.p6<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 3, 
                     verbose = TRUE, se.type='HC0')

#No clustering with day of week covariates/month covariates, cubic term (sharp, 3hrs/placebo)
RDest.p7<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3 + Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 3, 
                     verbose = TRUE, se.type='HC0')

#Clustering at facility type, day of week covariates/month covariates (sharp, 3hrs/placebo)
RDest.p8<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 3, 
                     verbose = TRUE, cluster=station.dat.4$facilityType, se.type='HC0')

#Cluster at location ID, day of week covariates/month covariates (sharp, 3hrs/placebo)
RDest.p9<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 3, 
                     verbose = TRUE, cluster=station.dat.4$locationId, se.type='HC0')

#Cluster at station ID, day of week covariates/month covariates (sharp, 3hrs/placebo)
RDest.p10<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri + monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov , data = station.dat.4, cutpoint = 3, 
                     verbose = TRUE, cluster=station.dat.4$stationId, se.type='HC0')



#----------------------------------
# FINAL: Display results
#----------------------------------

#Create matrix to display robustness table
table_s1<- matrix(nrow = 15,ncol = 6, dimnames = list(c("Cutoff 4 hours Sharp",
                                                       "(s.e. 1)",
                                                       "Cutoff 4 hours Sharp - Managers",
                                                       "(s.e. 2)",
                                                       "Cutoff 2 hours Sharp",
                                                       "(s.e. 3)",
                                                       "Cutoff 2 hours Sharp - Managers",
                                                       "(s.e. 4)",
                                                       "Cutoff 3 hours Sharp - Placebo",
                                                       "(s.e. 5)",
                                                       "Time dummies",
                                                       "Cubic polynomial",
                                                       "Clustering at facility type",
                                                       "Clustering at location ID",
                                                       "Clustering at station ID"
                                                       ),
                                                     c("(1)",
                                                       "(2)",
                                                       "(3)",
                                                       "(4)",
                                                       "(5)",
                                                       "(6)")))


table_s1[1, ] <- c(round(RDest4_1$est[1],4),round(RDest4_6$est[1],4),round(RDest4_7$est[1],4),round(RDest4_8$est[1],4),round(RDest4_9$est[1],4),round(RDest4_10$est[1],4))
table_s1[2, ] <- c(round(RDest4_1$se[1],4),round(RDest4_6$se[1],4),round(RDest4_7$se[1],4),round(RDest4_8$se[1],4),round(RDest4_9$se[1],4),round(RDest4_10$se[1],4))
table_s1[3, ] <- c(round(RDest_Man4_1$est[1],4),round(RDest_Man4_6$est[1],4),round(RDest_Man4_7$est[1],4),round(RDest_Man4_8$est[1],4),round(RDest_Man4_9$est[1],4),round(RDest_Man4_10$est[1],4))
table_s1[4, ] <- c(round(RDest_Man4_1$se[1],4),round(RDest_Man4_6$se[1],4),round(RDest_Man4_7$se[1],4),round(RDest_Man4_8$se[1],4),round(RDest_Man4_9$se[1],4),round(RDest_Man4_10$se[1],4))
table_s1[5, ] <- c(round(RDest2_1$est[1],4),round(RDest2_6$est[1],4),round(RDest2_7$est[1],4),round(RDest2_8$est[1],4),round(RDest2_9$est[1],4),round(RDest2_10$est[1],4))
table_s1[6, ] <- c(round(RDest2_1$se[1],4),round(RDest2_6$se[1],4),round(RDest2_7$se[1],4),round(RDest2_8$se[1],4),round(RDest2_9$se[1],4),round(RDest2_10$se[1],4))
table_s1[7, ] <- c(round(RDest_Man_1$est[1],4),round(RDest_Man_6$est[1],4),round(RDest_Man_7$est[1],4),round(RDest_Man_8$est[1],4),round(RDest_Man_9$est[1],4),round(RDest_Man_10$est[1],4))
table_s1[8, ] <- c(round(RDest_Man_1$se[1],4),round(RDest_Man_6$se[1],4),round(RDest_Man_7$se[1],4),round(RDest_Man_8$se[1],4),round(RDest_Man_9$se[1],4),round(RDest_Man_10$se[1],4))
table_s1[9, ] <- c(round(RDest.p1$est[1],4),round(RDest.p6$est[1],4),round(RDest.p7$est[1],4),round(RDest.p8$est[1],4),round(RDest.p9$est[1],4),round(RDest.p10$est[1],4))
table_s1[10, ] <- c(round(RDest.p1$se[1],4),round(RDest.p6$se[1],4),round(RDest.p7$se[1],4),round(RDest.p8$se[1],4),round(RDest.p9$se[1],4),round(RDest.p10$se[1],4))
table_s1[11, ] <- c("No","Yes","Yes","Yes","Yes","Yes")
table_s1[12, ] <- c("No","No","Yes","No","No","No")
table_s1[13, ] <- c("No","No","No","Yes","No","No")
table_s1[14, ] <- c("No","No","No","No","Yes","No")
table_s1[15, ] <- c("No","No","No","No","No","Yes")



table_s1

# _____________________________________________________________________________
# RD Results Table (3)
# _____________________________________________________________________________

table_3<- matrix(nrow = 5,ncol = 4, dimnames = list(c("Price effect, (4 hours), all users","Price effect (4 hours), managers","Behavioral effect (2 hours), all users"," Behavioral effect (2 hours), managers"," Placebo test (3 hours), all users"),
                                                    c("Optimal Bandwidth","RD Estimate", "Std Error", "Total sessions")))

table_3[1, ] <- c(round(obw.4,4), round(RDest4_3$est[1],4), round(RDest4_3$se[1],4), nrow(station.dat.4))
table_3[2, ] <- c(round(obw.4_Managers,4), round(RDest_Man4_3$est[1],4), round(RDest_Man4_3$se[1],4), nrow(station.dat.5))
table_3[3, ] <- c(round(obw.2,4), round(RDest2_3$est[1],4), round(RDest2_3$se[1],4), nrow(station.dat.4))
table_3[4, ] <- c(round(obw.2_Managers,4), round(RDest_Man_3$est[1],4), round(RDest_Man_3$se[1],4), nrow(station.dat.5))
table_3[5, ] <- c(round(obw.3,4), round(RDest.p3$est[1],4), round(RDest.p3$se[1],4), nrow(station.dat.4))

table_3

# _____________________________________________________________________________
# Figure (1): PRICE POLICY GRAPH
# _____________________________________________________________________________

#Specify hour limits and increments of 5 minutes (marginal cost incurred every 5 min after 4hrs)
x <- seq(0, 6, (1/12))

#Specify function corresponding to pricing scheme
fx <- (x > 0 & x <=4) *0+
  (x >4 & x < 4.5) * 0.5 +
  (x >= 4.5 & x < 10.5) * (x-4)

par(mar=c(5,5,1,1)+.1)
plot(x, fx, type="S", xlab="Duration of charge (hrs)",
     ylab="Cost (dollars)",cex.lab=1.5, cex.axis=1.3)
abline(v=4, lty=2)

figure_1 <- recordPlot()

figure_1

# _____________________________________________________________________________
# Figure (S1): HISTOGRAM OF TRANSACTION COUNT BY EV USER
# _____________________________________________________________________________

#Generate histogram using dataframe with count of transactions by user ID previously generated
figure_s3 <- ggplot(station.dat.agg2, aes(x=n))+
  geom_histogram(color="black", fill="orange2", binwidth = 10, alpha = 0.4)+
  labs(x="Number of transactions", y="Users")+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

figure_s3

# _____________________________________________________________________________
# Figure (4a): HISTOGRAM FOR PLUG-IN/PLUG-OUT (OF OBSERVATIONS CONSIDERED IN STUDY)
# _____________________________________________________________________________

figure_s4a <-ggplot(station.dat.4, aes(x = startTime, color=cols))+ 
  labs(title = "", x="Time of Day (hr)", y="Frequency of Transactions")+
  scale_x_discrete(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = 0:24)+
  scale_y_discrete(breaks = c(100,200,300,400,500), limits = c(0:500)) +
  geom_bar(aes(x=startTime,  fill="orange2"), color="orange3", alpha = 0.4, position = position_nudge(x = 0.5)) + 
  geom_bar(aes(x=endTime, fill="grey45"),color="grey30", alpha = 0.4, position = position_nudge(x = 0.5))+
  scale_fill_identity(name='',guide = 'legend',labels = c('Plug-out Times','Plug-in Times'))+
  theme_light()+
  theme(legend.position = c(0.05, 1),legend.justification = c(0, 1),legend.text = element_text(size=15))+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

figure_s4a

# _____________________________________________________________________________
# Figure (4b): HISTOGRAM FOR PLUG-IN/PLUG-OUT (OF EXTERNAL SAMPLE)
# _____________________________________________________________________________

figure_s4b <-ggplot(external.sample, aes(x = startTime, color=cols))+ 
  labs(title = "", x="Time of Day (hr)", y="Frequency of Transactions")+
  scale_x_discrete(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22,24),limits = 0:24)+
  scale_y_discrete(breaks = c(100,200), limits = c(0:200)) +
  geom_bar(aes(x=startTime,  fill="orange2"), color="orange3", alpha = 0.4, position = position_nudge(x = 0.5)) + 
  geom_bar(aes(x=endTime, fill="grey45"),color="grey30", alpha = 0.4, position = position_nudge(x = 0.5))+
  scale_fill_identity(name='',guide = 'legend',labels = c('Plug-out Times','Plug-in Times'))+
  theme_light()+ 
  theme(legend.position = c(.1, 1),legend.justification = c(0, 1),legend.text = element_text(size=15))+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

figure_s4b

figure_s4 <- grid.arrange(figure_s4a, figure_s4b, ncol=2)

# _____________________________________________________________________________
# Creation of generalizable material used in creation of figures 2 and 3.
# _____________________________________________________________________________

#Use n=6 for 30 days intervals, n=12 for 15 days intervals, n=25 for 7 days intervals
#One less than number of observations used for later calculations in matrix
n=169

#Use d=30 for 30 days intervals, d=15 for 15 days intervals, d=7 for 7 days intervals
d=1

# _____________________________________________________________________________
# Figure (2a): ESTIMATE OF MAIN SPECIFICATION USING DIFFERENT BANDWIDTHS (4hrs)
# _____________________________________________________________________________

#Generate row of matrix per month
r=8

#Matrix used to generate estimates for dynamic and main specification models with
#one column per metric needed in generating figures
mat_2a <- data.frame(matrix(0, ncol=6, nrow=r))
names(mat_2a)<-c("Month","Estimate","se", "p-generation", "lower ci", "upper ci")


#RDestimate
bndwdth<-c(0.25, 0.5,  0.75,  1,  1.5, 2, 2.5, 3)
for(i in 1:8){
  RRDD<-RDest.4<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs|monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov + Mon + Tues + Wed + Thurs + Fri, data = station.dat.4, bw = obw.4*(bndwdth[i]), cutpoint = 4, 
                            verbose = TRUE, cluster=station.dat.4$locationId, se.type='HC0')
  mat_2a[i,1]<-bndwdth[i]
  mat_2a[i,2]<-RRDD$ci[1]
  mat_2a[i,3]<-RRDD$ci[4]
  mat_2a[i,4]<-RRDD$est[1]
}

mat_2a[,1] <- mat_2a[,1]*100

figure_2a <- ggplot(mat_2a,  aes(mat_2a[,1])) + 
  labs(title = "", x="% of I-K Optimal Bandwidth", y=  expression("Estimate of RD Coefficient"))+
  geom_line(aes(y = mat_2a[,2], col = "CI_lower"), linetype= "dashed",colour="grey30", size = 1) + 
  geom_line(aes(y = mat_2a[,4], col = "estimate"),colour="black", size = 1)+
  geom_line(aes(y = mat_2a[,3], col = "CI_upper"), linetype="dashed",colour="grey30", size = 1)+
  scale_x_continuous(breaks = seq(0,300 , by = 50))+
  coord_cartesian(ylim=c(-0.4,0.1)) +
  geom_hline(yintercept=0, size = 1)+
  geom_ribbon(aes(ymin=mat_2a[,2], ymax=mat_2a[,3]), linetype=2, alpha=0.3, fill = "skyblue3")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=11, 
                                  face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans"))

figure_2a

# _____________________________________________________________________________
# Figure (2b): ESTIMATE OF MAIN SPECIFICATION USING DIFFERENT BANDWIDTHS (2hrs)
# _____________________________________________________________________________

#Generate row of matrix per month
r=8

date1<-as.Date("0014-11-18")
date2<-as.Date("0015-04-17")

#Matrix used to generate estimates for dynamic and main specification models with
#one column per metric needed in generating figures
mat_2b <- data.frame(matrix(0, ncol=6, nrow=r))
names(mat_2b)<-c("Month","Estimate","se", "p-generation", "lower ci", "upper ci")

#RDestimate
bndwdth<-c(0.25, 0.5,  0.75,  1,  1.5, 2, 2.5, 3)
for(i in 1:8){
  RRDD<-RDest.4<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs|monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov + Mon + Tues + Wed + Thurs + Fri, data = station.dat.4, bw = obw.2*(bndwdth[i]), cutpoint = 2, 
                            verbose = TRUE, cluster=station.dat.4$locationId, se.type='HC0')
  mat_2b[i,1]<-bndwdth[i]
  mat_2b[i,2]<-RRDD$ci[1]
  mat_2b[i,3]<-RRDD$ci[4]
  mat_2b[i,4]<-RRDD$est[1]
}

mat_2b[,1] <- mat_2b[,1]*100

figure_2b <- ggplot(mat_2b,  aes(mat_2b[,1])) + 
  labs(title = "", x="% of I-K Optimal Bandwidth", y=  expression("Estimate of RD Coefficient"))+
  geom_line(aes(y = mat_2b[,2], col = "CI_lower"), linetype= "dashed",colour="grey30", size = 1) + 
  geom_line(aes(y = mat_2b[,4], col = "estimate"),colour="black", size = 1)+
  geom_line(aes(y = mat_2b[,3], col = "CI_upper"), linetype="dashed",colour="grey30", size = 1)+
  scale_x_continuous(breaks = seq(0,300 , by = 50))+
  scale_y_continuous(breaks = seq(-0.8,.2 , by = 0.2))+
  geom_ribbon(aes(ymin=mat_2b[,2], ymax=mat_2b[,3]), linetype=2, alpha=0.3, fill = "salmon")+
  theme_bw()+
  coord_cartesian(ylim=c(-0.8,0.2)) +
  geom_hline(yintercept=0, size = 1)+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=11, 
                                  face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans"))

figure_2b


figure_2 <- grid.arrange(figure_2a, figure_2b, ncol=2)


# _____________________________________________________________________________
# Figure (3a): DYNAMIC ESTIMATES (4-HOUR CUTOFF)
# _____________________________________________________________________________
#Change according to the intervals. 7 rows for months, 13 for 15 days, 26 for 7 days, 182 for 1 day interval?
#Uses the last 170 days of observations to calculate estimates
r=170

#Matrix used to generate estimates for dynamic and main specification models with
#one column per metric needed in generating figures
mat_3a <- data.frame(matrix(0, ncol=6, nrow=r))
names(mat_3a)<-c("End date","Estimate","se", "p-generation", "lower ci", "upper ci")

date1<-as.Date("0014-11-18")
date2<-as.Date("0015-04-17")

subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]



ctpt<-4
#RDestimate with clustering by facility type and including day of week controls
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs |  monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov + Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

mat_3a[1,1]<-date2
mat_3a[1,2]<-RDest1$est[1]
mat_3a[1,3]<-RDest1$se[1]
mat_3a[1,4]<-RDest1$p[1]
mat_3a[1,5]<-RDest1$ci[1]
mat_3a[1,6]<-RDest1$ci[4]
mat_3a[1,7]<-150/7


for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs |  monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov + Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  mat_3a[1+i,1]<-date2
  mat_3a[1+i,2]<-RDest1$est[1]
  mat_3a[1+i,3]<-RDest1$se[1]
  mat_3a[1+i,4]<-RDest1$p[1]
  mat_3a[1+i,5]<-RDest1$ci[1]
  mat_3a[1+i,6]<-RDest1$ci[4]
  mat_3a[1+i,7]<-mat_3a[i,7]+1/7
}


figure_3a <- ggplot(mat_3a, aes(mat_3a[,7], mat_3a[,2])) + 
  labs(x="Weeks since start of program", y=  expression("Estimate of RD Coefficient"))+
  geom_line(aes(y = mat_3a[,5], col = "CI_lower"), linetype= "dashed",colour="grey30", size = 0.5) + 
  geom_line(aes(y = mat_3a[,2], col = "estimate"),colour="black", size = 0.5)+
  geom_line(aes(y = mat_3a[,6], col = "CI_upper"), linetype="dashed",colour="grey30", size = 0.5)+
  scale_x_continuous(breaks = seq(20,50, by = 2), limits=c(20,50),
                     sec.axis = sec_axis(~ . / 4, name = "Months since start of program", breaks = seq(5,12,1)))+
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1), labels=abbreviate)+
  geom_hline(yintercept=0, size = 1)+
  geom_ribbon(aes(ymin=mat_3a[,5], ymax=mat_3a[,6]), linetype=2, alpha=0.3, fill = "skyblue3")+
  coord_cartesian(xlim=c(22,44.75)) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=11, 
                                  face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans"))

figure_3a

# _____________________________________________________________________________
# Figure (3b): DYNAMIC ESTIMATES (2-HOUR CUTOFF)
# _____________________________________________________________________________

#Change according to the intervals. 7 rows for months, 13 for 15 days, 26 for 7 days, 182 for 1 day interval?
#Uses the last 170 days of observations to calculate estimates
r=170

#Matrix used to generate estimates for dynamic and main specification models with
#one column per metric needed in generating figures
mat_3b <- data.frame(matrix(0, ncol=6, nrow=r))
names(mat_3b)<-c("End date","Estimate","se", "p-value", "lower ci", "upper ci")

#First observation is on 11.18.2014
date1<-as.Date("0014-11-18")

#Last observations is on 10.04.2015
date2<-as.Date("0015-04-17")

as.Date("0015-04-17")-as.Date("0014-11-18")

#Subset only those observations ocurring between designated start and end dates
subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

#Set desired cutpoint to 2hrs
ctpt<<-2
#Regression discontinuity using delta.kwh.lag.ln and clustered by facility/day of week
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | monthJan + monthFeb + monthMar + monthJun + monthJul + monthAug + monthSep + monthOct + monthNov + Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

#initialize first row of matrix to reflect appropriate end date and relevant
#outputs from above RD
mat_3b[1,1]<-date2
mat_3b[1,2]<-RDest1$est[1]
mat_3b[1,3]<-RDest1$se[1]
mat_3b[1,4]<-RDest1$p[1]
mat_3b[1,5]<-RDest1$ci[1]
mat_3b[1,6]<-RDest1$ci[4]
mat_3b[1,7]<-150/7
for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | monthJan + monthFeb + monthMar + monthApr + monthMay+ monthJun + monthJul + monthAug + monthSep + monthOct + monthNov + Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  mat_3b[1+i,1]<-date2
  mat_3b[1+i,2]<-RDest1$est[1]
  mat_3b[1+i,3]<-RDest1$se[1]
  mat_3b[1+i,4]<-RDest1$p[1]
  mat_3b[1+i,5]<-RDest1$ci[1]
  mat_3b[1+i,6]<-RDest1$ci[4]
  mat_3b[1+i,7]<-mat_3b[i,7]+1/7
}

figure_3b <- ggplot(mat_3b, aes(mat_3b[,7], mat_3b[,2])) + 
  labs(x="Weeks since start of program", y=  expression("Estimate of RD Coefficient"))+
  geom_line(aes(y = mat_3b[,5], col = "CI_lower"), linetype= "dashed",colour="grey30",size = 0.5) + 
  geom_line(aes(y = mat_3b[,2], col = "estimate"),colour="black", size = 0.5)+
  geom_line(aes(y = mat_3b[,6], col = "CI_upper"), linetype="dashed",colour="grey30", size = 0.5)+
  scale_x_continuous(breaks = seq(20,50, by = 2), limits=c(20,50),
                     sec.axis = sec_axis(~ . / 4, name = "Months since start of program", breaks = seq(5,12,1)))+
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1), labels=abbreviate)+
  geom_hline(yintercept=0, size = 1)+
  geom_ribbon(aes(ymin=mat_3b[,5], ymax=mat_3b[,6]), linetype=2, alpha=0.2, fill = "salmon")+
  theme_bw()+
  coord_cartesian(xlim=c(22,44.75)) +
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=11, 
                                  face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans"))

figure_3b

figure_3 <- grid.arrange(figure_3a, figure_3b, ncol=2)

# _____________________________________________________________________________
# Figure S1: RD FIGURES
# _____________________________________________________________________________

station.dat.4$over4 = station.dat.4$chargeTimeHrs>4
figure_s1a <- ggplot(station.dat.4, aes(x = chargeTimeHrs, y = delta.kwh.lag.ln, color = over4)) + 
  geom_point(alpha = 0.4, stroke = 0.4, size=3) +
  geom_vline(xintercept=4, linetype="dashed") +
  coord_cartesian(ylim=c(-1,1),xlim=c(1, 5.5)) +
  stat_smooth(method="loess",formula = y~x, fill= "grey40", size = 1) +
  scale_x_continuous(breaks = seq(1, 5.5, by = 1)) +
  scale_colour_manual(values = c("grey5", "skyblue3")) +
  labs(title = "", x="Charge Time (hrs)", y=  expression("Change in Log of kWh Lag"), color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none", axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans"))

figure_s1a

station.dat.4$over2 = station.dat.4$chargeTimeHrs>2
figure_s1b <- ggplot(station.dat.4, aes(x = chargeTimeHrs, y = delta.kwh.lag.ln, color = over2)) + 
  geom_point(alpha = 0.4, stroke = 0.4, size=3) +
  geom_vline(xintercept=2, linetype="dashed") +
  coord_cartesian(ylim=c(-1,1),xlim=c(1, 5)) +
  stat_smooth(method="loess",formula = y~x, fill= "grey30", size = 1) +
  scale_colour_manual(values = c("grey5", "salmon")) +
  labs(title = "", x="Charge Time (hrs)", y=  expression("Change in Log of kWh Lag")) +
  scale_x_continuous(breaks = seq(1, 5, by = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none", axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(family = "sans"))

figure_s1b

figure_s1 <- grid.arrange(figure_s1a, figure_s1b, ncol=2)


#***************************************************************************
#Dynamic RD table
#***************************************************************************

table_4 <- matrix(nrow = 20,ncol = 6, dimnames = list(c("Month 4",
                                                        "SE4",
                                                        "Month 5",
                                                        "SE5",
                                                        "Month 6",
                                                        "SE6",
                                                        "Month 7", 
                                                        "SE7",
                                                        "Month 8", 
                                                        "SE8", 
                                                        "Month 9",
                                                        "SE9",
                                                        "Month 10",
                                                        "SE10",
                                                        "Month 11",
                                                        "SE11",
                                                        "Month 12",
                                                        "SE12",
                                                        "Day of the week dummies",
                                                        "Cube charge time"),
                                                      c("1", "2", "3", "4", "5","6")))


#MODEL (1)

#Testing period: first 3 months
date1<-as.Date("0014-11-18")
date2<-as.Date("0015-02-18")

#Creating charge time^3
station.dat.4$charge3<-station.dat.4$chargeTimeHrs^3
station.dat.4$charge2<-station.dat.4$chargeTimeHrs^2

subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

ctpt<-2

#RDestimate
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs , data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

summary(RDest1)

#9 months of the program excluding testing period
r=9

table_4[1,1]<-round(RDest1$est[1],3)
table_4[2,1]<-round(RDest1$se[1],3)


n=r-1

#Use d=30 for 30 days intervals
d=30

for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  summary(RDest1)
  
  table_4[1+2*i,1]<-round(RDest1$est[1],3)
  table_4[2+2*i,1]<-round(RDest1$se[1],3)
  
}





#MODEL (2)

#Testing period: first 3 months
date1<-as.Date("0014-11-18")
date2<-as.Date("0015-02-18")


subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

ctpt<-2

#RDestimate
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

summary(RDest1)

#9 months of the program excluding testing period
r=9

table_4[1,2]<-round(RDest1$est[1],3)
table_4[2,2]<-round(RDest1$se[1],3)


n=r-1

#Use d=30 for 30 days intervals
d=30

for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs| Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  table_4[1+2*i,2]<-round(RDest1$est[1],3)
  table_4[2+2*i,2]<-round(RDest1$se[1],3)
}



#MODEL (3)

#Testing period: first 3 months
date1<-as.Date("0014-11-18")
date2<-as.Date("0015-02-18")


subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

ctpt<-2

#RDestimate
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3+  Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

summary(RDest1)

#9 months of the program excluding testing period
r=9

table_4[1,3]<-round(RDest1$est[1],3)
table_4[2,3]<-round(RDest1$se[1],3)


n=r-1

#Use d=30 for 30 days intervals
d=30

for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs| charge3+Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  table_4[1+2*i,3]<-round(RDest1$est[1],3)
  table_4[2+2*i,3]<-round(RDest1$se[1],3)
  
  summary(RDest1)
}



#MODEL (4)

#Testing period: first 3 months
date1<-as.Date("0014-11-18")
date2<-as.Date("0015-02-18")

#Creating charge time^3
station.dat.4$charge3<-station.dat.4$chargeTimeHrs^3
station.dat.4$charge2<-station.dat.4$chargeTimeHrs^2

subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

ctpt<-4

#RDestimate
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs , data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

summary(RDest1)

#9 months of the program excluding testing period
r=9

table_4[1,4]<-round(RDest1$est[1],3)
table_4[2,4]<-round(RDest1$se[1],3)


n=r-1

#Use d=30 for 30 days intervals
d=30

for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  table_4[1+2*i,4]<-round(RDest1$est[1],3)
  table_4[2+2*i,4]<-round(RDest1$se[1],3)
}





#MODEL (5)

#Testing period: first 3 months
date1<-as.Date("0014-11-18")
date2<-as.Date("0015-02-18")


subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

ctpt<-4

#RDestimate
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

summary(RDest1)

#9 months of the program excluding testing period
r=9

table_4[1,5]<-round(RDest1$est[1],3)
table_4[2,5]<-round(RDest1$se[1],3)


n=r-1

#Use d=30 for 30 days intervals
d=30

for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs| Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  table_4[1+2*i,5]<-round(RDest1$est[1],3)
  table_4[2+2*i,5]<-round(RDest1$se[1],3)
}



#MODEL (6)

#Testing period: first 3 months
date1<-as.Date("0014-11-18")
date2<-as.Date("0015-02-18")


subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]

ctpt<-4

#RDestimate
RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs | charge3+  Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                   verbose = TRUE, cluster=subset1$locationId, se.type='HC0')

summary(RDest1)

#9 months of the program excluding testing period
r=9

table_4[1,6]<-round(RDest1$est[1],3)
table_4[2,6]<-round(RDest1$se[1],3)


n=r-1

#Use d=30 for 30 days intervals
d=30

for(i in 1:n)
{
  date2<-date2+d
  subset1<-station.dat.4[station.dat.4$created >= date1 & station.dat.4$created <= date2,]
  
  RDest1<-RDestimate(delta.kwh.lag.ln~chargeTimeHrs| charge3+Mon + Tues + Wed + Thurs + Fri, data = subset1, cutpoint = ctpt, 
                     verbose = TRUE, cluster=subset1$locationId, se.type='HC0')
  
  table_4[1+2*i,6]<-round(RDest1$est[1],3)
  table_4[2+2*i,6]<-round(RDest1$se[1],3)
  
  summary(RDest1)
}

#Input model specifications
table_4[19,] <- c("No", "Yes", "Yes", "No", "Yes", "Yes")
table_4[20,] <- c("No", "No", "Yes", "No", "No", "Yes")

