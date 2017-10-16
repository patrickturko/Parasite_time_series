
setwd("/media/FILES/Dropbox/Turko PhD Data/Turko et al Chapter 5 Greifensee/Time Series")
# Load packages and functions ---------------------------------------------


library(fANCOVA)
library(reshape2)
library(plyr)
library(R.oo)
library(R.utils)
library(ggplot2)
library(gridExtra)
library(scales)
source("webster.R")
source("heat.R")
heat_ci <- function(obj, lag = 0){
  a <- mean(obj$data$value[obj$data$variable == lag], na.rm = T)
  b <- sd(obj$data$value[obj$data$variable == lag], na.rm = T)
  qnorm(c(0.025,0.975), a, b)
  
}
# Load data ---------------------------------------------------------------

allo_mlg <- read.csv("inputs/allozyme_mlgs_date.csv")
allo_spp <- read.csv("inputs/allozyme_spp_date.csv", sep = ";")
msat_mlg <- read.csv("inputs/msat_mlg_time.csv", sep = ";")
msat_spp <- read.csv("inputs/msat_spp_time.csv", sep = ";")
  
d.demo <- read.csv("inputs/Demography.csv", sep =";")
algae <- read.csv("inputs/algae_biomass.csv", sep = " ")
caullerya <- read.csv("inputs/caullerya.txt", sep = "\t")
d.temperature <- read.csv("inputs/temperature_1942_2014.csv", sep = ";")
d.oxygen <- read.csv("inputs/O2.txt", sep = "\t")

allo_date_sums <- data.frame(date = allo_mlg$date, tot = rowSums(allo_mlg[, 2:ncol(allo_mlg)]))
# Fix dates of input files ------------------------------------------------

allo_mlg$date <- as.Date(allo_mlg$date)
allo_spp$date <- as.Date(allo_spp$date)
msat_mlg$date <- as.Date(msat_mlg$date)
msat_spp$date <- as.Date(msat_spp$date)
algae$date <- as.Date(algae$date)
# Set day zero ------------------------------------------------------------

min_date <- min(allo_mlg$date)
# Transform DEMOGRAPHY data  ------------------------------------------------

# Demography

d.demo <- within(d.demo, {
  d..day.1.[data_state!=""] <- NA
  b..day.1.[data_state!=""] <- NA
})

yrs <- as.Date(d.demo$Date, "%d.%m.%Y")
yrs <- data.frame(min=yrs[-length(yrs)], max=yrs[-1], diff=diff.Date(yrs))
yrs <- subset(yrs, diff>55)
yrs$n <- c(1:nrow(yrs))

#r
#-
d.r <- d.demo[,c("Date", "r..day.1.")]
d.r <- within(d.r, {date <- as.Date(Date, "%d.%m.%y")
                    value <- as.numeric(as.character(r..day.1.))
})
d.r <- d.r[!is.na(d.r$value),]

#d
#-
d.d <- d.demo[,c("Date", "d..day.1.")]
d.d <- within(d.d, {date <- as.Date(Date, "%d.%m.%y")
                    value <- as.numeric(as.character(d..day.1.))
})
d.d <- d.d[!is.na(d.d$value),]

#b
#-
d.b <- d.demo[,c("Date", "b..day.1.")]
d.b <- within(d.b, {date <- as.Date(Date, "%d.%m.%y")
                    value <- as.numeric(as.character(b..day.1.))
})
d.b <- d.b[!is.na(d.b$value),]

#n
#-
d.n <- d.demo[,c(1, 4:7)]
d.n$n <- rowSums(d.n[, 2:5])
d.n <- within(d.n, {date <- as.Date(Date, "%d.%m.%y")
                    value <- as.numeric(as.character(n))
})
d.n <- d.n[!is.na(d.n$value),]
# Transform ALGAE data ----------------------------------------------------


algae_rel <- dcast(algae, date ~ taxon2, value.var = "biomass")
algae_rel[, 2:ncol(algae_rel)] <- algae_rel[, 2:ncol(algae_rel)] / rowSums(algae_rel[, 2:ncol(algae_rel)])
algae_rel <- melt(algae_rel, value.name = "biomass", variable.name = "taxon2", id.vars = "date")
algae_rel[is.na(algae_rel)] <- 0
algae_sum <- ddply(algae, .(date), function(x) sum(x$biomass))

algae_matrix <- dcast(algae, date ~ taxon2)
algae_matrix[is.na(algae_matrix)] <- 0
algae_rel_matrix <- dcast(algae_rel, date ~ taxon2)
algae_rel_matrix <- algae_rel_matrix[rowSums(algae_rel_matrix[, 2:ncol(algae_rel_matrix)]) != 0, ]

algae_cyan <- algae_matrix[, c("date", "Cyanobacteria")]
algae_cyan$Cyanobacteria[is.na(algae_cyan$Cyanobacteria)] <- 0
# Transform CAULLERYA data ------------------------------------------------

caullerya$date <- as.Date(caullerya$date2, format = "%d.%m.%Y")
caullerya$value <- caullerya$caullerya.prevalence
# Transform TEMPERATURE data ----------------------------------------------

d.temperature <- within(d.temperature, date <- as.Date(date, "%d.%m.%Y"))
d.temperature <- subset(d.temperature, depth <=15 & date > min_date)
d.temp_aggr <- ddply(d.temperature, .(date), summarize, value = mean(temp, na.rm=T), n = length(temp))
temperature <- d.temp_aggr
# Transform OXYGEN data ---------------------------------------------------

names(d.oxygen) <- gsub(pattern="X",replacement="",x=names(d.oxygen))
d.oxygen <- melt(d.oxygen,na.rm=F)
d.oxygen <- within(d.oxygen, {date <- as.Date(Datum, "%d.%m.%Y")
                              variable <- as.numeric(as.character(variable))
})
d.oxygen <- subset(d.oxygen, variable <=15 & date > min_date)
oxygen <- ddply(d.oxygen, .(date), summarise, value = mean(value, na.rm=T), n = length(value))
# Turnover for DAPHNIA ----------------------------------------------------

window = 120

allo_mlg_turn <- webster(allo_mlg[, 2:ncol(allo_mlg)], 
                         dates = allo_mlg$date, 
                         window = window, dist = "bray",
                         window_type = "days", rarefy = TRUE)
allo_mlg_turn <- allo_mlg_turn[!is.na(allo_mlg_turn$dist), ]

allo_spp_turn <- webster(allo_spp[, 2:ncol(allo_spp)], 
                         dates = allo_spp$date, 
                         window = window, dist = "bray",
                         window_type = "days", rarefy = TRUE)
allo_spp_turn <- allo_spp_turn[!is.na(allo_spp_turn$dist), ]


msat_mlg_turn <- webster(msat_mlg[, 3:ncol(msat_mlg)], 
                         dates = msat_mlg$date, 
                         window = window, dist = "bray",
                         window_type = "days", rarefy = TRUE)
msat_mlg_turn <- msat_mlg_turn[!is.na(msat_mlg_turn$dist), ]


msat_spp_turn <- webster(msat_spp[, 2:ncol(msat_spp)], 
                         dates = msat_spp$date, 
                         window = window, dist = "bray",
                         window_type = "days", rarefy = TRUE)
msat_spp_turn <- msat_spp_turn[!is.na(msat_spp_turn$dist), ]
# Decompose ALLO MLG series -----------------------------------------------

mod1 <- lm(data = allo_mlg_turn, dist ~ date)
plot(allo_mlg_turn$date, allo_mlg_turn$dist)
abline(mod1)
summary(mod1)
anova(mod1)
allo_mlg_turn$dist2 <- residuals(mod1)

smp_gg <- ggplot(msat_mlg_turn, aes(x = date, y = all_sum))+
  geom_line() +
  geom_smooth(method = "lm")
smp_gg
# Turnover for ALGAE ------------------------------------------------------

window = 90

algae_turn <- webster(algae_matrix[, 2:ncol(algae_matrix)], 
                         dates = algae_matrix$date, 
                         window = window, dist = "bray",
                         window_type = "days", rarefy = FALSE)
algae_turn <- algae_turn[!is.na(algae_turn$dist), ]


algae_rel_turn <- webster(algae_rel_matrix[, 2:ncol(algae_rel_matrix)], 
                          dates = algae_rel_matrix$date, 
                          window = window, dist = "bray",
                          window_type = "days", rarefy = FALSE)
algae_rel_turn <- algae_rel_turn[!is.na(algae_rel_turn$dist), ]
# Set day zero to the same date in all series -----------------------------

allo_mlg_turn <- within(allo_mlg_turn, days <- difftime(date, min_date, units="days"))
allo_spp_turn <- within(allo_spp_turn, days <- difftime(date, min_date, units="days"))
msat_mlg_turn <- within(msat_mlg_turn, days <- difftime(date, min_date, units="days"))
msat_spp_turn <- within(msat_spp_turn, days <- difftime(date, min_date, units="days"))

d.r <- within(d.r, days <- difftime(date, min_date, units="days"))
d.b <- within(d.b, days <- difftime(date, min_date, units="days"))
d.d <- within(d.d, days <- difftime(date, min_date, units="days"))
d.n <- within(d.n, days <- difftime(date, min_date, units="days"))

algae_cyan <- within(algae_cyan, days <- difftime(date, min_date, units="days"))
algae_turn <- within(algae_turn, days <- difftime(date, min_date, units="days"))
algae_rel_turn <- within(algae_rel_turn, days <- difftime(date, min_date, units="days"))
algae_sum <- within(algae_sum, days <- difftime(date, min_date, units="days"))

caullerya <- within(caullerya, days <- difftime(date, min_date, units="days"))

temperature <- within(temperature, days <- difftime(date, min_date, units="days"))

oxygen <- within(oxygen, days <- difftime(date, min_date, units="days"))

lower_gap <- as.Date("2005-06-10")
upper_gap <- as.Date("2006-06-15")
lower_missing <- difftime(lower_gap, min_date, "days")
upper_missing <- difftime(upper_gap, min_date, "days")
# Smooth DAPHNIA TURNOVER series -----------------------------------------------------------

y_allo_m <- with(allo_mlg_turn, 
                 (loess.as(days, dist, plot = T, user.span = 0.03)))

y_allo_s <- with(allo_spp_turn, 
                 (loess.as(days, dist, plot = T, user.span = 0.03)))

y_msat_m <- with(msat_mlg_turn, 
                 (loess.as(days, dist, plot = T, user.span = 0.1)))

y_msat_s <- with(msat_spp_turn, 
                 (loess.as(days, dist, plot = T, user.span = 0.1)))


#  with(y_msat_m, plot(x, y, main = "Allozyme MLG Turnover", 
#                     las = 2, xlab = "days", type="p", cex = 0.1))
# with(y_msat_m, lines(x, fitted, lty="dashed"))
# Smooth DRIVER series ----------------------------------------------------

x_d.r <- with(d.r, 
            (loess.as(days, value, plot = T, user.span = 0.02)))
x_d.b <- with(d.b, 
            (loess.as(days, value, plot = T, user.span = 0.02)))
x_d.d <- with(d.d, 
            (loess.as(days, value, plot = T, user.span = 0.02)))
x_d.n <- with(d.n, 
            (loess.as(days, value, plot = T, user.span = 0.02)))

x_algae_cyan <- with(algae_cyan, 
                   (loess.as(days, Cyanobacteria, plot = T, user.span = 0.02)))
x_algae_turn <- with(algae_turn, 
                   (loess.as(days, dist, plot = T, user.span = 0.02)))
x_algae_rel_turn <- with(algae_rel_turn, 
                       (loess.as(days, dist, plot = T, user.span = 0.02)))
x_algae_sum <- with(algae_sum, 
                  (loess.as(days, V1, plot = T, user.span = 0.03)))

x_caullerya <- with(caullerya, 
                  (loess.as(days, caullerya.prevalence, plot = T, 
                            user.span = 0.04)))

x_temperatre <- with(temperature, 
                   (loess.as(days, value, plot = T, user.span = 0.04)))

x_oxygen <- with(oxygen, 
               (loess.as(days, value, plot = T, user.span = 0.04)))
# Driver graphs -----------------------------------------------------------

msat_min <- as.Date("2008-01-01")
msat_max <- max(msat_mlg_turn$date)

turngg1 <- ggplot(msat_mlg_turn, aes(x = date, y = dist)) +
  geom_line() +
  labs(title = expression(paste(italic(Daphnia), " microsatellite MLG turnover")), x = "Date", y = "Bray-curtis dissimilarity") +
  theme_bw() +
  scale_x_date(breaks = "3 month", labels = date_format("%Y - %b")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# turngg1
# ggsave("Results/3 Msat MLGs/msat_mlg_turn.pdf", turngg1)
# ggsave("Results/3 Msat MLGs/msat_mlg_turn.png", turngg1)


turngg2 <- ggplot(msat_spp_turn, aes(x = date, y = dist)) +
  geom_line()
# turngg2

turngg3 <- ggplot(allo_mlg_turn, aes(x = date, y = dist)) +
  geom_line() +
  labs(title = expression(paste(italic(Daphnia), " allozyme MLG turnover")), x = "Date", y = "Bray-curtis dissimilarity") +
  theme_bw() +
  stat_smooth(method = "lm") + 
  scale_x_date(breaks = "1 year", labels = date_format("%Y - %b")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# turngg3
# ggsave("Results/1 Allo MLGs/allo_mlg_turn.pdf", turngg3)
# ggsave("Results/1 Allo MLGs/allo_mlg_turn.png", turngg3)


turngg4 <- ggplot(allo_spp_turn, aes(x = date, y = dist)) +
  geom_line()
# turngg4


algaegg <- ggplot(algae, aes(x = date, y = biomass, fill = taxon2)) +
  geom_ribbon(position = "stack")
# algaegg
# ggsave("Results/5 drivers/Algae_total.pdf", algaegg,  width = 11, height = 8, units = "in")

turngg <- ggplot(algae_turn, aes(x = date, y = dist)) +
  geom_line()
# turngg
# ggsave("Results/5 drivers/Algae_turnover.pdf", turngg,  width = 11, height = 8, units = "in")

algaegg2 <- ggplot(algae_rel[algae_rel$biomass != 0,], 
                   aes(x = date, y = biomass, fill = taxon2)) +
  geom_ribbon(position = "stack")
# algaegg2
# ggsave("Results/5 drivers/Algae_relative.pdf", algaegg2,  width = 11, height = 8, units = "in")

turngg2 <- ggplot(algae_rel_turn, aes(x = date, y = dist)) +
  geom_line()
# turngg2
# ggsave("Results/5 drivers/Algae_relative_turnover.pdf", turngg,  width = 11, height = 8, units = "in")

cyangg <- ggplot(algae_cyan[algae_cyan$date > msat_min & algae_cyan$date < msat_max,], aes(x = date, y = Cyanobacteria)) +
  geom_line() +
  labs(title = "Cyanobacteria biomass", x = "Date", 
       y = expression(paste("Cyanobacteria biomass (g "^ -3, ")"))) +
  theme_bw() +
  scale_x_date(breaks = "3 month", labels = date_format("%Y - %b")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# cyangg
# ggsave("Results/5 drivers/cyano_biomass.pdf", cyangg)
# ggsave("Results/5 drivers/cyano_biomass.png", cyangg)


caulgg <- ggplot(caullerya[caullerya$date > msat_min & caullerya$date < msat_max, ] , aes(x = date, y = value)) +
  geom_line() +
  labs(title = expression(paste(italic(Caullerya), " prevalence in Greifensee")), 
       x = "Date", y = expression(paste("Fraction of ", italic(Daphnia), " infected"))) +
  theme_bw() +
  scale_x_date(breaks = "3 month", labels = date_format("%Y - %b")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# caulgg
# ggsave("Results/5 drivers/caullerya.pdf", caulgg)
# ggsave("Results/5 drivers/caullerya.png", caulgg)


demo_dgg <- ggplot(d.d, aes(x = date, y = value)) +
  geom_line()
# demo_dgg
# ggsave("Results/5 drivers/demo_death.pdf", demo_dgg,  width = 11, height = 8, units = "in")

demo_bgg <- ggplot(d.b, aes(x = date, y = value)) +
  geom_line()
# demo_bgg
# ggsave("Results/5 drivers/demo_birth.pdf", demo_bgg,  width = 11, height = 8, units = "in")

demo_rgg <- ggplot(d.r, aes(x = date, y = value)) +
  geom_line()
# demo_rgg
# ggsave("Results/5 drivers/demo_r.pdf", demo_rgg,  width = 11, height = 8, units = "in")

demo_ngg <- ggplot(d.n, aes(x = date, y = value)) +
  geom_line()
# demo_ngg
# ggsave("Results/5 drivers/demo_n.pdf", demo_ngg,  width = 11, height = 8, units = "in")

oxygg <- ggplot(oxygen, aes(x = date, y = value)) +
  geom_line()
# oxygg
# ggsave("Results/5 drivers/oxygen.pdf", oxygg,  width = 11, height = 8, units = "in")

tempgg <- ggplot(temperature[temperature$date > msat_min & temperature$date < msat_max,], aes(x = date, y = value)) +
  geom_line() +
  labs(title = "Temperature of Greifensee", x = "Date", 
       y = ("Integrated temperature (?C)")) +
  theme_bw() +
  scale_x_date(breaks = "3 month", labels = date_format("%Y - %b")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# tempgg
# ggsave("Results/5 drivers/temperature.pdf", tempgg)
# ggsave("Results/5 drivers/temperature.png", tempgg)
# One big figure with turnovers and drivers.  -----------------------------

allo_min <- as.Date("1998-01-01")
msat_max <- max(msat_mlg_turn$date)
rect <- data.frame(xmin = as.Date("2005-06-10"), xmax = as.Date("2006-06-15"), 
                   ymin = 0, ymax = 0.5, type = "Caullerya")


msat_to_bind <- data.frame(date = msat_mlg_turn$date, value = msat_mlg_turn$dist, type = "Microsat. Turn.", lab = "(a)", lab_y = 0.9, lab_x = allo_min)
allo_to_bind <- data.frame(date = allo_mlg_turn$date, value = allo_mlg_turn$dist, type = "Allozyme Turn.", lab = "(b)", lab_y = 0.6, lab_x = allo_min)
caul_to_bind <- data.frame(date = caullerya$date, value = caullerya$value, type = "Caullerya", lab = "(c)", lab_y = 0.4, lab_x = allo_min)
cyano_to_bind <- data.frame(date = algae_cyan$date, value = algae_cyan$Cyanobacteria, type = "Cyanobacteria", lab = "(d)", lab_y = 75, lab_x = allo_min)
temp_to_bind <- data.frame(date = temperature$date, value = temperature$value, type = "Temperature", lab = "(e)", lab_y = 15, lab_x = allo_min)

all_drivers <- rbind(msat_to_bind, allo_to_bind, caul_to_bind, cyano_to_bind, temp_to_bind)
all_drivers <- all_drivers[all_drivers$date < msat_max & all_drivers$date > allo_min, ]

all_gg <- ggplot(all_drivers, aes(x = date, y = value)) +
  geom_line() +
  facet_grid(type ~ ., scales = "free_y") +
  theme_bw() +
  scale_x_date(breaks = "1 year", limits = c(allo_min, msat_max), labels = date_format("%Y - %b")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
 # geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = F) +
  geom_smooth(data = subset(all_drivers, all_drivers$type == "Allozyme Turn."), method = "lm") +
  labs(x = "Date", y = "") +
  geom_text(aes(label = lab, y = lab_y, x = lab_x))
all_gg
ggsave("Results/6 Big graph/all_drivers.pdf", width = 11, height = 8.5, units = "in")
ggsave("Results/6 Big graph/all_drivers.png", width = 11, height = 8.5, units = "in")
# msat MLG heatmaps ------------------------------------------------------------
source("heat.R")

msat_m1 <- heat(y_msat_m, x_algae_cyan, sample = c(40, 60), lags = 5, 
                s_window = c(40, 60), l_window = c(-60, 60))
msat_m1
ggsave("Results/3 Msat MLGs/cyanobacteria.pdf", msat_m1)
heat_ci(msat_m1, 1)
heat_ci(msat_m1, 0)
heat_ci(msat_m1, -1)

msat_m2 <- heat(y_msat_m, x_algae_rel_turn, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m2
heat_ci(msat_m2, 1)
heat_ci(msat_m2, 0)
heat_ci(msat_m2, -1)

msat_m3 <- heat(y_msat_m, x_algae_sum, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m3
heat_ci(msat_m3, 1)
heat_ci(msat_m3, 0)
heat_ci(msat_m3, -1)


msat_m4 <- heat(y_msat_m, x_algae_turn, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m4
heat_ci(msat_m4, 1)
heat_ci(msat_m4, 0)
heat_ci(msat_m4, -2)


# caul
msat_m5 <- heat(y_msat_m, x_caullerya, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m5
heat_ci(msat_m5, 1)
heat_ci(msat_m5, 0)
heat_ci(msat_m5, -1)

msat_m6 <- heat(y_msat_m, x_d.b, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m6
heat_ci(msat_m6, 1)
heat_ci(msat_m6, 0)
heat_ci(msat_m6, -1)


msat_m7 <- heat(y_msat_m, x_d.d, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m7
heat_ci(msat_m7, 1)
heat_ci(msat_m7, 0)
heat_ci(msat_m7, -1)

msat_m8 <- heat(y_msat_m, x_d.n, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m8
heat_ci(msat_m8, 1)
heat_ci(msat_m8, 0)
heat_ci(msat_m8, -1)

# r
msat_m9 <- heat(y_msat_m, x_d.r, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m9
heat_ci(msat_m9, 1)
heat_ci(msat_m9, 0)
heat_ci(msat_m9, -1)

msat_m10 <- heat(y_msat_m, x_oxygen, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m10
heat_ci(msat_m10, 1)
heat_ci(msat_m10, 0)
heat_ci(msat_m10, -1)


msat_m11 <- heat(y_msat_m, x_temperatre, sample = c(40, 60), lags = 5, s_window = c(45, 55), l_window = c(-30, 30))
msat_m11
heat_ci(msat_m11, 1)
heat_ci(msat_m11, 0)
heat_ci(msat_m11, -1)
# allo MLG heatmaps ------------------------------------------------------------

allo_m1 <- heat(y_allo_m, x_algae_cyan, sample = c(15, 30), lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m1
heat_ci(allo_m1, 1)
heat_ci(allo_m1, 0)
heat_ci(allo_m1, -1)

allo_m2 <- heat(y_allo_m, x_algae_rel_turn, sample = c(15, 30), lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m2
heat_ci(allo_m2, 1)
heat_ci(allo_m2, 0)
heat_ci(allo_m2, -1)

# algae biomass
allo_m3 <- heat(y_allo_m, x_algae_sum, sample = c(15, 30), lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m3
heat_ci(allo_m3, 1)
heat_ci(allo_m3, 0)
heat_ci(allo_m3, -1)


allo_m4 <- heat(y_allo_m, x_algae_turn, sample = c(15, 30), blank_missing_days = T, lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m4
heat_ci(allo_m4, 1)
heat_ci(allo_m4, 0)
heat_ci(allo_m4, -1)

allo_m5 <- heat(y_allo_m, x_caullerya, sample = c(15, 30), blank_missing_days = T, lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m5
heat_ci(allo_m5, 1)
heat_ci(allo_m5, 0)
heat_ci(allo_m5, -1)

allo_m6 <- heat(y_allo_m, x_d.b, sample = c(15, 30), lags = 5, blank_missing_days = T, s_window = c(20, 25), l_window = c(-30, 30))
allo_m6
heat_ci(allo_m6, 1)
heat_ci(allo_m6, 0)
heat_ci(allo_m6, -1)

allo_m7 <- heat(y_allo_m, x_d.d, sample = c(15, 30), lags = 5, blank_missing_days = T, s_window = c(20, 25), l_window = c(-30, 30))
allo_m7
heat_ci(allo_m7, 1)
heat_ci(allo_m7, 0)
heat_ci(allo_m7, -1)

allo_m8 <- heat(y_allo_m, x_d.n, sample = c(15, 30), lags = 5, blank_missing_days = T, s_window = c(20, 25), l_window = c(-30, 30))
allo_m8
heat_ci(allo_m8, 1)
heat_ci(allo_m8, 0)
heat_ci(allo_m8, -1)

allo_m9 <- heat(y_allo_m, x_d.r, sample = c(15, 30), lags = 5, blank_missing_days = T, s_window = c(20, 25), l_window = c(-30, 30))
allo_m9
heat_ci(allo_m9, 1)
heat_ci(allo_m9, 0)
heat_ci(allo_m9, -1)

allo_m10 <- heat(y_allo_m, x_oxygen, sample = c(15, 30), blank_missing_days = T, lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m10
heat_ci(allo_m10, 1)
heat_ci(allo_m10, 0)
heat_ci(allo_m10, -1)

allo_m11 <- heat(y_allo_m, x_temperatre, sample = c(15, 30), blank_missing_days = T, lags = 5, s_window = c(20, 25), l_window = c(-30, 30))
allo_m11
heat_ci(allo_m11, 1)
heat_ci(allo_m11, 0)
heat_ci(allo_m11, -1)
# msat SPP heatmaps -------------------------------------------------------

msat_s1 <- heat(y_msat_s, x_algae_cyan)
msat_s1

msat_s2 <- heat(y_msat_s, x_algae_rel_turn)
msat_s2

msat_s3 <- heat(y_msat_s, x_algae_sum)
msat_s3

msat_s4 <- heat(y_msat_s, x_algae_turn)
msat_s4

msat_s5 <- heat(y_msat_s, x_caullerya)
msat_s5

msat_s6 <- heat(y_msat_s, x_d.b)
msat_s6

msat_s7 <- heat(y_msat_s, x_d.d)
msat_s7

# pop size
msat_s8 <- heat(y_msat_s, x_d.n)
msat_s8
ggsave("Results/4 Msat Spp/demo_n.pdf", msat_s8, width = 8, height = 9, units = "in")


msat_s9 <- heat(y_msat_s, x_d.r)
msat_s9

msat_s10 <- heat(y_msat_s, x_oxygen)
msat_s10

msat_s11 <- heat(y_msat_s, x_temperatre)
msat_s11
# allo SPP heatmaps -------------------------------------------------------

allo_s1 <- heat(y_allo_s, x_algae_cyan)
allo_s1

allo_s2 <- heat(y_allo_s, x_algae_rel_turn)
allo_s2

allo_s3 <- heat(y_allo_s, x_algae_sum)
allo_s3

allo_s4 <- heat(y_allo_s, x_algae_turn)
allo_s4

allo_s5 <- heat(y_allo_s, x_caullerya)
allo_s5

allo_s6 <- heat(y_allo_s, x_d.b)
allo_s6

allo_s7 <- heat(y_allo_s, x_d.d)
allo_s7

allo_s8 <- heat(y_allo_s, x_d.n)
allo_s8

allo_s9 <- heat(y_allo_s, x_d.r)
allo_s9

allo_s10 <- heat(y_allo_s, x_oxygen)
allo_s10

allo_s11 <- heat(y_allo_s, x_temperatre)
allo_s11
# DAPHNIA strip graphs ----------------------------------------------------
important_mlgs <- colSums(msat_mlg[, 3:ncol(msat_mlg)])
singles <- names(msat_mlg[3:length(names(msat_mlg))])[important_mlgs == 1]


msat_mlg_rel <- msat_mlg[, 3:ncol(msat_mlg)] / rowSums(msat_mlg[, 3:ncol(msat_mlg)])

max_prev <- sapply(msat_mlg_rel, max)
rares <- names(msat_mlg[3:length(names(msat_mlg))])[max_prev < 0.1]

msat_mlg_rel$date <- msat_mlg$date




msat_mlg_rel <- melt(msat_mlg_rel, id.vars = "date")
msat_mlg_rel$variable <- as.character(msat_mlg_rel$variable)

msat_mlg_rel$variable[msat_mlg_rel$variable %in% singles] <- "singles"
msat_mlg_rel$variable[msat_mlg_rel$variable %in% rares] <- "rare"


msat_mlg_rel <- ddply(msat_mlg_rel, .(date, variable), function(x) sum(x$value))
msat_mlg_rel$year <- as.numeric(format(msat_mlg_rel$date, "%Y"))

msat_mlg_rel_gg <- ggplot(msat_mlg_rel[msat_mlg_rel$year == 2010, ], aes(x = date, y = V1, fill = variable)) +
  geom_ribbon(position = "stack")
msat_mlg_rel_gg



important_a_mlgs <- colSums(allo_mlg[, 2:ncol(allo_mlg)])
singles_a <- names(allo_mlg[2:length(names(allo_mlg))])[important_a_mlgs == 1]


allo_mlg_rel <- allo_mlg[, 2:ncol(allo_mlg)] / rowSums(allo_mlg[, 2:ncol(allo_mlg)])

max_prev_a <- sapply(allo_mlg_rel, max)
rares_a <- names(allo_mlg[2:length(names(allo_mlg))])[max_prev_a < 0.1]

allo_mlg_rel$date <- allo_mlg$date
allo_mlg_rel <- melt(allo_mlg_rel, id.vars = "date")
allo_mlg_rel$variable <- as.character(allo_mlg_rel$variable)

allo_mlg_rel$variable[allo_mlg_rel$variable %in% singles_a] <- "singles"
allo_mlg_rel$variable[allo_mlg_rel$variable %in% rares_a] <- "rare"


allo_mlg_rel <- ddply(allo_mlg_rel, .(date, variable), function(x) sum(x$value))
allo_mlg_rel$year <- as.numeric(format(allo_mlg_rel$date, "%Y"))

allo_mlg_rel_gg <- ggplot(allo_mlg_rel, aes(x = date, y = V1, fill = variable)) +
  geom_ribbon(aes(ymax = V1), ymin = 0, alpha = 0.2)
allo_mlg_rel_gg

allo_mlg_rel_gg <- ggplot(allo_mlg_rel[allo_mlg_rel$year == 2001, ], aes(x = date, y = V1, fill = variable)) +
  geom_ribbon(position = "stack")
allo_mlg_rel_gg
# DAPHNIA diversity -------------------------------------------------------

div <- data.frame(date = msat_mlg$date, simpson = 0)
msat_mlg_rel2 <- msat_mlg[, 3:ncol(msat_mlg)] / rowSums(msat_mlg[, 3:ncol(msat_mlg)])

div$simpson <- diversity(msat_mlg_rel2, index = "simpson")

div$year <- as.factor(format(div$date, "%Y"))

div_plot <- ggplot(div, aes(x = date, y = simpson)) +
  geom_line() 
div_plot

div2 <- data.frame(date = allo_mlg$date, simpson = 0)
allo_mlg_rel2 <- allo_mlg[, 2:ncol(allo_mlg)] / rowSums(allo_mlg[, 2:ncol(allo_mlg)])

div2$simpson <- diversity(allo_mlg_rel2, index = "simpson")

div2$year <- as.factor(format(div2$date, "%Y"))

div_plot2 <- ggplot(div2, aes(x = date, y = simpson)) +
  geom_line() 
div_plot2
# contribution of MSAT MLGs to total -------------------------------------------

n <- data.frame(date = d.n$date, n = d.n$value)
n <- n[n$date %in% msat_mlg$date, ]

msat_mlg_rel3 <- msat_mlg[, 3:ncol(msat_mlg)] / rowSums(msat_mlg[, 3:ncol(msat_mlg)])
msat_mlg_rel3 <- msat_mlg_rel3[msat_mlg$date %in% n$date, ]
msat_mlg_rel4 <- msat_mlg_rel3 * n$n

important_mlgs <- colSums(msat_mlg[, 3:ncol(msat_mlg)])
#singles <- names(msat_mlg[3:length(names(msat_mlg))])[important_mlgs == 1]

max_prev <- sapply(msat_mlg_rel4, max)
rares <- names(msat_mlg_rel4)[max_prev < 0.1]

just_imp <- msat_mlg_rel4[!names(msat_mlg_rel4) %in% rares]
just_imp$date <- msat_mlg$date[msat_mlg$date %in% n$date]
just_imp$n <- n$n

just_imp <- melt(just_imp, id.vars = "date")
just_imp$type <- "clones"
just_imp$type[just_imp$variable == "n"] <- "total"

imp_gg <- ggplot(just_imp, aes(x = date, y = value, color = variable)) +
  geom_line() +
  facet_grid(type~., scales = "free_y")
imp_gg
# Contribution of ALLO MLGs to total --------------------------------------


n <- data.frame(date = d.n$date, n = d.n$value)
n <- n[n$date %in% allo_mlg$date, ]

allo_mlg_rel3 <- allo_mlg[, 2:ncol(allo_mlg)] / rowSums(allo_mlg[, 2:ncol(allo_mlg)])
allo_mlg_rel3 <- allo_mlg_rel3[allo_mlg$date %in% n$date, ]
allo_mlg_rel4 <- allo_mlg_rel3 * n$n

important_mlgs <- colSums(allo_mlg[, 2:ncol(allo_mlg)])
#singles <- names(allo_mlg[3:length(names(allo_mlg))])[important_mlgs == 1]

max_prev <- sapply(allo_mlg_rel4, max)
rares <- names(allo_mlg_rel4)[max_prev < 0.1]

just_imp <- allo_mlg_rel4[!names(allo_mlg_rel4) %in% rares]
just_imp$date <- allo_mlg$date[allo_mlg$date %in% n$date]
just_imp$n <- n$n

just_imp <- melt(just_imp, id.vars = "date")
just_imp$type <- "clones"
just_imp$type[just_imp$variable == "n"] <- "total"

imp_gg2 <- ggplot(just_imp, aes(x = date, y = value, color = variable)) +
  geom_line() +
  facet_grid(type~., scales = "free_y")
imp_gg2


deathgg <- ggplot(d.d, aes(x = date, y = value)) +
  geom_line() +
  scale_x_date(breaks = "year")
deathgg
# Correlation plots for msat MLG turnovers --------------------------------

source("heat.R")
msat_m1.1 <- heat(y_msat_m, x_algae_cyan, sample = c(40, 60), lags = 5, 
                s_window = c(40, 60), l_window = c(-60, 60), give_series = T)

msat_m2.1 <- heat(y_msat_m, x_algae_rel_turn, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m3.1 <- heat(y_msat_m, x_algae_sum, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m4.1 <- heat(y_msat_m, x_algae_turn, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)
# caul
msat_m5.1 <- heat(y_msat_m, x_caullerya, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), 
                give_series = T)

msat_m6.1 <- heat(y_msat_m, x_d.b, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m7.1 <- heat(y_msat_m, x_d.d, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m8.1 <- heat(y_msat_m, x_d.n, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m9.1 <- heat(y_msat_m, x_d.r, sample = c(40, 60), lags = 5, 
                s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m10.1 <- heat(y_msat_m, x_oxygen, sample = c(40, 60), lags = 5, 
                 s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m11.1 <- heat(y_msat_m, x_temperatre, sample = c(40, 60), lags = 5, 
                 s_window = c(45, 55), l_window = c(-30, 30), give_series = T)

msat_m1.1$index <- paste(msat_m1.1$interval, msat_m1.1$start_push, 
                         sep = "-")

selection <- sample(unique(msat_m1.1$index), size = 56, replace = F)
to_plot <- msat_m1.1[msat_m1.1$index %in% selection, ]

corrs <- ddply(to_plot, .(index), function(x){
  corr <- x[1, "lag_0"]
  names(corr) <- "corr"
  corr})

corrs$corr <- round(corrs$corr, 3)

ggmsat_m1.1 <- ggplot(to_plot, aes(x = series2, y = series1)) +
  geom_point() +
  scale_color_discrete() +
  facet_wrap(~ index, ncol = 7) +
  theme_bw() +
  stat_smooth(method = "lm") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  geom_text(data = corrs, aes(label = corr), x = 12, y = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 # scale_x_continuous(breaks = seq(0, 0.15, 0.05))+
  labs(title = "Cyanobacteria Biomass", 
       y = expression(paste(italic(Daphnia), " microsat. Turnover")), 
       x = expression(paste("Cyanobacteria biomass (g m"^ -3, ")")))
ggmsat_m1.1


msat_m5.1$index <- paste(msat_m5.1$interval, msat_m5.1$start_push, 
                         sep = "-")

selection <- sample(unique(msat_m5.1$index), size = 56, replace = F)
to_plot <- msat_m5.1[msat_m5.1$index %in% selection, ]

corrs <- ddply(to_plot, .(index), function(x){
  corr <- x[1, "lag_0"]
  names(corr) <- "corr"
  corr})

corrs$corr <- round(corrs$corr, 3)

ggmsat_m5.1 <- ggplot(to_plot, aes(x = series2, y = series1)) +
  geom_point() +
  scale_color_discrete() +
  facet_wrap(~ index, ncol = 7) +
  theme_bw() +
  stat_smooth(method = "lm") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  geom_text(data = corrs, aes(label = corr), x = 0.15, y = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(0, 0.15, 0.05))+
  labs(title = "Parasite Prevalence", 
       y = expression(paste(italic(Daphnia), " microsat. Turnover")), 
       x = expression(paste("Fraction of ", italic(Daphnia), " infected by ", 
                            italic(Caullerya))))
ggmsat_m5.1




msat_m11.1$index <- paste(msat_m11.1$interval, msat_m11.1$start_push, 
                         sep = "-")

selection <- sample(unique(msat_m11.1$index), size = 56, replace = F)
to_plot <- msat_m11.1[msat_m11.1$index %in% selection, ]

corrs <- ddply(to_plot, .(index), function(x){
  corr <- x[1, "lag_0"]
  names(corr) <- "corr"
  corr})

corrs$corr <- round(corrs$corr, 3)

ggmsat_m11.1 <- ggplot(to_plot, aes(x = series2, y = series1)) +
  geom_point() +
  scale_color_discrete() +
  facet_wrap(~ index, ncol = 7) +
  theme_bw() +
  stat_smooth(method = "lm") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  geom_text(data = corrs, aes(label = corr), x = 12, y = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_x_continuous(breaks = seq(0, 0.15, 0.05))+
  labs(title = "Temperature", 
       y = expression(paste(italic(Daphnia), " microsat. Turnover")), 
       x = "Intgrated temperature (?C)")
ggmsat_m11.1
  
ggsave("Results/Msat MLG Turn Correlations/Cyanos.pdf", ggmsat_m1.1, 
       width = 8, height = 9, units = "in")

ggsave("Results/Msat MLG Turn Correlations/Cyanos.png", ggmsat_m1.1, 
       width = 8, height = 9, units = "in")

ggsave("Results/Msat MLG Turn Correlations/Caul.pdf", ggmsat_m5.1, 
       width = 8, height = 9, units = "in")

ggsave("Results/Msat MLG Turn Correlations/Caul.png", ggmsat_m5.1, 
       width = 8, height = 9, units = "in")

ggsave("Results/Msat MLG Turn Correlations/Temp.pdf", ggmsat_m11.1, 
       width = 8, height = 9, units = "in")

ggsave("Results/Msat MLG Turn Correlations/Temp.png", ggmsat_m11.1, 
       width = 8, height = 9, units = "in")