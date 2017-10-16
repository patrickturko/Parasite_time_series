# 
#  series1 <- y_allo_m; series2 <- x_algae_sum; blank_missing_days <- F ; 
#  sample = c(40, 80) 
#  lags = 5; s_window = c(45, 55); l_window = c(-24, 24)

# sample = c(15, 30); lags = 5; s_window = c(20, 25); l_window = c(-30, 30)

heat <- function(series1, series2, blank_missing_days = F, sample = c(1, 60), 
                 lags = 23, s_window = c(10, 20), l_window = c(-24, 24), 
                 give_series = F, log_s1 = F, log_s2 = F){
  
  lower_gap <- as.Date("2005-06-10")
  upper_gap <- as.Date("2006-06-15")
  lower_missing <- difftime(lower_gap,min_date, "days")
  upper_missing <- difftime(upper_gap,min_date, "days")
  
  interval <- c(s_window[1]:s_window[2], s_window[2], s_window[2]:s_window[1], s_window[1])
  int_length <- length(interval)
  l_days <- c(rep(l_window[2], (int_length / 2 -1)), 
              rep(l_window[1], int_length / 2), l_window[2])
  variable<-l_days/interval
  value<-rep(0,int_length)
  coords<-data.frame(interval,variable,value,l_days)
  
m.x <- matrix(nrow = sample[2], ncol = sample[2])

max_min <- as.numeric(max(min(series1$x), min(series2$x)))
min_max <- as.numeric(min(max(series1$x), max(series2$x)))

all_series <- data.frame(series1 = NA, series2 = NA,  interval = NA, 
                         start_push = NA)
all_series <- all_series[-1, ]

corr_names <- paste("lag", seq(-lags, lags, 1), sep = "_")
corr_data <- as.data.frame(matrix(ncol = lags * 2 + 1))
names(corr_data) <- corr_names
corr_data <- corr_data[-1, ]

all_series <- cbind(all_series, corr_data)

for(u in sample[1]: sample[2]){
  n <- 30
  m.d <- matrix(nrow = n, ncol = lags * 2 + 1)
  
  
  for(k in 1:n){
    m.d[k, 2] <- j <- k
    lag <- u
    fit.dates <- seq(max_min + j, min_max, lag)
    
    d.fit1 <- predict(series1, fit.dates, se = T)
   if(log_s1 == TRUE) d.fit1$fit <- log(d.fit1$fit)
    d.fit2 <- predict(series2, fit.dates, se = T)
   if(log_s2 == TRUE) d.fit2$fit <- log(d.fit2$fit)
   
    if(log_s1 == F)
    d.fit1 <- within(d.fit1, fit[fit < 0] <- 0)
    if(log_s2 == F)
    d.fit2 <- within(d.fit2, fit[fit < 0] <- 0)
    
    if(blank_missing_days == TRUE){
      d.fit1 <- within(d.fit1, fit[(fit.dates > lower_missing & 
                                      fit.dates < upper_missing)] <- NA)
      d.fit2 <- within(d.fit2, fit[(fit.dates > lower_missing & 
                                      fit.dates < upper_missing)] <- NA)
    } 
    
    d.fit1 <- ts(d.fit1$fit)
    d.fit2 <- ts(d.fit2$fit)
    
    d.cc <- ccf(d.fit1, d.fit2, na.action = na.pass, plot = F, lag.max = lags)
    m.d[k, 1:ncol(m.d)] <- d.cc$acf
    
    
    check <- data.frame(series1 = as.numeric(d.fit1), 
                        series2 = as.numeric(d.fit2),  
                        interval = u, start_push = k)
    check2 <- as.data.frame(matrix(ncol = lags * 2 + 1, nrow = length(as.numeric(d.fit1))))
    names(check2) <- corr_names
    check2[1, ] <- d.cc$acf
    check3 <- cbind(check, check2)
    all_series <- rbind(all_series, check3)
    
  }
  m.x[u,1] <- u
  m.x[u,2:(ncol(m.d) + 1)] <- apply(m.d, 2, mean)
}

m.x <- as.data.frame(m.x)
names(m.x) <- c("interval",c(-lags:lags))
# write.table(m.x, "../results/cross-correlations/cross_correlations_temperature.csv", row.names=F, sep=";", quote=F)

d.contour <- melt(as.data.frame(m.x), id=1)
d.contour <- within(d.contour, {variable <- as.numeric(as.character(variable))
                                l_days <- interval*variable})


p8 <- ggplot(d.contour, aes(x = interval, y = variable)) +
  geom_tile(aes(fill = value)) + 
  scale_fill_gradientn(colours = rainbow(8)) +
  stat_contour(aes(z = value)) +
  scale_x_continuous(breaks = seq(0, 60, 2))  + 
  scale_y_continuous(breaks = seq(-22, 22, 2)) +
  labs(x = "sampling interval [d]", y = "lag [d]") +
  geom_path(data = coords, col = "grey", size = 1)

if(give_series == T) return(all_series)

if(give_series == F) return(p8)

}