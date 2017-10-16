 #data = allo_mlg[, 3:ncol(allo_mlg)]; dates = allo_mlg$date; window = 30
 #window_type = "days"; dist = "bray"

webster <- function(data, dates, window, window_type, dist, rarefy = F){
  library(vegan)
  mids <- dates[2:length(dates)] - dates[1:(length(dates) -1 )]
  mids <- as.Date(dates[1:(length(dates) -1 )] + (mids / 2))
  
  answer <- data.frame(date = mids, dist = NA, left_smps = 0, right_smps = 0, 
                       left_days = 0, right_days = 0, left = as.Date("1900-01-01"), 
                       right = as.Date("1900-01-01"), left_sum = 0, right_sum = 0, 
                       all_sum = 0)
  
  if(window_type == "sample"){
    start <- window + 1
    end <- nrow(answer) - window
  }
  
  if(window_type == "days"){
    start <- which(mids - (dates[1] + window) >= 0)[1]
    ends <- which(mids - (dates[length(dates)] - window) <=0)
    end <- ends[length(ends)]
  }
  
  
  for(i in start:end){
    left <- dates[dates < answer$date[i]]
    
    right <- dates[dates > answer$date[i]]
    
    if(window_type == "sample"){
      left <- left[(length(left) - window + 1): length(left)]
      right <- right[1: window]
      # right <- right[(length(right) - window): length(right)]
    }
    if(window_type == "days"){
      left <- left[left >= (answer$date[i] - window)]
      right <- right[right <= (answer$date[i] + window)]
    }
    left_index <- which(dates %in% left)
    right_index <- which(dates %in% right)
    
    comps <- data.frame(rbind(colSums(data[left_index, ]), colSums(data[right_index, ])))
    left_sum <- sum(comps[1,])
    right_sum <- sum(comps[2,])
    all_sum <- left_sum + right_sum
    
    if(rarefy == "TRUE"){
      min <- min(rowSums(comps))
      comps[1, ] <- rrarefy(comps[1, ], sample = min)
      comps[2, ] <- rrarefy(comps[2, ], sample = min)
    }
     
    answer$dist[i] <- vegdist(comps, method = dist)
    answer$left_smps[i] <- length(left)
    answer$right_smps[i] <- length(right)
    answer$left_days[i] <- as.numeric(max(left) -min(left))
    answer$right_days[i] <- as.numeric(max(right) -min(right))
    answer$left[i] <- left[1]
    answer$right[i] <- right[length(right)]
    answer$window <- window
    answer$left_sum[i] <- left_sum
    answer$right_sum[i] <- right_sum
    answer$all_sum[i] <- all_sum
    
  
  }
  answer
}