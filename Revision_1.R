setwd("/media/FILES/Dropbox/Turko PhD Data/Turko et al Chapter 5 Greifensee/Time Series")
setwd("E:/Dropbox/Turko PhD Data/Turko et al Chapter 5 Greifensee/Time Series")

library(plyr)
library(ggplot2)
library(ggExtra)
library(reshape2)

# Load, clean, and prepare data -------------------------------------------

# Load and clean data
mlgs <- read.table(file = "inputs/confirmed_mlgs.txt", header = TRUE, sep = "\t")
mlgs <- mlgs[, c("mlg", "pop", "X2")]
names(mlgs) <- c("mlg", "date", "type")
mlgs$date <- as.Date(mlgs$date)
mlgs$mlg <- as.factor(mlgs$mlg)
mlgs$date[mlgs$date == "2012-09-20" & mlgs$type == "Random"] <- "2012-09-06"


# Extract only those dates which have both infected and random sample
matched <- ddply(mlgs, .variables = "date", .fun = function(x){
  if(length(unique(x$type)) == 2)
    x
})

# Count how many individuals of each MLG are present in each sample at each date
counts <- dcast(matched, date + type ~ mlg, fun.aggregate = length)
counts <- melt(counts, id.vars = c("date", "type"))
names(counts) <- c("date", "type", "mlg", "count")

#Prepare a sample size table for the paper
ss <- ddply(counts, .(date, type), function(x) sum(x$count))
ss <- dcast(ss, date ~ type)  

  
#Express counts as percentages instead.
counts <- ddply(counts, .(date, type), function(x){
  total <- sum(x$count)
  x$percent <- x$count / total
  x
})

# Rename the mlgs so that each date has it's own set
# For the percent data
percent_by_date <- ddply(counts, .(date), function(x){
  tmp <- dcast(x, mlg ~ type, value.var = "percent")
  tmp <- tmp[tmp$Infected != 0 | tmp$Random !=0, ]
  tmp$mlg2 <- paste(x$date[1], tmp$mlg, sep = "_")
  tmp <- tmp[order(tmp$Random), ]
  tmp$mlg2 <- factor(tmp$mlg2, levels = tmp$mlg2)
  tmp
})

#And for the count data
counts_by_date <- ddply(counts, .(date), function(x){
  tmp <- dcast(x, mlg ~ type, value.var = "count")
  tmp <- tmp[tmp$Infected != 0 | tmp$Random !=0, ]
  tmp$mlg2 <- paste(x$date[1], tmp$mlg, sep = "_")
  tmp <- tmp[order(tmp$Random), ]
  tmp$mlg2 <- factor(tmp$mlg2, levels = tmp$mlg2)
  tmp
})

# Refactor variables to get the desired order in the plot
percent_by_date <- melt(percent_by_date, id.vars = c("date", "mlg", "mlg2"))
names(percent_by_date)[4:5] <- c("type", "percent")
percent_by_date$type <- factor(percent_by_date$type, levels = c("Random", "Infected"))
percent_by_date[percent_by_date == 0] <- NA

counts_by_date <- melt(counts_by_date, id.vars = c("date", "mlg", "mlg2"))
names(counts_by_date)[4:5] <- c("type", "count")
counts_by_date$type <- factor(counts_by_date$type, levels = c("Random", "Infected"))
counts_by_date[counts_by_date == 0] <- NA






# Create Graphs -----------------------------------------------------------


#Create the graph (for counts)
pop.graph <- ggplot(counts_by_date, 
                    aes(x = mlg2, y = count)) +
  facet_grid(date~type, scales = "free") +
  geom_bar(stat = "identity", width = 0.8, color = "lightgrey") +
  theme_bw() +
  removeGrid() +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab(label = "Count") +
  xlab(label = "MLG") +
  theme(aspect.ratio = 1) +
  coord_flip()
pop.graph

# ggsave("Rand_vs_Inf.pdf", pop.graph, height = 10, width = 4, units = "in", dpi = 600)

#for percents
pop.graph2 <- ggplot(percent_by_date, 
                    aes(x = mlg2, y = percent * 100)) +
  facet_grid(date~type, scales = "free_y") +
  geom_bar(stat = "identity", width = 0.8, fill = "lightgrey") +
  theme_bw() +
  removeGrid() +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab(label = "Percent of Sample") +
  xlab(label = "MLG") +
  theme(aspect.ratio = 1) +
  coord_flip() 
pop.graph2

# ggsave("Rand_vs_Inf_percent.pdf", pop.graph2, height = 10, width = 4, units = "in", dpi = 600)





# Stats 1: Fisher test on whole pop ------------------------------------------------------------

#Do a Fisher's exact test to check whether infected and random samples differ at each date
dates <- unique(counts$date)
pop.fish <- data.frame(date = dates, p_value = 0)

for(i in 1:length(dates)){
  tmp <- counts[counts$date == dates[i], 1:4] 
  tmp <- dcast(tmp, mlg ~ type)
  tmp <- tmp[tmp$Infected != 0 | tmp$Random !=0, ]
  fish <- fisher.test(as.matrix(tmp[, 2:3]), workspace = 1e+09)
  pop.fish[i, "p_value"] <- fish$p.value
}
# Apply Holm-bonferroni correction

pop.fish <- pop.fish[order(pop.fish$p_value), ]
pop.fish$threshhold <- 0.05 / (nrow(pop.fish) - 1:nrow(pop.fish) + 1)
pop.fish$sig <- "no"
pop.fish$sig[pop.fish$p_value <= pop.fish$threshhold] <- "yes"

#write to file
write.table(pop.fish, "whole_population.csv", row.names = F)





# Stats 2: Fisher test on common clones (top5 %)-----------------------------------

# Take only those over 5%, and reshape to suit test
common <- ddply(counts, .(date, mlg), function(x){
  if(max(x$percent) >= 0.05) 
    x$common <- "yes"
  else
    x$common <- "no"
  x
})

#Make holder & do the test
dates2 <- unique(common[common$common == "yes", c("date", "mlg")])
dates2$p_value = 1

for(i in 1:nrow(dates2)){
  date <- dates2$date[i]
  mlg <- dates2$mlg[i]
  com <- common[common$date == date & common$mlg == mlg, ]
  com_rand <- com$count[com$type == "Random"]
  com_inf <- com$count[com$type == "Infected"]
  other <- common[common$date == date & common$mlg != mlg, ]
  other_rand <- sum(other$count[other$type == "Random"])
  other_inf <- sum(other$count[other$type == "Infected"])
  to_test <- data.frame(random = c(com_rand, other_rand),
                        infected = c(com_inf, other_inf))
  tmp_p <- fisher.test(to_test)$p.value
  dates2[i, "p_value"] <- tmp_p
  }

# Apply Holm-bonferroni correction

dates2 <- dates2[order(dates2$p_value), ]
dates2$threshhold <- 0.05 / (nrow(dates2) - 1:nrow(dates2) + 1)
dates2$sig <- "no"
dates2$sig[dates2$p_value <= dates2$threshhold] <- "yes"

#write to file
write.table(dates2, "over_under_inf.csv", row.names = F)




# Stats 2.5 Fisher test on two most common clones -----------------------------------
#Rank all the mlgs in each sample type . 
counts2 <- ddply(counts, .(date, type), function(x){
  x <- x[order(x$count, decreasing = TRUE), ]
  x$rank <- 0
  for(i in 1:nrow(x)){
    if(x$count[i] != 0 & i == 1)
      x$rank[i] <- 1
    if(i > 1){
      if(x$count[i] != 0 & x$count[i] == x$count[i - 1])
        x$rank[i] <- x$rank[i - 1] 
      if(x$count[i] != 0 & x$count[i] < x$count[i - 1])
        x$rank[i] <- x$rank[i - 1] + 1
    }
  }
  x
})


# Take only the top clones in random (which are also over 5%), and reshape to suit test
common2 <- ddply(counts2, .(date, mlg), function(x){
  if(x$rank[2] %in% c(1,2) & x$percent[2] > 0.052)
    x$common <- "yes"
  else
    x$common <- "no"
  x
})




#Make holder & do the test
dates3 <- unique(common2[common2$common == "yes", c("date", "mlg")])
dates3$p_value = 1

for(i in 1:nrow(dates3)){
  date <- dates3$date[i]
  mlg <- dates3$mlg[i]
  com <- common2[common2$date == date & common2$mlg == mlg, ]
  com_rand <- com$count[com$type == "Random"]
  com_inf <- com$count[com$type == "Infected"]
  other <- common2[common2$date == date & common2$mlg != mlg, ]
  other_rand <- sum(other$count[other$type == "Random"])
  other_inf <- sum(other$count[other$type == "Infected"])
  to_test <- data.frame(random = c(com_rand, other_rand),
                        infected = c(com_inf, other_inf))
  tmp_p <- fisher.test(to_test)$p.value
  dates3[i, "p_value"] <- tmp_p
}

# Git rid of second-class clone where there are two first-rank clones
dates3 <- dates3[dates3$mlg != 32, ]

# Apply Holm-bonferroni correction

dates3 <- dates3[order(dates3$p_value), ]
dates3$threshhold <- 0.05 / (nrow(dates3) - 1:nrow(dates3) + 1)
dates3$sig <- "no"
dates3$sig[dates3$p_value <= dates3$threshhold] <- "yes"

#write to file
#write.table(dates3, "over_under_inf_top_two.csv", row.names = F)

#Subset common clones to report to Justyna

common3 <- common2[common2$type == "Random" & common2$common == "yes", ]
common3 <- common3[with(common3, order(date, rank)), ]
#write.table(common3, "list_of_common_clones.csv", quote = F, row.names = F)
 
# calculate the probability that this many significant results could 
# be found by chance alone Moran 2003 "Arguments for rejecting the sequential
# Bonferroni in ecological studies"

alpha <- 0.05
N <- nrow(dates3)
K <- nrow(dates3[dates3$p_value < alpha, ])
NmK <- N - K

first_term <- factorial(N) / (factorial(NmK) * factorial(K))
second_term <- alpha^K * (1 - alpha)^NmK
probability <- first_term * second_term



# Redo graph in light of stats --------------------------------------------
# Make list of clones tested for mis-infection
mis_inf_clones <- paste(dates3$date, dates3$mlg, sep = "_")

#Filter percent data to just these clones 
percent_by_date2 <- percent_by_date[percent_by_date$mlg2 %in% mis_inf_clones, ]

#Add info on significance
percent_by_date2$mlg2 <- paste(percent_by_date2$date, percent_by_date2$mlg, sep = "_")
dates3$mlg2 <- mis_inf_clones
percent_by_date2 <- merge(percent_by_date2, dates3)

#Mark over vs under-infection
percent_by_date2 <- ddply(percent_by_date2, .(mlg2), function(x){
  
  x[is.na(x)] <- 0
  if(x[x$type == "Infected", "percent"] < x[x$type == "Random", "percent"])
    x$Infection <- "Under"
  if(x[x$type == "Infected", "percent"] > x[x$type == "Random", "percent"])
    x$Infection <- "Over"
  x
})

#Mark non-signifcant differences 
percent_by_date2[percent_by_date2$p_value > 0.05, "Infection"] <- "Proportionate"

#Add colored bars to original ggplot

pop.graph2.red <- pop.graph2 +
  geom_bar(data = percent_by_date2, aes(x = mlg2, y = percent * 100, 
                                        fill = Infection), 
           stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("red", "black", "blue"))
pop.graph2.red  

#Add p-values to colored plot
#Crate the annotation data frame
ps <- pop.fish
ps$p_value[ps$p_value > 0.00001] <- round(ps$p_value[ps$p_value > 0.00001], 3)
ps$p_value[ps$p_value < 0.00001] <- "< 0.00001"
ps$type <- factor("Infected", levels = c("Random", "Infected"))
ps$face <- "plain"
ps$face[ps$sig == "yes"] <- "italic"

#Calculate the label positions 
for(i in 1:nrow(ps)){
  ps$pos[i] <- nrow(percent_by_date[percent_by_date$date == ps[i, "date"], ]) * 0.4
}

#add the annotations to the graph
final.graph <- pop.graph2.red +
  geom_text(data = ps, aes(label = p_value, x = pos, fontface = face), 
            y = 45, hjust = "inward", vjust = "inward", 
            size = 3) +
  theme(legend.position="none")
final.graph

#save it
ggsave("Rand_vs_Inf_percent_annotated3.pdf", final.graph, height = 10, width = 6, units = "in", dpi = 600)

