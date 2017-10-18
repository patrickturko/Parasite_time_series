
library(plyr)
library(ggplot2)
library(ggExtra)
library(reshape2)

#Load and clean data
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

#Express counts as percentages instead.
counts <- ddply(counts, .(date, type), function(x){
  total <- sum(x$count)
  x$percent <- x$count / total
  x
})

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

####Make a graph showing how sample types differ at each date
# First, rename the mlgs so that each date has it's own set
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

#Second, create the graph (for counts)
pop.graph <- ggplot(counts_by_date, 
                    aes(x = mlg2, y = count)) +
  facet_grid(date~type, scales = "free") +
  geom_bar(stat = "identity", width = 0.8) +
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
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  removeGrid() +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab(label = "Percent of Sample") +
  xlab(label = "MLG") +
  theme(aspect.ratio = 1) +
  coord_flip() +
  geom_hline(yintercept = 5, size = 0.1)
pop.graph2

ggsave("Rand_vs_Inf_percent.pdf", pop.graph2, height = 10, width = 4, units = "in", dpi = 600)

