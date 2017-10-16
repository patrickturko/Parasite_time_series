
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

#Do a Fisher's exact test to check whether infected and random samples differ at each date
dates <- unique(counts$date)
pop.fish <- data.frame(date = dates, p_value = 0)

for(i in 1:length(dates)){
  tmp <- counts[counts$date == dates[i],] 
  tmp <- dcast(tmp, mlg ~ type)
  tmp <- tmp[tmp$Infected != 0 | tmp$Random !=0, ]
  fish <- fisher.test(as.matrix(tmp[, 2:3]), workspace = 1e+09)
  pop.fish[i, "p_value"] <- fish$p.value
}

####Make a graph showing how sample types differ at each date
# First, rename the mlgs so that each date has it's own set

mlgs_by_date <- ddply(counts, .(date), function(x){
  tmp <- dcast(x, mlg ~ type)
  tmp <- tmp[tmp$Infected != 0 | tmp$Random !=0, ]
  tmp$mlg2 <- paste(x$date[1], tmp$mlg, sep = "_")
  tmp <- tmp[order(tmp$Random), ]
  tmp$mlg2 <- factor(tmp$mlg2, levels = tmp$mlg2)
  tmp
})

mlgs_by_date <- melt(mlgs_by_date, id.vars = c("date", "mlg", "mlg2"))
mlgs_by_date$type <- factor(mlgs_by_date$type, levels = c("Random", "Infected"))
names(mlgs_by_date)[4:5] <- c("type", "count")

#Second, create the graph 
pop.graph <- ggplot(mlgs_by_date, 
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

ggsave("test.pdf", pop.graph, height = 10, width = 4, units = "in", dpi = 600)

       