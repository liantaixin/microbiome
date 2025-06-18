# Clean up all objects in the work environment
rm(list = ls())


library(devtools)
library(ggplot2)
library(ggimage)
library(plyr)

library(MicEco)

# Read data
read_table <- read.csv(file = "D:/phD/16s/16s_files/NCM/gut/high/all.csv", row.names = 1)

res <- neutral.fit(t(read_table))


m <- res[[1]][1]
N <- res[[1]][4]
Nm<- N*m
r2 <- res[[1]][3]
out <- res[[2]]

# Processing data
out$group <- with(out, ifelse(freq < Lower, "#D8BFD8", 
                              ifelse(freq > Upper, "#CD5C5C", "#4682B4")))


# Plot the model results
p1 <- ggplot(data = out) +
  geom_line(aes(x = log(p), y = freq.pred), size = 1.5, linetype = 1) +
  geom_line(aes(x = log(p), y = Lower), size = 1.5, linetype = 8) +
  geom_line(aes(x = log(p), y = Upper), size = 1.5, linetype = 8) +
  geom_point(aes(x = log(p), y = freq, color = group), size = 4) +
  xlab("Log10(mean relative abundance)") +
  ylab("Occurrence frequency") +
  scale_colour_manual(values = c("#4682B4", "#CD5C5C", "#29A7A5")) +  # Color Scheme
  annotate("text", x = -10, y = 0.99, label = paste("Rsqr = ", round(r2, 2)), size = 5.5) +
  annotate("text", x = -10, y = 0.94, label = paste("Nm = ", round(Nm, 2)), size = 5.5) +
  annotate("text", x = -10, y = 0.89, label = paste("N = ", round(N, 2)), size = 5.5) +
  annotate("text", x = -10, y = 0.84, label = paste("m = ", round(m, 3)), size = 5.5) + 
    theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(size = 0.2, colour = "black"),
    axis.line.y = element_line(size = 0.2, colour = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 15),
    text = element_text(family = "sans", size = 15),
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )


p1

###########The proportion of each category can be calculated by adding a pie chart to the graph########


#low = nrow(out[out[,6]== "#2E8B57",])
#med = nrow(out[out[,6]== "#4682B4",])
#high = nrow(out[out[,6]== "#CD5C5C",])
#type <- c('med','high','low')
#nums <- c(med,high,low)
#df <- data.frame(type = type, nums = nums)

data_summary <- as.data.frame(table(out$group))
colnames(data_summary) <- c("group", "nums")
data_summary$type<-(c("Med","High","Low"))
data_summary$percentage <- round(data_summary$nums / sum(data_summary$nums) * 100, 1)
data_summary$label <- paste(data_summary$type, paste(data_summary$percentage, "%", sep = ''))

# pie chart
p2 <- ggplot(data_summary, aes(x = "", y = nums, fill = group)) +
  geom_bar(stat = "identity", width = 3) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c("#4682B4","#CD5C5C","#29A7A5"),  
    labels = data_summary$label
  ) +
  theme_void() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 12)
  )
p2

# Embed the pie chart into the main image
p_final <- p1 + geom_subview(subview = p2, x = -3.8, y = 0.15, w = 5, h = 5)

print(p_final)

