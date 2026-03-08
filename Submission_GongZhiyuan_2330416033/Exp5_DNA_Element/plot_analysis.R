# 定义文件
file_name <- "score_distance.txt"

# 读取数据
data <- read.table(file_name, header = TRUE, sep = "\t")

# 1. 散点图 (Score vs Distance) [cite: 311]
png(file = "plot_scatter.png", width=800, height=800)
plot(data$Distance, data$Score, 
     type="p", pch=20, col=rgb(0,0,1,0.5),
     main="TATA-box Score vs Distance to TSS",
     xlab="Distance (bp)", ylab="PWM Score")
# 添加趋势线
abline(lm(data$Score ~ data$Distance), col="red", lwd=2)
dev.off()

# 2. 距离频数分布直方图 [cite: 372]
png(file = "plot_hist_distance.png", width=800, height=600)
hist(data$Distance, 
     breaks=100, 
     col="lightblue", 
     main="Distribution of Distances to Nearest Gene",
     xlab="Distance (bp)")
# 标注中位数
abline(v=median(data$Distance), col="red", lwd=2, lty=2)
legend("topright", legend=c("Median"), col=c("red"), lty=2, lwd=2)
dev.off()

# 3. 得分频数分布直方图 [cite: 338]
png(file = "plot_hist_score.png", width=800, height=600)
hist(data$Score, 
     breaks=50, 
     col="lightgreen", 
     main="Distribution of TATA-box Scores",
     xlab="PWM Score")
dev.off()

# 4. 计算相关性 [cite: 295]
cor_result <- cor.test(data$Distance, data$Score, method = "pearson")
print(cor_result)
