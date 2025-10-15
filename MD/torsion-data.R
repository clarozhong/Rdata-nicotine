# 扭转角分析 - 多图可视化并保存PDF
library(bio3d)

# 创建输出目录
output_dir <- "~/ZYJ/torsion/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 读取轨迹文件和PDB文件
dcd <- read.dcd("~/ZYJ/torsion/equ1-1_skip20.dcd")
pdb <- read.pdb("~/ZYJ/torsion/CaNi.pdb")

# 查看轨迹信息
cat("轨迹帧数:", nrow(dcd), "\n")
cat("总原子数:", nrow(pdb$atom), "\n")

# 使用您指定的原子索引
torsion_atoms <- c(3826, 3827, 3831, 3832)
cat("计算原子", torsion_atoms, "的扭转角\n")
cat("原子名称:", pdb$atom$elety[torsion_atoms], "\n")

# 提取坐标索引
xyz_inds <- atom2xyz(torsion_atoms)
cat("XYZ坐标索引:", xyz_inds, "\n")

# 手动实现向量叉积函数
vector_cross <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

# 自定义二面角计算函数
calculate_dihedral <- function(p1, p2, p3, p4) {
  # 计算向量
  b1 <- p2 - p1
  b2 <- p3 - p2
  b3 <- p4 - p3
  
  # 计算法向量（使用手动叉积）
  n1 <- vector_cross(b1, b2)
  n2 <- vector_cross(b2, b3)
  
  # 单位化法向量
  norm_n1 <- sqrt(sum(n1^2))
  norm_n2 <- sqrt(sum(n2^2))
  n1 <- n1 / norm_n1
  n2 <- n2 / norm_n2
  
  # 计算夹角
  cos_angle <- sum(n1 * n2)
  cos_angle <- max(-1, min(1, cos_angle))  # 确保在有效范围内
  
  # 计算角度（弧度）
  angle_rad <- acos(cos_angle)
  
  # 确定符号
  sign_check <- sum(b1 * vector_cross(n1, n2))
  if (sign_check < 0) {
    angle_rad <- -angle_rad
  }
  
  # 转换为角度
  angle_deg <- angle_rad * (180 / pi)
  return(angle_deg)
}

# 计算扭转角
cat("正在计算扭转角...\n")
torsion_angles <- numeric(nrow(dcd))
for(i in 1:nrow(dcd)) {
  coords <- matrix(dcd[i, xyz_inds], ncol = 3, byrow = TRUE)
  torsion_angles[i] <- calculate_dihedral(coords[1,], coords[2,], coords[3,], coords[4,])
}

# 时间轴 (1250帧对应50ns)
time_ns <- (1:nrow(dcd)) * (50 / nrow(dcd))

# 基本统计
cat("\n=== 基本统计 ===\n")
cat("平均值:", round(mean(torsion_angles), 2), "度\n")
cat("标准差:", round(sd(torsion_angles), 2), "度\n")
cat("最小值:", round(min(torsion_angles), 2), "度\n")
cat("最大值:", round(max(torsion_angles), 2), "度\n")
cat("范围:", round(max(torsion_angles) - min(torsion_angles), 2), "度\n")

# 1. 时间序列图
pdf(file.path(output_dir, "torsion_time_series.pdf"), width = 10, height = 6)
plot(time_ns, torsion_angles, type = "l", 
     xlab = "Time (ns)", ylab = "Torsion Angle (degrees)",
     main = "Torsion Angle Variation over 50ns MD Simulation",
     col = "blue", lwd = 2, ylim = c(-180, 180))

# 添加参考线
abline(h = mean(torsion_angles), col = "red", lty = 2, lwd = 2)
abline(h = 0, col = "gray", lty = 3)
abline(h = c(-180, -90, 90, 180), col = "lightgray", lty = 3)

legend("topright", 
       legend = c(paste("Mean:", round(mean(torsion_angles), 2), "°"),
                  paste("SD:", round(sd(torsion_angles), 2), "°"),
                  paste("Range: [", round(min(torsion_angles), 1), ",", 
                        round(max(torsion_angles), 1), "]")),
       bty = "n", cex = 0.8)
dev.off()

# 2. 分布直方图
pdf(file.path(output_dir, "torsion_histogram.pdf"), width = 8, height = 6)
hist_data <- hist(torsion_angles, breaks = 30, plot = FALSE)
hist(torsion_angles, breaks = 30, 
     main = "Torsion Angle Distribution",
     xlab = "Torsion Angle (degrees)", 
     ylab = "Frequency",
     col = "lightblue", border = "black",
     xlim = c(-180, 180))

abline(v = mean(torsion_angles), col = "red", lwd = 2, lty = 2)
abline(v = mean(torsion_angles) + sd(torsion_angles), col = "orange", lwd = 1, lty = 3)
abline(v = mean(torsion_angles) - sd(torsion_angles), col = "orange", lwd = 1, lty = 3)

legend("topright", 
       legend = c(paste("Mean:", round(mean(torsion_angles), 2)),
                  paste("SD:", round(sd(torsion_angles), 2)),
                  "Mean ± SD"),
       col = c("red", "orange", "orange"),
       lty = c(2, 3, 3), lwd = c(2, 1, 1),
       bty = "n")
dev.off()

# 3. 密度图
pdf(file.path(output_dir, "torsion_density.pdf"), width = 8, height = 6)
density_data <- density(torsion_angles)
plot(density_data, 
     main = "Torsion Angle Density Distribution",
     xlab = "Torsion Angle (degrees)",
     ylab = "Density",
     col = "darkblue", lwd = 2,
     xlim = c(-180, 180))

polygon(density_data, col = rgb(0, 0, 1, 0.3), border = "darkblue")
abline(v = mean(torsion_angles), col = "red", lwd = 2, lty = 2)

legend("topright", 
       legend = c("Density", paste("Mean:", round(mean(torsion_angles), 2))),
       col = c("darkblue", "red"),
       lty = c(1, 2), lwd = c(2, 2),
       bty = "n")
dev.off()

# 4. 箱线图
pdf(file.path(output_dir, "torsion_boxplot.pdf"), width = 6, height = 8)
boxplot(torsion_angles, 
        main = "Torsion Angle Distribution",
        ylab = "Torsion Angle (degrees)",
        col = "lightgreen",
        ylim = c(-180, 180))

# 添加统计信息
text(1, max(torsion_angles) - 10, 
     paste("Mean:", round(mean(torsion_angles), 2)), pos = 3, cex = 0.8)
text(1, min(torsion_angles) + 10, 
     paste("SD:", round(sd(torsion_angles), 2)), pos = 1, cex = 0.8)
dev.off()

# 5. 多图组合
pdf(file.path(output_dir, "torsion_comprehensive.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# 5.1 时间序列
plot(time_ns, torsion_angles, type = "l", 
     xlab = "Time (ns)", ylab = "Torsion Angle (degrees)",
     main = "Time Series", col = "blue", lwd = 1)
abline(h = mean(torsion_angles), col = "red", lty = 2)

# 5.2 直方图
hist(torsion_angles, breaks = 30, 
     main = "Distribution",
     xlab = "Torsion Angle (degrees)", 
     col = "lightblue", border = "black")
abline(v = mean(torsion_angles), col = "red", lwd = 2, lty = 2)

# 5.3 密度图
plot(density(torsion_angles), 
     main = "Density",
     xlab = "Torsion Angle (degrees)",
     col = "darkgreen", lwd = 2)
polygon(density(torsion_angles), col = rgb(0, 1, 0, 0.3))

# 5.4 Q-Q图（正态性检验）
qqnorm(torsion_angles, main = "Q-Q Plot (Normality Check)")
qqline(torsion_angles, col = "red")

mtext("Comprehensive Torsion Angle Analysis", outer = TRUE, line = -2, cex = 1.5)
dev.off()

# 6. 累积分布函数图
pdf(file.path(output_dir, "torsion_ecdf.pdf"), width = 8, height = 6)
ecdf_data <- ecdf(torsion_angles)
plot(ecdf_data, 
     main = "Empirical Cumulative Distribution Function",
     xlab = "Torsion Angle (degrees)",
     ylab = "Cumulative Probability",
     col = "purple", lwd = 2)

# 添加分位数线
q50 <- quantile(torsion_angles, 0.5)
q25 <- quantile(torsion_angles, 0.25)
q75 <- quantile(torsion_angles, 0.75)

abline(v = q50, col = "red", lty = 2, lwd = 2)
abline(v = c(q25, q75), col = "orange", lty = 3, lwd = 1)

legend("topleft", 
       legend = c("ECDF", 
                  paste("Median:", round(q50, 2)),
                  paste("Q1-Q3: [", round(q25, 1), ",", round(q75, 1), "]")),
       col = c("purple", "red", "orange"),
       lty = c(1, 2, 3), lwd = c(2, 2, 1),
       bty = "n")
dev.off()

# 保存详细统计结果
stats_summary <- data.frame(
  Statistic = c("Frames", "Mean", "SD", "Min", "Max", "Range", 
                "Q1", "Median", "Q3", "IQR"),
  Value = c(length(torsion_angles),
            round(mean(torsion_angles), 3),
            round(sd(torsion_angles), 3),
            round(min(torsion_angles), 3),
            round(max(torsion_angles), 3),
            round(max(torsion_angles) - min(torsion_angles), 3),
            round(quantile(torsion_angles, 0.25), 3),
            round(quantile(torsion_angles, 0.50), 3),
            round(quantile(torsion_angles, 0.75), 3),
            round(IQR(torsion_angles), 3))
)

write.csv(stats_summary, file.path(output_dir, "torsion_statistics.csv"), row.names = FALSE)
write.csv(data.frame(Frame = 1:nrow(dcd), Time_ns = time_ns, Torsion_Angle = torsion_angles),
          file.path(output_dir, "torsion_angle_data.csv"), row.names = FALSE)

# 输出生成的文件列表
cat("\n=== 生成的文件 ===\n")
generated_files <- list.files(output_dir, full.names = TRUE)
for(file in generated_files) {
  file_size <- file.info(file)$size / 1024  # KB
  cat(basename(file), sprintf("(%.1f KB)\n", file_size))
}

cat("\n所有图表已保存到:", output_dir, "\n")
cat("包括:\n")
cat("- torsion_time_series.pdf: 时间序列图\n")
cat("- torsion_histogram.pdf: 分布直方图\n")
cat("- torsion_density.pdf: 密度分布图\n")
cat("- torsion_boxplot.pdf: 箱线图\n")
cat("- torsion_comprehensive.pdf: 综合多图\n")
cat("- torsion_ecdf.pdf: 累积分布函数图\n")
cat("- torsion_statistics.csv: 统计摘要\n")
cat("- torsion_angle_data.csv: 原始数据\n")

