# 加载所需R包
library(ggplot2)
library(tidyr)

# 1. 录入数据
temperature <- c(40, 46, 52, 58, 64, 70)
# 第二列数据为 TPL 组
tpl_data <- c(1, 0.92612348, 0.899072402, 0.649285426, 0.383236885, 0.159939283)
# 第三列数据为 DMSO 组
dmso_data <- c(1, 0.6685332, 0.4411345, 0.3414461, 0.2091509, 0.0811387)

# 构建数据框
cetsa_data <- data.frame(
  Temperature = temperature,
  TPL = tpl_data,
  DMSO = dmso_data   
)

# 将宽数据转换为长数据，适配 ggplot2 格式
data_long <- pivot_longer(cetsa_data,
                          cols = c("DMSO", "TPL"),
                          names_to = "Group",
                          values_to = "Relative_Level")

# 2. 绘制CETSA图
p <- ggplot(data_long, aes(x = Temperature, y = Relative_Level, color = Group, shape = Group)) +
  # 添加折线和数据点
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  
  # === 关键修改：替换为你指定的高级配色 ===
  # 注意：记得在色号前面加上 "#"
  scale_color_manual(values = c("DMSO" = "#8ebcdb", "TPL" = "#ace0cf")) +
  
  # 设置X轴和Y轴刻度
  scale_x_continuous(breaks = temperature) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.2, by = 0.2)) +
  
  # 设置坐标轴标签 (去除大标题)
  labs(
    x = "Temperature (°C)",
    y = "Relative S100A8 Level"
  ) +
  
  # 使用经典主题 (白底无网格线，适合期刊发表)
  theme_classic() +
  
  # 细节美化
  theme(
    axis.text = element_text(size = 12, color = "black"),       # 坐标轴数字
    axis.title = element_text(size = 14, face = "bold"),        # 坐标轴标题
    axis.line = element_line(linewidth = 0.8, color = "black"), # 坐标轴线
    axis.ticks = element_line(linewidth = 0.8, color = "black"),# 刻度线
    legend.position = c(0.85, 0.85),                            # 图例放在右上角
    legend.title = element_blank(),                             # 去掉图例标题
    legend.text = element_text(size = 12),                      # 图例文字大小
    legend.background = element_rect(fill = "transparent")      # 图例背景透明
  )

# 3. 保存为 PDF 文件
save_path <- "D:/桌面/sepsis/实验/S100A8_CETSA_CustomColor.pdf"

ggsave(filename = save_path, plot = p, width = 6, height = 4, device = "pdf")

print("绘图完成！新颜色的PDF已保存至指定路径。")
