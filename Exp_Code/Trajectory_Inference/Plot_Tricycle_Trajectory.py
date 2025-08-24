
library(ggplot2)
library(dplyr)
library(patchwork) 

file_path <- "/root/Cycle/HarmoCycle_v1/Fig3/Projection_U5.csv"
name_tag = 'Tricycle_U5'
project_df <- read.csv(file_path)


project_df$true <- project_df$stage
# 把 true 为 0 的换为 G1
project_df$true[project_df$true == 0] <- "G1"

# 把 true 为 1 的换为 S
project_df$true[project_df$true == 1] <- "S"

# 把 true 为 2 的换为 G2M
project_df$true[project_df$true == 2] <- "G2M"

# 查看结果
project_df$pca_adjust_angle <- project_df$angles



# 确保 'stage' 列被视为分类变量 (factor)
project_df$stage <- as.factor(project_df$true)


# --- 步骤 2: 定义通用颜色和主题 ---

# 定义和参考代码一致的离散颜色方案
all_stage_colors <- c(
  "G1"  = "#C85B43", 
  "S"   = "#4F85A2", 
  "G2M" = "#E3A897",
  # 以下为备用，如果您的数据包含这些分期
  "G1.S" = "#E69F00", 
  "G2" = "#56B4E9", 
  "G2.M" = "#009E73", # G2M 和 G2.M 使用相同颜色
  "M.G1" = "#CC79A7", 
  "Unclassified" = "grey60"
)

# 创建一个可重用的自定义主题，以匹配参考风格
# 这能让代码更整洁
theme_custom <- function() {
  theme_classic() +
  theme(
    # --- 标题样式 (模拟分面标题) ---
    plot.title = element_text(
      size = 16, 
      face = "bold", 
      hjust = 0.5, # 标题居中
      margin = margin(b = 10) # 标题和图形之间的距离
    ),
    
    # --- 坐标轴样式 ---
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
    
    # --- 图例样式 ---
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )
}


# --- 步骤 3: 创建第一张图 (按角度着色) ---

p_angle <- ggplot(project_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pca_adjust_angle), size = 5, alpha = 0.7) +
  # 使用 'viridis' 色谱，它在视觉上比默认色谱更优秀
  scale_color_viridis_c(option = "C") + 
  labs(
      title = paste0("Angle (", name_tag, " et al)"),
      x = "PC 1",
      y = "PC 2",
      color = "Angle" # 设置连续色标的图例标题
  ) +
  theme_custom()


# --- 步骤 4: 创建第二张图 (按真实Stage着色) ---

p_stage <- ggplot(project_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = stage), size = 5, alpha = 0.7) +
  # 应用我们预定义的颜色
  scale_color_manual(
    name = "Cell Cycle Stage",
    values = all_stage_colors
  ) +
  labs(
      title = paste0("True Stage (", name_tag, " et al)"),
      x = "PC 1",
      y = "PC 2"
  ) +
  theme_custom() +
  # 调整图例中点的大小，使其更清晰
  guides(
    color = guide_legend(
      override.aes = list(size = 5), # 对应 override.aes
      title.position = "top", 
      title.hjust = 0.5
    )
  )


# --- 步骤 5: 使用 patchwork 组合两张图 ---

# 使用 '+' 操作符将两张图并排排列
final_plot <- p_angle + p_stage

# 显示最终的组合图
print(final_plot)

# 如果想保存图片True Stage (
ggsave(paste0("/root/Cycle/HarmoCycle_v1/Fig3/", name_tag, ".pdf"), final_plot, width = 12, height = 6, dpi = 300)
print(project_df)
