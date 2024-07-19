library(ggplot2)
library(reshape2)
library(scales)

heatmap_values = read.table("spreadsheets/normalized_majorclass_Z_expression.csv", sep = ',', header = TRUE)
heatmap_values_long = melt(heatmap_values, variable.name = "gene", value.name = 'Z')

zheat = ggplot(heatmap_values_long, aes(gene, majorclass, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", limits = c(-2.5,2.5), oob = squish) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black")
  )

ggsave(filename = "figures/zheat_majorclass.png", zheat)

zheat = ggplot(heatmap_values_long, aes(gene, majorclass, fill = Z > 1.5)) +
  geom_tile() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black")
  )

ggsave(filename = "figures/zheat_majorclass_threshold1.5.png", zheat)

make_threshold_heatmap <- function(Zthres) {
  zheat = ggplot(heatmap_values_long, aes(gene, majorclass, fill = Z > Zthres)) +
    geom_tile() +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.text.y = element_text(color = "black")
    ) +
    labs(fill = paste0("Z > ", Zthres))

  ggsave(filename = paste0("figures/zheat_majorclass_threshold.", Zthres, ".png"), zheat)
}
make_threshold_heatmap(2.5)

library(data.table)

metadata = fread("spreadsheets/obs.csv")

raw_counts_celltype = ggplot(metadata, aes(x = majorclass, y = nCount_RNA)) +
  geom_boxplot(outlier.size = 0) + theme_bw()
ggsave("figures/raw_counts_celltype_boxplot.png", raw_counts_celltype)
