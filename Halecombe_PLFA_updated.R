#=========================
# Halecombe Quarry PCA & Statistical Analysis
#=========================

# 1. Load libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)

# 2. Load Halecombe data
plfa_halecombe <- read_csv("PLFA Stats_Halecombe(transposed) (3).csv") %>%
  filter(!is.na(ID)) %>%
  mutate(Group = case_when(
    grepl("95", ID) ~ "1995",
    grepl("05", ID) ~ "2005",
    grepl("12", ID) ~ "2012",
    grepl("EX", ID) ~ "Reference",
    TRUE             ~ "Unknown"
  ))

# 3. Filter valid rows
valid_rows <- plfa_halecombe %>%
  select(-ID, -Group) %>%
  mutate(row_ok = if_all(everything(), ~ is.finite(.))) %>%
  pull(row_ok)

plfa_halecombe <- plfa_halecombe[valid_rows, ]
pca_input <- plfa_halecombe %>% select(-ID, -Group)

# 4. Run PCA
pca_model <- rda(pca_input, scale = TRUE)
eigvals <- eigenvals(pca_model)
var_pct <- eigvals / sum(eigvals) * 100

# 5. Scree plot with variance labels
scree_df <- tibble(PC = factor(1:length(eigvals)), Eigenvalue = eigvals, Variance = round(var_pct, 2))
scree_plot <- ggplot(scree_df, aes(x = PC, y = Eigenvalue)) +
  geom_col(fill = NA, color = "steelblue", linewidth = NA) +
  geom_line(aes(group = 1), color = "firebrick") +
  geom_point(color = "firebrick") +
  geom_text(aes(label = paste0(Variance, "%")), vjust = -0.3, size = 2) +
  theme_minimal(base_size = 14) +
  labs( x = "Principal Component", y = "Eigenvalue")
print(scree_plot)

# 6. Site & species scores (PC1 & PC2 only)
site_scores <- scores(pca_model, scaling = 2, display = "sites") %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  mutate(Sample = plfa_halecombe$ID, Group = plfa_halecombe$Group)

species_scores <- scores(pca_model, scaling = 2, display = "species") %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column("Biomarker")

# 7. PCA plot (PC1 vs PC2 with sample labels)
p1 <- ggplot(site_scores, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 3.2, max.overlaps = 25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  labs(x = sprintf("PC1 (%.2f%%)", var_pct[1]),
    y = sprintf("PC2 (%.2f%%)", var_pct[2])
  )
print(p1)

# 8. Correlation circle (PC1 vs PC2)
circle_df <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 200)),
  y = sin(seq(0, 2 * pi, length.out = 200))
)
corr_plot <- ggplot(species_scores, aes(PC1, PC2)) +
  geom_path(data = circle_df, aes(x, y), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.02, "npc")), color = "black") +
  geom_text_repel(aes(label = Biomarker), size = 3) +
  theme_minimal(base_size = 14) +
  labs(x = sprintf("PC1 (%.2f%%)", var_pct[1]),
       y = sprintf("PC2 (%.2f%%)", var_pct[2])) +
  coord_fixed()
print(corr_plot)

# 9. ANOVA + Tukey for PC1 and PC2
site_scores_stats <- site_scores %>% select(any_of(c("Group", "PC1", "PC2")))

for (pc in c("PC1", "PC2")) {
  aov_model <- aov(as.formula(paste(pc, "~ Group")), data = site_scores_stats)
  aov_summary <- summary(aov_model)
  capture.output(aov_summary, file = paste0("halecombe_anova_", tolower(pc), ".txt"))
  
  tuk <- TukeyHSD(aov_model)
  tuk_df <- as.data.frame(tuk$Group, rownames = "Comparison")
  write_csv(tuk_df, paste0("halecombe_tukey_", tolower(pc), ".csv"))
}

# 10. MANOVA (PC1 and PC2)
pcs <- intersect(c("PC1", "PC2"), colnames(site_scores))
man <- manova(as.formula(paste("cbind(", paste(pcs, collapse = ","), ") ~ Group")), data = site_scores_stats)
capture.output(summary(man), file = "halecombe_manova.txt")

# 11. ANOVA on fungal biomarker 18:2w6,9
fungal_col <- names(plfa_halecombe)[str_detect(names(plfa_halecombe), "18[:punct:]?2.*w6.*9")][1]

if (!is.na(fungal_col)) {
  fungal_aov <- aov(reformulate("Group", response = fungal_col), data = plfa_halecombe)
  capture.output(summary(fungal_aov), file = "halecombe_anova_fungal_18_2w6_9.txt")
  tuk_fungal <- TukeyHSD(fungal_aov)
  write_csv(as.data.frame(tuk_fungal$Group, rownames = "Comparison"), "halecombe_tukey_fungal_18_2w6_9.csv")
} else {
  message("Fungal marker 18:2w6,9 not found in column names.")
}
