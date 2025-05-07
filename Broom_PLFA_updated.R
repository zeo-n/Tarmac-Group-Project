#=========================
# Broom Quarry PCA & ANOVA (Styled like br_plfa.R)
#=========================

# Load libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(ggpubr)

# 1. Load Broom data
plfa_broom <- read_csv("PLFA Stats_Broom(transposed) (2).csv") %>%
  filter(!is.na(ID)) %>%
  mutate(Group = case_when(
    grepl("96", ID) ~ "1996",
    grepl("04", ID) ~ "2004",
    grepl("11", ID) ~ "2011",
    grepl("17", ID) ~ "2017",
    grepl("EX", ID) ~ "Reference",
    TRUE             ~ "Unknown"
  ))

# 2. Prepare PCA input
valid_rows <- plfa_broom %>%
  select(-ID, -Group) %>%
  mutate(row_ok = if_all(everything(), ~ is.finite(.))) %>%
  pull(row_ok)

plfa_broom <- plfa_broom[valid_rows, ]
pca_input <- plfa_broom %>% select(-ID, -Group)

# 3. Run PCA
pca_model <- rda(pca_input, scale = TRUE)
eigvals <- eigenvals(pca_model)
var_pct <- eigvals / sum(eigvals) * 100

# 4. Scree plot (Eigenvalues)
var_pct <- eigvals / sum(eigvals) * 100
scree_df <- tibble(PC = factor(1:length(eigvals)), Eigenvalue = eigvals, Variance = round(var_pct, 2))

scree_plot <- ggplot(scree_df, aes(x = PC, y = Eigenvalue)) +
  geom_col(fill = NA, color = "steelblue", linewidth = NA) +
  geom_line(aes(group = 1), color = "firebrick") +
  geom_point(color = "firebrick") +
  theme_minimal(base_size = 14) +
  geom_text(aes(label = paste0(Variance, "%")), vjust = -0.3, size = 2.0) +
  labs( x = "Principal Component", y = "Eigenvalue")

print(scree_plot)

# 5. PCA scores
site_scores <- scores(pca_model, scaling = 2, display = "sites") %>%
  as.data.frame() %>%
  mutate(Sample = plfa_broom$ID, Group = plfa_broom$Group)
species_scores <- scores(pca_model, scaling = 2, display = "species") %>%
  as.data.frame() %>%
  rownames_to_column("Biomarker")

# 6. PCA Plot (PC1 vs PC2)
p1 <- ggplot(site_scores, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 3.2, max.overlaps = 25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  labs( x = sprintf("PC1 (%.2f%%)", var_pct[1]),
       y = sprintf("PC2 (%.2f%%)", var_pct[2]))

print(p1)

# 7. Correlation Circle
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

# 8. ANOVA + Tukey (PC1–PC2 only)
for (pc in c("PC1", "PC2")) {
  aov_model <- aov(as.formula(paste(pc, "~ Group")), data = site_scores_stats)
  aov_summary <- summary(aov_model)
  capture.output(aov_summary, file = paste0("broom_anova_", tolower(pc), ".txt"))
  
  tuk <- TukeyHSD(aov_model)
  tuk_df <- as.data.frame(tuk$Group, rownames = "Comparison")
  write_csv(tuk_df, paste0("broom_tukey_", tolower(pc), ".csv"))
}

# 9. MANOVA
man <- manova(as.formula(paste("cbind(", paste(pcs, collapse = ","), ") ~ Group")), data = site_scores_stats)
capture.output(summary(man), file = "broom_manova.txt")

# 10. ANOVA on fungal biomarker 18:2w6,9
# Ensure correct column name exists (may vary depending on import)
fungal_col <- names(plfa_broom)[str_detect(names(plfa_broom), "18[:punct:]?2.*w6.*9")][1]

if (!is.na(fungal_col)) {
  fungal_aov <- aov(reformulate("Group", response = fungal_col), data = plfa_broom)
  capture.output(summary(fungal_aov), file = "Broom_ANOVA_fungal_18_2w6_9.txt")
  tuk_fungal <- TukeyHSD(fungal_aov)
  write_csv(as.data.frame(tuk_fungal$Group, rownames = "Comparison"), "Broom_TukeyHSD_fungal_18_2w6_9.csv")
} else {
  message("Fungal marker 18:2w6,9 not found in column names.")
}
