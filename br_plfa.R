# =========================================================
# PLFA BROOM QUARRY – PCA + ONE-WAY ANOVA + STYLED SCREE PLOT
# =========================================================

library(tidyverse)
library(vegan)
library(ggplot2)

# Step 1: Load Broom data
plfa_broom <- read.csv("PLFA_Broom_Analysis_Input.csv")

# Step 2: Create group (from Year Established)
plfa_broom <- plfa_broom %>%
  mutate(Group = as.character(`Year.Established`))

# Step 3: Prepare input for PCA
pca_input <- plfa_broom %>%
  select(-Sample.Number, -Sample, -Year.Established, -Quarry, -Group)

# Step 4: Run PCA (scaled)
plfa_pca <- rda(pca_input, scale = TRUE)

# Step 5: Extract PCA scores (site scores = sample coordinates)
pca_scores <- as.data.frame(scores(plfa_pca, scaling = 2)$sites)
pca_scores$Sample <- plfa_broom$Sample
pca_scores$Group <- plfa_broom$Group

# Step 6: Styled Scree Plot (Eigenvalues like your reference)
eigenvalues <- eigenvals(plfa_pca)
scree_df <- data.frame(PC = seq_along(eigenvalues), Eigenvalue = eigenvalues)

ggplot(scree_df, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "orange", linewidth = 1) +
  geom_point(color = "orange", size = 2) +
  scale_x_continuous(breaks = scree_df$PC) +
  theme_minimal(base_size = 14) +
  labs(title = "Scree Plot of PCA Eigenvalues (BROOM)",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"))

# Step 7: PCA plot (PC1 vs PC2) with sample labels, no origin lines
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  theme_minimal() +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Broom Quarry – PCA (PC1 vs PC2)", x = "PC1", y = "PC2")

# Step 8: Total PLFA per sample
plfa_broom$Total_PLFA <- rowSums(pca_input, na.rm = TRUE)

# Step 9: One-way ANOVA on Total PLFA by Group
anova_result <- aov(Total_PLFA ~ Group, data = plfa_broom)
summary(anova_result)

# Step 10: Export results
write.csv(pca_scores, "Broom_PCA_Sample_Scores.csv", row.names = FALSE)
capture.output(summary(anova_result), file = "Broom_PCA_ANOVA_Result.txt")
