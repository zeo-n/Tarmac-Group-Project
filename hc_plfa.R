# =========================================================
# PLFA HALECOMBE QUARRY – PCA + ONE-WAY ANOVA + STYLED SCREE PLOT
# =========================================================

library(tidyverse)
library(vegan)
library(ggplot2)

# Step 1: Load Halecombe data
plfa_halecombe <- read.csv("PLFA_Halecombe_Analysis_Input.csv")

# Step 2: Create Group column from Year Established
plfa_halecombe <- plfa_halecombe %>%
  mutate(Group = as.character(Year.Established))  # use correct column name here

# Step 3: Prepare PCA input (remove non-compound metadata)
pca_input <- plfa_halecombe %>%
  select(-Sample.Number, -Sample, -Year.Established, -Quarry, -Group)

# Step 4: Run PCA
plfa_pca <- rda(pca_input, scale = TRUE)

# Step 5: Extract PCA site scores (sample positions)
pca_scores <- as.data.frame(scores(plfa_pca, scaling = 2)$sites)
pca_scores$Sample <- plfa_halecombe$Sample
pca_scores$Group <- plfa_halecombe$Group

# Step 6: Scree Plot – Line + Point style (styled like example)
eigenvalues <- eigenvals(plfa_pca)
scree_df <- data.frame(PC = seq_along(eigenvalues), Eigenvalue = eigenvalues)

ggplot(scree_df, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "orange", linewidth = 1) +
  geom_point(color = "orange", size = 2) +
  scale_x_continuous(breaks = scree_df$PC) +
  theme_minimal(base_size = 14) +
  labs(title = "Scree Plot of PCA Eigenvalues – Halecombe",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"))

# Step 7: PCA plot – with sample names and bold black origin lines
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Halecombe Quarry – PCA (PC1 vs PC2)",
       x = "PC1", y = "PC2")

# Step 8: Calculate Total PLFA for each sample
plfa_halecombe$Total_PLFA <- rowSums(pca_input, na.rm = TRUE)

# Step 9: One-way ANOVA on Total PLFA by Group
anova_result <- aov(Total_PLFA ~ Group, data = plfa_halecombe)
summary(anova_result)

# Step 10: Export results
write.csv(pca_scores, "Halecombe_PCA_Sample_Scores.csv", row.names = FALSE)
capture.output(summary(anova_result), file = "Halecombe_PCA_ANOVA_Result.txt")
