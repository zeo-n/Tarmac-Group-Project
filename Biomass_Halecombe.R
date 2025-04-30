# =========================================================
# TARMAC GROUP PROJECT - BIOMASS ANALYSIS (ONLY {HC} HALECOMBE SHEET GROUPED)
# =========================================================

# Loading necessary libraries
library(tidyverse)

# Step 1: Load HC CSV file
hc_data <- read_csv("C Biomass and N Biomass (1)(HC).csv")

# Step 2: Select important columns and create grouping variable
hc_clean <- hc_data %>%
  select(`SAMPLE CODE`, Site, `Year Established`, `Biomass carbon (ug/g)`) %>%
  mutate(Group = case_when(
    `Year Established` == "1980"      ~ "1980",
    `Year Established` == "1995"      ~ "1995",
    `Year Established` == "2005"      ~ "2005",
    `Year Established` == "2012"      ~ "2012",
    `Year Established` == "Existing"  ~ "Existing",
    TRUE                              ~ "Other"    # fallback for any unexpected entries
  ))

# View first few rows
head(hc_clean)

# Step 3: Basic Summary Statistics by Group
hc_summary <- hc_clean %>%
  group_by(Group) %>%
  summarise(
    Mean_Biomass = mean(`Biomass carbon (ug/g)`, na.rm = TRUE),
    SD_Biomass   = sd(`Biomass carbon (ug/g)`, na.rm = TRUE),
    n            = n()
  )

print(hc_summary)

# Step 4: ANOVA to compare Biomass Carbon across Groups
hc_anova <- aov(`Biomass carbon (ug/g)` ~ Group, data = hc_clean)
summary(hc_anova)

# Optional: If the ANOVA is significant, run a Tukey HSD post-hoc test
# TukeyHSD(hc_anova)

# Step 5: Boxplot of Biomass Carbon by Group
ggplot(hc_clean, aes(x = Group, y = `Biomass carbon (ug/g)`, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Biomass Carbon (µg/g) by Year Group (HC Samples)",
    x     = "Group (Year Established or Existing)",
    y     = "Biomass Carbon (µg/g)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 6: Histogram of Biomass Carbon (optional)
ggplot(hc_clean, aes(x = `Biomass carbon (ug/g)`, fill = Group)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  theme_minimal() +
  labs(
    title = "Distribution of Biomass Carbon (HC Samples)",
    x     = "Biomass Carbon (µg/g)",
    y     = "Frequency"
  )

# =========================================================
# BAR PLOT: Mean Biomass Carbon (µg/g) per Group WITH ERROR BARS
# =========================================================

# Compute error‐bar limits
hc_summary <- hc_summary %>%
  mutate(
    ymin = Mean_Biomass - SD_Biomass,
    ymax = Mean_Biomass + SD_Biomass
  )

# Create bar plot with error bars
ggplot(hc_summary, aes(x = Group, y = Mean_Biomass, fill = Group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  theme_minimal() +
  labs(
    title = "Mean Biomass Carbon (µg/g) by Restoration Year / Existing Woodland (HC)",
    x     = "Group",
    y     = "Mean Biomass Carbon (µg/g)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =========================================================
# END OF SCRIPT
# =========================================================