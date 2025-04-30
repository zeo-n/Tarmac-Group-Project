# =========================================================
# TARMAC GROUP PROJECT – BIOMASS ANALYSIS FOR BROOM (BR)  
# =========================================================

# Loading necessary libraries
library(tidyverse)
library(readxl)

# Step 1: Load only BR sheet
br_data <- read_excel("C Biomass and N Biomass (1) (1).xlsx", sheet = "BR")

# Step 2: Select important columns and create grouping variable
br_clean <- br_data %>%
  select(`SAMPLE CODE`, Site, `Year Established`, `Biomass carbon (ug/g)`) %>%
  mutate(Group = case_when(
    `Year Established` == 1996       ~ "1996",
    `Year Established` == 2004       ~ "2004",
    `Year Established` == 2011       ~ "2011",
    `Year Established` == 2017       ~ "2017",
    `Year Established` == "Existing" ~ "Existing",
    TRUE                             ~ "Other"    # fallback for unexpected entries
  ))

# View the first few cleaned rows
head(br_clean)

# Step 3: Basic Summary Statistics by Group with Standard Error
br_summary <- br_clean %>%
  group_by(Group) %>%
  summarise(
    Mean_Biomass = mean(`Biomass carbon (ug/g)`, na.rm = TRUE),
    SE_Biomass   = sd(`Biomass carbon (ug/g)`, na.rm = TRUE) / sqrt(n()),
    n            = n()
  )

print(br_summary)

# Step 4: ANOVA to compare Biomass Carbon across Groups
br_anova <- aov(`Biomass carbon (ug/g)` ~ Group, data = br_clean)
summary(br_anova)

# Optional: If the ANOVA is significant, run a Tukey HSD post-hoc test
# TukeyHSD(br_anova)

# Step 5: Boxplot of Biomass Carbon by Group
ggplot(br_clean, aes(x = Group, y = `Biomass carbon (ug/g)`, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Biomass Carbon (µg/g) by Year Group (BR Samples)",
    x     = "Group (Restoration Year or Existing)",
    y     = "Biomass Carbon (µg/g)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 6: Histogram of Biomass Carbon (optional)
ggplot(br_clean, aes(x = `Biomass carbon (ug/g)`, fill = Group)) +
  geom_histogram(alpha = 0.6, bins = 15, position = "identity") +
  theme_minimal() +
  labs(
    title = "Distribution of Biomass Carbon (BR Samples)",
    x     = "Biomass Carbon (µg/g)",
    y     = "Frequency"
  )

# =========================================================
# BAR PLOT: Mean Biomass Carbon (µg/g) per Group WITH STANDARD ERROR BARS
# =========================================================

# Compute error-bar limits (Mean ± SE)
br_summary <- br_summary %>%
  mutate(
    ymin = Mean_Biomass - SE_Biomass,
    ymax = Mean_Biomass + SE_Biomass
  )

# Create bar plot with SE error bars
ggplot(br_summary, aes(x = Group, y = Mean_Biomass, fill = Group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  theme_minimal() +
  labs(
    title = "Mean Biomass Carbon (µg/g) by Restoration Year / Existing Woodland (BR)",
    x     = "Group",
    y     = "Mean Biomass Carbon (µg/g)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =========================================================
# END OF SCRIPT
# =========================================================