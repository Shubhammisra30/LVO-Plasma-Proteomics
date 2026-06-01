##########################################################################################################################################
#All samples analysis (GMH + YNHH)
##########################################################################################################################################

#########################################################################
# Scatter plot for concordance-discordance analysis (LVO vs. non-LVO AIS)
#########################################################################

library(ggplot2)
library(ggrepel)
library(readr)

# Load the CSV file
common_proteins <- read_csv("Concordant_DEPs_GMH_YNHH.csv", show_col_types = FALSE)
View(common_proteins)

# Define color mapping for 'Direction'
direction_colors <- c(
  "Concordant UP" = "darkgreen",
  "Discordant" = "red"
)

# Plot with fixed axis limits (0 to 4) and 0.5 breaks
ggplot(common_proteins, aes(
  x = GMH_FC, 
  y = YNHH_FC, 
  label = Protein_ID, 
  colour = Direction
)) +
  geom_point() +
  scale_colour_manual(values = direction_colors) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text_repel(max.overlaps = Inf) +
  labs(
    x = "Fold Change (GMH)",
    y = "Fold Change (YNHH)",
    color = "Direction",
    title = "Fold Changes: Concordant and Discordant DEPs in GMH vs. YNHH cohorts for classifying LVO"
  ) +
  theme_minimal()

#################################################################################################################################################

#####################################################################################################################################
# Sub-analysis of Samples within 30 hours from Yale cohort (GMH + YNHH_30 hrs)
#####################################################################################################################################

##########################################################################
# Scatter plot for concordance-discordance analysis (LVO vs. non-LVO AIS)
##########################################################################

library(ggplot2)
library(ggrepel)
library(readr)

# Load the CSV file
common_proteins <- read_csv("Concordant_DEPs_GMH_YNHH_30hrs.csv", show_col_types = FALSE)

# Define color mapping for 'Direction'
direction_colors <- c(
  "Concordant UP" = "darkgreen",
  "Concordant DOWN" = "blue",
  "Discordant" = "red"
)

# Plot with label repel settings to reduce overlap
ggplot(common_proteins, aes(
  x = GMH_FC, 
  y = YNHH_30hrs_FC, 
  label = Protein_ID, 
  colour = Direction
)) +
  geom_point() +
  scale_colour_manual(values = direction_colors) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text_repel(
    max.overlaps = Inf,
    box.padding = 0.5,         # Space around labels
    point.padding = 0.3,       # Space between point and label
    force = 1.5,               # Force to repel overlapping labels
    force_pull = 0.1,          # Pull labels closer to origin
    segment.size = 0.2,        # Thinner lines
    segment.color = "gray70",  # Lighter color for clarity
    min.segment.length = 0     # Draw all leader lines
  ) +
  labs(
    x = "Fold Change (GMH)",
    y = "Fold Change (YNHH_30hrs)",
    color = "Direction",
    title = "Fold Changes: Cross-validated Proteins (GMH vs. Sub-YNHH within 30 hours)"
  ) +
  theme_minimal()
