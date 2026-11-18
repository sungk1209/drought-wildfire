library(tidyverse)
library(patchwork)

# Figure 3 (same design as figs/Figure3.png).
# Working directory: project root.
# Fix: scale_x_discrete(limits = bins) in severity_plot(). Without limits, the
# reference bin drawn in a separate layer was placed at the end of the x axis
# while the labels kept the original order, so panel (c) was shifted by one bin.

dir.create("figs/revision_final", showWarnings = FALSE, recursive = TRUE)
red <- "#B2182B"
index_colours <- c("SPI-3" = red, "SPEI-3" = "#E08214")
bins <- c("> 0.5", "-0.5 to 0.5", "-1.0 to -0.5", "-1.5 to -1.0", "<= -1.5")
bin_labels <- c("> 0.5", "−0.5 to 0.5", "−1.0 to −0.5", "−1.5 to −1.0", "≤ −1.5")

figure_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.tag = element_text(face = "bold"))

# a. Binary drought comparison
binary <- read_csv("data_out/tables/Table2_main_regression.csv", show_col_types = FALSE) %>%
  filter(str_detect(model, "^M[13] ")) %>%
  mutate(adjustment = if_else(str_starts(model, "M1"), "Unadjusted", "Fully adjusted"),
         outcome_name = recode(outcome, "log(burned area)" = "Burned area",
                               "log(containment duration)" = "Suppression duration"),
         outcome_label = paste0(outcome_name, "\n(n = ", scales::comma(n), ")"),
         outcome_label = fct_reorder(outcome_label, as.integer(outcome_name == "Burned area")),
         across(c(multiplier, ci_lo, ci_hi), ~ 100 * (.x - 1)),
         y_pos = as.numeric(outcome_label) + if_else(adjustment == "Unadjusted", 0.07, -0.07))

p_binary <- ggplot(binary, aes(multiplier, y_pos, colour = adjustment)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.06) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Unadjusted" = "grey55", "Fully adjusted" = red),
                      breaks = c("Unadjusted", "Fully adjusted"),
                      labels = c("Unadjusted", "Adjusted")) +
  scale_y_continuous(breaks = 1:2, labels = levels(binary$outcome_label), limits = c(0.5, 2.5)) +
  scale_x_continuous(breaks = seq(0, 50, 10), limits = c(0, 52)) +
  labs(x = "Difference from non-drought (%)", y = NULL, tag = "a") + figure_theme +
  theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
        legend.direction = "horizontal")

# b–c. Estimates by drought severity. The open diamond marks the reference bin.
effects <- read_csv("data_out/tables/revision/Figure3_index_bin_effects_rev.csv", show_col_types = FALSE) %>%
  filter(sens == "main") %>%
  mutate(bin = factor(bin, levels = bins), reference = bin == "-0.5 to 0.5",
         index = factor(index, levels = names(index_colours)))

severity_plot <- function(source_row, outcome, y_title, x_title = "Drought Severity") {
  d <- effects %>% filter(row == source_row, y == outcome)
  compare_indices <- source_row == "a"
  dodge <- position_dodge(width = if (compare_indices) 0.45 else 0)

  p <- ggplot(d, aes(bin, effect_pct, colour = index, group = index)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4)

  if (compare_indices) p <- p + geom_line(position = dodge, linewidth = 0.45)
  p +
    geom_errorbar(data = filter(d, !reference), aes(ymin = lo_pct, ymax = hi_pct),
                  width = 0.14, position = dodge) +
    geom_point(data = filter(d, !reference), aes(shape = index), position = dodge, size = 2) +
    geom_point(data = filter(d, reference), shape = 23, fill = "white",
               position = dodge, size = 2, show.legend = FALSE) +
    scale_colour_manual(values = index_colours) +
    scale_shape_manual(values = c("SPI-3" = 16, "SPEI-3" = 15)) +
    scale_x_discrete(limits = bins, labels = bin_labels) +
    coord_cartesian(ylim = if (compare_indices) c(-12, 52) else c(-29, 34), expand = FALSE,
                    xlim = c(0.5, 5.5)) +
    labs(x = x_title, y = y_title) +
    figure_theme + theme(axis.text.x = element_text(angle = 50, hjust = 1),
                         legend.position = "none")
}

p_area <- severity_plot("a", "log_burned", "Burned area\ndifference (%)") +
  labs(tag = "b") +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1),
        legend.direction = "horizontal")
p_duration <- severity_plot("a", "log_dur", "Suppression duration\ndifference (%)")
p_heli <- severity_plot("b", "log_heli", "Helicopter deployment\ndifference (%)", NULL) + labs(tag = "c")
p_crew <- severity_plot("b", "log_crew", "Personnel deployment\ndifference (%)")
p_resource_duration <- severity_plot("b", "log_dur", "Suppression duration\ndifference (%)", NULL)

# wrap_elements(): panel (a) is not aligned with (b) and (c), so its long outcome
# labels do not push the y-axis titles of (b) and (c) away from their panels.
figure3 <- wrap_elements(full = p_binary) / (p_area | p_duration) /
  (p_heli | p_crew | p_resource_duration) +
  plot_layout(heights = c(0.75, 1, 1))

print(figure3)
ggsave("figs/revision_final/Figure3.png", figure3, width = 7, height = 7, units = "in", dpi = 600, bg = "white")
ggsave("figs/revision_final/Figure3.pdf", figure3, width = 7, height = 7, units = "in", bg = "white")
