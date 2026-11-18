library(tidyverse)
library(patchwork)
dir.create("figs/revision_final", showWarnings = FALSE, recursive = TRUE)

# Figure 5 (same design as the figure in the manuscript).
# Working directory: project root.
# (a) incidents in the binary comparison with complete covariates
#     (non-drought n = 8,134; drought-associated n = 2,512)
# (b) quantile-regression ratios from Figure5_quantile_effects_rev.csv

red  <- "#B2182B"
grey <- "grey55"

figure_theme <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold"))

# ---------------------------------------------------------------- (a) distributions
incidents <- read_csv("data_out/incident_panel.csv", show_col_types = FALSE) %>%
  filter(!is.na(compound)) %>%
  drop_na(log_burned, rhmin, wsmax, tmean, dry_spell, wet_days_mo, cause, sido, season, year) %>%
  mutate(group = if_else(compound == 1, "Drought-associated", "Non-drought"),
         group = factor(group, levels = c("Non-drought", "Drought-associated")),
         log10_ha = log10(burned_ha))

n_group <- incidents %>% count(group)

x_breaks <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
x_labels <- c("0.01", "0.1", "1", "10", "100", "1,000", "10,000")

box_y <- c("Non-drought" = -0.12, "Drought-associated" = -0.25)
box_df <- incidents %>%
  group_by(group) %>%
  summarise(q1 = quantile(log10_ha, 0.25), med = median(log10_ha), q3 = quantile(log10_ha, 0.75),
            lo = max(min(log10_ha), q1 - 1.5 * (q3 - q1)),
            hi = min(max(log10_ha), q3 + 1.5 * (q3 - q1))) %>%
  mutate(y = box_y[as.character(group)])

label_df <- tibble(group = factor(c("Non-drought", "Drought-associated"), levels = levels(incidents$group)),
                   x = log10(25000), y = c(0.62, 0.42)) %>%
  left_join(n_group, by = "group") %>%
  mutate(label = paste0(group, "\n(n = ", scales::comma(n), ")"))

p_a <- ggplot(incidents, aes(log10_ha, fill = group, colour = group)) +
  geom_density(alpha = 0.35, linewidth = 0.5, bw = "nrd0") +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_segment(data = box_df, aes(x = lo, xend = hi, y = y, yend = y), inherit.aes = FALSE,
               colour = c(grey, red), linewidth = 0.4) +
  geom_segment(data = box_df, aes(x = lo, xend = lo, y = y - 0.015, yend = y + 0.015), inherit.aes = FALSE,
               colour = c(grey, red), linewidth = 0.4) +
  geom_segment(data = box_df, aes(x = hi, xend = hi, y = y - 0.015, yend = y + 0.015), inherit.aes = FALSE,
               colour = c(grey, red), linewidth = 0.4) +
  geom_rect(data = box_df, aes(xmin = q1, xmax = q3, ymin = y - 0.04, ymax = y + 0.04, fill = group),
            inherit.aes = FALSE, colour = c(grey, red), alpha = 0.8, linewidth = 0.4) +
  geom_segment(data = box_df, aes(x = med, xend = med, y = y - 0.04, yend = y + 0.04), inherit.aes = FALSE,
               colour = "black", linewidth = 0.7) +
  geom_text(data = label_df, aes(x = x, y = y, label = label, colour = group), inherit.aes = FALSE,
            hjust = 1, vjust = 1, size = 3.5, lineheight = 0.95) +
  scale_fill_manual(values = c("Non-drought" = grey, "Drought-associated" = red)) +
  scale_colour_manual(values = c("Non-drought" = grey, "Drought-associated" = red)) +
  scale_x_continuous(breaks = log10(x_breaks), labels = x_labels) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.2), labels = sprintf("%.1f", seq(0, 0.6, 0.2))) +
  coord_cartesian(xlim = c(-2, log10(70000)), ylim = c(-0.32, 0.7), expand = FALSE) +
  labs(x = "Burned area (ha; log scale)", y = "Density", tag = "a") +
  figure_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# ---------------------------------------------------------------- (b) quantile ratios
qr <- read_csv("data_out/tables/revision/Figure5_quantile_effects_rev.csv", show_col_types = FALSE) %>%
  mutate(q_label = paste0(round(tau * 100), "th"),
         q_label = factor(q_label, levels = q_label))

p_b <- ggplot(qr, aes(q_label, ratio, group = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(colour = red, linewidth = 0.5) +
  geom_errorbar(data = filter(qr, tau > 0.1), aes(ymin = lo_ratio, ymax = hi_ratio),
                width = 0.12, colour = red) +
  geom_point(colour = red, size = 2) +
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 4), labels = c("×1", "×1.5", "×2", "×3", "×4")) +
  coord_cartesian(ylim = c(0.93, 4.4)) +
  labs(x = "Quantile", y = "Ratio (drought / non-drought)", tag = "b") +
  figure_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

figure5 <- p_a | p_b

print(figure5)
ggsave("figs/revision_final/Figure5.png", figure5, width = 7.2, height = 3.3, units = "in", dpi = 600, bg = "white")
ggsave("figs/revision_final/Figure5.pdf", figure5, width = 7.2, height = 3.3, units = "in", bg = "white")
