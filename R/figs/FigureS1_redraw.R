library(tidyverse)
dir.create("figs/revision_final", showWarnings = FALSE, recursive = TRUE)

# Figure S1 (same design as figures_v2/FigureS1.png).
# Working directory: project root.
# Fix: group = family in geom_point. The points were dodged by the fill groups
# (family x significance) while the error bars were dodged by family, so at
# 6-month, -1.0 the SPI and SPEI points were drawn at each other's positions.

grid_df <- read_csv("data_out/tables/revision/FigureS1_index_sensitivity_rev.csv",
                    show_col_types = FALSE) %>%
  mutate(family    = factor(family, levels = c("SPI", "SPEI")),
         scale_lab = factor(scale_lab, levels = c("1-month", "3-month", "6-month", "12-month")),
         thr_lab   = factor(sprintf("%.1f", thr), levels = c("-0.8", "-1.0", "-1.5")),
         sig       = if_else(p < 0.05, "sig", "ns"))

index_cols <- c("SPI" = "#B45309", "SPEI" = "#2171B5")
dodge <- position_dodge(width = 0.5)

pS1 <- ggplot(grid_df, aes(x = thr_lab, y = effect_pct, colour = family, group = family)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo_pct, ymax = hi_pct),
                width = 0.15, linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = interaction(family, sig)), shape = 21, size = 2.8, stroke = 0.9,
             position = dodge) +
  scale_colour_manual(values = index_cols) +
  scale_fill_manual(values = c("SPI.sig" = "#B45309", "SPEI.sig" = "#2171B5",
                               "SPI.ns" = "white", "SPEI.ns" = "white"), guide = "none") +
  facet_wrap(~ scale_lab, nrow = 1) +
  labs(x = "Drought threshold", y = "Difference in burned area (%)", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

print(pS1)
ggsave("figs/revision_final/FigureS1.png", pS1, width = 10, height = 4.2, dpi = 300, bg = "white")
ggsave("figs/revision_final/FigureS1.pdf", pS1, width = 10, height = 4.2, bg = "white")
