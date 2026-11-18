############################################################
# 17_figures_revision_rev.r

############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(sandwich)
  library(lmtest)
  library(quantreg)
})

select <- dplyr::select
filter <- dplyr::filter

project_root <- if (file.exists(file.path("data_out", "incident_panel.csv"))) {
  "."
} else if (file.exists(file.path("..", "data_out", "incident_panel.csv"))) {
  ".."
} else {
  stop("Run this script from the project root or the R directory.")
}

panel_path  <- file.path(project_root, "data_out", "incident_panel.csv")
ledger_path <- file.path(project_root, "data_raw", "kfs_fires.csv")
grid_path   <- file.path(project_root, "data_out", "sensitivity_grid.csv")
table_path  <- file.path(project_root, "data_out", "tables", "revision")
fig_path    <- file.path(project_root, "figs", "revision")

dir.create(table_path, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

## -------------------------------------------------------------------------
## Shared Nature / npj-style settings
## -------------------------------------------------------------------------

font_family <- "Arial"

colours <- list(
  drought     = "#2166AC",
  non_drought = "#BDBDBD",
  buffer      = "#A6CEE3",
  spi         = "#9ECAE1",
  spei        = "#2171B5",
  ink         = "#202020",
  grid        = "#D9D9D9"
)

hazard_colours <- c(
  "Single-hazard"      = colours$non_drought,
  "Compound (drought)" = colours$drought
)

theme_npj <- theme_classic(base_size = 12.5, base_family = font_family) +
  theme(
    axis.title       = element_text(size = 13, colour = colours$ink),
    axis.text        = element_text(size = 11.2, colour = colours$ink),
    axis.line        = element_line(linewidth = 0.55, colour = colours$ink),
    axis.ticks       = element_line(linewidth = 0.5, colour = colours$ink),
    axis.ticks.length = unit(2.2, "mm"),
    strip.background = element_blank(),
    strip.text       = element_text(size = 12.2, face = "bold", colour = colours$ink),
    legend.position  = "top",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 11.2),
    legend.key.height = unit(4.2, "mm"),
    legend.key.width  = unit(7.0, "mm"),
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    plot.caption     = element_blank(),
    plot.tag         = element_text(size = 14, face = "bold", family = font_family),
    plot.margin      = margin(8, 10, 8, 8),
    panel.spacing    = unit(8, "mm")
  )

save_figure <- function(plot, stem, width, height) {
  png_file <- file.path(fig_path, paste0(stem, ".png"))
  pdf_file <- file.path(fig_path, paste0(stem, ".pdf"))
  png_tmp <- tempfile(pattern = paste0(stem, "_"), tmpdir = fig_path, fileext = ".png")
  pdf_tmp <- tempfile(pattern = paste0(stem, "_"), tmpdir = fig_path, fileext = ".pdf")

  ggsave(
    png_tmp, plot, width = width, height = height, units = "in",
    dpi = 600, bg = "white", limitsize = FALSE
  )
  if (!file.copy(png_tmp, png_file, overwrite = TRUE)) {
    stop("Could not write PNG: ", png_file)
  }
  unlink(png_tmp)

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(
    pdf_tmp, plot, width = width, height = height, units = "in",
    device = pdf_device, bg = "white", limitsize = FALSE
  )
  if (!file.copy(pdf_tmp, pdf_file, overwrite = TRUE)) {
    stop("Could not write PDF: ", pdf_file)
  }
  unlink(pdf_tmp)
}

ctrl_vars <- c("rhmin", "wsmax", "tmean", "dry_spell", "wet_days_mo")
fe_terms  <- "factor(season) + factor(sido) + factor(year) + factor(cause)"

panel <- read_csv(panel_path, show_col_types = FALSE) %>%
  mutate(cause = if_else(is.na(cause) | cause == "", "Unknown", cause)) %>%
  add_count(cause, name = "cause_n") %>%
  mutate(cause = if_else(cause_n < 30, "Other", cause)) %>%
  select(-cause_n) %>%
  mutate(
    cluster = paste0(sigungu_key, "_", year),
    comp_all = as.integer(spi3 <= -1),
    hazard = factor(
      if_else(compound == 1, "Compound (drought)", "Single-hazard"),
      levels = names(hazard_colours)
    )
  )

cl_row <- function(fit, cluster, term = "compound") {
  ct <- coeftest(fit, vcov = vcovCL(fit, cluster = cluster))
  tibble(
    estimate = ct[term, 1],
    se = ct[term, 2],
    p_value = ct[term, 4]
  ) %>%
    mutate(
      mult = exp(estimate),
      lo = exp(estimate - 1.96 * se),
      hi = exp(estimate + 1.96 * se)
    )
}

## =========================================================================
## FIGURE 2
## Annual burned area and concentration of burned area under drought
## =========================================================================

thr_fig2    <- -1.0
buffer_thr  <- -0.5
mega_years  <- c(2000, 2022, 2025)
label_years <- mega_years

share_df <- panel %>%
  mutate(comp_all = as.integer(spi3 <= thr_fig2)) %>%
  filter(!is.na(comp_all)) %>%
  group_by(year) %>%
  summarise(n_panel = n(), share = 100 * mean(comp_all), .groups = "drop")

area_df <- read_csv(ledger_path, show_col_types = FALSE) %>%
  group_by(year) %>%
  summarise(n_all = n(), burned = sum(burned_ha, na.rm = TRUE), .groups = "drop")

annual_df <- left_join(area_df, share_df, by = "year") %>%
  filter(is.finite(burned), burned > 0, is.finite(share)) %>%
  arrange(year)

# Map the drought share to the logarithmic burned-area scale. The inverse
# transformation is supplied to sec_axis, so both axes remain interpretable.
area_floor <- 50
area_top   <- 100000
share_top  <- 50
log_span   <- log10(area_top) - log10(area_floor)

share_to_area <- function(x) {
  10^(log10(area_floor) + (x / share_top) * log_span)
}

area_to_share <- function(y) {
  share_top * (log10(y) - log10(area_floor)) / log_span
}

p2a <- ggplot(annual_df, aes(x = year)) +
  geom_rect(
    aes(xmin = year - 0.36, xmax = year + 0.36, ymin = area_floor, ymax = burned),
    fill = "#D9D9D9", colour = "white", linewidth = 0.25
  ) +
  geom_line(
    aes(y = share_to_area(share)),
    colour = colours$drought, linewidth = 1.05
  ) +
  geom_point(
    aes(y = share_to_area(share)),
    colour = colours$drought, fill = "white", shape = 21,
    size = 2.4, stroke = 0.75
  ) +
  geom_text(
    data = filter(annual_df, year %in% label_years),
    aes(y = burned, label = comma(round(burned))),
    vjust = -0.65, size = 3.55, family = font_family,
    fontface = "bold", colour = colours$ink
  ) +
  scale_y_log10(
    name = "Annual burned area (ha; log scale)",
    limits = c(area_floor, area_top),
    breaks = c(100, 1000, 10000, 100000),
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.02)),
    sec.axis = sec_axis(
      transform = area_to_share,
      name = "Drought-associated wildfires (%)",
      breaks = seq(0, share_top, 10),
      labels = label_number(suffix = "%")
    )
  ) +
  scale_x_continuous(
    name = NULL,
    breaks = seq(2000, 2025, 5),
    minor_breaks = NULL,
    expand = expansion(add = 0.7)
  ) +
  theme_npj +
  theme(
    axis.title.y.right = element_text(colour = colours$drought),
    axis.text.y.right  = element_text(colour = colours$drought),
    axis.line.y.right  = element_line(colour = colours$drought),
    axis.ticks.y.right = element_line(colour = colours$drought)
  )

group_levels <- c(
  sprintf("Drought (SPI-3 <= %.1f)", thr_fig2),
  sprintf("Buffer (%.1f < SPI-3 <= %.1f)", thr_fig2, buffer_thr),
  sprintf("Non-drought (SPI-3 > %.1f)", buffer_thr)
)

group_colours <- setNames(
  c(colours$drought, colours$buffer, colours$non_drought),
  group_levels
)

conc_base <- panel %>%
  filter(!is.na(spi3), !is.na(burned_ha)) %>%
  mutate(
    grp = case_when(
      spi3 <= thr_fig2   ~ group_levels[1],
      spi3 <= buffer_thr ~ group_levels[2],
      TRUE               ~ group_levels[3]
    ),
    grp = factor(grp, levels = group_levels)
  )

conc_df <- bind_rows(
  conc_base %>% mutate(period = "All years"),
  conc_base %>%
    filter(!year %in% mega_years) %>%
    mutate(period = "Excluding 2000, 2022 and 2025")
) %>%
  group_by(period, grp) %>%
  summarise(
    n_fires = n(),
    area_ha = sum(burned_ha),
    .groups = "drop_last"
  ) %>%
  mutate(
    `Wildfire incidents` = 100 * n_fires / sum(n_fires),
    `Burned area` = 100 * area_ha / sum(area_ha)
  ) %>%
  ungroup() %>%
  pivot_longer(
    c(`Wildfire incidents`, `Burned area`),
    names_to = "metric", values_to = "pct"
  ) %>%
  mutate(
    period = factor(period, levels = c("All years", "Excluding 2000, 2022 and 2025")),
    metric = factor(metric, levels = c("Burned area", "Wildfire incidents"))
  )

p2b <- ggplot(conc_df, aes(x = pct, y = metric, fill = grp)) +
  geom_col(
    width = 0.66, colour = "white", linewidth = 0.45,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(
      group = grp,
      label = if_else(pct >= 6, sprintf("%.0f%%", pct), ""),
      colour = if_else(grp == group_levels[2], colours$ink, "white")
    ),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3.75, family = font_family, fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_manual(values = group_colours) +
  scale_x_continuous(
    name = "Share of total (%)",
    breaks = seq(0, 100, 20),
    labels = label_number(suffix = "%"),
    expand = expansion(mult = c(0, 0.01))
  ) +
  facet_wrap(~ period, ncol = 1) +
  labs(y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_npj +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = colours$grid, linewidth = 0.35),
    strip.text = element_text(hjust = 0),
    legend.position = "bottom"
  )

figure2 <- p2a / p2b +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

write_csv(annual_df, file.path(table_path, "Figure2a_annual_burned_area_rev.csv"))
write_csv(conc_df, file.path(table_path, "Figure2b_burned_area_concentration_rev.csv"))
save_figure(figure2, "Figure2_annual_burned_area_rev", width = 7.4, height = 8.8)

## =========================================================================
## FIGURE 3
## Adjusted outcome differences by SPI-3 and SPEI-3 category
## =========================================================================

bin_breaks <- c(-Inf, -1.5, -1.0, -0.5, 0.5, Inf)
bin_labels <- c("<= -1.5", "-1.5 to -1.0", "-1.0 to -0.5", "-0.5 to 0.5", "> 0.5")
ref_label  <- "-0.5 to 0.5"
x_order    <- rev(bin_labels)
bin_terms  <- c(
  b_s15 = bin_labels[1], b_s10 = bin_labels[2],
  b_s05 = bin_labels[3], b_wet = bin_labels[5]
)

base_df <- panel %>%
  filter(
    !is.na(spi3), !is.na(rhmin), !is.na(wsmax), !is.na(tmean),
    !is.na(dry_spell), !is.na(wet_days_mo)
  ) %>%
  add_count(date, name = "n_concurrent") %>%
  mutate(
    log_heli = log1p(heli_n),
    log_crew = log1p(crew_n),
    log_equip = log1p(equip_n)
  )

bin_spec <- tribble(
  ~idx,    ~y,           ~extra,                                           ~row, ~label,
  "spi3",  "log_burned", "",                                               "a",  "Burned area",
  "spi3",  "log_dur",    "",                                               "a",  "Suppression duration",
  "spei3", "log_burned", "",                                               "a",  "Burned area",
  "spei3", "log_dur",    "",                                               "a",  "Suppression duration",
  "spi3",  "log_heli",   "+ log_burned",                                   "b",  "Helicopters deployed",
  "spi3",  "log_crew",   "+ log_burned",                                   "b",  "Personnel deployed",
  "spi3",  "log_dur",    "+ log_burned + log_heli + log_crew + log_equip", "b",  "Suppression duration\n(same size and resources)"
)

bin_res <- tibble()

for (sens in c("main", "concurrent")) {
  for (i in seq_len(nrow(bin_spec))) {
    tmp <- base_df %>%
      filter(!is.na(.data[[bin_spec$idx[i]]])) %>%
      mutate(
        bin = cut(
          .data[[bin_spec$idx[i]]], bin_breaks,
          labels = bin_labels, right = TRUE
        ),
        b_s15 = as.integer(bin == bin_labels[1]),
        b_s10 = as.integer(bin == bin_labels[2]),
        b_s05 = as.integer(bin == bin_labels[3]),
        b_wet = as.integer(bin == bin_labels[5])
      ) %>%
      filter(!is.na(.data[[bin_spec$y[i]]]))

    if (bin_spec$row[i] == "b") {
      tmp <- tmp %>% filter(!is.na(heli_n), !is.na(crew_n), !is.na(equip_n))
    }

    f <- paste(
      bin_spec$y[i], "~ b_s15 + b_s10 + b_s05 + b_wet +",
      paste(ctrl_vars, collapse = " + "), bin_spec$extra[i], "+", fe_terms
    )
    if (sens == "concurrent") f <- paste(f, "+ log(n_concurrent)")

    fit <- lm(as.formula(f), data = tmp)

    one <- map_dfr(names(bin_terms), function(k) {
      cl_row(fit, tmp$cluster, term = k) %>%
        mutate(bin = bin_terms[[k]], n_bin = sum(tmp[[k]]))
    })

    one <- bind_rows(
      one,
      tibble(
        estimate = 0, se = NA_real_, p_value = NA_real_, mult = 1,
        lo = 1, hi = 1, bin = ref_label,
        n_bin = sum(tmp$bin == ref_label)
      )
    )

    bin_res <- bind_rows(
      bin_res,
      one %>%
        mutate(
          idx = bin_spec$idx[i], y = bin_spec$y[i], row = bin_spec$row[i],
          label = bin_spec$label[i], sens = sens, n = nobs(fit)
        )
    )
  }
}

bin_res <- bin_res %>%
  mutate(
    effect_pct = 100 * (mult - 1),
    lo_pct = 100 * (lo - 1),
    hi_pct = 100 * (hi - 1),
    index = factor(
      if_else(idx == "spi3", "SPI-3", "SPEI-3"),
      levels = c("SPI-3", "SPEI-3")
    ),
    bin = factor(bin, levels = x_order),
    label = factor(
      label,
      levels = c(
        "Burned area", "Suppression duration", "Helicopters deployed",
        "Personnel deployed", "Suppression duration\n(same size and resources)"
      )
    )
  )

plot_res <- bin_res %>% filter(sens == "main")
index_colours <- c("SPI-3" = colours$spi, "SPEI-3" = colours$spei)
dodge <- position_dodge(width = 0.72)

p3a <- ggplot(
  filter(plot_res, row == "a"),
  aes(
    x = bin, y = effect_pct, ymin = lo_pct, ymax = hi_pct,
    fill = index, group = index
  )
) +
  geom_hline(yintercept = 0, colour = colours$ink, linewidth = 0.55) +
  geom_crossbar(
    width = 0.62, linewidth = 0.45, colour = colours$ink,
    position = dodge
  ) +
  geom_point(
    data = function(d) filter(d, bin == ref_label),
    shape = 23, size = 3.0, stroke = 0.65,
    colour = colours$ink, position = dodge, show.legend = FALSE
  ) +
  scale_fill_manual(values = index_colours) +
  facet_wrap(~ label, nrow = 1) +
  labs(
    x = "SPI-3",
    y = "Difference relative to reference (%)"
  ) +
  theme_npj +
  theme(
    axis.text.x = element_text(size = 10.4, angle = 30, hjust = 1),
    panel.grid.major.y = element_line(colour = colours$grid, linewidth = 0.35)
  )

p3b <- ggplot(
  filter(plot_res, row == "b"),
  aes(x = bin, y = effect_pct, ymin = lo_pct, ymax = hi_pct)
) +
  geom_hline(yintercept = 0, colour = colours$ink, linewidth = 0.55) +
  geom_crossbar(
    width = 0.58, linewidth = 0.45,
    colour = colours$ink, fill = colours$spi
  ) +
  geom_point(
    data = function(d) filter(d, bin == ref_label),
    shape = 23, size = 3.0, stroke = 0.65,
    fill = colours$spi, colour = colours$ink
  ) +
  facet_wrap(~ label, nrow = 1) +
  labs(
    x = "SPI-3",
    y = "Difference relative to reference (%)"
  ) +
  theme_npj +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10.4, angle = 30, hjust = 1),
    panel.grid.major.y = element_line(colour = colours$grid, linewidth = 0.35)
  )

figure3 <- (p3a + plot_spacer() + plot_layout(widths = c(2, 1))) / p3b +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

write_csv(bin_res, file.path(table_path, "Figure3_index_bin_effects_rev.csv"))
save_figure(figure3, "Figure3_drought_categories_rev", width = 7.8, height = 8.7)

## =========================================================================
## FIGURE 4
## Interaction between antecedent drought and concurrent fire weather
## =========================================================================

wx_vars <- c("rhmin", "wsmax")
wx_labels <- c(
  rhmin = "Minimum relative humidity\n(high to low)",
  wsmax = "Maximum wind speed\n(weak to strong)"
)

inter_df <- panel %>%
  filter(
    !is.na(compound), !is.na(log_burned),
    !is.na(rhmin), !is.na(wsmax), !is.na(tmean),
    !is.na(dry_spell), !is.na(wet_days_mo)
  )

strata_df <- tibble()

for (v in wx_vars) {
  tmp <- inter_df %>%
    mutate(
      wx_t = ntile(.data[[v]], 3),
      wx_t = if (v == "rhmin") 4 - wx_t else wx_t,
      comp_t1 = compound * (wx_t == 1),
      comp_t2 = compound * (wx_t == 2),
      comp_t3 = compound * (wx_t == 3)
    )

  f <- as.formula(paste(
    "log_burned ~ comp_t1 + comp_t2 + comp_t3 + factor(wx_t) +",
    paste(ctrl_vars, collapse = " + "), "+", fe_terms
  ))
  fit <- lm(f, data = tmp)

  for (k in 1:3) {
    strata_df <- bind_rows(
      strata_df,
      cl_row(fit, tmp$cluster, term = paste0("comp_t", k)) %>%
        mutate(
          wx = v, tercile = k,
          n_drought = sum(tmp$compound == 1 & tmp$wx_t == k),
          n_all = sum(tmp$wx_t == k)
        )
    )
  }
}

strata_df <- strata_df %>%
  mutate(
    panel_lab = factor(wx_labels[wx], levels = wx_labels),
    condition = factor(
      recode(as.character(tercile),
             `1` = "Favourable", `2` = "Intermediate", `3` = "Adverse"),
      levels = c("Favourable", "Intermediate", "Adverse")
    ),
    effect_pct = 100 * (mult - 1),
    lo_pct = 100 * (lo - 1),
    hi_pct = 100 * (hi - 1),
    sig = if_else(p_value < 0.05, "p < 0.05", "p >= 0.05")
  )

p4a <- ggplot(
  strata_df,
  aes(x = condition, y = effect_pct, group = panel_lab)
) +
  geom_hline(yintercept = 0, colour = colours$ink, linewidth = 0.55) +
  geom_errorbar(
    aes(ymin = lo_pct, ymax = hi_pct),
    width = 0.10, linewidth = 0.75, colour = colours$drought
  ) +
  geom_line(colour = colours$drought, linewidth = 0.9) +
  geom_point(
    aes(fill = sig), shape = 21, size = 3.6, stroke = 0.95,
    colour = colours$drought
  ) +
  scale_fill_manual(
    values = c("p < 0.05" = colours$drought, "p >= 0.05" = "white")
  ) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  facet_wrap(~ panel_lab, nrow = 1) +
  labs(
    x = "Fire weather condition",
    y = "Difference in burned area (%)"
  ) +
  theme_npj +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10.8),
    panel.grid.major.y = element_line(colour = colours$grid, linewidth = 0.35)
  )

large_df <- tibble()

for (v in wx_vars) {
  tmp <- inter_df %>%
    filter(!is.na(large_fire)) %>%
    mutate(
      adverse = if (v == "rhmin") {
        .data[[v]] < median(.data[[v]])
      } else {
        .data[[v]] > median(.data[[v]])
      }
    )

  large_df <- bind_rows(
    large_df,
    tmp %>%
      group_by(hazard, adverse) %>%
      summarise(n = n(), n_large = sum(large_fire), .groups = "drop") %>%
      mutate(wx = v)
  )
}

large_df <- large_df %>%
  mutate(
    pct = 100 * n_large / n,
    panel_lab = factor(wx_labels[wx], levels = wx_labels),
    wx_cond = factor(
      if_else(adverse, "Adverse", "Favourable"),
      levels = c("Favourable", "Adverse")
    )
  )

p4b <- ggplot(large_df, aes(x = wx_cond, y = pct, fill = hazard)) +
  geom_col(
    position = position_dodge(width = 0.75), width = 0.66,
    colour = colours$ink, linewidth = 0.35
  ) +
  geom_text(
    aes(label = sprintf("%.2f", pct)),
    position = position_dodge(width = 0.75),
    vjust = -0.55, size = 3.65, family = font_family, colour = colours$ink
  ) +
  scale_fill_manual(
    values = hazard_colours,
    labels = c(
      "Single-hazard" = "Non-drought",
      "Compound (drought)" = "Drought-associated"
    )
  ) +
  scale_y_continuous(
    name = "Fires >= 100 ha (%)",
    expand = expansion(mult = c(0, 0.15))
  ) +
  facet_wrap(~ panel_lab, nrow = 1) +
  labs(x = "Fire weather condition") +
  theme_npj +
  theme(
    panel.grid.major.y = element_line(colour = colours$grid, linewidth = 0.35)
  )

inter_test <- tibble()

for (v in wx_vars) {
  for (y in c("log_burned", "log_dur", "large_fire")) {
    tmp <- inter_df %>% filter(!is.na(.data[[y]]))
    if (y == "large_fire") tmp <- tmp %>% mutate(large_fire = 100 * large_fire)

    f <- as.formula(paste(
      y, "~ compound *", v, "+",
      paste(setdiff(ctrl_vars, v), collapse = " + "), "+", fe_terms
    ))
    fit <- lm(f, data = tmp)

    inter_test <- bind_rows(
      inter_test,
      cl_row(fit, tmp$cluster, term = paste0("compound:", v)) %>%
        select(estimate, se, p_value) %>%
        mutate(outcome = y, wx = v, n = nobs(fit))
    )
  }
}

figure4 <- p4a / p4b +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

write_csv(strata_df, file.path(table_path, "Figure4a_drought_effect_by_fireweather_rev.csv"))
write_csv(large_df, file.path(table_path, "Figure4b_large_fire_by_fireweather_rev.csv"))
write_csv(inter_test, file.path(table_path, "Table_interaction_tests_rev.csv"))
save_figure(figure4, "Figure4_fire_weather_interaction_rev", width = 7.4, height = 8.0)

## =========================================================================
## FIGURE 5
## Intuitive distribution view plus adjusted quantile ratios
## =========================================================================

taus <- c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
n_boot <- 200

qr_panel <- panel %>%
  filter(
    !is.na(compound), is.finite(log_burned),
    if_all(all_of(ctrl_vars), is.finite)
  )

qr_fml <- as.formula(paste(
  "log_burned ~ compound +",
  paste(ctrl_vars, collapse = " + ")
))

quantile_cache <- file.path(project_root, "data_out", "tables", "Fig4_quantile.csv")
recompute_quantiles <- identical(
  tolower(Sys.getenv("RECOMPUTE_QUANTILES", "false")),
  "true"
)

if (file.exists(quantile_cache) && !recompute_quantiles) {
  
  qr_df <- read_csv(quantile_cache, show_col_types = FALSE) %>%
    select(tau, estimate, se)
} else {
  set.seed(20261101)
  qr_df <- map_dfr(taus, function(tau) {
    fit <- rq(qr_fml, tau = tau, data = qr_panel)
    s <- summary(fit, se = "boot", R = n_boot)$coefficients

    tibble(
      tau = tau,
      estimate = s["compound", 1],
      se = s["compound", 2]
    )
  })
}

qr_df <- qr_df %>%
  mutate(
    lo = estimate - 1.96 * se,
    hi = estimate + 1.96 * se,
    ratio = exp(estimate),
    lo_ratio = exp(lo),
    hi_ratio = exp(hi),
    effect_pct = 100 * (ratio - 1),
    lo_pct = 100 * (lo_ratio - 1),
    hi_pct = 100 * (hi_ratio - 1),
    quantile_label = factor(
      case_when(
        tau == 0.10 ~ "10th",
        tau == 0.25 ~ "25th",
        tau == 0.50 ~ "Median",
        tau == 0.75 ~ "75th",
        tau == 0.90 ~ "90th",
        tau == 0.95 ~ "95th",
        tau == 0.99 ~ "99th"
      ),
      levels = c("10th", "25th", "Median", "75th", "90th", "95th", "99th")
    )
  )

dist_panel <- qr_panel %>%
  transmute(
    group = factor(
      if_else(compound == 1, "Drought-associated", "Non-drought"),
      levels = c("Drought-associated", "Non-drought")
    ),
    burned_ha = pmax(burned_ha, 0.01),
    x = log10(pmax(burned_ha, 0.01))
  )

group_positions <- c("Drought-associated" = 1.0, "Non-drought" = 2.35)
dist_colours <- c(
  "Drought-associated" = colours$drought,
  "Non-drought" = "#B8C0CC"
)

x_limits <- c(-2, 5)
density_df <- map_dfr(levels(dist_panel$group), function(g) {
  x <- filter(dist_panel, group == g)$x
  den <- density(x, from = x_limits[1], to = x_limits[2], n = 512, adjust = 1.05)
  tibble(
    group = factor(g, levels = levels(dist_panel$group)),
    x = den$x,
    density = 0.56 * den$y / max(den$y),
    baseline = unname(group_positions[g])
  )
})

box_df <- dist_panel %>%
  group_by(group) %>%
  summarise(
    n = n(),
    q1 = quantile(x, 0.25),
    median = median(x),
    q3 = quantile(x, 0.75),
    iqr = IQR(x),
    lower = max(min(x), q1 - 1.5 * iqr),
    upper = min(max(x), q3 + 1.5 * iqr),
    .groups = "drop"
  ) %>%
  mutate(
    baseline = unname(group_positions[as.character(group)]),
    y_box = baseline - 0.20
  )

set.seed(20261101)
point_df <- dist_panel %>%
  group_by(group) %>%
  group_modify(~ slice_sample(.x, n = min(260, nrow(.x)))) %>%
  ungroup() %>%
  mutate(
    baseline = unname(group_positions[as.character(group)]),
    y = baseline - 0.43 + runif(n(), -0.065, 0.065)
  )

group_counts <- box_df %>%
  mutate(label = paste0(as.character(group), "\n(n = ", comma(n), ")"))

log_break_values <- c(0.01, 0.1, 1, 10, 100, 1000, 10000, 100000)

p5a <- ggplot() +
  geom_ribbon(
    data = density_df,
    aes(x = x, ymin = baseline, ymax = baseline + density, fill = group),
    alpha = 0.62, colour = NA
  ) +
  geom_line(
    data = density_df,
    aes(x = x, y = baseline + density, colour = group),
    linewidth = 0.75
  ) +
  geom_segment(
    data = density_df %>% distinct(group, baseline),
    aes(x = x_limits[1], xend = x_limits[2], y = baseline, yend = baseline),
    colour = colours$ink, linewidth = 0.4
  ) +
  geom_point(
    data = point_df,
    aes(x = x, y = y, colour = group),
    size = 0.95, alpha = 0.42
  ) +
  geom_segment(
    data = box_df,
    aes(x = lower, xend = upper, y = y_box, yend = y_box),
    colour = colours$ink, linewidth = 0.65
  ) +
  geom_segment(
    data = box_df,
    aes(x = lower, xend = lower, y = y_box - 0.065, yend = y_box + 0.065),
    colour = colours$ink, linewidth = 0.65
  ) +
  geom_segment(
    data = box_df,
    aes(x = upper, xend = upper, y = y_box - 0.065, yend = y_box + 0.065),
    colour = colours$ink, linewidth = 0.65
  ) +
  geom_rect(
    data = box_df,
    aes(xmin = q1, xmax = q3, ymin = y_box - 0.11, ymax = y_box + 0.11, fill = group),
    colour = colours$ink, linewidth = 0.65
  ) +
  geom_segment(
    data = box_df,
    aes(x = median, xend = median, y = y_box - 0.11, yend = y_box + 0.11),
    colour = colours$ink, linewidth = 1.0
  ) +
  scale_fill_manual(values = dist_colours, guide = "none") +
  scale_colour_manual(values = dist_colours, guide = "none") +
  scale_x_continuous(
    name = "Burned area (ha; log scale)",
    limits = x_limits,
    breaks = log10(log_break_values),
    labels = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000", "100,000"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name = NULL,
    limits = c(0.38, 3.20),
    breaks = group_positions,
    labels = group_counts$label[match(names(group_positions), as.character(group_counts$group))],
    expand = expansion(mult = c(0, 0))
  ) +
  theme_npj +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11.5, hjust = 1, lineheight = 1.05),
    panel.grid.major.x = element_line(colour = "#ECECEC", linewidth = 0.35)
  )

p5b <- ggplot(qr_df, aes(x = quantile_label, y = ratio, group = 1)) +
  geom_hline(yintercept = 1, colour = colours$ink, linewidth = 0.55) +
  geom_line(colour = colours$spei, linewidth = 0.9) +
  geom_errorbar(
    aes(ymin = lo_ratio, ymax = hi_ratio),
    width = 0.12, linewidth = 0.7, colour = colours$spei
  ) +
  geom_point(
    aes(fill = tau >= 0.90), shape = 21, size = 3.6,
    stroke = 0.85, colour = colours$spei
  ) +
  geom_text(
    aes(label = paste0("\u00d7", sprintf("%.2f", ratio))),
    nudge_y = 0.17, family = font_family, size = 3.35,
    colour = colours$ink, vjust = 0
  ) +
  scale_fill_manual(
    values = c(`FALSE` = colours$spi, `TRUE` = colours$spei),
    guide = "none"
  ) +
  scale_y_continuous(
    name = "Drought / non-drought ratio",
    labels = function(x) paste0("x", sprintf("%.1f", x)),
    expand = expansion(mult = c(0.06, 0.12))
  ) +
  labs(x = "Burned-area quantile") +
  theme_npj +
  theme(
    panel.grid.major.y = element_line(colour = colours$grid, linewidth = 0.35)
  )

figure5 <- p5a / p5b +
  plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

write_csv(qr_df, file.path(table_path, "Figure5_quantile_effects_rev.csv"))
save_figure(figure5, "Figure5_quantile_distribution_rev", width = 7.6, height = 8.8)

## =========================================================================
## FIGURE S1

## =========================================================================

if (file.exists(grid_path)) {
  grid_df <- read_csv(grid_path, show_col_types = FALSE) %>%
    mutate(
      family = factor(
        if_else(str_detect(index, "spei"), "SPEI", "SPI"),
        levels = c("SPI", "SPEI")
      ),
      scale = as.integer(str_extract(index, "[0-9]+")),
      scale_lab = factor(
        paste0(scale, "-month"),
        levels = c("1-month", "3-month", "6-month", "12-month")
      ),
      thr_lab = factor(sprintf("%.1f", thr), levels = c("-0.8", "-1.0", "-1.5")),
      effect_pct = 100 * (exp(beta) - 1),
      lo_pct = 100 * (exp(beta - 1.96 * SE) - 1),
      hi_pct = 100 * (exp(beta + 1.96 * SE) - 1),
      sig = if_else(p < 0.05, "sig", "ns")
    )

  index_colours_s1 <- c("SPI" = colours$spi, "SPEI" = colours$spei)

  pS1 <- ggplot(grid_df, aes(x = thr_lab, y = effect_pct, colour = family)) +
    geom_hline(yintercept = 0, colour = colours$ink, linewidth = 0.55) +
    geom_errorbar(
      aes(ymin = lo_pct, ymax = hi_pct),
      width = 0.12, linewidth = 0.7,
      position = position_dodge(width = 0.5)
    ) +
    geom_point(
      aes(fill = interaction(family, sig)),
      shape = 21, size = 3.4, stroke = 0.9,
      position = position_dodge(width = 0.5)
    ) +
    scale_colour_manual(values = index_colours_s1) +
    scale_fill_manual(
      values = c(
        "SPI.sig" = colours$spi, "SPEI.sig" = colours$spei,
        "SPI.ns" = "white", "SPEI.ns" = "white"
      ),
      guide = "none"
    ) +
    facet_wrap(~ scale_lab, nrow = 1) +
    labs(
      x = "Drought threshold",
      y = "Difference in burned area (%)"
    ) +
    theme_npj +
    theme(
      panel.grid.major.y = element_line(colour = colours$grid, linewidth = 0.35)
    )

  write_csv(grid_df, file.path(table_path, "FigureS1_index_sensitivity_rev.csv"))
  save_figure(pS1, "FigureS1_index_sensitivity_rev", width = 7.8, height = 3.8)
}

message("Revision figures written to: ", normalizePath(fig_path))
