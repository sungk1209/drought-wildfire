# ============================================================================
# 16_make_tables_figures.R — 원고용 Table / Figure 생성

#   source("R/16_make_tables_figures.R")
#   make_all_outputs()
# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages({
  library(data.table); library(sandwich); library(lmtest); library(quantreg)
})

TDIR <- "data_out/tables"; FDIR <- "figs"
CTRL <- c("rhmin", "wsmax", "tmean", "dry_spell", "wet_days_mo")
FE   <- "factor(season) + factor(sido) + factor(year) + factor(cause)"
fml  <- function(y, rhs) as.formula(paste(y, "~", rhs))
cl   <- function(m, g) coeftest(m, vcov = vcovCL(m, cluster = g))

# 흑백 인쇄에서도 구분되는 팔레트
COL <- c(control = "#9CA3AF", compound = "#B45309",
         line = "#1F2937", accent = "#0F766E", light = "#E5E7EB")

.load <- function() {
  d <- fread("data_out/incident_panel.csv")
  d <- d[!is.na(compound)]
  need <- intersect(c("log_burned","compound", CTRL), names(d))
  d <- d[complete.cases(d[, ..need])]
  numv <- need[vapply(need, function(v) is.numeric(d[[v]]), logical(1))]
  d <- d[Reduce(`&`, lapply(numv, function(v) is.finite(d[[v]])))]
  d[, cause := fifelse(is.na(cause) | !nzchar(cause), "미상", cause)]
  tb <- d[, .N, by = cause]; d[cause %in% tb[N < 30]$cause, cause := "기타"]
  d[, cluster := paste0(sigungu_key, "_", year)]
  d[]
}

# PDF + PNG 동시 출력. expr 은 호출자 환경에서 평가해야 한다(env 인자).
.plot_both <- function(name, expr, w = 7, h = 5, env = parent.frame()) {
  for (dev_i in 1:2) {
    if (dev_i == 1) pdf(file.path(FDIR, paste0(name, ".pdf")), width = w, height = h,
                        useDingbats = FALSE)
    else png(file.path(FDIR, paste0(name, ".png")), width = round(w * 300),
             height = round(h * 300), res = 300)
    par(mar = c(4.2, 4.4, 2.2, 1.2), mgp = c(2.6, 0.7, 0), las = 1,
        cex.axis = 0.9, cex.lab = 1.0)
    ok <- try(eval(expr, envir = env), silent = TRUE)
    dev.off()
    if (inherits(ok, "try-error")) {
      message("  !! ", name, " 실패: ", conditionMessage(attr(ok, "condition")))
      return(invisible(FALSE))
    }
  }
  message("  figs/", name, ".pdf / .png")
  invisible(TRUE)
}

# ===========================================================================
# TABLE 1 — 기술통계: 복합 vs 단독 
# ===========================================================================
table1 <- function(d) {
  f <- function(x, lab) {
    a <- x[compound == 1]; b <- x[compound == 0]
    data.table(sample = lab,
               n_control = nrow(b), n_compound = nrow(a),
               mean_control = round(mean(b$burned_ha), 2),
               mean_compound = round(mean(a$burned_ha), 2),
               ratio_mean = round(mean(a$burned_ha) / mean(b$burned_ha), 2),
               med_control = round(median(b$burned_ha), 2),
               med_compound = round(median(a$burned_ha), 2),
               p90_control = round(quantile(b$burned_ha, .9), 2),
               p90_compound = round(quantile(a$burned_ha, .9), 2),
               pct100_control = round(100 * mean(b$large_fire), 2),
               pct100_compound = round(100 * mean(a$large_fire), 2),
               ratio_pct100 = round(mean(a$large_fire) / mean(b$large_fire), 2),
               dur_med_control = round(median(b$dur_hr, na.rm = TRUE), 2),
               dur_med_compound = round(median(a$dur_hr, na.rm = TRUE), 2))
  }
  t <- rbind(
    f(d, "All fires"),
    f(d[burned_ha <= quantile(burned_ha, .999)], "Excl. top 0.1%"),
    f(d[burned_ha <= quantile(burned_ha, .99)],  "Excl. top 1%"))
  fwrite(t, file.path(TDIR, "Table1_descriptives.csv"))
  message("  Table1_descriptives.csv")
  t
}

# ===========================================================================
# TABLE 2 — 주회귀: 피해면적 + 진화시간 
# ===========================================================================
table2 <- function(d) {
  specs <- list(
    "M1 Unadjusted"            = "compound",
    "M2 + fire weather"        = paste("compound", paste(CTRL, collapse = " + "), sep = " + "),
    "M3 + season/region/year/cause" =
      paste("compound", paste(CTRL, collapse = " + "), FE, sep = " + "))
  one <- function(y, lab) rbindlist(lapply(names(specs), function(nm) {
    x <- d[is.finite(get(y))]
    ct <- cl(lm(fml(y, specs[[nm]]), x), x$cluster)
    data.table(outcome = lab, model = nm,
               beta = round(ct["compound",1], 4), SE = round(ct["compound",2], 4),
               p = signif(ct["compound",4], 3),
               multiplier = round(exp(ct["compound",1]), 3),
               ci_lo = round(exp(ct["compound",1] - 1.96*ct["compound",2]), 3),
               ci_hi = round(exp(ct["compound",1] + 1.96*ct["compound",2]), 3),
               n = nrow(x))
  }))
  t <- rbind(one("log_burned", "log(burned area)"),
             one("log_dur",    "log(containment duration)"))
  fwrite(t, file.path(TDIR, "Table2_main_regression.csv"))
  message("  Table2_main_regression.csv")
  t
}

# ===========================================================================
# TABLE 3 — 진화자원 통제 
# ===========================================================================
table3 <- function(d) {
  x <- d[is.finite(heli_n) & is.finite(log_dur)]
  if (nrow(x) < 200) { message("  (Table3 생략: 진화자원 표본 부족)"); return(NULL) }
  x[, `:=`(l_heli = log1p(heli_n), l_crew = log1p(crew_n), l_equip = log1p(equip_n))]
  cov <- paste(paste(CTRL, collapse = " + "), "+", FE)
  specs <- list(
    "R1 no resource control" = paste("compound +", cov),
    "R2 + helicopters"       = paste("compound + l_heli +", cov),
    "R3 + heli/crew/equip"   = paste("compound + l_heli + l_crew + l_equip +", cov),
    "R4 + burned area"       = paste("compound + l_heli + l_crew + l_equip + log_burned +", cov))
  t <- rbindlist(lapply(names(specs), function(nm) {
    ct <- cl(lm(fml("log_dur", specs[[nm]]), x), x$cluster)
    data.table(model = nm, beta = round(ct["compound",1], 4),
               SE = round(ct["compound",2], 4), p = signif(ct["compound",4], 3),
               multiplier = round(exp(ct["compound",1]), 3), n = nrow(x))
  }))
  # 자원 투입량 자체에 대한 효과 (규모 통제)
  t2 <- rbindlist(lapply(c("l_heli","l_crew","l_equip"), function(y) {
    ct <- cl(lm(fml(y, paste("compound + log_burned +", cov)), x), x$cluster)
    data.table(outcome = y, beta = round(ct["compound",1], 4),
               SE = round(ct["compound",2], 4), p = signif(ct["compound",4], 3),
               pct_change = round(100 * (exp(ct["compound",1]) - 1), 1))
  }))
  fwrite(t,  file.path(TDIR, "Table3a_resource_controlled_duration.csv"))
  fwrite(t2, file.path(TDIR, "Table3b_resource_deployment.csv"))
  message("  Table3a_resource_controlled_duration.csv / Table3b_resource_deployment.csv")
  list(dur = t, dep = t2)
}

# ===========================================================================
# TABLE 4 — 시간규모: SPI-3 vs SPI-12
# ===========================================================================
table4 <- function(d) {
  ix <- intersect(c("spi1","spi3","spi6","spi12"), names(d))
  outs <- c(log_burned = "log(burned area)", log_dur = "log(duration)")
  # (a) 개별 투입
  ta <- rbindlist(lapply(names(outs), function(y) rbindlist(lapply(ix, function(v) {
    x <- d[is.finite(get(v)) & is.finite(get(y))]
    x[, z := scale(get(v))[,1]]
    ct <- cl(lm(fml(y, paste("z +", paste(CTRL, collapse = " + "), "+", FE)), x), x$cluster)
    data.table(outcome = outs[[y]], index = toupper(v),
               beta_std = round(ct["z",1], 4), SE = round(ct["z",2], 4),
               p = signif(ct["z",4], 3))
  }))))
  # (b) 동시 투입
  tb <- rbindlist(lapply(names(outs), function(y) {
    x <- d[is.finite(spi3) & is.finite(spi12) & is.finite(get(y))]
    x[, `:=`(z3 = scale(spi3)[,1], z12 = scale(spi12)[,1])]
    ct <- cl(lm(fml(y, paste("z3 + z12 +", paste(CTRL, collapse = " + "), "+", FE)), x),
             x$cluster)
    data.table(outcome = outs[[y]],
               spi3_beta = round(ct["z3",1], 4),  spi3_p = signif(ct["z3",4], 3),
               spi12_beta = round(ct["z12",1], 4), spi12_p = signif(ct["z12",4], 3))
  }))
  # (c) 조합별
  tc <- NULL
  if ("drought_type" %in% names(d)) {
    x <- d[!is.na(drought_type)]
    x[, drought_type := relevel(factor(drought_type), ref = "none")]
    desc <- x[, .(n = .N, med_ha = round(median(burned_ha), 2),
                  med_dur = round(median(dur_hr, na.rm = TRUE), 2),
                  pct100 = round(100 * mean(large_fire), 2)), by = drought_type]
    tc <- rbindlist(lapply(names(outs), function(y) {
      x2 <- x[is.finite(get(y))]
      ct <- cl(lm(fml(y, paste("drought_type +", paste(CTRL, collapse = " + "),
                               "+", FE)), x2), x2$cluster)
      k <- grep("^drought_type", rownames(ct), value = TRUE)
      data.table(outcome = outs[[y]], group = sub("drought_type", "", k),
                 beta = round(ct[k,1], 4), SE = round(ct[k,2], 4),
                 p = signif(ct[k,4], 3), multiplier = round(exp(ct[k,1]), 3))
    }))
    tc <- merge(tc, desc, by.x = "group", by.y = "drought_type", all.x = TRUE)
  }
  fwrite(ta, file.path(TDIR, "Table4a_timescale_separate.csv"))
  fwrite(tb, file.path(TDIR, "Table4b_timescale_joint.csv"))
  if (!is.null(tc)) fwrite(tc, file.path(TDIR, "Table4c_timescale_groups.csv"))
  message("  Table4a/b/c_timescale*.csv")
  list(sep = ta, joint = tb, grp = tc)
}


tables_SI <- function(d) {

  s1 <- dcast(d[, .N, by = .(cause, compound)], cause ~ compound,
              value.var = "N", fill = 0)
  setnames(s1, c("0","1"), c("control","compound"), skip_absent = TRUE)
  s1[, `:=`(pct_control = round(100*control/sum(control), 1),
            pct_compound = round(100*compound/sum(compound), 1))]
  s1[, diff_pp := round(pct_compound - pct_control, 1)]
  setorder(s1, -pct_compound)
  ct <- suppressWarnings(chisq.test(table(d$cause, d$compound)))
  s1 <- rbind(s1, data.table(cause = sprintf("chi2=%.1f, df=%d, p=%.3g",
              ct$statistic, ct$parameter, ct$p.value)), fill = TRUE)
  fwrite(s1, file.path(TDIR, "TableS1_cause_composition.csv"))

  # S2 균형표
  vs <- intersect(c(CTRL, "precip_prev7d", "stn_km"), names(d))
  s2 <- rbindlist(lapply(vs, function(v) {
    a <- d[compound==1][[v]]; b <- d[compound==0][[v]]
    sdp <- sqrt((var(a, na.rm=TRUE) + var(b, na.rm=TRUE))/2)
    data.table(variable = v, compound = round(mean(a, na.rm=TRUE), 3),
               control = round(mean(b, na.rm=TRUE), 3),
               SMD = round((mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE))/sdp, 3))
  }))
  fwrite(s2, file.path(TDIR, "TableS2_balance.csv"))

  # S3 민감도
  if (file.exists("data_out/sensitivity_grid.csv"))
    file.copy("data_out/sensitivity_grid.csv",
              file.path(TDIR, "TableS3_sensitivity_grid.csv"), overwrite = TRUE)

  # S4 거리
  if (file.exists("data_out/sigungu_station_map.csv")) {
    s4 <- fread("data_out/sigungu_station_map.csv")
    fwrite(s4, file.path(TDIR, "TableS4_station_distance.csv"))
  }
  message("  TableS1-S4_*.csv")
  invisible(NULL)
}

# ===========================================================================
# FIGURE 1 
# ===========================================================================
fig1_quantile <- function(d, R = 200) {
  rhs <- paste("compound", paste(CTRL, collapse = " + "), sep = " + ")
  taus <- c(.10,.25,.50,.75,.90,.95,.99)
  res <- rbindlist(lapply(taus, function(q) {
    s <- summary(rq(fml("log_burned", rhs), tau = q, data = d),
                 se = "boot", R = R)$coefficients
    data.table(tau = q, b = s["compound",1], se = s["compound",2])
  }))
  res[, `:=`(lo = b - 1.96*se, hi = b + 1.96*se)]
  fwrite(res, file.path(TDIR, "Fig1_quantile_data.csv"))

  .plot_both("Fig1_quantile_effects", quote({
    ylim <- range(c(res$lo, res$hi, 0)) + c(-0.05, 0.08)
    plot(res$tau, res$b, type = "n", ylim = ylim, xlim = c(0.05, 1.02),
         xlab = expression("Quantile of burned area ("*tau*")"),
         ylab = "Effect of compound drought on log(burned area)",
         main = "Drought amplification concentrates in the upper tail")
    abline(h = 0, col = COL[["line"]], lty = 2)
    polygon(c(res$tau, rev(res$tau)), c(res$lo, rev(res$hi)),
            col = rgb(180/255, 83/255, 9/255, 0.18), border = NA)
    lines(res$tau, res$b, col = COL[["compound"]], lwd = 2.4)
    points(res$tau, res$b, pch = 19, col = COL[["compound"]], cex = 1.1)
    # 배수 축
    mm <- round(exp(res$b), 2)
    text(res$tau, res$hi + 0.03, sprintf("%.2fx", mm), cex = 0.72,
         col = COL[["line"]])
    legend("topleft", bty = "n", cex = 0.85,
           legend = c("Point estimate", "95% CI (bootstrap)"),
           col = c(COL[["compound"]], NA), lwd = c(2.4, NA),
           fill = c(NA, rgb(180/255, 83/255, 9/255, 0.18)), border = NA)
  }), w = 7.2, h = 5)
  res
}

# ===========================================================================
# FIGURE 2 — SPI3 
# ===========================================================================
fig2_threshold <- function(d) {
  x <- d[is.finite(spi3)]
  brk <- c(-Inf, -2, -1.5, -1, -0.5, 0, Inf)
  lab <- c("< -2", "-2 to -1.5", "-1.5 to -1", "-1 to -0.5", "-0.5 to 0", "> 0")
  x[, bin := cut(spi3, brk, labels = lab)]
  s <- x[, .(n = .N, p90 = quantile(burned_ha, .9),
             pct100 = 100*mean(large_fire),
             med_dur = median(dur_hr, na.rm = TRUE)), by = bin][order(bin)]
  fwrite(s, file.path(TDIR, "Fig2_threshold_data.csv"))

  .plot_both("Fig2_spi3_bins", quote({
    par(mfrow = c(1, 2), mar = c(6.5, 4.4, 2.4, 1))
    bp <- barplot(s$pct100, names.arg = s$bin, las = 2,
                  col = COL[["compound"]], border = NA,
                  ylab = "Fires >= 100 ha (%)",
                  main = "(a) Large-fire probability")
    text(bp, s$pct100, sprintf("%.2f", s$pct100), pos = 3, cex = 0.75, xpd = NA)
    mtext(expression(SPI-3~"(preceding month)"), side = 1, line = 5, cex = 0.85)
    bp2 <- barplot(s$p90, names.arg = s$bin, las = 2,
                   col = COL[["accent"]], border = NA,
                   ylab = "90th percentile burned area (ha)",
                   main = "(b) Upper-tail burned area")
    text(bp2, s$p90, sprintf("%.2f", s$p90), pos = 3, cex = 0.75, xpd = NA)
    mtext(expression(SPI-3~"(preceding month)"), side = 1, line = 5, cex = 0.85)
  }), w = 9, h = 5)
  s
}

# ===========================================================================
# FIGURE 3 
# ===========================================================================
fig3_synergy <- function(d) {
  x <- copy(d)[, hi_ws := as.integer(wsmax > median(wsmax))]
  cel <- x[, .(m = mean(burned_ha), p90 = quantile(burned_ha, .9),
               pct100 = 100*mean(large_fire), n = .N), by = .(compound, hi_ws)]
  setorder(cel, compound, hi_ws)
  fwrite(cel, file.path(TDIR, "Fig3_synergy_data.csv"))
  g <- function(cc, ww, col) cel[compound == cc & hi_ws == ww][[col]]
  base <- g(0,0,"pct100"); dO <- g(1,0,"pct100"); wO <- g(0,1,"pct100")
  obs  <- g(1,1,"pct100"); expd <- (dO - base) + (wO - base) + base

  .plot_both("Fig3_nonadditivity", quote({
    v <- c(base, dO, wO, expd, obs)
    nm <- c("Neither", "Drought\nonly", "High wind\nonly",
            "Additive\nexpectation", "Observed\njoint")
    bp <- barplot(v, names.arg = nm, border = NA,
                  col = c(COL[["control"]], COL[["compound"]], COL[["accent"]],
                          COL[["light"]], COL[["compound"]]),
                  ylab = "Fires >= 100 ha (%)",
                  main = "Non-additivity of drought and high wind",
                  ylim = c(0, max(v)*1.25))
    text(bp, v, sprintf("%.2f", v), pos = 3, cex = 0.8, xpd = NA)
    segments(bp[4], expd, bp[5], expd, lty = 2, col = COL[["line"]])
    arrows(bp[5], expd, bp[5], obs, length = .07, code = 3, col = COL[["line"]])
    text(bp[5], (expd+obs)/2, sprintf("  +%.2f pp", obs-expd), pos = 4,
         cex = 0.8, xpd = NA)
  }), w = 7.5, h = 5)
  cel
}

# ===========================================================================
# FIGURE 4 
# ===========================================================================
fig4_resource <- function(d) {
  t <- table3(d)
  if (is.null(t)) return(NULL)
  r <- t$dur
  .plot_both("Fig4_resource_controlled", quote({
    par(mar = c(4.4, 12, 2.4, 1.5))
    y <- rev(seq_len(nrow(r)))
    lo <- r$beta - 1.96*r$SE; hi <- r$beta + 1.96*r$SE
    plot(r$beta, y, xlim = range(c(lo, hi, 0)) + c(-0.03, 0.03),
         ylim = c(0.5, nrow(r)+0.5), yaxt = "n", pch = 19, cex = 1.2,
         col = COL[["compound"]], xlab = "Effect on log(containment duration)",
         ylab = "", main = "Duration remains longer after controlling resources")
    abline(v = 0, lty = 2, col = COL[["line"]])
    segments(lo, y, hi, y, col = COL[["compound"]], lwd = 2.2)
    axis(2, at = y, labels = r$model, las = 1, cex.axis = 0.85, tick = FALSE)
    text(hi, y, sprintf("  %.2fx (p=%.3g)", r$multiplier, r$p),
         pos = 4, cex = 0.72, xpd = NA)
  }), w = 8.5, h = 4.2)
  r
}

# ===========================================================================
# FIGURE 5 
# ===========================================================================
fig5_distance <- function() {
  p <- "data_out/sigungu_station_map.csv"
  if (!file.exists(p)) return(NULL)
  m <- fread(p)
  .plot_both("FigS1_station_distance", quote({
    hist(m$stn_km, breaks = 24, col = COL[["accent"]], border = "white",
         xlab = "Distance from municipality centroid to ASOS station (km)",
         ylab = "Number of municipalities",
         main = "Spatial matching quality")
    abline(v = median(m$stn_km), lty = 2, lwd = 2, col = COL[["compound"]])
    legend("topright", bty = "n", cex = 0.85,
           legend = sprintf("median = %.1f km\n95th pct = %.1f km\nmax = %.1f km",
                            median(m$stn_km), quantile(m$stn_km, .95),
                            max(m$stn_km)))
  }), w = 7, h = 4.6)
  invisible(m)
}

# ===========================================================================
make_all_outputs <- function() {
  dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(FDIR, showWarnings = FALSE)
  d <- .load()
  message("분석 표본: ", nrow(d), "건 / compound ", sum(d$compound))
  message("\n[Tables]")
  table1(d); table2(d); table4(d); tables_SI(d)
  message("\n[Figures] 본문")
  fig1_dose_response(d)     # Fig 1: 용량-반응 (강도 x 기간)
  fig2_cascade(d)           # Fig 2: 인과 사슬 총효과
  fig1_quantile(d)          # Fig 3: 분위회귀
  fig3_synergy(d)           # Fig 4: 비가산성
  message("\n[Figures] SI")
  fig2_threshold(d); fig5_water_constraint(d); fig4_resource(d); fig5_distance()
  message("\n완료: ", TDIR, " / ", FDIR)
}

# ===========================================================================
# FIGURE 5 — 소방용수 제약 
#
# ===========================================================================
.smooth_by <- function(x, y, f = 0.6, n_out = 80) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 40) return(NULL)
  lw <- lowess(x[ok], y[ok], f = f)
  xo <- seq(quantile(x[ok], .02), quantile(x[ok], .98), length.out = n_out)
  data.table(x = xo, y = approx(lw$x, lw$y, xout = xo, rule = 2)$y)
}

fig5_water_constraint <- function(d, span = 0.6) {
  x <- d[is.finite(heli_n) & is.finite(dur_hr) & is.finite(burned_ha)]
  if (nrow(x) < 150) { message("  (Fig5 생략: 진화자원 표본 부족)"); return(NULL) }
  x[, lb := log10(burned_ha)]

  cur <- function(yv) {
    a <- .smooth_by(x[compound == 0]$lb, x[compound == 0][[yv]], span)
    b <- .smooth_by(x[compound == 1]$lb, x[compound == 1][[yv]], span)
    list(ctl = a, cmp = b)
  }
  c_heli <- cur("heli_n")
  x[, ldur := log10(pmax(dur_hr, 0.25))]
  c_dur  <- cur("ldur")

  # (c) 자원 수준별 진화시간 — 표본 30건 미만 셀 제외
  min_cell <- 30
  x[, heli_bin := cut(heli_n, c(-1, 0, 2, 5, Inf),
                      labels = c("0", "1-2", "3-5", ">5"))]
  hc <- x[!is.na(heli_bin), .(n = .N, dur = median(dur_hr),
                              lo = quantile(dur_hr, .25),
                              hi = quantile(dur_hr, .75)),
          by = .(heli_bin, compound)][order(heli_bin, compound)]
  ok_bins <- hc[, .(mn = min(n)), by = heli_bin][mn >= min_cell]$heli_bin
  hc <- hc[heli_bin %in% ok_bins]
  fwrite(hc, file.path(TDIR, "Fig5c_duration_by_heli_data.csv"))

  # 모형 기반 보정 효과(규모 연속 통제) — 캡션에 넣을 수치
  cov <- paste(paste(CTRL, collapse = " + "), "+", FE)
  adj <- rbindlist(lapply(c("heli_n","dur_hr"), function(y) {
    x2 <- copy(x); x2[, yy := log1p(get(y))]
    ct <- cl(lm(fml("yy", paste("compound + log_burned +", cov)), x2), x2$cluster)
    data.table(outcome = y, beta = round(ct["compound",1], 4),
               p = signif(ct["compound",4], 3),
               pct = round(100*(exp(ct["compound",1])-1), 1))
  }))
  fwrite(adj, file.path(TDIR, "Fig5_adjusted_effects.csv"))
  print(adj)

  .plot_both("Fig5_water_constraint", quote({
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(1, 1.05))
    par(mar = c(4.6, 4.6, 3.0, 1.0), mgp = c(2.7, 0.7, 0), las = 1)
    axlab <- function() {
      at <- c(-2,-1,0,1,2,3,4)
      axis(1, at = at, labels = c("0.01","0.1","1","10","100","1,000","10,000"))
      mtext("Burned area (ha, log scale)", side = 1, line = 2.5, cex = 0.85)
    }
    pan <- function(cv, ylab, main, tag, ylog_lab = NULL) {
      yr <- range(c(cv$ctl$y, cv$cmp$y), na.rm = TRUE)
      plot(cv$ctl$x, cv$ctl$y, type = "n", xaxt = "n",
           yaxt = if (is.null(ylog_lab)) "s" else "n",
           ylim = yr + c(-0.05, 0.12)*diff(yr),
           xlab = "", ylab = ylab, main = main)
      axlab()
      if (!is.null(ylog_lab)) axis(2, at = ylog_lab$at, labels = ylog_lab$lab)
      polygon(c(cv$cmp$x, rev(cv$ctl$x)), c(cv$cmp$y, rev(cv$ctl$y)),
              col = rgb(180/255, 83/255, 9/255, 0.13), border = NA)
      lines(cv$ctl$x, cv$ctl$y, col = COL[["control"]], lwd = 2.8)
      lines(cv$cmp$x, cv$cmp$y, col = COL[["compound"]], lwd = 2.8, lty = 1)
      mtext(tag, side = 3, line = 1.4, adj = -0.20, font = 2, cex = 1.0)
    }
    pan(c_heli, "Helicopters deployed", "Fewer helicopters at the same fire size", "(a)")
    legend("topleft", bty = "n", cex = 0.85, lwd = 2.8,
           col = c(COL[["control"]], COL[["compound"]]),
           legend = c("Single-hazard fire", "Compound (drought) fire"))
    pan(c_dur, "Suppression duration (h)", "Longer suppression at the same fire size",
        "(b)", ylog_lab = list(at = log10(c(1,2,5,10,50,200)),
                               lab = c("1","2","5","10","50","200")))

    # (c)
    par(mar = c(6.2, 4.6, 3.2, 1.0))
    w  <- dcast(hc, heli_bin ~ compound, value.var = "dur")
    wl <- dcast(hc, heli_bin ~ compound, value.var = "lo")
    wh <- dcast(hc, heli_bin ~ compound, value.var = "hi")
    for (z in list(w, wl, wh)) setnames(z, c("0","1"), c("ctl","cmp"), skip_absent = TRUE)
    xx <- seq_len(nrow(w))
    yr <- range(c(wl$ctl, wl$cmp, wh$ctl, wh$cmp), na.rm = TRUE)
    plot(xx, w$ctl, type = "n", xaxt = "n", log = "y",
         ylim = c(max(yr[1]*0.75, 0.2), yr[2]*1.9), xlim = c(0.8, nrow(w)+0.25),
         xlab = "Helicopters deployed",
         ylab = "Suppression duration (h, log scale)",
         main = "At the same level of air support, compound fires still take longer")
    axis(1, at = xx, labels = as.character(w$heli_bin))
    off <- 0.07
    for (g in c("ctl","cmp")) {
      col <- if (g == "ctl") COL[["control"]] else COL[["compound"]]
      dx  <- if (g == "ctl") -off else off
      segments(xx+dx, wl[[g]], xx+dx, wh[[g]], col = col, lwd = 1.6)
      lines(xx+dx, w[[g]], col = col, lwd = 2.8)
      points(xx+dx, w[[g]], pch = if (g=="ctl") 19 else 17, col = col, cex = 1.4)
    }
    text(xx, pmax(w$cmp, w$ctl)*1.35, sprintf("%.2fx", w$cmp/w$ctl), cex = 0.8,
         col = COL[["line"]], xpd = NA)
    mtext("(c)", side = 3, line = 1.8, adj = -0.08, font = 2, cex = 1.0)
    legend("topleft", bty = "n", cex = 0.9,
           legend = c("Single-hazard fire", "Compound (drought) fire"),
           col = c(COL[["control"]], COL[["compound"]]), pch = c(19, 17), lwd = 2.8)
    nl <- dcast(hc, heli_bin ~ compound, value.var = "n")
    setnames(nl, c("0","1"), c("ctl","cmp"), skip_absent = TRUE)
    mtext(sprintf("n (single / compound): %s",
                  paste(sprintf("%d/%d", nl$ctl, nl$cmp), collapse = "   ")),
          side = 1, line = 3.5, cex = 0.75, col = COL[["line"]])
    mtext("Curves: LOWESS. Bars: interquartile range. Cells with n < 30 omitted.",
          side = 1, line = 4.4, cex = 0.7, col = COL[["line"]])
  }), w = 10, h = 9.0)
  list(heli = c_heli, dur = c_dur, byheli = hc, adjusted = adj)
}

# ===========================================================================
# FIGURE 
#  
# ===========================================================================
fig1_dose_response <- function(d) {
  has_dur <- "dry_months_1" %in% names(d) && sum(is.finite(d$dry_months_1)) > 500

  # (좌) 가뭄 강도: SPI-3 구간
  b1 <- c(-Inf, -2, -1.5, -1, -0.5, 0, 0.5, Inf)
  l1 <- c("< -2", "-2\nto -1.5", "-1.5\nto -1", "-1\nto -0.5",
          "-0.5\nto 0", "0\nto 0.5", "> 0.5")
  x <- d[is.finite(spi3)]
  x[, sbin := cut(spi3, b1, labels = l1)]
  s1 <- x[!is.na(sbin), .(n = .N,
                          med = median(burned_ha),
                          p90 = quantile(burned_ha, .9),
                          pct100 = 100*mean(large_fire),
                          dur = median(dur_hr, na.rm = TRUE)), by = sbin][order(sbin)]
  fwrite(s1, file.path(TDIR, "Fig1a_dose_intensity.csv"))

  # (우) 가뭄 지속기간: 연속 개월수
  s2 <- NULL
  if (has_dur) {
    y <- d[is.finite(dry_months_1)]
    y[, dbin := cut(dry_months_1, c(-1, 0, 1, 2, 3, 6, Inf),
                    labels = c("0", "1", "2", "3", "4-6", ">6"))]
    s2 <- y[!is.na(dbin), .(n = .N,
                            med = median(burned_ha),
                            p90 = quantile(burned_ha, .9),
                            pct100 = 100*mean(large_fire),
                            dur = median(dur_hr, na.rm = TRUE)), by = dbin][order(dbin)]
    s2 <- s2[n >= 30]
    fwrite(s2, file.path(TDIR, "Fig1b_dose_duration.csv"))
  }

  .plot_both("Fig1_dose_response", quote({
    par(mfrow = c(2, if (has_dur) 2 else 1), mar = c(5.2, 4.8, 3.0, 1.2),
        mgp = c(2.9, 0.7, 0), las = 1)

    bars <- function(s, xl, yv, ylab, main, tag, col, fmt = "%.2f") {
      v <- s[[yv]]
      bp <- barplot(v, names.arg = as.character(s[[1]]), col = col, border = NA,
                    ylab = ylab, main = main, ylim = c(0, max(v)*1.28),
                    cex.names = 0.82)
      text(bp, v, sprintf(fmt, v), pos = 3, cex = 0.78, xpd = NA)
      mtext(sprintf("n = %s", paste(s$n, collapse = " / ")), side = 1,
            line = 3.4, cex = 0.68, col = COL[["line"]])
      mtext(xl, side = 1, line = 2.3, cex = 0.86)
      mtext(tag, side = 3, line = 1.4, adj = -0.20, font = 2, cex = 1.05)
      invisible(bp)
    }
    bars(s1, "SPI-3 (month before ignition)", "pct100", "Fires >= 100 ha (%)",
         "Deeper drought, higher large-fire risk", "(a)", COL[["compound"]])
    if (has_dur)
      bars(s2, "Consecutive months with SPI-3 <= -1", "pct100",
           "Fires >= 100 ha (%)", "Longer drought, higher large-fire risk",
           "(b)", COL[["compound"]])
    bars(s1, "SPI-3 (month before ignition)", "p90",
         "90th pct burned area (ha)", "Deeper drought, larger fires",
         "(c)", COL[["accent"]])
    if (has_dur)
      bars(s2, "Consecutive months with SPI-3 <= -1", "p90",
           "90th pct burned area (ha)", "Longer drought, larger fires",
           "(d)", COL[["accent"]])
  }), w = if (has_dur) 10 else 6, h = 9)
  list(intensity = s1, duration = s2)
}

# ===========================================================================
# FIGURE 2
# ===========================================================================
fig2_cascade <- function(d) {
  cov <- paste(paste(CTRL, collapse = " + "), "+", FE)
  outs <- list(
    list(v = "log_burned", lab = "Burned area",          tf = TRUE),
    list(v = "log_dur",    lab = "Suppression duration", tf = TRUE),
    list(v = "heli_n",     lab = "Helicopters deployed", tf = "log1p"),
    list(v = "crew_n",     lab = "Personnel deployed",   tf = "log1p"),
    list(v = "equip_n",    lab = "Equipment deployed",   tf = "log1p"))
  res <- rbindlist(lapply(outs, function(o) {
    if (!o$v %in% names(d)) return(NULL)
    x <- copy(d)
    if (identical(o$tf, "log1p")) x[, yy := log1p(get(o$v))] else x[, yy := get(o$v)]
    x <- x[is.finite(yy)]
    if (nrow(x) < 200) return(NULL)
    ct <- cl(lm(fml("yy", paste("compound +", cov)), x), x$cluster)
    data.table(outcome = o$lab, beta = ct["compound",1], se = ct["compound",2],
               p = ct["compound",4], n = nrow(x))
  }))
  res[, `:=`(mult = exp(beta), lo = exp(beta - 1.96*se), hi = exp(beta + 1.96*se))]
  fwrite(res, file.path(TDIR, "Fig2_cascade_total_effects.csv"))
  print(res[, .(outcome, mult = round(mult,3), lo = round(lo,3),
                hi = round(hi,3), p = signif(p,3), n)])

  .plot_both("Fig2_cascade_total_effects", quote({
    par(mar = c(5.0, 11.5, 3.2, 9.5), mgp = c(2.8, 0.7, 0), las = 1)
    y <- rev(seq_len(nrow(res)))
    xr <- range(c(res$lo, res$hi, 1))
    plot(res$mult, y, xlim = c(min(xr)*0.96, max(xr)*1.08), ylim = c(0.4, nrow(res)+0.6),
         log = "x", yaxt = "n", pch = 19, cex = 1.5, col = COL[["compound"]],
         xlab = "Multiplicative effect of antecedent drought (total effect)",
         ylab = "", main = "One cascade: bigger fires, more resources, longer suppression")
    abline(v = 1, lty = 2, col = COL[["line"]], lwd = 1.6)
    segments(res$lo, y, res$hi, y, col = COL[["compound"]], lwd = 2.6)
    axis(2, at = y, labels = res$outcome, las = 1, cex.axis = 0.95, tick = FALSE)
    pf <- function(p) if (p < 1e-300) "p < 1e-300" else
      sprintf("p = %s", format(signif(p, 2), scientific = TRUE))
    text(res$hi, y, sprintf("  %.2fx   %s", res$mult, vapply(res$p, pf, character(1))),
         pos = 4, cex = 0.78, xpd = NA)
    mtext("Total effects. Burned area is not controlled: it is a mediator of the cascade.",
          side = 1, line = 3.4, cex = 0.74, col = COL[["line"]])
  }), w = 9.5, h = 4.8)
  res
}

# ===========================================================================
# FIGURE — 연도별 피해면적과 가뭄 연관 산불 비율
# 
# ===========================================================================
fig_annual_burned_share <- function(label_years = c(2000, 2022, 2025),
                                    log_bars = FALSE,
                                    area_source = c("ledger", "panel"),
                                    thr = -1.0) {
  area_source <- match.arg(area_source)
  p <- fread("data_out/incident_panel.csv")
  p[, comp_all := as.integer(spi3 <= thr)]      # 완충구간 없이 전수 판정

  share <- p[is.finite(comp_all), .(n_panel = .N,
                                    share = 100 * mean(comp_all)), by = year]

  if (area_source == "ledger" && file.exists("data_raw/kfs_fires.csv")) {
    l <- fread("data_raw/kfs_fires.csv")
    area <- l[, .(n_all = .N, burned = sum(burned_ha, na.rm = TRUE)), by = year]
  } else {
    area <- p[, .(n_all = .N, burned = sum(burned_ha, na.rm = TRUE)), by = year]
  }
  s <- merge(area, share, by = "year", all = TRUE)[order(year)]
  s <- s[is.finite(burned) & is.finite(share)]
  fwrite(s, file.path(TDIR, "Fig_annual_burned_share_data.csv"))
  print(s)

  .plot_both("Fig_annual_burned_share", quote({
    par(mar = c(4.4, 5.4, 3.0, 5.2), mgp = c(3.2, 0.7, 0), las = 1)
    yb <- s$burned
    if (log_bars) yb <- pmax(yb, 1)
    ymax <- max(yb) * 1.22

    bp <- barplot(yb, names.arg = NA, col = COL[["light"]], border = NA,
                  ylim = c(0, ymax), log = if (log_bars) "y" else "",
                  ylab = "", axes = FALSE)
    at <- pretty(c(0, max(yb)), 5)
    axis(2, at = at, labels = format(at, big.mark = ",", scientific = FALSE),
         cex.axis = 0.88)
    mtext("Annual burned area (ha)", side = 2, line = 4.0, las = 0, cex = 0.98)
    axis(1, at = bp, labels = s$year, cex.axis = 0.78, tick = FALSE, line = -0.4)

    # 극단 연도 값 라벨
    sel <- which(s$year %in% label_years)
    if (length(sel))
      text(bp[sel], yb[sel], format(round(s$burned[sel]), big.mark = ","),
           pos = 3, cex = 0.78, font = 2, col = COL[["line"]], xpd = NA)

    # 우측축: 가뭄 연관 비율
    par(new = TRUE)
    plot(bp, s$share, type = "n", axes = FALSE, xlab = "", ylab = "",
         ylim = c(0, max(s$share) * 1.35), xlim = range(bp) + c(-0.6, 0.6))
    lines(bp, s$share, col = COL[["compound"]], lwd = 2.8)
    points(bp, s$share, pch = 19, col = COL[["compound"]], cex = 1.0)
    axis(4, cex.axis = 0.88, col = COL[["compound"]],
         col.axis = COL[["compound"]])
    mtext("Drought-associated wildfires (%)", side = 4, line = 3.0, las = 0,
          cex = 0.98, col = COL[["compound"]])

    title(main = "Annual burned area and drought-associated wildfire share")
    legend("topleft", bty = "n", cex = 0.85, inset = c(0, -0.02),
           legend = c("Burned area (left axis)",
                      sprintf("Drought-associated share (right axis, SPI-3 <= %.1f)", thr)),
           fill = c(COL[["light"]], NA), border = NA,
           col = c(NA, COL[["compound"]]), lwd = c(NA, 2.8), pch = c(NA, 19),
           merge = TRUE)
  }), w = 10, h = 5.2)
  s
}
