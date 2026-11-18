# ============================================================================
# 14_analysis_q2.R 
# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages({
  library(data.table); library(sandwich); library(lmtest)
  library(quantreg);   library(MASS)
})

RES <- character(0)
say  <- function(...) { s <- paste0(...); cat(s, "\n"); RES <<- c(RES, s) }
pr   <- function(x) { o <- capture.output(print(x)); cat(o, sep = "\n"); RES <<- c(RES, o) }
head_ <- function(t) { say(""); say(strrep("=", 74)); say(t); say(strrep("=", 74)) }

# 사전처치 공변량만. wet_days_mo 는 행동 경로 분리를 위한 단기 날씨 통제.
CTRL <- c("rhmin", "wsmax", "tmean", "dry_spell", "wet_days_mo")
FE   <- "factor(season) + factor(sido) + factor(year) + factor(cause)"
fml  <- function(y, rhs) as.formula(paste(y, "~", rhs))
cl   <- function(m, g) coeftest(m, vcov = vcovCL(m, cluster = g))

load_panel <- function(path = "data_out/incident_panel.csv") {
  d <- fread(path)
  d <- d[!is.na(compound)]
  need <- c("log_burned", "compound", CTRL)
  need <- intersect(need, names(d))
  d <- d[complete.cases(d[, ..need])]
  # complete.cases 는 NA 만 잡는다. Inf/-Inf 가 남으면 glm 이 'NA/NaN/Inf in x' 로 죽는다.
  numv <- need[vapply(need, function(v) is.numeric(d[[v]]), logical(1))]
  if (length(numv)) {
    keepr <- Reduce(`&`, lapply(numv, function(v) is.finite(d[[v]])))
    nbad <- sum(!keepr)
    if (nbad) message("비유한값(Inf/NaN) 행 제외: ", nbad)
    d <- d[keepr]
  }
  d[, cause := fifelse(is.na(cause) | !nzchar(cause), "미상", cause)]
  # 표본이 적은 원인 범주는 '기타' 로 묶어 고정효과 과적합을 막는다
  tb <- d[, .N, by = cause]
  d[cause %in% tb[N < 30]$cause, cause := "기타"]
  d[]
}

# ------------------------------------------------------- A0 원인구성 점검
cause_balance <- function(d) {
  head_("A0. 발생원인 구성 점검 — 비교가능한 산불끼리 비교했는가")
  t <- dcast(d[, .N, by = .(cause, compound)], cause ~ compound,
             value.var = "N", fill = 0)
  setnames(t, c("0","1"), c("control","compound"), skip_absent = TRUE)
  t[, `:=`(pct_control  = round(100 * control / sum(control), 1),
           pct_compound = round(100 * compound / sum(compound), 1))]
  t[, diff_pp := round(pct_compound - pct_control, 1)]
  pr(t[order(-pct_compound)])
  ct <- suppressWarnings(chisq.test(table(d$cause, d$compound)))
  say(sprintf("\n  카이제곱 검정: X2=%.1f, df=%d, p=%.3g",
              ct$statistic, ct$parameter, ct$p.value))
  say("  해석: 구성비 차이가 작으면 '유사한 종류의 산불끼리 비교' 근거가 된다.")
  say("        차이가 크면 원인 고정효과 통제가 필수")

}

# ------------------------------------------------------- A 균형표
balance <- function(d) {
  head_("A. 처치/대조 균형표 (Standardized Mean Difference)")
  vs <- intersect(c(CTRL, "precip_prev7d", "res_anom", "stn_km"), names(d))
  t <- rbindlist(lapply(vs, function(v) {
    a <- d[compound == 1][[v]]; b <- d[compound == 0][[v]]
    s <- sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
    data.table(variable = v, compound = round(mean(a, na.rm = TRUE), 3),
               control = round(mean(b, na.rm = TRUE), 3),
               SMD = round((mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / s, 3))
  }))
  pr(t)
  say("\n  |SMD| > 0.10 은 회귀에서 반드시 통제 (또는 G 의 매칭으로 균형화)")
}

# ------------------------------------------------------- B, C 주회귀
main_reg <- function(d) {
  head_("B. 주회귀 — log(피해면적) ~ 가뭄 + 통제")
  specs <- list(
    "M1 무통제"            = "compound",
    "M2 +화재기상"         = paste("compound", paste(CTRL, collapse = " + "), sep = " + "),
    "M3 +계절/지역/연도/원인" = paste("compound", paste(CTRL, collapse = " + "), FE, sep = " + "))
  t <- rbindlist(lapply(names(specs), function(nm) {
    m <- lm(fml("log_burned", specs[[nm]]), d); ct <- cl(m, d$cluster)
    data.table(model = nm, beta = round(ct["compound",1], 4),
               SE = round(ct["compound",2], 4), p = signif(ct["compound",4], 3),
               mult = round(exp(ct["compound",1]), 3),
               n = nobs(m), R2 = round(summary(m)$r.squared, 3))
  }))
  pr(t)
  say("\n  mult = 배수효과. M1->M3 에서 계수가 유지되면 가뭄 효과가")
  say("  화재기상/계절/원인 교란으로 설명되지 않는다는 뜻 (R1-4.1 대응).")

  head_("C. Gamma GLM (log link) — 분포 가정 강건성")
  gfit <- try(glm(fml("burned_ha", specs[["M3 +계절/지역/연도/원인"]]), d,
                  family = Gamma(link = "log"),
                  control = glm.control(maxit = 100)), silent = TRUE)
  if (inherits(gfit, "try-error")) {
    say("  Gamma GLM 적합 실패: ", conditionMessage(attr(gfit, "condition")))
    say("  대체로 피해면적 분포의 극단 왜도(최대 5만 ha) 때문이다.")

    # 대안 1: 준포아송 (과대산포 허용, 로그링크)
    q <- try(glm(fml("burned_ha", specs[["M3 +계절/지역/연도/원인"]]), d,
                 family = quasipoisson(link = "log"),
                 control = glm.control(maxit = 100)), silent = TRUE)
    if (!inherits(q, "try-error")) {
      ct <- cl(q, d$cluster)
      say(sprintf("  [대안] 준포아송(log link): beta=%+.4f (SE %.4f, p=%.3g, 배수 %.2f)",
                  ct["compound",1], ct["compound",2], ct["compound",4],
                  exp(ct["compound",1])))
    }
    # 대안 2: 상위 0.1% 절단 후 Gamma (극단값 영향 확인)
    dd <- d[burned_ha <= quantile(burned_ha, 0.999)]
    g2 <- try(glm(fml("burned_ha", specs[["M3 +계절/지역/연도/원인"]]), dd,
                  family = Gamma(link = "log"),
                  control = glm.control(maxit = 100)), silent = TRUE)
    if (!inherits(g2, "try-error")) {
      ct <- cl(g2, dd$cluster)
      say(sprintf("  [대안] Gamma, 상위 0.1%% 절단: beta=%+.4f (p=%.3g, 배수 %.2f, n=%d)",
                  ct["compound",1], ct["compound",4], exp(ct["compound",1]), nrow(dd)))
    }
  } else {
    ct <- cl(gfit, d$cluster)
    say(sprintf("  beta=%+.4f (SE %.4f, p=%.3g, 배수 %.2f)",
                ct["compound",1], ct["compound",2], ct["compound",4],
                exp(ct["compound",1])))
  }

  # 진화시간도 같은 사양으로
  head_("B2. 주회귀 — log(진화시간)")
  dd <- d[is.finite(log_dur)]
  m <- lm(fml("log_dur", specs[["M3 +계절/지역/연도/원인"]]), dd)
  ct <- cl(m, dd$cluster)
  say(sprintf("  beta=%+.4f (SE %.4f, p=%.3g, 배수 %.2f, n=%d)",
              ct["compound",1], ct["compound",2], ct["compound",4],
              exp(ct["compound",1]), nobs(m)))
}

# ------------------------------------------------------- D 분위회귀
quantiles <- function(d, R = 100) {
  head_("D. 분위회귀 — 효과가 어느 분위에서 나타나는가")
  rhs <- paste("compound", paste(CTRL, collapse = " + "), sep = " + ")
  t <- rbindlist(lapply(c(.10,.25,.50,.75,.90,.95,.99), function(q) {
    m <- rq(fml("log_burned", rhs), tau = q, data = d)
    s <- summary(m, se = "boot", R = R)$coefficients
    data.table(tau = q, beta = round(s["compound",1], 4),
               SE = round(s["compound",2], 4), p = signif(s["compound",4], 3),
               mult = round(exp(s["compound",1]), 3))
  }))
  pr(t)
  say("\n  beta 가 tau 와 함께 커지면 '가뭄은 대형화재 야기 '")
  say("  -> 평균 비교보다 훨씬 강한 증폭 근거 (R1-2 대응). SE 는 부트스트랩.")
}

# ------------------------------------------------------- E 상호작용
synergy <- function(d) {
  head_("E. 비가산성(synergy) 검정 — R1-2 '벤다이어그램' 비판 대응")
  x <- copy(d)
  x[, ws_c := wsmax - mean(wsmax)][, rh_c := rhmin - mean(rhmin)]
  for (v in c("ws_c", "rh_c")) {
    other <- setdiff(CTRL, if (v == "ws_c") "wsmax" else "rhmin")
    m <- lm(fml("log_burned", paste0("compound * ", v, " + ",
                paste(other, collapse = " + "), " + ", FE)), x)
    ct <- cl(m, x$cluster); k <- grep(":", rownames(ct), value = TRUE)[1]
    say(sprintf("  %-18s beta=%+.4f  SE=%.4f  p=%.3g",
                if (v == "ws_c") "가뭄 x 최대풍속" else "가뭄 x 최저습도",
                ct[k,1], ct[k,2], ct[k,4]))
  }
  say("\n  상호작용항이 유의한 양수면 '결합효과 > 각각의 합' = 정량화된 synergy.")

  x[, hi_ws := as.integer(wsmax > median(wsmax))]
  cel <- dcast(x[, .(m = mean(burned_ha)), by = .(compound, hi_ws)],
               compound ~ hi_ws, value.var = "m")
  b <- cel[compound == 0][["0"]]; a1 <- cel[compound == 1][["0"]]
  a2 <- cel[compound == 0][["1"]]; ob <- cel[compound == 1][["1"]]
  ex <- (a1 - b) + (a2 - b) + b
  say(sprintf("\n  2x2 초가산성(평균 ha): 기저 %.2f / 가뭄만 %.2f / 강풍만 %.2f",
              b, a1, a2))
  say(sprintf("  가산 기대 %.2f vs 실제 동시발생 %.2f -> 초과 %+.2f ha (%.2f배)",
              ex, ob, ob - ex, ob / ex))
}

# ------------------------------------------------------- F 임계값
threshold <- function(d) {
  head_("F. 가뭄 임계값 탐색 — 비선형 급변점")
  x <- d[is.finite(spi3)]
  best <- NULL
  scan <- rbindlist(lapply(seq(-2.0, 0.5, by = 0.05), function(cc) {
    xx <- copy(x)[, seg := pmin(spi3 - cc, 0)]
    m <- lm(fml("log_burned", paste("spi3 + seg +", paste(CTRL, collapse = " + "))), xx)
    ss <- sum(residuals(m)^2)
    if (is.null(best) || ss < best$ss) best <<- list(c = cc, ss = ss, m = m)
    data.table(c = cc, ssr = ss, beta_seg = coef(m)[["seg"]])
  }))
  s <- summary(best$m)$coefficients
  say(sprintf("  최적 급변점 SPI3 = %.2f", best$c))
  say(sprintf("  임계값 이하 추가 기울기 = %+.4f (p=%.3g)", s["seg",1], s["seg",4]))
  say("\n  구간별 피해면적")
  x[, bin := cut(spi3, c(-Inf,-2,-1.5,-1,-0.5,0,Inf))]
  pr(x[, .(n = .N, mean_ha = round(mean(burned_ha), 2),
           med_ha = round(median(burned_ha), 3),
           p90 = round(quantile(burned_ha, .9), 2),
           large_pct = round(mean(large_fire), 4)), by = bin][order(bin)])
  fwrite(scan, "data_out/threshold_scan.csv")
}

# ------------------------------------------------------- G 성향점수 매칭
psm <- function(d) {
  head_("G. 성향점수 매칭 — 사전처치 공변량만 사용")
  f <- paste("compound ~", paste(CTRL, collapse = " + "),
             "+ factor(season) + factor(sido)")
  p <- glm(as.formula(f), d, family = binomial())
  x <- copy(d)[, ps := predict(p, type = "response")]
  lg <- function(z) { z <- pmin(pmax(z, 1e-6), 1-1e-6); log(z/(1-z)) }
  x[, lps := lg(ps)]
  cal <- 0.2 * sd(x$lps)
  tr <- which(x$compound == 1); ctl <- which(x$compound == 0)
  cv <- x$lps[ctl]; o <- order(cv); cs <- cv[o]; ci <- ctl[o]
  used <- logical(length(ci)); pt <- integer(0); pc <- integer(0)
  for (i in tr) {
    v <- x$lps[i]; j <- findInterval(v, cs)
    k <- max(1, j-15):min(length(cs), j+15)
    k <- k[!used[k] & abs(cs[k]-v) <= cal]
    if (!length(k)) next
    kk <- k[which.min(abs(cs[k]-v))]
    used[kk] <- TRUE; pt <- c(pt, i); pc <- c(pc, ci[kk])
  }
  if (length(pt) < 30) { say("  매칭쌍 부족"); return(invisible()) }
  df <- x$log_burned[pt] - x$log_burned[pc]
  se <- sd(df)/sqrt(length(df))
  say(sprintf("  매칭쌍 %d (캘리퍼 %.3f)", length(pt), cal))
  say(sprintf("  ATT(log 면적)=%+.4f SE=%.4f t=%.2f 배수=%.2f",
              mean(df), se, mean(df)/se, exp(mean(df))))
  say(sprintf("  매칭 후 평균: compound %.2f ha vs control %.2f ha",
              mean(x$burned_ha[pt]), mean(x$burned_ha[pc])))
  say("\n  매칭 후 균형(SMD):")
  for (v in CTRL) {
    a <- x[[v]][pt]; b <- x[[v]][pc]
    say(sprintf("    %-14s %+.3f", v, (mean(a)-mean(b))/sqrt((var(a)+var(b))/2)))
  }
}

# ------------------------------------------------------- H 매개분석
mediation <- function(d, B = 500) {   # 논문용으로 B=2000 까지 올려도 됨
  head_("H. 매개분석 — 가뭄 -> 진화시간 -> 피해면적")
  x <- d[is.finite(log_dur)]
  cov <- paste(paste(CTRL, collapse = " + "),
               "+ factor(season) + factor(sido) + factor(year) + factor(cause)")
  am <- lm(fml("log_dur",    paste("compound +", cov)), x)
  bm <- lm(fml("log_burned", paste("compound + log_dur +", cov)), x)
  tm <- lm(fml("log_burned", paste("compound +", cov)), x)
  a <- coef(am)[["compound"]]; b <- coef(bm)[["log_dur"]]
  cp <- coef(bm)[["compound"]]; cc <- coef(tm)[["compound"]]
  sa <- summary(am)$coefficients; sb <- summary(bm)$coefficients
  say(sprintf("  a (가뭄->log진화시간)   = %+.4f  p=%.3g", a, sa["compound",4]))
  say(sprintf("  b (log진화시간->log면적)= %+.4f  p=%.3g", b, sb["log_dur",4]))
  say(sprintf("  c  총효과=%+.4f   c' 직접효과=%+.4f", cc, cp))
  say(sprintf("  a*b 간접효과=%+.4f   매개비율=%.1f%%", a*b, 100*a*b/cc))

  # 부트스트랩: 모형행렬을 한 번만 만들고 .lm.fit 으로 재적합한다
  # (lm() 을 매 반복 호출하면 formula 처리 비용이 지배적이라 10배 이상 느리다)
  Xa <- model.matrix(am); Xb <- model.matrix(bm)
  ya <- x$log_dur; yb <- x$log_burned
  ja <- match("compound", colnames(Xa)); jb <- match("log_dur", colnames(Xb))
  n <- nrow(Xa); bs <- numeric(B)
  for (i in seq_len(B)) {
    id <- sample.int(n, n, TRUE)
    fa <- .lm.fit(Xa[id, , drop = FALSE], ya[id])
    fb <- .lm.fit(Xb[id, , drop = FALSE], yb[id])
    bs[i] <- fa$coefficients[ja] * fb$coefficients[jb]
  }
  bs <- bs[is.finite(bs)]
  q <- quantile(bs, c(.025,.975), na.rm = TRUE)
  say(sprintf("  간접효과 부트스트랩 95%% CI = [%+.4f, %+.4f] (B=%d)", q[1], q[2], length(bs)))
  say("\n  간접효과가 유의하면 '가뭄이 진화를 지연시켜 피해를 키운다'는 사슬")
  say("  물리 데이터로 확인")
}

# ------------------------------------------------------- K 시간규모 분화
timescale <- function(d) {
  head_("K. 가뭄 시간규모별 효과 — SPI-3(단기) vs SPI-12(장기) 누적 강수부족")
  ix <- intersect(c("spi1","spi3","spi6","spi12"), names(d))
  if (length(ix) < 2) { say("  SPI12 없음 — 건너뜀"); return(invisible()) }
  outs <- c(log_burned = "log(피해면적)", log_dur = "log(진화시간)")
  if ("res_anom" %in% names(d) && sum(is.finite(d$res_anom)) > 100)
    outs <- c(outs, res_anom = "저수율 편차")

  t <- rbindlist(lapply(names(outs), function(y) rbindlist(lapply(ix, function(v) {
    x <- d[is.finite(get(v)) & is.finite(get(y))]
    if (nrow(x) < 100) return(NULL)
    x[, z := scale(get(v))[,1]]
    rhs <- if (y == "res_anom") "z + factor(season) + factor(sido)"
           else paste("z +", paste(CTRL, collapse = " + "), "+", FE)
    ct <- cl(lm(fml(y, rhs), x), x$cluster)
    data.table(outcome = outs[[y]], index = v, beta_std = round(ct["z",1], 4),
               SE = round(ct["z",2], 4), p = signif(ct["z",4], 3))
  }))))
  pr(t)
  say("\n  결과변수별로 |beta_std| 가 가장 큰 누적기간이 그 결과와 가장 강하게")
  say("  연관된 시간규모 SPI-12 는 강수 누적만 쓰는 기상학적 지수이므로")
  say("  수문학적 가뭄(SRI/SWSI/저수율)과 동일시하지 말 것. 저수율과의 상관을")
  say("  별도로 제시해야 장기 누적이 저류를 반영한다는 주장성립")

  if (all(c("spi3","spi12") %in% names(d))) {
    say("\n  [두 규모 동시 투입]")
    for (y in names(outs)) {
      x <- d[is.finite(spi3) & is.finite(spi12) & is.finite(get(y))]
      x[, `:=`(z3 = scale(spi3)[,1], z12 = scale(spi12)[,1])]
      rhs <- if (y == "res_anom") "z3 + z12 + factor(season) + factor(sido)"
             else paste("z3 + z12 +", paste(CTRL, collapse = " + "), "+", FE)
      ct <- cl(lm(fml(y, rhs), x), x$cluster)
      say(sprintf("    %-14s SPI3 %+.4f (p=%.3g) | SPI12 %+.4f (p=%.3g)",
                  outs[[y]], ct["z3",1], ct["z3",4], ct["z12",1], ct["z12",4]))
    }
  }
  if ("drought_type" %in% names(d)) {
    say("\n  [시간규모 조합별 — 기준: none / short=SPI3<=-1, long=SPI12<=-1]")
    x <- d[!is.na(drought_type)]
    x[, drought_type := relevel(factor(drought_type), ref = "none")]
    pr(x[, .(n = .N, burned_med = round(median(burned_ha), 3),
             dur_med = round(median(dur_hr, na.rm = TRUE), 2),
             large_pct = round(mean(large_fire), 4)), by = drought_type][order(-n)])
    for (y in c("log_burned","log_dur")) {
      x2 <- x[is.finite(get(y))]
      ct <- cl(lm(fml(y, paste("drought_type +", paste(CTRL, collapse = " + "),
                               "+", FE)), x2), x2$cluster)
      say(sprintf("\n    %s:", outs[[y]]))
      for (k in grep("^drought_type", rownames(ct), value = TRUE))
        say(sprintf("      %-16s %+.4f (SE %.4f, p=%.3g, 배수 %.2f)",
                    sub("drought_type","",k), ct[k,1], ct[k,2], ct[k,4], exp(ct[k,1])))
    }
  }
}

# ------------------------------------------------------- R 진화자원 통제
resource_control <- function(d) {
  head_("R. 진화자원 통제 — 소방용수 논증의 결정적 검정")
  x <- d[is.finite(heli_n) & is.finite(log_dur)]
  if (nrow(x) < 200) { say("  진화자원 표본 부족 (n=", nrow(x), ") — 건너뜀"); return(invisible()) }
  say("  표본: ", nrow(x), "건 (", min(x$year), "-", max(x$year), ")")
  x[, `:=`(l_heli = log1p(heli_n), l_crew = log1p(crew_n), l_equip = log1p(equip_n))]
  cov <- paste(paste(CTRL, collapse = " + "), "+", FE)

  specs <- list(
    "R1 자원 미통제"        = paste("compound +", cov),
    "R2 +헬기"              = paste("compound + l_heli +", cov),
    "R3 +헬기/인력/장비"    = paste("compound + l_heli + l_crew + l_equip +", cov),
    "R4 +규모까지 통제"     = paste("compound + l_heli + l_crew + l_equip + log_burned +", cov))
  t <- rbindlist(lapply(names(specs), function(nm) {
    m <- lm(fml("log_dur", specs[[nm]]), x); ct <- cl(m, x$cluster)
    data.table(model = nm, beta = round(ct["compound",1], 4),
               SE = round(ct["compound",2], 4), p = signif(ct["compound",4], 3),
               mult = round(exp(ct["compound",1]), 3), n = nobs(m))
  }))
  pr(t)
  say("\n  핵심 해석: R1 -> R3 에서 계수가 유지되면 진화 지연이 '자원 부족'으로")
  say("  설명되지 않는다는 뜻이다. R4 에서 피해면적까지 통제한 뒤에도 남으면")
  say("  '같은 규모, 같은 자원인데도 더 오래 걸린다' -> 용수 제약 논거 성립")
 

  say("\n  [가뭄이 투입 자원 자체를 늘리는가]")
  for (y in c("l_heli","l_crew","l_equip")) {
    ct <- cl(lm(fml(y, paste("compound + log_burned +", cov)), x), x$cluster)
    say(sprintf("    %-8s beta=%+.4f (p=%.3g)", y, ct["compound",1], ct["compound",4]))
  }
  say("    규모를 통제한 뒤에도 자원 투입이 늘면 '같은 크기 산불에 더 많은 자원이")
  say("    필요했다'는 뜻 -> 진화 난이도 증가 증거.")

  if ("stage" %in% names(x) && uniqueN(x$stage) > 1) {
    say("\n  [산불단계 분포]")
    pr(dcast(x[, .N, by = .(stage, compound)], stage ~ compound,
             value.var = "N", fill = 0))
  }
}

# ------------------------------------------------------- J 민감도
sensitivity <- function(d) {
  head_("J. 민감도 격자")
  src <- fread("data_out/incident_panel.csv")
  src[, cause := fifelse(is.na(cause) | !nzchar(cause), "미상", cause)]
  rhs <- paste("compound", paste(CTRL, collapse = " + "), FE, sep = " + ")
  t <- rbindlist(lapply(intersect(c("spi1","spi3","spi6","spi12","spei3","spei6","spei12"),
                                  names(src)), function(ix)
    rbindlist(lapply(c(-0.8,-1.0,-1.5), function(thr) {
      x <- src[is.finite(get(ix))]
      x <- x[complete.cases(x[, c("log_burned", CTRL), with = FALSE])]
      x[, compound := as.integer(get(ix) <= thr)]
      if (sum(x$compound) < 50) return(NULL)
      ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
      data.table(index = ix, thr = thr, n_comp = sum(x$compound),
                 beta = round(ct["compound",1], 4), SE = round(ct["compound",2], 4),
                 p = signif(ct["compound",4], 3), mult = round(exp(ct["compound",1]), 3))
    }))))
  pr(t); fwrite(t, "data_out/sensitivity_grid.csv")
  say("\n  모든 조합에서 부호/유의성이 유지되는지 확인")

  say("\n  [추가 민감도]")
  # 관측소 거리 먼 건 제외
  if ("far_flag" %in% names(d)) {
    x <- d[far_flag == 0]
    ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
    say(sprintf("    관측소 60km 초과 제외: beta=%+.4f (p=%.3g, n=%d)",
                ct["compound",1], ct["compound",4], nrow(x)))
  }
  # 2025년(초대형) 제외
  x <- d[year != 2025]
  ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
  say(sprintf("    2025년 제외: beta=%+.4f (p=%.3g)  <- 특정 연도 주도 제거거",
              ct["compound",1], ct["compound",4]))
  # 봄철(산불 계절)만
  x <- d[month %in% 2:5]
  ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
  say(sprintf("    봄철(2-5월)만: beta=%+.4f (p=%.3g, n=%d)",
              ct["compound",1], ct["compound",4], nrow(x)))
  # 중복 복합체 제외 (dup_flag 가 있으면)
  if ("dup_flag" %in% names(d)) {
    x <- d[dup_flag == 0]
    ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
    say(sprintf("    중복 복합체 제외: beta=%+.4f (p=%.3g, n=%d)",
                ct["compound",1], ct["compound",4], nrow(x)))
  }
  # 상위 0.1% 절단 (초대형 사건 영향)
  x <- d[burned_ha <= quantile(burned_ha, .999)]
  ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
  say(sprintf("    피해면적 상위 0.1%% 제외: beta=%+.4f (p=%.3g)",
              ct["compound",1], ct["compound",4]))
  # 상위 1% 절단
  x <- d[burned_ha <= quantile(burned_ha, .99)]
  ct <- cl(lm(fml("log_burned", rhs), x), x$cluster)
  say(sprintf("    피해면적 상위 1%% 제외: beta=%+.4f (p=%.3g)",
              ct["compound",1], ct["compound",4]))
}

# ------------------------------------------------------- main
run_analysis <- function(boot_mediation = 500, boot_quantile = 100) {
  d <- load_panel()
  head_("데이터 개요")
  say(sprintf("  분석 표본: %d 사건 (%d-%d, %d 시군구)", nrow(d),
              min(d$year), max(d$year), uniqueN(d$sigungu_key)))
  say(sprintf("  compound %d / control %d", sum(d$compound), sum(1-d$compound)))
  say(sprintf("  진화시간 확보: %d건 / 진화자원 확보: %d건",
              sum(is.finite(d$log_dur)), sum(is.finite(d$heli_n))))
  say("\n  추정 대상: 발화 조건부 증폭 (amplification conditional on ignition)")

  cause_balance(d); balance(d); main_reg(d)
  quantiles(d, R = boot_quantile); synergy(d); threshold(d); psm(d)
  mediation(d, B = boot_mediation); timescale(d); resource_control(d)
  sensitivity(d)

  writeLines(RES, "data_out/analysis_report.txt")
  say("\n[출력] data_out/analysis_report.txt")
}

if (sys.nframe() == 0L) run_analysis()
