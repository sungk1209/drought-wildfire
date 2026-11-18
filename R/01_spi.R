# ============================================================================
# 01_spi.R — SPI / SPEI 

# ============================================================================

suppressPackageStartupMessages({
  library(MASS)
})

# ---------------------------------------------------------------- PET
pet_hargreaves <- function(tmean, tmax, tmin, month, lat_deg) {
  # Hargreaves-Samani 월 PET (mm). tmax/tmin 이 있으면 Thornthwaite 보다 권장.
  lat <- lat_deg * pi / 180
  doy_mid <- c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  dec <- 0.409 * sin(2 * pi * doy_mid / 365 - 1.39)
  dr  <- 1 + 0.033 * cos(2 * pi * doy_mid / 365)
  ws  <- acos(pmin(pmax(-tan(lat) * tan(dec), -1), 1))
  Ra  <- (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(lat) * sin(dec) + cos(lat) * cos(dec) * sin(ws))   # MJ m-2 d-1
  ndays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  tr <- pmax(tmax - tmin, 0)
  pet_day <- 0.0023 * (tmean + 17.8) * sqrt(tr) * Ra[month] * 0.408
  pmax(pet_day, 0) * ndays[month]
}

pet_thornthwaite <- function(tmean, month, year, lat_deg) {
  t <- pmax(tmean, 0)
  i <- (t / 5)^1.514
  I <- ave(i, year, FUN = function(z) sum(z, na.rm = TRUE))
  a <- 6.75e-7 * I^3 - 7.71e-5 * I^2 + 1.792e-2 * I + 0.49239
  pet_un <- 16 * (10 * t / ifelse(I == 0, NA, I))^a
  lat <- lat_deg * pi / 180
  doy_mid <- c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  dec <- 0.409 * sin(2 * pi * doy_mid / 365 - 1.39)
  ws  <- acos(pmin(pmax(-tan(lat) * tan(dec), -1), 1))
  daylight <- 24 / pi * ws
  ndays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  out <- pet_un * (daylight[month] / 12) * (ndays[month] / 30)
  ifelse(is.finite(out), out, 0)
}

# ---------------------------------------------------------------- 분포 적합
.fit_gamma_month <- function(x) {
  # 무강수 확률 q 를 분리한 혼합 감마. floc=0 고정(2모수)은 SPI 표준 관행.
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  q  <- mean(x <= 0)
  xp <- x[x > 0]
  if (length(xp) < 8) return(NULL)
  fit <- try(suppressWarnings(MASS::fitdistr(xp, "gamma")), silent = TRUE)
  if (inherits(fit, "try-error")) {
    # 적률법 대체
    m <- mean(xp); v <- var(xp)
    if (v <= 0) return(NULL)
    return(list(q = q, shape = m^2 / v, rate = m / v))
  }
  list(q = q, shape = fit$estimate[["shape"]], rate = fit$estimate[["rate"]])
}

.dllog <- function(x, a, b) {          # 2모수 log-logistic 밀도
  (b / a) * (x / a)^(b - 1) / (1 + (x / a)^b)^2
}
.pllog <- function(x, a, b) 1 / (1 + (a / x)^b)

.fit_llog_month <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  shift <- if (min(x) <= 0) -min(x) + 1 else 0
  y <- x + shift
  nll <- function(p) {
    a <- exp(p[1]); b <- exp(p[2])
    d <- .dllog(y, a, b)
    if (any(!is.finite(d)) || any(d <= 0)) return(1e10)
    -sum(log(d))
  }
  op <- try(optim(c(log(median(y)), log(2)), nll, method = "Nelder-Mead"),
            silent = TRUE)
  if (inherits(op, "try-error")) return(NULL)
  list(a = exp(op$par[1]), b = exp(op$par[2]), shift = shift)
}

.to_z <- function(p) qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))

# ---------------------------------------------------------------- 공개 함수
#' @param x     월 시계열 값 (강수 합계 또는 P-PET), 시간순 정렬
#' @param ym    x 에 대응하는 월 시작일 (Date)
#' @param scale 누적 개월수 (SPI3 -> 3)
compute_spi <- function(x, ym, scale = 3, base = c(1991, 2020)) {
  n <- length(x)
  acc <- rep(NA_real_, n)
  if (n >= scale) {
    cs <- c(0, cumsum(ifelse(is.na(x), 0, x)))
    ok <- c(0, cumsum(as.integer(!is.na(x))))
    for (i in scale:n) {
      if (ok[i + 1] - ok[i + 1 - scale] == scale) acc[i] <- cs[i + 1] - cs[i + 1 - scale]
    }
  }
  mo <- as.integer(format(ym, "%m")); yr <- as.integer(format(ym, "%Y"))
  inbase <- yr >= base[1] & yr <= base[2]
  out <- rep(NA_real_, n)
  for (m in 1:12) {
    pm <- .fit_gamma_month(acc[mo == m & inbase])
    if (is.null(pm)) next
    idx <- which(mo == m & is.finite(acc))
    v <- acc[idx]
    p <- ifelse(v <= 0, pm$q / 2,
                pm$q + (1 - pm$q) * pgamma(v, shape = pm$shape, rate = pm$rate))
    out[idx] <- .to_z(p)
  }
  out
}

compute_spei <- function(d, ym, scale = 3, base = c(1991, 2020)) {
  # d = P - PET 월 시계열
  n <- length(d)
  acc <- rep(NA_real_, n)
  if (n >= scale) {
    cs <- c(0, cumsum(ifelse(is.na(d), 0, d)))
    ok <- c(0, cumsum(as.integer(!is.na(d))))
    for (i in scale:n) {
      if (ok[i + 1] - ok[i + 1 - scale] == scale) acc[i] <- cs[i + 1] - cs[i + 1 - scale]
    }
  }
  mo <- as.integer(format(ym, "%m")); yr <- as.integer(format(ym, "%Y"))
  inbase <- yr >= base[1] & yr <= base[2]
  out <- rep(NA_real_, n)
  for (m in 1:12) {
    pm <- .fit_llog_month(acc[mo == m & inbase])
    if (is.null(pm)) next
    idx <- which(mo == m & is.finite(acc))
    out[idx] <- .to_z(.pllog(acc[idx] + pm$shift, pm$a, pm$b))
  }
  out
}

#' 실효습도  EH(t) = (1-r) * sum_{i=0}^{nd-1} r^i * RH(t-i)
#' 한국 기상청/산림청 관행 r=0.7, 5일. 논문 SI에 식을 그대로 제시할 것.
effective_humidity <- function(rh, r = 0.7, nd = 5) {
  w <- (1 - r) * r^(0:(nd - 1))
  n <- length(rh)
  out <- rep(NA_real_, n)
  for (i in nd:n) {
    v <- rh[i:(i - nd + 1)]
    if (!any(is.na(v))) out[i] <- sum(v * w)
  }
  out
}

#' 직전 무강수 지속일수
dry_spell_len <- function(precip, thr = 1.0) {
  n <- length(precip); out <- integer(n); run <- 0L
  for (i in seq_len(n)) {
    if (!is.na(precip[i]) && precip[i] > thr) run <- 0L else run <- run + 1L
    out[i] <- run
  }
  out
}
