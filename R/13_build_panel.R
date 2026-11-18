# ============================================================================
# 13_build_panel.R — 사건 단위 분석 테이블 구축  [산불피해대장 2000-2025 기준]
#

# 입력
#   data_raw/kfs_fires.csv            11_parse_damage_ledger.R 산출
#   data_raw/kma_asos_daily.csv       00_download.R 산출 (1990-01-01 부터)
#   data_raw/station_meta.csv
#   data_raw/sigungu_centroid.csv
#   data_out/sigungu_station_map.csv   12_station_remap.R 산출
#   data_raw/진화자원이력.xlsx         (선택) 헬기/인력/장비
#   data_raw/reservoir_daily.csv       (선택) 저수율
#
# 출력
#   data_out/incident_panel.csv
#   data_out/spi_panel.csv
#   data_out/build_report.txt
# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages(library(data.table))
source("R/01_spi.R")

CFG <- list(
  baseline      = c(1991, 2020),      
  spi_scales    = c(1, 3, 6, 12),       
  spei_scales   = c(3, 6, 12),
  compound_index = "spi3",
  compound_thr  = -1.0,                 
  control_thr   = -0.5,                
  min_burned_ha = 0.01
)

LOG <- character(0)
say <- function(...) { s <- paste0(...); cat(s, "\n"); LOG <<- c(LOG, s) }

# ---------------------------------------------------------------- 1. SPI 패널
build_spi_panel <- function(asos, stn) {
  b <- CFG$baseline
  asos <- as.data.table(asos)
  asos[, ym := as.Date(format(date, "%Y-%m-01"))]
  mon <- asos[, .(precip = if (sum(!is.na(precip)) >= 20) sum(precip, na.rm = TRUE)
                           else NA_real_,
                  tmean = mean(tmean, na.rm = TRUE),
                  tmax  = mean(tmax,  na.rm = TRUE),
                  tmin  = mean(tmin,  na.rm = TRUE)),
              by = .(stn_id, ym)][order(stn_id, ym)]

  out <- rbindlist(lapply(split(mon, mon$stn_id), function(g) {
    g <- g[order(ym)]
    r <- data.table(stn_id = g$stn_id[1], ym = g$ym)
    for (sc in CFG$spi_scales)
      r[[paste0("spi", sc)]] <- compute_spi(g$precip, g$ym, sc, b)
    lat <- stn$lat[match(g$stn_id[1], stn$stn_id)]
    if (any(is.finite(g$tmax)) && any(is.finite(g$tmin)) && is.finite(lat)) {
      pet <- pet_hargreaves(g$tmean, g$tmax, g$tmin,
                            as.integer(format(g$ym, "%m")), lat)
      for (sc in CFG$spei_scales)
        r[[paste0("spei", sc)]] <- compute_spei(g$precip - pet, g$ym, sc, b)
    }
    r
  }), fill = TRUE)
  say("[1] SPI/SPEI 패널: ", uniqueN(out$stn_id), "지점 x ", uniqueN(out$ym),
      "개월 (기준기간 ", b[1], "-", b[2], ")")
  out
}

# ---------------------------------------------------------------- 2. 기상 파생
weather_derived <- function(asos) {
  a <- as.data.table(asos)[order(stn_id, date)]
  a[, eh := effective_humidity(rhmean), by = stn_id]
  a[, dry_spell := dry_spell_len(precip), by = stn_id]
  a[, precip_prev7d := {
      s <- shift(precip, 1, type = "lag")
      frollsum(fifelse(is.na(s), 0, s), 7, na.rm = TRUE)
    }, by = stn_id]
  # 당월 강수일수: 단기 날씨(=행동 경로)를 통제하기 위한 변수
  a[, ym := as.Date(format(date, "%Y-%m-01"))]
  wetd <- a[, .(wet_days_mo = sum(precip >= 1, na.rm = TRUE)), by = .(stn_id, ym)]
  say("[2] 기상 파생변수 생성: ", nrow(a), " 지점-일")
  list(daily = a, wetd = wetd)
}

# ---------------------------------------------------------------- 3. 처치 배정
assign_drought <- function(fires, spi_panel) {

  f <- copy(fires)
  f[, ym_prev := as.Date(format(as.Date(date) - as.numeric(format(as.Date(date), "%d")),
                                "%Y-%m-01"))]
  sp <- copy(spi_panel); setnames(sp, "ym", "ym_prev")
  m <- merge(f, sp, by = c("stn_id", "ym_prev"), all.x = TRUE)

  ix <- CFG$compound_index
  m[, compound := fifelse(get(ix) <= CFG$compound_thr, 1L,
                   fifelse(get(ix) >  CFG$control_thr, 0L, NA_integer_))]


  if (all(c("spi3", "spi12") %in% names(m))) {
    thr <- CFG$compound_thr
    m[, short_dry := as.integer(spi3  <= thr)]
    m[, long_dry  := as.integer(spi12 <= thr)]
    m[, drought_type := fifelse(is.na(short_dry) | is.na(long_dry), NA_character_,
        fifelse(short_dry == 1 & long_dry == 1, "both",
        fifelse(short_dry == 1 & long_dry == 0, "short_only",
        fifelse(short_dry == 0 & long_dry == 1, "long_only", "none"))))]
    tb <- m[!is.na(drought_type), .N, by = drought_type][order(-N)]
    say("[3b] 가뭄 유형: ", paste(sprintf("%s %d", tb$drought_type, tb$N),
                                  collapse = " / "))
  }
  say("[3] 처치 배정 (", ix, " <= ", CFG$compound_thr, "): compound ",
      sum(m$compound == 1, na.rm = TRUE), " / control ",
      sum(m$compound == 0, na.rm = TRUE), " / 완충구간 제외 ",
      sum(is.na(m$compound) & !is.na(m[[ix]])), " / 지수결측 ",
      sum(is.na(m[[ix]])))
  m
}

# ---------------------------------------------------------------- 4. 진화자원
attach_resources <- function(d, path = "data_raw/진화자원이력.xlsx") {
  if (!file.exists(path) || !requireNamespace("readxl", quietly = TRUE)) {
    say("[4] 진화자원이력 없음 — heli_n/crew_n/equip_n 생략")
    d[, `:=`(heli_n = NA_real_, crew_n = NA_real_, equip_n = NA_real_,
             stage = NA_character_)]
    return(d)
  }
  source("R/10_sigungu.R", local = TRUE)
  r <- as.data.table(readxl::read_excel(path, col_types = "text"))
  setnames(r, 1:7, c("loc","ig","ct","heli_n","crew_n","equip_n","stage"))
  tk <- strsplit(trimws(r$loc), "[[:space:]]+")
  sd <- vapply(tk, function(z) if (length(z)) z[1] else NA_character_, character(1))
  sg <- vapply(tk, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  r[, sido := split_combined_sido(sd, sg)]
  r[, sigungu_key := resolve_sigungu(sido, sg)]
  r[, date := as.Date(substr(ig, 1, 10))]
  for (v in c("heli_n","crew_n","equip_n"))
    r[[v]] <- suppressWarnings(as.numeric(r[[v]]))
  r <- unique(r[!is.na(date) & !is.na(sigungu_key)],
              by = c("date","sido","sigungu_key"))

  out <- merge(d, r[, .(date, sido, sigungu_key, heli_n, crew_n, equip_n, stage)],
               by = c("date","sido","sigungu_key"), all.x = TRUE)
  out <- unique(out, by = "fire_id")
  sub <- out[year >= 2022]
  say(sprintf("[4] 진화자원 결합: 전체 %.1f%% / 2022년 이후 %.1f%%",
              100 * mean(!is.na(out$heli_n)), 100 * mean(!is.na(sub$heli_n))))
  out
}

# ---------------------------------------------------------------- 5. 저수율
attach_reservoir <- function(d, path = "data_raw/reservoir_daily.csv") {
  if (!file.exists(path)) {
    say("[5] 저수율 없음 — 확보 권장 (소방용수 주장의 물리 근거)")
    d[, `:=`(res_rate = NA_real_, res_anom = NA_real_)]
    return(d)
  }
  r <- fread(path); r[, date := as.Date(date)]
  if (!"sigungu_key" %in% names(r)) {
    source("R/10_sigungu.R", local = TRUE)
    r[, sigungu_key := resolve_sigungu(NA_character_, sigungu)]
  }
  r[, md := format(date, "%m-%d")]
  base <- r[year(date) %between% CFG$baseline,
            .(res_norm = mean(res_rate, na.rm = TRUE)), by = .(sigungu_key, md)]
  r <- merge(r, base, by = c("sigungu_key","md"), all.x = TRUE)
  r[, res_anom := res_rate - res_norm]
  out <- merge(d, unique(r[, .(sigungu_key, date, res_rate, res_anom)],
                         by = c("sigungu_key","date")),
               by = c("sigungu_key","date"), all.x = TRUE)
  say(sprintf("[5] 저수율 결합: 매칭률 %.1f%%", 100 * mean(!is.na(out$res_rate))))
  out
}

# ---------------------------------------------------------------- main
build_panel <- function() {
  dir.create("data_out", showWarnings = FALSE)
  say(strrep("=", 70))
  say("사건 단위 분석 테이블 구축 (산불피해대장 2000-2025)")
  say(strrep("=", 70))

  f <- fread("data_raw/kfs_fires.csv")
  # date 열 복원: fread 가 빈 열을 logical 로 읽거나, 파서 버전에 따라
  # date 가 없을 수 있다. ignition_dt 에서 안전하게 재생성한다.
  if (!"date" %in% names(f) || !inherits(f$date, c("Date","IDate")) ||
      all(is.na(f$date))) {
    if (!"ignition_dt" %in% names(f))
      stop("kfs_fires.csv 에 ignition_dt 가 없습니다. ",
           "R/11_parse_damage_ledger.R 로 다시 파싱하세요.")
    f[, date := as.Date(substr(as.character(ignition_dt), 1, 10))]
    message("[0] date 열을 ignition_dt 에서 재생성했습니다.")
  } else {
    f[, date := as.Date(date)]
  }
  if (all(is.na(f$date))) stop("date 를 만들 수 없습니다. ignition_dt 형식을 확인하세요.")
  for (v in c("year","month")) if (!v %in% names(f))
    f[[v]] <- as.integer(format(f$date, if (v == "year") "%Y" else "%m"))
  if (!"season" %in% names(f))
    f[, season := c("DJF","DJF","MAM","MAM","MAM","JJA","JJA","JJA",
                    "SON","SON","SON","DJF")[month]]
  say("[0] 산불 ", nrow(f), "건 (", min(f$year), "-", max(f$year), ")")
  n0 <- nrow(f); f <- f[!is.na(sigungu_key)]
  if (n0 != nrow(f)) say("    시군구 미해결 제외: ", n0 - nrow(f), "건")

  stn  <- fread("data_raw/station_meta.csv")
  asos <- fread("data_raw/kma_asos_daily.csv"); asos[, date := as.Date(date)]
  say("[0] ASOS ", nrow(asos), "행 / ", uniqueN(asos$stn_id), "지점 (",
      format(min(asos$date), "%Y-%m"), " ~ ", format(max(asos$date), "%Y-%m"), ")")

  smap <- fread("data_out/sigungu_station_map.csv")
  f <- merge(f, smap[, .(sido, sigungu_key, stn_id, stn_km, far_flag, remapped)],
             by = c("sido","sigungu_key"), all.x = TRUE)
  n1 <- nrow(f); f <- f[!is.na(stn_id)]
  if (n1 != nrow(f)) say("[0] 관측소 매핑 실패 제외: ", n1 - nrow(f), "건")

  spi <- build_spi_panel(asos, stn); fwrite(spi, "data_out/spi_panel.csv")
  wx  <- weather_derived(asos)

  d <- assign_drought(f, spi)
  wcols <- intersect(c("stn_id","date","rhmin","rhmean","eh","wsmax","wsmean",
                       "tmean","dry_spell","precip_prev7d"), names(wx$daily))
  d <- merge(d, wx$daily[, ..wcols], by = c("stn_id","date"), all.x = TRUE)
  d[, ym := as.Date(format(date, "%Y-%m-01"))]
  d <- merge(d, wx$wetd, by = c("stn_id","ym"), all.x = TRUE)
  d <- attach_resources(d)
  d <- attach_reservoir(d)

  d[, log_dur := log(pmax(dur_hr, 0.25))]
  d[, cluster := paste0(sigungu_key, "_", year)]

  keep <- intersect(c("fire_id","date","year","month","season","sido","sigungu",
                      "sigungu_key","emd","cause","stn_id","stn_km","far_flag",
                      "remapped","burned_ha","log_burned","large_fire","dur_hr",
                      "log_dur","spi1","spi3","spi6","spi12","spei3","spei6",
                      "spei12","compound","short_dry","long_dry","drought_type",
                      "rhmin","rhmean","eh","wsmax","wsmean","tmean","dry_spell",
                      "precip_prev7d","wet_days_mo","heli_n","crew_n","equip_n",
                      "stage","res_rate","res_anom","cluster"), names(d))
  fwrite(d[, ..keep], "data_out/incident_panel.csv")
  say(strrep("-", 70))
  say("[출력] data_out/incident_panel.csv (", nrow(d), "행 x ", length(keep), "열)")

  sub <- d[!is.na(compound)]
  say("\n기술통계 (compound vs control)")
  desc <- sub[, .(n = .N,
                  burned_mean = round(mean(burned_ha), 2),
                  burned_med  = round(median(burned_ha), 3),
                  p90_ha      = round(quantile(burned_ha, .9), 2),
                  large_pct   = round(mean(large_fire), 4),
                  dur_med     = round(median(dur_hr, na.rm = TRUE), 2)),
              by = compound][order(compound)]
  print(desc); LOG <<- c(LOG, capture.output(print(desc)))
  writeLines(LOG, "data_out/build_report.txt")
  invisible(d)
}

if (sys.nframe() == 0L) build_panel()
