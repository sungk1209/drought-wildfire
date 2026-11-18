# ============================================================================
# 12_station_remap.R — 기록 길이를 고려한 ASOS 지점 재매핑

# 사용:
#   source("R/12_station_remap.R")
#   aud <- audit_stations()                  # data_out/station_audit.csv
#   map <- remap_sigungu(min_years = 28)     # data_out/sigungu_station_map.csv
# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages(library(data.table))

.hav <- function(la1, lo1, la2, lo2) {
  R <- 6371; p1 <- la1 * pi/180; p2 <- la2 * pi/180
  a <- sin((p2 - p1)/2)^2 + cos(p1) * cos(p2) * sin((lo2 - lo1) * pi/360)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

#' 지점별 실제 관측 기간과 결측률
#' @param base SPI 기준기간. 
audit_stations <- function(asos_path = "data_raw/kma_asos_daily.csv",
                           base = c(1991, 2020),
                           out = "data_out/station_audit.csv") {
  stopifnot(file.exists(asos_path))
  a <- fread(asos_path)
  a[, date := as.Date(date)]
  a[, yr := year(date)]

  # 월 단위로 강수 관측일수를 세고, 20일 이상인 달만 필터링링
  mo <- a[, .(ndays = sum(!is.na(precip))), by = .(stn_id, yr, mo = month(date))]
  mo[, valid := as.integer(ndays >= 20)]

  aud <- mo[, .(
    yr_min       = min(yr),
    yr_max       = max(yr),
    n_valid_mo   = sum(valid),
    base_yr_full = uniqueN(yr[valid == 1 & yr >= base[1] & yr <= base[2]]),
    base_mo_full = sum(valid == 1 & yr >= base[1] & yr <= base[2])
  ), by = stn_id]
  aud[, base_coverage := round(base_mo_full / (12 * (base[2] - base[1] + 1)), 3)]
  setorder(aud, -base_yr_full, stn_id)

  dir.create("data_out", showWarnings = FALSE)
  fwrite(aud, out)
  message("지점 ", nrow(aud), "개 감사 완료 -> ", out)
  cat("\n기준기간(", base[1], "-", base[2], ") 확보 연수 분포\n", sep = "")
  print(aud[, .N, by = .(bin = cut(base_yr_full,
        c(-1, 0, 5, 10, 20, 25, 28, 30),
        labels = c("0", "1-5", "6-10", "11-20", "21-25", "26-28", "29-30")))][order(bin)])
  cat("\n기준기간 20년 미만 지점 (SPI-12 부적합)\n")
  print(aud[base_yr_full < 20, .(stn_id, yr_min, yr_max, base_yr_full)])
  invisible(aud)
}

#' 기록이 충분한 지점 중 최근접으로 시군구 재배정
#' @param min_years 
#' @param max_km    
remap_sigungu <- function(centroid_path = "data_raw/sigungu_centroid.csv",
                          meta_path = "data_raw/station_meta.csv",
                          audit_path = "data_out/station_audit.csv",
                          min_years = 28, max_km = 60,
                          out = "data_out/sigungu_station_map.csv") {
  cen <- fread(centroid_path)
  stn <- fread(meta_path)
  aud <- fread(audit_path)

  # 합성자료 
  if (!"sigungu_key" %in% names(cen)) {
    if ("sigungu" %in% names(cen) && exists("strip_sgg")) {
      cen[, sigungu_key := strip_sgg(sigungu)]
      message("sigungu_key 열이 없어 sigungu 에서 생성했습니다.")
    } else {
      stop("중심좌표 파일에 sigungu_key 열이 없습니다: ", centroid_path,
           "\n  -> 실제 sigungu_centroid.csv (230행, 열: sido/sigungu/sigungu_key/lat/lon)",
           " 로 교체하세요.")
    }
  }
  if (nrow(cen) < 100)
    warning("중심좌표가 ", nrow(cen), "행뿐입니다. 합성자료(12행)로 덮여 있을 수 ",
            "있습니다. 실제 파일은 230행입니다.")

  stn <- merge(stn, aud[, .(stn_id, base_yr_full, yr_min, yr_max)],
               by = "stn_id", all.x = TRUE)
  stn <- stn[is.finite(lat) & is.finite(lon)]
  ok <- stn[!is.na(base_yr_full) & base_yr_full >= min_years]
  message(sprintf("지점 %d개 중 기준 통과(>=%d년) %d개",
                  nrow(stn), min_years, nrow(ok)))
  if (nrow(ok) < 10) stop("기준 통과 지점이 너무 적습니다. min_years 를 낮추세요.")

  res <- rbindlist(lapply(seq_len(nrow(cen)), function(i) {
    d_all <- .hav(cen$lat[i], cen$lon[i], stn$lat, stn$lon)
    j_all <- which.min(d_all)                      # 기록 무시 최근접 (비교용)
    d_ok  <- .hav(cen$lat[i], cen$lon[i], ok$lat, ok$lon)
    j_ok  <- which.min(d_ok)                       # 기준 통과 중 최근접
    data.table(
      sido = cen$sido[i], sigungu_key = cen$sigungu_key[i],
      stn_id = ok$stn_id[j_ok], stn_km = round(d_ok[j_ok], 1),
      stn_base_yr = ok$base_yr_full[j_ok],
      nearest_any = stn$stn_id[j_all],
      nearest_km  = round(d_all[j_all], 1))
  }))
  res[, remapped := as.integer(stn_id != nearest_any)]
  res[, km_penalty := round(stn_km - nearest_km, 1)]
  res[, far_flag := as.integer(stn_km > max_km)]

  fwrite(res, out)
  message(sprintf("재매핑: %d개 시군구 중 %d개 변경 (기록 부족 지점 회피)",
                  nrow(res), sum(res$remapped)))
  message(sprintf("  거리 %dkm 초과: %d개 (민감도 분석에서 제외 테스트)",
                  max_km, sum(res$far_flag)))
  cat("\n거리 분포(km)\n"); print(round(quantile(res$stn_km, c(.5,.75,.9,.95,1)), 1))
  if (sum(res$remapped)) {
    cat("\n재매핑으로 거리가 가장 많이 늘어난 시군구\n")
    print(head(res[remapped == 1][order(-km_penalty),
                .(sido, sigungu_key, stn_id, stn_km, nearest_any, nearest_km,
                  km_penalty)], 12))
  }
  invisible(res)
}
