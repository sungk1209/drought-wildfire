# ============================================================================
# 11_parse_damage_ledger.R 
#
#
# 사용:
#   source("R/11_parse_damage_ledger.R")
#   parse_damage_ledger("data_raw/산불피해대장_2000년-2025년.xlsx")
# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages(library(data.table))

parse_damage_ledger <- function(path, sheet = "합본", skip = 4,
                                out = "data_raw/kfs_fires.csv") {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("readxl 이 필요합니다: install.packages(\"readxl\")")
  d <- readxl::read_excel(path, sheet = sheet, skip = skip,
                          col_names = FALSE, col_types = "text")
  setDT(d)
  if (ncol(d) < 14) stop("열이 14개 미만입니다. skip/sheet 를 확인하세요.")
  setnames(d, 1:14, c("iy","im","id","it","cy","cm","cd","ct",
                      "sido","sigungu","emd","ri","cause","area"))
  d <- d[!is.na(iy) & grepl("^[0-9]{4}$", trimws(iy))]
  message("원자료 ", nrow(d), "행")

  # 시각 파싱. 세 형태를 모두 처리
.
  .time_hm <- function(T) {
    tt <- trimws(as.character(T))
    hh <- suppressWarnings(as.integer(sub("^(\\d{1,2}):.*$", "\\1", tt)))
    mi <- suppressWarnings(as.integer(sub("^\\d{1,2}:(\\d{2}).*$", "\\1", tt)))
    frac_like <- !grepl(":", tt) & grepl("^0?\\.[0-9]+$|^[01](\\.[0-9]+)?$", tt)
    if (any(frac_like, na.rm = TRUE)) {
      fv <- suppressWarnings(as.numeric(tt[frac_like]))
      secs <- round(fv * 86400)
      hh[frac_like] <- (secs %/% 3600) %% 24
      mi[frac_like] <- (secs %% 3600) %/% 60
    }
    list(h = hh, m = mi)
  }
  .mk <- function(Y, M, D, T) {
    y  <- suppressWarnings(as.integer(Y)); mo <- suppressWarnings(as.integer(M))
    dd <- suppressWarnings(as.integer(D))
    th <- .time_hm(T); hh <- th$h; mi <- th$m
    hh[is.na(hh) | hh > 23] <- 12L; mi[is.na(mi) | mi > 59] <- 0L
    ok <- !is.na(y) & !is.na(mo) & !is.na(dd) & mo >= 1 & mo <= 12 & dd >= 1 & dd <= 31
    o <- rep(as.POSIXct(NA), length(y))
    o[ok] <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00",
                                y[ok], mo[ok], dd[ok], hh[ok], mi[ok]),
                        tz = "Asia/Seoul")
    o
  }

  f <- data.table(
    ignition_dt = .mk(d$iy, d$im, d$id, d$it),
    contain_dt  = .mk(d$cy, d$cm, d$cd, d$ct),
    sido        = trimws(d$sido),
    sigungu     = trimws(d$sigungu),
    emd         = trimws(d$emd),
    ri          = trimws(d$ri),
    cause       = trimws(d$cause),
    burned_ha   = suppressWarnings(as.numeric(d$area)))

  n0 <- nrow(f)
  f <- f[!is.na(ignition_dt) & !is.na(burned_ha) & burned_ha > 0 &
           !is.na(sigungu) & nzchar(sigungu)]
  message("필수항목 결측/0면적 제외: ", n0, " -> ", nrow(f))

  f[, dur_hr := as.numeric(difftime(contain_dt, ignition_dt, units = "hours"))]
  bad <- !is.na(f$dur_hr) & (f$dur_hr < 0 | f$dur_hr > 24 * 30)
  f[bad, dur_hr := NA_real_]
  message("진화시간 이상치 결측화: ", sum(bad), "건 / 확보율 ",
          sprintf("%.1f%%", 100 * mean(!is.na(f$dur_hr))))

  f[, cause_raw := cause]
  f[, cause := gsub("[[:space:]]", "", cause)]
  f[, cause := sub("추정$", "", cause)]
  f[, cause := fcase(
      grepl("입산자|등산객", cause),            "입산자실화",
      grepl("담뱃불|담배불", cause),            "담뱃불실화",
      grepl("논두렁|밭두렁|논밭두렁", cause),   "논밭두렁소각",
      grepl("쓰레기", cause),                   "쓰레기소각",
      grepl("농산", cause),                     "농산부산물소각",
      grepl("성묘", cause),                     "성묘객실화",
      grepl("어린이|불장난", cause),            "어린이불장난",
      grepl("방화", cause),                     "방화",
      grepl("건축물|주택|화재비화", cause),     "건축물실화",
      grepl("산업|용접|공사", cause),           "산업현장실화",
      grepl("낙뢰|기상", cause),                "기상요인",
      grepl("군사|훈련", cause),                "군사훈련",
      grepl("미상|조사중|불명", cause),         "원인미상",
      default = "기타")]
  say_cause <- f[, .N, by = cause][order(-N)]
  message("발생원인 정규화: ", uniqueN(f$cause_raw), "개 -> ",
          uniqueN(f$cause), "개 범주")

  f[, `:=`(date  = as.Date(ignition_dt),
           year  = as.integer(format(ignition_dt, "%Y")),
           month = as.integer(format(ignition_dt, "%m")))]
  f[, season := c("DJF","DJF","MAM","MAM","MAM","JJA","JJA","JJA",
                  "SON","SON","SON","DJF")[month]]
  f[, log_burned := log(burned_ha + 0.01)]
  f[, large_fire := as.integer(burned_ha >= 100)]
  setorder(f, ignition_dt)
  f[, fire_id := sprintf("F%06d", .I)]

  if (file.exists("R/10_sigungu.R")) {
    source("R/10_sigungu.R", local = TRUE)
    f[, sigungu_key := resolve_sigungu(sido, sigungu)]
    message("시군구 표준화 실패: ", sum(is.na(f$sigungu_key)), "건")
  }

  keep <- intersect(c("fire_id","ignition_dt","contain_dt","date","year","month",
                      "season","sido","sigungu","sigungu_key","emd","ri","cause",
                      "cause_raw","burned_ha","log_burned","large_fire","dur_hr"),
                    names(f))
  fo <- copy(f[, ..keep])
  fo[, ignition_dt := format(ignition_dt, "%Y-%m-%d %H:%M:%S", tz = "Asia/Seoul")]
  fo[, contain_dt  := format(contain_dt,  "%Y-%m-%d %H:%M:%S", tz = "Asia/Seoul")]
  fwrite(fo, out)
  message("저장: ", out, " (", nrow(fo), "건, ",
          min(f$year), "-", max(f$year), ")")

  cat("\n연도별 요약\n")
  print(f[, .(n = .N,
              dur_ok = round(mean(!is.na(dur_hr)), 3),
              burned_med = round(median(burned_ha), 2),
              burned_sum = round(sum(burned_ha), 1),
              n30  = sum(burned_ha >= 30),
              n100 = sum(burned_ha >= 100)), by = year][order(year)])
  invisible(f)
}

#' 진화자원이력(2022-2026) 결합: 헬기/인력/장비/산불단계
#' 발생일 + 시도 + 시군구_key 로 매칭한다.
attach_suppression <- function(fires_path = "data_raw/kfs_fires.csv",
                               res_path = "data_raw/진화자원이력.xlsx",
                               out = "data_raw/kfs_fires_with_resources.csv") {
  if (!file.exists(res_path)) { message("진화자원이력 없음 — 생략"); return(invisible(NULL)) }
  source("R/10_sigungu.R", local = TRUE)
  f <- fread(fires_path)
  r <- as.data.table(readxl::read_excel(res_path, col_types = "text"))
  setnames(r, 1:7, c("loc","ig","ct","heli_n","crew_n","equip_n","stage"))

  tk <- strsplit(trimws(r$loc), "[[:space:]]+")
  sd <- vapply(tk, function(z) if (length(z)) z[1] else NA_character_, character(1))
  sg <- vapply(tk, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  r[, sido := split_combined_sido(sd, sg)]
  r[, sigungu_key := resolve_sigungu(sido, sg)]
  r[, date := as.Date(substr(ig, 1, 10))]
  for (v in c("heli_n","crew_n","equip_n"))
    r[[v]] <- suppressWarnings(as.numeric(r[[v]]))

  f[, date := as.Date(date)]
  m <- merge(f, r[, .(date, sido, sigungu_key, heli_n, crew_n, equip_n, stage)],
             by = c("date", "sido", "sigungu_key"), all.x = TRUE)
  m <- unique(m, by = "fire_id")
  message(sprintf("진화자원 결합: %d건 중 %d건 매칭 (%.1f%%)",
                  nrow(m), sum(!is.na(m$heli_n)), 100 * mean(!is.na(m$heli_n))))
  sub <- m[year >= 2022]
  message(sprintf("  2022년 이후만: %d건 중 %d건 (%.1f%%)",
                  nrow(sub), sum(!is.na(sub$heli_n)), 100 * mean(!is.na(sub$heli_n))))
  fwrite(m, out)
  message("저장: ", out)
  invisible(m)
}


# ============================================================================
# dedupe_complex() — 동일 산불 복합체의 중복 등재 병합
#
# 원자료는 하나의 대형 산불이 여러 읍면으로 확산된 경우를 각각 등재하기도 한다.
#   2025-03-22 의성: 안평 46,575 ha / 안계 52,707 ha  (진화 종료 시각 동일)
#   2000-04-15 삼척: 동해안 산불 복합체 13,343 / 4,053 ha
# 그대로 두면 피해면적이 이중 계상되어 상위 꼬리가 왜곡된다.
#
# 판정: 같은 시군구 + 같은 진화종료시각 + 해당 그룹에 min_ha 이상이 2건 이상
# 병합: 발화는 최초 시각, 피해면적은 method 에 따라
#   "max"  (기본) 가장 큰 기록을 복합체 피해면적으로 본다 — 이중계상 방지
#   "sum"  합산. 읍면별 독립 집계가 확실할 때만
#   "flag" 병합하지 않고 dup_flag 만 남긴다 (민감도 분석용)
#
# min_ha 를 낮추면 실제로 별개인 중소형 산불까지 묶일 수 있다. 출력된 후보를
# 눈으로 확인한 뒤 결정할 것 (기본 1000 ha: 초대형 복합체만 대상)
# ============================================================================
dedupe_complex <- function(path = "data_raw/kfs_fires.csv",
                           min_ha = 1000, method = c("max","sum","flag"),
                           out = path) {
  method <- match.arg(method)
  f <- fread(path)
  f[, contain_dt := as.character(contain_dt)]
  f[, .gkey := paste(sigungu_key, contain_dt)]

  # 같은 그룹에 min_ha 이상 기록이 2건 이상인 경우만 복합체로 본다
  cand <- f[burned_ha >= min_ha & !is.na(contain_dt) & nzchar(contain_dt),
            .(n_big = .N), by = .gkey][n_big > 1]
  f[, dup_flag := as.integer(.gkey %in% cand$.gkey)]

  if (!nrow(cand)) {
    message("중복 복합체 없음 (>= ", format(min_ha, big.mark = ","), " ha 2건 이상)")
    f[, .gkey := NULL][, n_merged := 1L]
    fwrite(f, out); return(invisible(f))
  }

  message("중복 의심 복합체 ", nrow(cand), "건:")
  for (k in cand$.gkey) {
    sub <- f[.gkey == k][order(-burned_ha)]
    message(sprintf("  %-8s 진화 %s : %s ha (%d건)",
                    sub$sigungu_key[1], sub$contain_dt[1],
                    paste(format(sub$burned_ha, big.mark = ","), collapse = " + "),
                    nrow(sub)))
  }

  if (method == "flag") {
    message("병합하지 않고 dup_flag 만 표시했습니다 (", sum(f$dup_flag), "행).")
    f[, .gkey := NULL][, n_merged := 1L]
    fwrite(f, out); return(invisible(f))
  }

  keep <- f[dup_flag == 0L]
  keep[, n_merged := 1L]
  grp  <- f[dup_flag == 1L][order(ignition_dt)]
  agg  <- grp[, {
    ba <- if (method == "max") max(burned_ha) else sum(burned_ha)
    ig <- ignition_dt[1]; ct <- contain_dt[1]
    dh <- as.numeric(difftime(as.POSIXct(ct, tz = "Asia/Seoul"),
                              as.POSIXct(ig, tz = "Asia/Seoul"), units = "hours"))
    c(.SD[1, ], list(burned_ha = ba, log_burned = log(ba + 0.01),
                     large_fire = as.integer(ba >= 100), dur_hr = dh,
                     emd = paste(unique(emd), collapse = "+"),
                     n_merged = .N))
  }, by = .gkey, .SDcols = setdiff(names(f),
       c(".gkey","burned_ha","log_burned","large_fire","dur_hr","emd"))]
  agg[, .gkey := NULL]

  common <- intersect(names(keep), names(agg))
  f2 <- rbind(keep[, ..common], agg[, ..common])
  setorder(f2, ignition_dt)
  f2[, fire_id := sprintf("F%06d", .I)]
  message(sprintf("병합(%s): %d건 -> %d건 (%d건 감소)",
                  method, nrow(f), nrow(f2), nrow(f) - nrow(f2)))
  fwrite(f2, out)
  message("저장: ", out)
  invisible(f2)
}
