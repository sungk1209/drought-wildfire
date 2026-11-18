# ============================================================================
# 00_download.R — 실제 공공데이터 자동 수집
#
# ※ 인증키를 먼저 발급받아야 한다 (아래 README_DATA.md 참조).
#    Sys.setenv(KMA_KEY = "...", DATAGO_KEY = "...")  또는 .Renviron 에 저장
#
# 주의: 이 스크립트는 인터넷 접속이 되는 사용자 PC에서 실행해야 한다.
#       KFS 산불통계는 파일데이터(CSV)라 브라우저로 1회 내려받는 편이 빠르다.
# ============================================================================

dir.create("data_raw", showWarnings = FALSE)

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))

kma_key    <- Sys.getenv("KMA_KEY")      # https://apihub.kma.go.kr 발급
datago_key <- Sys.getenv("DATAGO_KEY")   # https://www.data.go.kr 발급(Decoding키)

# ---------------------------------------------------------------------------
# 1. 기상청 ASOS 일자료  (apihub, kma_sfcdd3.php = 일자료 기간조회)
#
# 응답 구조
#   START7777
#   #  1. TM     : 관측시각 (KST)        <- help=1 일 때만. 주석은 EUC-KR
#   #  5. WD_MAX : 최대풍향
#   #  6. WS_MAX : 최대풍속 (m/s)
#   ...
#   20250301 136 1.2 ...                 <- 공백구분 데이터
#   7777END
#
# 열 인덱스를 하드코딩하지 않는다. help=1 로 한 번 받아 헤더를 파싱해
# 이름 -> 열번호 매핑을 만든 뒤 재사용한다.
# 이유: WD_MAX/WS_MAX, WD_INS/WS_INS 처럼 풍향과 풍속이 바로 붙어 있어
#       인덱스가 한 칸만 틀려도 풍속 자리에 36방위 풍향값이 조용히 들어간다.
# ---------------------------------------------------------------------------

#' apihub 응답 읽기 (주석이 EUC-KR 이므로 변환)
.read_apihub <- function(u) {
  txt <- try(suppressWarnings(readLines(u, warn = FALSE)), silent = TRUE)
  if (inherits(txt, "try-error")) return(NULL)
  iconv(txt, from = "EUC-KR", to = "UTF-8", sub = "")
}

#' 응답을 헤더줄/데이터줄로 분리 (START7777, 7777END 마커 제거)
.split_response <- function(txt) {
  list(header = grep("^#\\s*[0-9]+\\.\\s+[A-Za-z]", txt, value = TRUE),
       body   = txt[!grepl("^#|^START7777|^7777END", txt) & nzchar(trimws(txt))])
}

#' 헤더줄 -> 이름/열번호 매핑 (네트워크 없이 단위검증 가능)
.parse_header <- function(header_lines) {
  if (!length(header_lines)) return(stats::setNames(integer(0), character(0)))
  stats::setNames(
    as.integer(sub("^#\\s*([0-9]+)\\..*$", "\\1", header_lines)),
    sub("^#\\s*[0-9]+\\.\\s+([A-Za-z0-9_]+).*$", "\\1", header_lines))
}

#' 헤더에서 이름 -> 열번호 매핑 산출
kma_column_map <- function(key = kma_key, stn = 108, tm = "20250301") {
  stopifnot(nzchar(key))
  u <- sprintf(paste0("https://apihub.kma.go.kr/api/typ01/url/kma_sfcdd3.php",
                      "?tm1=%s&tm2=%s&stn=%s&disp=0&help=1&authKey=%s"),
               tm, tm, stn, key)
  txt <- .read_apihub(u)
  if (is.null(txt)) stop("API 응답 실패. 인증키와 네트워크를 확인하세요.")
  h <- .split_response(txt)$header
  if (!length(h))
    stop("헤더를 찾지 못했습니다. 응답 원문 앞부분:\n",
         paste(utils::head(txt, 6), collapse = "\n"))
  map <- .parse_header(h)
  need <- c("TM","STN","RN_DAY","TA_AVG","TA_MAX","TA_MIN",
            "HM_AVG","HM_MIN","WS_AVG","WS_INS")
  miss <- setdiff(need, names(map))
  if (length(miss)) warning("헤더에 없는 항목: ", paste(miss, collapse = ", "))
  message("열 매핑: WS_INS=", map[["WS_INS"]], "(최대순간풍속) ",
          "RN_DAY=", map[["RN_DAY"]], "(일강수량) ",
          "HM_MIN=", map[["HM_MIN"]], "(최저습도)")
  map
}

kma_asos_daily <- function(stn, from = "19900101", to = "20251231",
                           key = kma_key, chunk_years = 5, pause = 0.4) {
  stopifnot(nzchar(key))
  yr1 <- as.integer(substr(from, 1, 4)); yr2 <- as.integer(substr(to, 1, 4))
  brk <- unique(c(seq(yr1, yr2, by = chunk_years), yr2 + 1))
  res <- list()
  for (i in seq_len(length(brk) - 1)) {
    t1 <- sprintf("%04d0101", brk[i])
    t2 <- sprintf("%04d1231", min(brk[i + 1] - 1, yr2))
    u <- sprintf(paste0("https://apihub.kma.go.kr/api/typ01/url/kma_sfcdd3.php",
                        "?tm1=%s&tm2=%s&stn=%s&disp=0&help=0&authKey=%s"),
                 t1, t2, stn, key)
    txt <- .read_apihub(u)
    if (is.null(txt)) { warning("실패: 지점 ", stn, " ", t1); next }
    body <- .split_response(txt)$body
    if (!length(body)) next
    d <- try(read.table(text = body, header = FALSE, fill = TRUE,
                        stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(d, "try-error")) { warning("파싱 실패: 지점 ", stn, " ", t1); next }
    res[[length(res) + 1]] <- d
    Sys.sleep(pause)
  }
  if (!length(res)) return(NULL)
  nc <- max(vapply(res, ncol, integer(1)))
  res <- lapply(res, function(x) { while (ncol(x) < nc) x[[ncol(x) + 1]] <- NA; x })
  do.call(rbind, res)
}

#' 여러 지점을 순회 수집해 표준 스키마로 저장
#' @param stns 생략 시 산불 시군구의 최근접 지점만 자동 선택
#' @param col_map 생략 시 kma_column_map() 으로 헤더에서 자동 산출 (권장)
fetch_all_asos <- function(stns = NULL, from = "19730101", to = "20251231",
                           col_map = NULL, key = kma_key,
                           cache_dir = "data_raw/asos_cache",
                           refetch = FALSE,
                           out = "data_raw/kma_asos_daily.csv") {
  if (is.null(stns)) {
    stns <- all_asos_stations()
    if (length(stns) < 50)
      warning("지점이 ", length(stns), "개뿐입니다. station_meta.csv 가 ",
              "합성자료이거나 불완전할 수 있습니다 (실제 ASOS 는 약 96개).")
    message("지점 ", length(stns), "개 수집")
  }
  if (as.integer(substr(from, 1, 4)) > 1990)
    warning("SPI-12 의 1991-2020 기준기간을 채우려면 from 을 19900101 이전으로 두세요.")
  if (is.null(col_map)) col_map <- kma_column_map(key = key, stn = stns[1])

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  g_of <- function(d, k) {
    if (!k %in% names(col_map)) return(rep(NA, nrow(d)))
    j <- col_map[[k]]
    if (j > ncol(d)) return(rep(NA, nrow(d)))
    d[[j]]
  }

  for (s in stns) {
    cf <- file.path(cache_dir, sprintf("asos_%s.csv", s))
    # 이어받기: 이미 받은 지점은 건너뛴다 (refetch=TRUE 로 강제 재수집)
    if (!refetch && file.exists(cf)) {
      message(sprintf("  지점 %s (%d/%d): 캐시 사용", s, match(s, stns), length(stns)))
      next
    }
    message(sprintf("  지점 %s (%d/%d) 수집 중...", s, match(s, stns), length(stns)))
    d <- kma_asos_daily(s, from, to, key = key)
    if (is.null(d) || !nrow(d)) { warning("수집 실패: 지점 ", s); next }
    one <- data.frame(
      stn_id = as.integer(g_of(d, "STN")),
      date   = as.Date(as.character(g_of(d, "TM")), format = "%Y%m%d"),
      precip = suppressWarnings(as.numeric(g_of(d, "RN_DAY"))),
      tmean  = suppressWarnings(as.numeric(g_of(d, "TA_AVG"))),
      tmax   = suppressWarnings(as.numeric(g_of(d, "TA_MAX"))),
      tmin   = suppressWarnings(as.numeric(g_of(d, "TA_MIN"))),
      rhmean = suppressWarnings(as.numeric(g_of(d, "HM_AVG"))),
      rhmin  = suppressWarnings(as.numeric(g_of(d, "HM_MIN"))),
      wsmean = suppressWarnings(as.numeric(g_of(d, "WS_AVG"))),
      wsmax  = suppressWarnings(as.numeric(g_of(d, "WS_INS"))),  # 최대순간풍속
      stringsAsFactors = FALSE)
    one <- one[!is.na(one$date), ]
    write.csv(one, cf, row.names = FALSE)
    message(sprintf("      %d행 저장 (%s ~ %s)", nrow(one),
                    format(min(one$date), "%Y-%m"), format(max(one$date), "%Y-%m")))
  }

  # 캐시를 모아 하나로 합친다
  files <- list.files(cache_dir, pattern = "^asos_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("수집된 자료가 없습니다.")
  d <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  d$date <- as.Date(d$date)
  for (v in c("precip","tmean","tmax","tmin","rhmean","rhmin","wsmean","wsmax"))
    d[[v]][!is.na(d[[v]]) & d[[v]] <= -9] <- NA     # 기상청 결측코드
  d$precip[is.na(d$precip)] <- 0                    # 무강수는 공란인 경우가 있음
  d <- d[!is.na(d$date) & !is.na(d$stn_id), ]
  d <- d[order(d$stn_id, d$date), ]
  write.csv(d, out, row.names = FALSE)
  message("저장: ", out, " (", nrow(d), "행, ", length(unique(d$stn_id)), "지점, ",
          format(min(d$date), "%Y-%m"), " ~ ", format(max(d$date), "%Y-%m"), ")")

  cat("\n지점별 관측기간 (앞 10개)\n")
  tb <- do.call(rbind, lapply(split(d, d$stn_id), function(g) data.frame(
    stn_id = g$stn_id[1], from = format(min(g$date), "%Y-%m"),
    to = format(max(g$date), "%Y-%m"), n = nrow(g),
    precip_ok = round(mean(!is.na(g$precip)), 3))))
  print(head(tb[order(tb$from), ], 10), row.names = FALSE)
  invisible(d)
}

# ---------------------------------------------------------------------------
# 2. ASOS 지점 메타 (지점번호/위경도)
#    apihub: 지상관측 지점정보, 또는 지상기상월보 getSfcStnLstTbl
# ---------------------------------------------------------------------------
fetch_station_meta <- function(key = kma_key) {
  url <- sprintf(
    "https://apihub.kma.go.kr/api/typ01/url/stn_inf.php?inf=SFC&stn=&tm=20250101&authKey=%s", key)
  txt <- .read_apihub(url)
  if (is.null(txt)) {
    message("자동 수집 실패. 기상자료개방포털에서 지점정보를 내려받아 ",
            "data_raw/station_meta.csv (stn_id,name,lat,lon) 로 저장하세요.")
    return(invisible(NULL))
  }
  body <- .split_response(txt)$body
  d <- read.table(text = body, header = FALSE, stringsAsFactors = FALSE,
                  fill = TRUE)
  meta <- data.frame(stn_id = as.integer(d[[1]]),
                     lon = as.numeric(d[[2]]), lat = as.numeric(d[[3]]))
  meta <- meta[is.finite(meta$lat) & is.finite(meta$lon), ]
  write.csv(meta, "data_raw/station_meta.csv", row.names = FALSE)
  invisible(meta)
}

# ---------------------------------------------------------------------------
# 2b. 수집할 지점 목록 결정
#
# 지점을 손으로 적지 말 것. 아래 둘 중 하나로 '유도'한다.
#   all_asos_stations()      : station_meta.csv 의 전 지점 (권장)
#   stations_for_fires()     : 산불 데이터에 등장하는 시군구의 최근접 지점만
#
# ASOS 는 96개 지점뿐이라 전 지점을 받아도 30~60분이면 끝난다.
# 전 지점을 받아두면 나중에 대상 시군구를 넓힐 때 다시 받지 않아도 된다.
# ---------------------------------------------------------------------------
#' 지점 메타 읽기 (컬럼명 차이를 흡수한다)
#' 기상자료개방포털 내려받기 파일은 station_id / 위도 / 경도 등으로 나오므로
#' stn_id / lat / lon 으로 표준화한 뒤 사용한다.
read_station_meta <- function(meta = "data_raw/station_meta.csv") {
  stopifnot(file.exists(meta))
  s <- try(read.csv(meta, check.names = FALSE), silent = TRUE)
  if (inherits(s, "try-error"))
    s <- read.csv(meta, check.names = FALSE, fileEncoding = "CP949")
  nm <- names(s)
  pick <- function(cands, pat) {
    hit <- nm[tolower(nm) %in% tolower(cands)]
    if (length(hit)) return(hit[1])
    hit <- nm[grepl(pat, nm, ignore.case = TRUE)]
    if (length(hit)) hit[1] else NA_character_
  }
  c_id  <- pick(c("stn_id","station_id","stn","지점","지점번호"), "stn|station|지점")
  c_lat <- pick(c("lat","latitude","위도"), "^lat|위도")
  c_lon <- pick(c("lon","lng","longitude","경도"), "^lon|^lng|경도")
  if (any(is.na(c(c_id, c_lat, c_lon))))
    stop("지점 메타에서 지점번호/위도/경도 열을 찾지 못했습니다. 열 이름: ",
         paste(nm, collapse = ", "))
  out <- data.frame(
    stn_id = suppressWarnings(as.integer(s[[c_id]])),
    lat    = suppressWarnings(as.numeric(s[[c_lat]])),
    lon    = suppressWarnings(as.numeric(s[[c_lon]])))
  out <- out[is.finite(out$stn_id) & is.finite(out$lat) & is.finite(out$lon), ]
  # 한국 영역 밖 좌표(오기/결측코드) 제거
  out <- out[out$lat > 32 & out$lat < 40 & out$lon > 124 & out$lon < 133, ]
  out <- out[!duplicated(out$stn_id), ]
  out[order(out$stn_id), ]
}

all_asos_stations <- function(meta = "data_raw/station_meta.csv") {
  m <- read_station_meta(meta)
  message("지점 메타: ", nrow(m), "개 (", min(m$stn_id), "-", max(m$stn_id), ")")
  m$stn_id
}

#' 표준 스키마로 다시 저장해 이후 단계(12_station_remap 등)가 바로 읽게 한다.
normalize_station_meta <- function(meta = "data_raw/station_meta.csv") {
  m <- read_station_meta(meta)
  write.csv(m, meta, row.names = FALSE)
  message("표준화 저장: ", meta, " (", nrow(m), "개 지점, 열: stn_id, lat, lon)")
  invisible(m)
}

#' 산불 데이터의 시군구 -> 최근접 ASOS 지점 집합
#' @return list(stn_ids, mapping)  mapping 은 거리(km) 포함, 점검용
stations_for_fires <- function(fires = "data_raw/kfs_fires.csv",
                               centroid = "data_raw/sigungu_centroid.csv",
                               meta = "data_raw/station_meta.csv",
                               max_km = 50) {
  f   <- read.csv(fires, stringsAsFactors = FALSE)
  cen <- read.csv(centroid, stringsAsFactors = FALSE)
  stn <- read.csv(meta)

  hav <- function(la1, lo1, la2, lo2) {
    R <- 6371; p1 <- la1 * pi/180; p2 <- la2 * pi/180
    a <- sin((p2 - p1)/2)^2 + cos(p1) * cos(p2) * sin((lo2 - lo1) * pi/360)^2
    2 * R * asin(pmin(1, sqrt(a)))
  }

  want <- unique(f$sigungu[!is.na(f$sigungu)])
  cen  <- cen[cen$sigungu %in% want, ]
  miss <- setdiff(want, cen$sigungu)
  if (length(miss))
    warning("중심좌표 없는 시군구 ", length(miss), "개 (행정구역 crosswalk 점검): ",
            paste(utils::head(miss, 10), collapse = ", "))

  map <- do.call(rbind, lapply(seq_len(nrow(cen)), function(i) {
    d <- hav(cen$lat[i], cen$lon[i], stn$lat, stn$lon)
    j <- which.min(d)
    data.frame(sido = cen$sido[i], sigungu = cen$sigungu[i],
               stn_id = stn$stn_id[j], stn_km = round(d[j], 1))
  }))
  far <- sum(map$stn_km > max_km)
  message("대상 시군구 ", nrow(map), "개 -> 필요한 ASOS 지점 ",
          length(unique(map$stn_id)), "개")
  if (far > 0) message("  ", max_km, "km 초과 ", far,
                       "건 — 민감도 분석에서 제외 테스트 대상")
  write.csv(map, "data_raw/sigungu_station_map.csv", row.names = FALSE)
  list(stn_ids = sort(unique(map$stn_id)), mapping = map)
}

# ---------------------------------------------------------------------------
# 3. 한국농어촌공사 전국 저수지 일별 저수율 (시군구별)
#    파일데이터: https://www.data.go.kr/data/15113475/fileData.do
#    오픈API   : https://www.data.go.kr/data/15099919/openapi.do (county 파라미터)
# ---------------------------------------------------------------------------
krc_reservoir <- function(county, from, to, key = datago_key) {
  # 저수지 위치(시군구) 조회는 최대 31일 -> 월 단위로 끊어서 호출
  stopifnot(nzchar(key))
  url <- sprintf(paste0("https://apis.data.go.kr/B552149/reservoirlevel/",
                        "getReservoirLevelList?serviceKey=%s&county=%s",
                        "&stdt=%s&eddt=%s&numOfRows=1000&pageNo=1&type=json"),
                 utils::URLencode(key, reserved = TRUE),
                 utils::URLencode(county, reserved = TRUE), from, to)
  try(readLines(url, warn = FALSE), silent = TRUE)
}

# ---------------------------------------------------------------------------
# 4. 산림청 산불통계 (파일데이터 — 브라우저 다운로드 권장)
# ---------------------------------------------------------------------------
kfs_note <- function() {
  cat("
산림청 산불통계데이터 (파일데이터, 로그인 없이 다운로드 가능)
  https://www.data.go.kr/data/15121380/fileData.do
  포함 컬럼: 발생일시(년/월/일/시간/요일), 진화종료시간(년/월/일/시간),
             발생장소(관서/시도/시군구/읍면/동리), 발생원인(구분/세부/기타),
             피해면적(합계)
  -> data_raw/kfs_fires_raw.csv 로 저장한 뒤 normalize_kfs() 실행

대안(오픈API): https://www.data.go.kr/data/3070842/openapi.do
")
}

#' 산림청 원자료 -> 표준 스키마 변환
normalize_kfs <- function(path = "data_raw/kfs_fires_raw.csv") {
  d <- read.csv(path, fileEncoding = "EUC-KR", stringsAsFactors = FALSE)
  nm <- names(d)
  pick <- function(pat) nm[grepl(pat, nm)][1]
  gi <- function(pat) suppressWarnings(as.integer(d[[pick(pat)]]))

  y <- gi("발생.*년"); m <- gi("발생.*월"); dd <- gi("발생.*일$|발생.*일[^시]")
  hh <- gi("발생.*시")
  cy <- gi("진화.*년"); cm <- gi("진화.*월"); cd <- gi("진화.*일")
  ch <- gi("진화.*시")

  mk <- function(Y, M, D, H) as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00",
                                                Y, M, D, ifelse(is.na(H), 12, H)),
                                        tz = "Asia/Seoul")
  out <- data.frame(
    ignition_dt = mk(y, m, dd, hh),
    contain_dt  = mk(cy, cm, cd, ch),
    sido        = d[[pick("시도")]],
    sigungu     = d[[pick("시군구")]],
    cause       = d[[pick("발생원인.*구분|원인구분")]],
    burned_ha   = suppressWarnings(as.numeric(d[[pick("피해면적")]])),
    stringsAsFactors = FALSE)
  out$fire_id <- sprintf("F%06d", seq_len(nrow(out)))
  write.csv(out, "data_raw/kfs_fires.csv", row.names = FALSE)
  message("저장: data_raw/kfs_fires.csv (", nrow(out), "건)")
  invisible(out)
}
