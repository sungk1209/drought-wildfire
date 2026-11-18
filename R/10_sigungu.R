# ============================================================================
# 10_sigungu.R 
#

# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))
suppressPackageStartupMessages(library(data.table))

CENTROID_PATH <- "data_raw/sigungu_centroid.csv"
GWANGJU_GU <- c("동", "서", "남", "북", "광산")
INCHEON_ALIAS <- c("검단" = "서", "영종" = "중", "서해" = "옹진")

# 과거 행정구역명 + 원자료 오타 (산불피해대장 2000~2025 에서 확인된 것)
#   통합/폐지: 마산·진해 -> 창원(2010), 연기 -> 세종(2012),
#              북제주·남제주 -> 제주·서귀포(2006), 청원 -> 청주(2014, 좌표에 별도 존재)
#   오타     : 체천 -> 제천, 굥주 -> 공주, 광화 -> 강화
# "시도|시군구" 형태를 키로 쓴다. 시도 무관이면 "*|시군구".
HIST_ALIAS <- c(
  "경남|마산"   = "창원",
  "경남|진해"   = "창원",
  "경남|마산합포" = "창원",
  "경남|마산회원" = "창원",
  "부산|진"     = "부산진",
  "충남|연기"   = "세종",
  "제주|북제주" = "제주",
  "제주|남제주" = "서귀포",
  "제주|화북"   = "제주",
  "인천|광화"   = "강화",
  "충북|체천"   = "제천",
  "충남|굥주"   = "공주",
  "경북|화북"   = "상주"
)

.centroid <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      if (!file.exists(CENTROID_PATH))
        stop("중심좌표 파일이 없습니다: ", CENTROID_PATH)
      cache <<- fread(CENTROID_PATH)
    }
    cache
  }
})

#' "전남광주" 합칭을 광주/전남으로 분리
#' @param sido_tok 첫 토큰, @param sgg_tok 둘째 토큰
split_combined_sido <- function(sido_tok, sgg_tok) {
  s <- trimws(as.character(sido_tok))
  fifelse(s == "전남광주",
          fifelse(trimws(as.character(sgg_tok)) %in% GWANGJU_GU, "광주", "전남"),
          s)
}

#' 시군구명 표준화 -> sigungu_key (접미사 없는 표준명). 실패 시 NA.
resolve_sigungu <- function(sido, sigungu) {
  cen <- .centroid()
  allk <- unique(cen$sigungu_key)
  by_sido <- split(cen$sigungu_key, cen$sido)

  sd <- trimws(as.character(sido))
  x  <- gsub("\\s+", "", as.character(sigungu))
  x[is.na(x) | x %in% c("", "undefined", "NA", "-", "nan")] <- NA_character_

  strip1 <- function(z) sub("(시|군|구)$", "", z)

  vapply(seq_along(x), function(i) {
    if (is.na(x[i])) return(NA_character_)
    if (!is.na(sd[i]) && sd[i] == "세종") return("세종")
    # 0) 과거 행정구역명/오타 별칭을 가장 먼저 적용
    ak <- paste0(ifelse(is.na(sd[i]), "*", sd[i]), "|", x[i])
    if (ak %in% names(HIST_ALIAS)) return(unname(HIST_ALIAS[ak]))
    ak2 <- paste0("*|", x[i])
    if (ak2 %in% names(HIST_ALIAS)) return(unname(HIST_ALIAS[ak2]))
    cand <- if (!is.na(sd[i]) && !is.null(by_sido[[sd[i]]])) by_sido[[sd[i]]] else allk
    # 1) 원본 그대로, 2) 접미사 1개 제거  (양구 -> 양 으로 잘리는 것 방지)
    for (t in c(x[i], strip1(x[i]))) if (t %in% cand) return(t)
    # 3) 인천 지역 별칭
    if (!is.na(sd[i]) && sd[i] == "인천" && x[i] %in% names(INCHEON_ALIAS))
      return(unname(INCHEON_ALIAS[x[i]]))
    # 4) 최장 접두 일치 (용인처인 -> 용인, 천안동남 -> 천안)
    for (t in c(x[i], strip1(x[i]))) {
      pref <- cand[nchar(cand) >= 2 & startsWith(t, cand)]
      if (length(pref)) return(pref[which.max(nchar(pref))])
    }
    # 5) 시도 무시하고 전국 대조
    for (t in c(x[i], strip1(x[i]))) if (t %in% allk) return(t)
    NA_character_
  }, character(1))
}

#' 중심좌표(lat/lon) 결합. sigungu_key 가 없으면 먼저 만든다.
attach_centroid <- function(d, sido_col = "sido", sgg_col = "sigungu") {
  d <- as.data.table(d)
  if (!"sigungu_key" %in% names(d))
    d[, sigungu_key := resolve_sigungu(get(sido_col), get(sgg_col))]
  cen <- .centroid()[, .(sido, sigungu_key, sigungu_std = sigungu, lat, lon)]
  # 시도+키로 먼저 붙이고, 실패분은 키만으로 재시도
  out <- merge(d, cen, by.x = c(sido_col, "sigungu_key"),
               by.y = c("sido", "sigungu_key"), all.x = TRUE, sort = FALSE)
  if (any(is.na(out$lat))) {
    cen2 <- unique(cen[, .(sigungu_key, lat2 = lat, lon2 = lon)], by = "sigungu_key")
    out <- merge(out, cen2, by = "sigungu_key", all.x = TRUE, sort = FALSE)
    out[is.na(lat), `:=`(lat = lat2, lon = lon2)]
    out[, c("lat2", "lon2") := NULL]
  }
  miss <- sum(is.na(out$lat))
  message(sprintf("중심좌표 결합: %d건 중 실패 %d건 (%.2f%%)",
                  nrow(out), miss, 100 * miss / nrow(out)))
  if (miss > 0) {
    bad <- out[is.na(lat), .N, by = .(get(sgg_col))][order(-N)]
    message("  실패 상위: ",
            paste(utils::head(sprintf("%s(%d)", bad[[1]], bad$N), 8), collapse = ", "))
  }
  out[]
}
