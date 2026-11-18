# ============================================================================
# run_all.R — 복합 가뭄-산불 분석 파이프라인
#
# 작업 디렉터리: 저장소 최상위 폴더
#   Rscript R/run_all.R
#
# 사전 준비: data_raw/ 에 산불피해대장, 진화자원이력, ASOS 일자료, 지점정보
#   (ASOS 는 R/00_download.R 의 fetch_all_asos() 로 수집, README 참조)
# ============================================================================

invisible(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")))

need <- c("data.table", "sandwich", "lmtest", "quantreg", "MASS", "readxl",
          "tidyverse", "patchwork", "scales")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("패키지 설치 필요: install.packages(c(",
                       paste0('"', miss, '"', collapse = ", "), "))")

for (p in c("data_out", "figs")) dir.create(p, showWarnings = FALSE)
step <- function(x) cat("\n\n##########  ", x, "  ##########\n\n")

# ---------------------------------------------------------------- 1. 산불
step("1. 산불피해대장 파싱 (2000-2025) + 동일 복합체 중복 병합")
source("R/11_parse_damage_ledger.R")
LEDGER <- Sys.glob("data_raw/*산불피해대장*.xls*")
if (!length(LEDGER)) stop("산불피해대장 엑셀이 없습니다: data_raw/ 에 넣어주세요")
invisible(parse_damage_ledger(LEDGER[1]))
invisible(dedupe_complex(min_ha = 1000, method = "max"))

# ---------------------------------------------------------------- 2. 기상
step("2. ASOS 자료 확인")
req <- c("data_raw/kma_asos_daily.csv", "data_raw/station_meta.csv")
no <- req[!file.exists(req)]
if (length(no)) stop("ASOS 자료 없음:\n  ", paste(no, collapse = "\n  "),
                     "\n  R/00_download.R 의 fetch_all_asos() 로 먼저 수집하세요")
if (!file.exists("data_raw/sigungu_centroid.csv"))
  file.copy("sigungu_centroid.csv", "data_raw/sigungu_centroid.csv")

# ---------------------------------------------------------------- 3. 지점
step("3. 지점 감사 + 기록 길이 기반 재매핑")
source("R/12_station_remap.R")
invisible(audit_stations())
invisible(remap_sigungu(min_years = 28))

# ---------------------------------------------------------------- 4. 결합
step("4. 사건 단위 분석 테이블 구축")
source("R/13_build_panel.R")
invisible(build_panel())

# ---------------------------------------------------------------- 5. 분석
step("5. 발화 조건부 증폭 분석")
source("R/14_analysis_q2.R")
run_analysis(boot_mediation = 2000, boot_quantile = 300)

# ---------------------------------------------------------------- 6. 표·그림
step("6. 표 및 분위회귀 결과 (data_out/tables)")
source("R/16_make_tables_figures.R")
make_all_outputs()

step("7. 리비전 추정치 (data_out/tables/revision)")
source("R/17_figures_revision_rev.r")

step("8. 최종 그림 (figs/revision_final)")
source("R/figs/Figure3_redraw.R")
source("R/figs/Figure5_redraw.R")
source("R/figs/FigureS1_redraw.R")

cat("\n완료. data_out/analysis_report.txt 와 figs/revision_final/ 을 확인하세요.\n")
