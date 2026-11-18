source("R/00_download.R")
kma_key <- "XX"


kma_probe()
fetch_all_asos(stns = all_asos_stations(), from = "19730101", to = "20251231")
file.copy("sigungu_centroid.csv", "data_raw/sigungu_centroid.csv", overwrite = TRUE)
source("R/12_station_remap.R")
remap_sigungu(min_years = 28)



map <- data.table::fread("data_out/sigungu_station_map.csv")
map[stn_km > 60, far_flag := 1L]
data.table::fwrite(map, "data_out/sigungu_station_map.csv")

a<- Sys.glob("data_raw/*산불피해대장*")     

source("R/11_parse_damage_ledger.R")
parse_damage_ledger(a)

source("R/13_build_panel.R"); build_panel()
source("R/14_analysis_q2.R"); run_analysis(boot_mediation = 500, boot_quantile = 100)


# 초대형 산불 sample size 
d <- data.table::fread("data_out/incident_panel.csv")
d[!is.na(compound)][order(-burned_ha)][1:12,
                                       .(year, sido, sigungu, burned_ha, compound, spi3, spi12)]
d[!is.na(compound), .(n = .N, mean_ha = mean(burned_ha),
                      med = median(burned_ha),
                      p99 = quantile(burned_ha, .99),
                      large_pct = mean(large_fire)), by = compound]

# average except for the top 5 evetns
s <- d[!is.na(compound)][burned_ha < quantile(burned_ha, 0.999)]
s[, .(n = .N, mean_ha = mean(burned_ha)), by = compound]

f <- data.table::fread("data_raw/kfs_fires.csv")

# 중복 의심 사례 전수 점검 (2000-2025 전체)
dup <- f[burned_ha >= 100][, .N, by = .(sigungu_key, contain_dt)][N > 1]
print(dup)


#=========================================
  
  source("R/11_parse_damage_ledger.R")
parse_damage_ledger("data_raw/산불피해대장(2000년-2025년) (1).xlsx")
dedupe_complex(min_ha = 1000, method = "max")

source("R/13_build_panel.R"); build_panel()
source("R/14_analysis_q2.R"); run_analysis(boot_mediation = 2000, boot_quantile = 300)
source("R/16_make_tables_figures.R"); make_all_outputs()
