## ----setup, include=FALSE-----------------------------------------------------
library(formatR)
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", fig.show = "hold", time_it = TRUE, dpi = 100)
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 60), tidy = T, collapse = TRUE)
knitr::opts_knit$set(root.dir = '/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')
library(rgl)
library(terra)
library(lme4)
library(MASS)
library(lmerTest)
library(MuMIn)
library(terra)
library(sf)
library(ggplot2)
library(ggeffects)
library(merTools)
library(glmnet)
library(stats)
library(ggcorrplot)
library(RColorBrewer)
library(webshot)
library(kableExtra)
library(formatR)
library(dplyr)
library(stringr)
library(flextable)

library(tmap)

knitr::knit_hooks$set(webgl = hook_webgl)
rgl::setupKnitr(autoprint = TRUE)

terraOptions(
    memfrac = 0.1
)

setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

tmap_options(max.raster = c(plot = 1e7, view = 1e6))


## -----------------------------------------------------------------------------
#| include: false
#| echo: false
#| message: false
#| warning: false
#| cache: true
knitr::purl("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/SOIL CARBON/All_WA/analysis/WA_SOC_controls/All_WA_SOC_Prediction.qmd", output = "/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/SOIL CARBON/All_WA/analysis/WA_SOC_controls/All_WA_SOC_Prediction.R")

#source("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/SOIL CARBON/All_WA/analysis/WA_SOC_controls/All_WA_SOC_Prediction.R")
source("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/SOIL CARBON/All_WA/analysis/WA_SOC_controls/All_WA_SOC_Modeling_Analysis.R")


## -----------------------------------------------------------------------------
wa_dat_format <- wa_dat |> mutate(landscapeClassWetUpl = case_when(WIP >= 0.50 ~ "Wetland",
                                            WIP < 0.50 ~ "Upland"),
                                  landscapeClassWetUplMes = case_when(WIP <= 0.25 | WIP >= 0.75 ~ "Mesic",
                                               .default = landscapeClassWetUpl)) 

wetupltbl <- wa_dat_format |> 
    group_by(site, landscapeClassWetUpl, lower_depth) |> 
    summarise(meanSOC = mean(SOC_stock_spline), sdSOC = sd(SOC_stock_spline), n= n()) |>
    mutate(meanSOC = as.character(round(meanSOC, 2)),
           sdSOC = as.character(round(sdSOC, 2)),
           site = case_match(site, "COL" ~ "Colville",
                             "MAS" ~ "Mashel",
                             "HOH" ~ "Hoh"),
           SOC = paste0(meanSOC, " ± ", sdSOC)) |> 
    select(site, landscapeClassWetUpl, SOC, everything(), -meanSOC, -sdSOC) |> 
    rename("Study Area" = site,
           "Landscape Class" = landscapeClassWetUpl,
           "Depth (cm)" = lower_depth,
           #"Mean SOC (g cm$^{-1}$)" = meanSOC,
           "SOC (g cm$^{-1}$)" = SOC) |> 
    flextable() |> merge_v(j = c("Study Area", "Landscape Class")) |> colformat_num(big.mark = ",") |> 
    set_table_properties(layout = "autofit", align = "center") |>
    align(align = "center", part = "all") |>
    mk_par(
    part = "header", j = 3,
    value = as_paragraph("Mean SOC Stock (Mg ha", as_sup("-1"), ")")
    )

wetuplmestbl <- wa_dat_format |> group_by(site, landscapeClassWetUplMes, lower_depth) |> 
    summarise(meanSOC = mean(SOC_stock_spline), sdSOC = sd(SOC_stock_spline), n = n()) |>
    mutate(meanSOC = as.character(round(meanSOC, 2)),
           sdSOC = as.character(round(sdSOC, 2)),
           site = case_match(site, "COL" ~ "Colville",
                             "MAS" ~ "Mashel",
                             "HOH" ~ "Hoh"),
           SOC = paste0(meanSOC, " ± ", sdSOC)) |> 
    select(site, landscapeClassWetUplMes, SOC, everything(), -meanSOC, -sdSOC) |> 
    rename("Study Area" = site,
           "Landscape Class" = landscapeClassWetUplMes,
           "Depth (cm)" = lower_depth,
           #"Mean SOC (g cm$^{-1}$)" = meanSOC,
           "SOC (g cm$^{-1}$)" = SOC) |> 
    flextable() |>merge_v(j = c("Study Area", "Landscape Class")) |> colformat_num(big.mark = ",") |> 
    set_table_properties(layout = "autofit", align = "center") |>
    align(align = "center", part = "all") |>
    mk_par(
    part = "header", j = 3,
    value = as_paragraph("Mean SOC Stock (Mg ha", as_sup("-1"), ")")
    )

wetupltbl
wetuplmestbl


## -----------------------------------------------------------------------------
pcdf <- read.csv("SOIL CARBON/All_WA/data/dataframes/All_WA_PartialCorrelations.csv") |> select(-X) |>
    mutate(num = rep(1:5, each = 5))

pcdf |>
    filter(str_detect(var, "NULL")) |>
    select(var, est) |>
    mutate(WIP = NA,
           DTM = NA,
           CHM = NA,
           HLI = NA,
           lower_depth = NA)
pcdf |>
    filter(!str_detect(var, "NULL") &
           num == 2)  |> select(est) |>
    t() |> as.data.frame() |> tibble::remove_rownames()

names <- pcdf |>
    filter(!str_detect(var, "NULL")) 

names <- c("stat", names[,1])

pcdf |>
    filter(!str_detect(var, "NULL")) |> select(num, est) |> 
    t() |>  as.data.frame() |> tibble::rownames_to_column(var = "stat") |> setNames(names)



## -----------------------------------------------------------------------------
SOC_df <- read.csv("SOIL CARBON/All_WA/data/dataframes/All_WA_MapSOC_wetuplmid_100mask.csv")

SOC_df |> mutate(LandscapeClass = case_when(
             str_detect(Name, "wet") ~ "Wetland",
             str_detect(Name, "upl") ~ "Upland",
             str_detect(Name, "mid") ~ "Mesic",
             .default = "All"),
             StudyArea = case_when(str_detect(Name, "hoh_") ~ "Hoh",
                                   str_detect(Name, "mas_") ~ "Mashel",
                                   str_detect(Name, "col_") ~ "Colville"),
             Total_area = round(Total_area, digits = 0),
             AverageSOC_Mgha = as.character(signif(AverageSOC_Mgha, digits = 3)),
             total_Carbon_Tg = round(total_Carbon_Tg, digits = 1)) |> 
    select(StudyArea, LandscapeClass, everything(), -Name) |> 
        rename("Study Area" = StudyArea,
               "Landscape Class" = LandscapeClass,
                 "Total Area (ha)" = Total_area,
                 "Mean SOC Stock (Mg ha$^{-1}$)" = AverageSOC_Mgha,
                 "Total SOC Stock (Tg)" = total_Carbon_Tg) |> 
    flextable() |> merge_v("Study Area") |> colformat_num(big.mark = ",") |> 
    mk_par(
    part = "header", j = 4,
    value = as_paragraph("Mean SOC Stock (Mg ha", as_sup("-1"), ")")
    )
    # 
    # kbl(align = "rcccc", escape = F, format.args = list(big.mark = ",")) |> 
    # collapse_rows(columns = 1) |> 
    # kable_paper(full_width = F, position = "center", 
    #             bootstrap_options = c("striped", "hover"), fixed_thead = TRUE, html_font = "Aptos")


## -----------------------------------------------------------------------------
mod7graph


## -----------------------------------------------------------------------------
hoh_SOC_sum100mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/hohSOCsum100mask.tif") |> 
    aggregate(fact = 4)
mas_SOC_sum100mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/masSOCsum100mask.tif")|> 
    aggregate(fact = 4)
col_SOC_sum100mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/colSOCsum100mask.tif")|> 
    aggregate(fact = 4)

hoh_SOC_sum200mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/hohSOCsum200.tif")|> 
    aggregate(fact = 4)
mas_SOC_sum200mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/masSOCsum200.tif")|> 
    aggregate(fact = 4)
col_SOC_sum200mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/colSOCsum200.tif")|> 
    aggregate(fact = 4)

hoh_poly <- st_read("SOIL CARBON/All_WA/data/Vectors/HOH_POLYGON_711.gpkg") |> st_transform("EPSG:2855")
mas_poly <- st_read("SOIL CARBON/All_WA/data/Vectors/Mashel_HUC12_2856.shp")
col_poly <- st_read("SOIL CARBON/All_WA/data/Vectors/ColvilleHUC_2855.shp")

hoh_hs <- rast("SOIL CARBON/All_WA/data/Rasters/Hoh_hillshade.tif")|> 
    aggregate(fact = 4)
mas_hs <- rast("SOIL CARBON/All_WA/data/Rasters/Mas_hillshade.tif")|> 
    aggregate(fact = 4)
col_hs <- rast("SOIL CARBON/All_WA/data/Rasters/Col_hillshade.tif")|> 
    aggregate(fact = 4)

hoh_wip <- rast("SOIL CARBON/All_WA/data/Rasters/Hoh_WIP_Final_2855.tif")|> 
    aggregate(fact = 4)
mas_wip <- rast("SOIL CARBON/All_WA/data/Rasters/Mashel_WIP_Final_2856.tif", lyrs = "WET")|> 
    aggregate(fact = 4)
col_wip <- rast("SOIL CARBON/All_WA/data/Rasters/Colville_WIP_Final_2855.tif", lyrs = "WET")|> 
    aggregate(fact = 4)


## -----------------------------------------------------------------------------
tm_func <- function(rast, hs, type){
    rastname <- tools::toTitleCase(str_extract(deparse(substitute(rast)), "hoh|mas|col"))
    
    if(type == "SOC"){
        tm_shape(rast) +
        tm_raster("sum", palette = viridis::magma(n = 20),
                  title = paste0(rastname, " SOC"), style = "cont", legend.reverse = TRUE) +
        tm_shape(hs) +
            tm_raster(col = "hillshade", palette = grey(0:100/100),
                     title = "", style = "cont", alpha = 0.1, legend.show = FALSE) +
            tm_layout(legend.position = c("left", "bottom"), frame = FALSE) + 
            tm_scale_bar(position = c("right", "bottom"))
    } else if(type == "WIP") {
        tm_shape(rast) +
        tm_raster("WET", palette = brewer.pal(n = 9, "YlGnBu"),
                  title = paste0(rastname, " WIP"), style = "cont", legend.reverse = TRUE) +
        tm_shape(hs) +
            tm_raster(col = "hillshade", palette = grey(0:100/100),
                     title = "", style = "cont", alpha = 0.1, legend.show = FALSE) +
            tm_layout(legend.position = c("left", "bottom"), frame = FALSE) + 
            tm_scale_bar(position = c("right", "bottom"))
    }
}

hoh_soc_tm <- tm_func(hoh_SOC_sum100mask, hoh_hs, type = "SOC")
mas_soc_tm <- tm_func(mas_SOC_sum100mask, mas_hs, type = "SOC")
col_soc_tm <- tm_func(col_SOC_sum100mask, col_hs, type = "SOC")

hoh_wip_tm <- tm_func(hoh_wip, hoh_hs, type = "WIP")
mas_wip_tm <- tm_func(mas_wip, mas_hs, type = "WIP")
col_wip_tm <- tm_func(col_wip, col_hs, type = "WIP")

tmap_arrange(hoh_wip_tm, mas_wip_tm, col_wip_tm,
             hoh_soc_tm, mas_soc_tm, col_soc_tm, 
             ncol = 3, nrow = 2)

