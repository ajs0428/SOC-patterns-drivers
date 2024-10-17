library(lme4)
library(lmerTest)
library(terra)
library(tidyverse)
library(ggcorrplot)
library(MuMIn)
library(mgcv)


r.sq <- function(y,y.fitted){
    res <- y-y.fitted
    1-sum(res^2)/sum((y-mean(y))^2)
}

setwd('/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')

wa_hor <- read.csv("SOIL CARBON/All_WA/data/dataframes/All_WA_horizons.csv")

stocks_func <- function(horizons){
    stopifnot("sample_ID" %in% names(horizons))
    stopifnot("depth_cm" %in% names(horizons))
    stopifnot("carbon_perc" %in% names(horizons))
    stopifnot("BD_g_cm3" %in% names(horizons))
    stopifnot("rock_perc" %in% names(horizons))
    
    stocks <- horizons |>  dplyr::group_by(sample_ID) |>
        mutate(top = case_when(is.na(depth_cm - lag(depth_cm)) ~ 0,
                               .default = lag(depth_cm)), 
               bottom = depth_cm, 
               center = abs(top - (top - bottom)/2)) |>
        mutate(thick30 = case_when(is.na(depth_cm - lag(depth_cm)) ~ bottom,
                                   bottom > 30 & top <= 30 ~ bottom - (bottom -30)-lag(bottom),
                                   bottom > 30 & top >= 30 ~ 0,
                                   bottom < 30 & bottom - lag(bottom) < 0 ~ bottom,
                                   bottom < 30 & bottom - lag(bottom) > 0 ~ bottom - lag(bottom),
                                   bottom == 30 ~ bottom - lag(bottom),
                                   .default = 0),
               SOC_hor_stock30 = ((carbon_perc/100)*BD_g_cm3*thick30*(1-(rock_perc/100)))) |>
        # clay_hor_30 = (((clay/100)*thick30)/30*(1-(rock_perc/100))),
        # silt_hor_30 = (((silt/100)*thick30)/30*(1-(rock_perc/100))),
        # sand_hor_30 = (((sand/100)*thick30)/30*(1-(rock_perc/100)))) |>
        mutate(thick60 = case_when(is.na(depth_cm - lag(depth_cm)) ~ bottom,
                                   bottom > 60 & top <= 60 ~ bottom - (bottom -60)-lag(bottom),
                                   bottom > 60 & top >= 60 ~ 0,
                                   bottom < 60 & bottom - lag(bottom) < 0 ~ bottom,
                                   bottom < 60 & bottom - lag(bottom) > 0 ~ bottom - lag(bottom),
                                   bottom == 60 ~ bottom - lag(bottom),
                                   .default = 0),
               SOC_hor_stock60 = ((carbon_perc/100)*BD_g_cm3*thick60*(1-(rock_perc/100)))) |>
        # clay_hor_60 = (((clay/100)*thick60)/60*(1-(rock_perc/100))),
        # silt_hor_60 = (((silt/100)*thick60)/60*(1-(rock_perc/100))),
        # sand_hor_60 = (((sand/100)*thick60)/60*(1-(rock_perc/100)))) |>
        mutate(thick100 = case_when(is.na(depth_cm - lag(depth_cm)) ~ bottom,
                                    bottom > 100 & top < 100 ~ bottom - (bottom -100)-lag(bottom),
                                    bottom > 100 & top == 0 ~ bottom - (bottom -100)-0,
                                    bottom > 100 & top > 100 ~ 0,
                                    bottom < 100 & bottom - lag(bottom) < 0 ~ bottom,
                                    bottom < 100 & bottom - lag(bottom) > 0 ~ bottom - lag(bottom),
                                    bottom == 100 ~ bottom - lag(bottom),
                                    .default = 0),
               SOC_hor_stock100 = ((carbon_perc/100)*BD_g_cm3*thick100*(1-(rock_perc/100)))) |>
        # clay_hor_100 = (((clay/100)*thick100)/100*(1-(rock_perc/100))),
        # silt_hor_100 = (((silt/100)*thick100)/100*(1-(rock_perc/100))),
        # sand_hor_100 = (((sand/100)*thick100)/100*(1-(rock_perc/100))) ) |>
        mutate(thick120 = case_when(is.na(depth_cm - lag(depth_cm)) ~ bottom,
                                    bottom > 120 & top <= 120 ~ bottom - (bottom -120)-lag(bottom),
                                    bottom > 120 & top >= 120 ~ 0,
                                    bottom < 120 & bottom - lag(bottom) < 0 ~ bottom,
                                    bottom < 120 & bottom - lag(bottom) > 0 ~ bottom - lag(bottom),
                                    bottom == 120 ~ bottom - lag(bottom),
                                    .default = 0),
               SOC_hor_stock120 = ((carbon_perc/100)*BD_g_cm3*thick120*(1-(rock_perc/100)))) |>
        # clay_hor_120 = (((clay/100)*thick120)/120*(1-(rock_perc/100))),
        # silt_hor_120 = (((silt/100)*thick120)/120*(1-(rock_perc/100))),
        # sand_hor_120 = (((sand/100)*thick120)/120*(1-(rock_perc/100))) ) |>
        mutate(thickfull = case_when(is.na(depth_cm - lag(depth_cm)) ~ bottom,
                                     bottom > 999 & top <= 999 ~ bottom - (bottom -999)-lag(bottom),
                                     bottom > 999 & top >= 999 ~ 0,
                                     bottom < 999 & bottom - lag(bottom) < 0 ~ bottom,
                                     bottom < 999 & bottom - lag(bottom) > 0 ~ bottom - lag(bottom),
                                     bottom == 999 ~ bottom - lag(bottom),
                                     .default = 0),
               SOC_hor_stockfull = ((carbon_perc/100)*BD_g_cm3*thickfull*(1-(rock_perc/100)))) |>
        # clay_hor_full = (((clay/100))*(1-(rock_perc/100))),
        # silt_hor_full = (((silt/100))*(1-(rock_perc/100))),
        # sand_hor_full = (((sand/100))*(1-(rock_perc/100))) ) |>
        dplyr::summarise(SOC_stock_100 = sum(SOC_hor_stock100)*100,
                         SOC_stock_30 = sum(SOC_hor_stock30)*100,
                         SOC_stock_60 = sum(SOC_hor_stock60)*100,
                         SOC_stock_120 = sum(SOC_hor_stock120)*100,
                         SOC_stock_full = sum(SOC_hor_stockfull)*100,
                         # clay_30 = sum(clay_hor_30),
                         # silt_30 = sum(silt_hor_30),
                         # sand_30 = sum(sand_hor_30),
                         # clay_60 = sum(clay_hor_60),
                         # silt_60 = sum(silt_hor_60),
                         # sand_60 = sum(sand_hor_60),
                         # clay_100 = sum(clay_hor_100),
                         # silt_100 = sum(silt_hor_100),
                         # sand_100 = sum(sand_hor_100),
                         # clay_120 = sum(clay_hor_120),
                         # silt_120 = sum(silt_hor_120),
                         # sand_120 = sum(sand_hor_120),
                         # clay_full = sum(clay_hor_full),
                         # silt_full = sum(silt_hor_full),
                         # sand_full = sum(sand_hor_full),
                         #stock_check = mean(SOC120),
                         #WIP = mean(WIP),
                         site = first(site),
                         maxdepth = max(depth_cm),
                         n = n(),
                         lat = mean(lat, na.rm = T),
                         lon = mean(lon, na.rm = T),
                         # sand = mean(sand, na.rm = T),
                         # silt = mean(silt, na.rm = T),
                         # clay = mean(clay, na.rm = T),
                         # ph = mean(ph, na.rm = T),
                         # hlf = mean(hlf_perc, na.rm = T)
                         # .by = sample_ID
        ) |> 
        # dplyr::summarise(across(where(is.numeric) & starts_with("SOC"), sum),
        #                  across(where(is.character), first),
        #                  across(c("depth_cm"), max),
        #                  across(c("lat", "lon"), ~ mean(.x, na.rm = T)),
        #                  across(c("cn", "ph", "sand", "silt", "clay"), ~ mean(.x, na.rm= T))) |>
        dplyr::arrange(sample_ID, desc(lat), .by_group = TRUE) |> 
        mutate(SOC_stock_0to30 = SOC_stock_30,
               SOC_stock_30to60 = SOC_stock_60 - SOC_stock_30,
               SOC_stock_60to100 = SOC_stock_100 - SOC_stock_60,
               SOC_stock_100to120 = SOC_stock_120 - SOC_stock_100) |> 
        as.data.frame()
    return(stocks)
}


wa_stocks <- stocks_func(wa_hor)
hoh_stocks_pts <- wa_stocks |> 
    filter(site == "HOH") |> 
    vect( geom = c("lon", "lat")) |> 
    terra::set.crs("EPSG:4326") |> 
    terra::project("EPSG:2855")
mas_stocks_pts <- wa_stocks |> 
    filter(site == "MAS") |> 
    vect( geom = c("lon", "lat")) |> 
    terra::set.crs("EPSG:4326") |> 
    terra::project("EPSG:2856")
col_stocks_pts <- wa_stocks |> 
    filter(site == "COL") |> 
    vect( geom = c("lon", "lat")) |> 
    terra::set.crs("EPSG:4326") |> 
    terra::project("EPSG:2855")


hoh_stack <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/Hoh_PredictorStack_Class.tif") |> 
    tidyterra::rename(EVI = EVI_median)
hoh_slp <- terra::terrain(hoh_stack$DTM, v = "slope", neighbors = 8)
hoh_dep50 <- rast("WIP/Hoh_data/TopoIndices/dev_50.tif") |> project(hoh_stack$DTM)
hoh_dep300 <- rast("WIP/Hoh_data/TopoIndices/dev_300.tif") |> project(hoh_stack$DTM)
hoh_dep1000 <- rast("WIP/Hoh_data/TopoIndices/dev_1000.tif") |> project(hoh_stack$DTM)

hoh_stack_plus <- c(hoh_stack, hoh_slp, hoh_dep50, hoh_dep300, hoh_dep1000)

mas_stack <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/Mas_PredictorStack_Class.tif") |> 
    tidyterra::rename(EVI = EVI_median)
mas_slp <- terra::terrain(mas_stack$DTM, v = "slope", neighbors = 8)
mas_dep50 <- rast("Mashel/Mashel_WIP_2024/DataExport/dev_50.tif") |> project(mas_stack$DTM)
mas_dep300 <- rast("Mashel/Mashel_WIP_2024/DataExport/dev_300.tif") |> project(mas_stack$DTM)
mas_dep1000 <- rast("Mashel/Mashel_WIP_2024/DataExport/dev_1000.tif") |> project(mas_stack$DTM)

mas_stack_plus <- c(mas_stack, mas_slp, mas_dep50, mas_dep300, mas_dep1000)

col_stack <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/Col_PredictorStack_Class.tif") |> 
    tidyterra::rename(EVI = EVI_median)
col_slp <- terra::terrain(col_stack$DTM, v = "slope", neighbors = 8)
col_dep50 <- rast("Colville/Colville_WIP_2024/DataExport/dev_50.tif") |> project(col_stack$DTM)
col_dep300 <- rast("Colville/Colville_WIP_2024/DataExport/dev_300.tif") |> project(col_stack$DTM)
col_dep1000 <- rast("Colville/Colville_WIP_2024/DataExport/dev_1000.tif") |> project(col_stack$DTM)

col_stack_plus <- c(col_stack, col_slp, col_dep50, col_dep300, col_dep1000)

hoh_stack_plus$GEO <- subst(hoh_stack_plus$GEO, from = c(0,1,2,3), c("Eocene", "Miocene", "Pleistocene", "Quaternary"))
mas_stack_plus$GEO <- subst(mas_stack_plus$GEO, from = c(0,1,2,3), c("Eocene", "Miocene", "Pleistocene", "Quaternary")) 
col_stack_plus$GEO <- subst(col_stack_plus$GEO, from = c(0,1,2,3), c("Eocene", "Pleistocene", "PreTertiary", "Quaternary")) 

hoh_stocks_ext <- terra::extract(hoh_stack_plus, hoh_stocks_pts, method = "simple", bind = TRUE)
mas_stocks_ext <- terra::extract(mas_stack_plus, mas_stocks_pts, method = "simple", bind = TRUE)
col_stocks_ext <- terra::extract(col_stack_plus, col_stocks_pts, method = "simple", bind = TRUE)

wa_stocks_dat <- rbind(data.frame(hoh_stocks_ext),data.frame(mas_stocks_ext), data.frame(col_stocks_ext)) |> 
    mutate(
        GEO = as.factor(GEO),
        EVT = as.factor(EVT),
        geomorphons = as.factor(geomorphons),
        site = as.factor(site),
        site = fct_reorder(site, SOC_stock_100, .fun = "median")) |>
    dplyr::rename_with(~gsub("_median", "", .x, fixed = TRUE))

columns_to_exclude <- c("SOC_stock_100") 

wa_stocks_scale <- wa_stocks_dat |> 
    dplyr::select(sample_ID, SOC_stock_100, 
                  site, CHM, DTM, geomorphons, GEO, WIP, EVT,
                  EVI, NDVI, NDYI, VPD,
                  MNDWI, HLI, MAT, MAP) |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude),
        ~dplyr::case_when(TRUE ~ scale(.))),
        site = fct_reorder(site, SOC_stock_100, .fun = "median")) 

wa_stocks_wet <- filter(wa_stocks_dat, WIP >=0.5)
wa_stocks_upl <- filter(wa_stocks_dat, WIP <0.5)

wa_stocks_num <- wa_stocks_wet |> 
    dplyr::select(where(is.numeric)) |> as.matrix()

corplot <- ggcorrplot(cor(wa_stocks_num), method = "square", type = "full", lab = T, lab_size = 3)

corplot


mod7 <- lm((SOC_stock_100) ~ 
               EVI+ HLI + MAT + MAP +(GEO)+WIP*geomorphons,
             #family = Gamma(link = "log"),
             data = wa_stocks_dat
            )
summary(mod7)
anova(mod7)

ggplot(wa_stocks_scale, aes(x = SOC_stock_100, y = fitted(mod7))) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0)

modg <- glm((SOC_stock_100) ~ 
                EVI+ HLI + MAT + MAP +(GEO)+WIP,
            family = Gamma(link = "log"),
            data = wa_stocks_wet)
summary(modg)

gmod <- glm((SOC_stock_100) ~ 
                  CHM + HLI + GEO + MAP+EVI + DTM + VPD + WIP + slope + dev_300+dev_1000, 
            family = Gamma(link = "log"),
              data = wa_stocks_dat, na.action = 'na.fail')

gmodw <- glm((SOC_stock_100) ~ 
                CHM + HLI + GEO + MAP+EVI + DTM + VPD + WIP + slope + dev_300+dev_1000, 
            family = Gamma(link = "log"),
            data = wa_stocks_wet, na.action = 'na.fail')

gmodu <- glm((SOC_stock_100) ~ 
                CHM + HLI + GEO + MAP+EVI + DTM + VPD + WIP + slope + dev_300+dev_1000, 
            family = Gamma(link = "log"),
            data = wa_stocks_upl, na.action = 'na.fail')

modd <- dredge(gmod, rank = "AIC", extra = "R^2", m.lim = c(1,15))
moddw <- dredge(gmodw, rank = "AIC", extra = "R^2", m.lim = c(1,15))
moddu <- dredge(gmodu, rank = "AIC", extra = "R^2", m.lim = c(1,15))

head(modd, 10) 
modd1 <- get.models(modd, 1)[[1]] # fewest terms and lowest AIC
modd1w <- get.models(moddw, 1)[[1]] # fewest terms and lowest AIC
modd1u <- get.models(moddu, 1)[[1]] # fewest terms and lowest AIC

summary(modd1)
summary(modd1w)
summary(modd1u)

anova(modd1, test = "F")
anova(modd1w, test = "F")
anova(modd1u, test = "F")

r2mod1 <- 1-(modd1$deviance/modd1$null.deviance)
r2mod1w <- 1-(modd1w$deviance/modd1w$null.deviance)
r2mod1u <- 1-(modd1u$deviance/modd1u$null.deviance)

performance::r2_nagelkerke(modd1)
performance::r2_nagelkerke(modd1w)
performance::r2_nagelkerke(modd1u)

r2glmm::r2beta(modd1, method = "lm")
r2glmm::r2beta(modd1w, method = "lm")
r2glmm::r2beta(modd1u, method = "lm")

confint(modd1)
confint(modd1w)
confint(modd1u)

ggplot(wa_stocks_upl, aes(x = WIP, y = SOC_stock_100, colour = WIP)) +
    geom_point() + 
    geom_smooth(span = 0.5)

gammod <- gam((SOC_stock_100) ~ EVI+ HLI + MAT + s(MAP, k = 12) + (GEO) + s(WIP, k = 8) + (slope),
              family = Gamma(link = "log"),
              data = wa_stocks_dat)
summary(gammod)

ggplot(wa_stocks_dat, aes(x = SOC_stock_100, y = predict(modd1, type = "response"))) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "loess")
ggplot(wa_stocks_wet, aes(x = SOC_stock_100, y = predict(modd1w, type = "response"))) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "loess")
ggplot(wa_stocks_upl, aes(x = SOC_stock_100, y = predict(modd1u, type = "response"))) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "loess")


wa_stocks_dat_plot <- wa_stocks_dat |> mutate(wetupl = case_when(WIP >= 0.5 ~ "WET",
                                                                 .default = "UPL"))
ggplot() +
    geom_point(aes(x = SOC_stock_100, y = fitted(modd1w)), data = wa_stocks_wet, colour = "aquamarine3") +
    geom_point(aes(x = SOC_stock_100, y = fitted(modd1u)), data = wa_stocks_upl, colour = "tan") +
    geom_smooth(aes(x = SOC_stock_100, y = fitted(modd1w)), data = wa_stocks_wet, method = "lm", colour = "aquamarine3", se = F) +
    geom_smooth(aes(x = SOC_stock_100, y = fitted(modd1u)), data = wa_stocks_upl, method = "lm", colour = "tan", se = F) +
    geom_smooth(aes(x = SOC_stock_100, y = fitted(modd1)), data = wa_stocks_dat, method = "lm", se = F) +
    geom_abline(slope = 1, intercept = 0) +
    ylim(0, 700) + 
    xlim(0,700)
   

r.sq(wa_stocks_wet$SOC_stock_100, predict(modd1w, newdata = wa_stocks_wet, type = "response"))
r.sq(wa_stocks_upl$SOC_stock_100, predict(modd1u, newdata = wa_stocks_upl, type = "response"))
r.sq(wa_stocks_dat$SOC_stock_100, predict(modd1, newdata = wa_stocks_dat, type = "response"))

library(segmented)

sg <- segmented(modd1, ~WIP,  data = wa_stocks_dat)
summary(sg)
slope(sg)
plot(sg, term = "WIP")
plot(wa_stocks_dat$SOC_stock_100, predict(sg, type = "response")) + abline(0,1)

##############################################################################################################################
predict(modg, type = "response")

hoh_stack_plus_wet <- terra::mask(hoh_stack_plus, hoh_stack_plus$WIP >= 0.5, maskvalues = FALSE, updatevalue = NA)
hoh_stack_plus_upl <- terra::mask(hoh_stack_plus, hoh_stack_plus$WIP >= 0.5, maskvalues = TRUE, updatevalue = NA)

testmap <- predict(hoh_stack, modg, type = "response")
plot((testmap))
testmapw <- predict(hoh_stack_plus_wet, modd1w, type = "response")
testmapw_mask <- terra::mask(testmapw, hoh_stack_plus_wet$MNDWI_median > 0.3, maskvalues = TRUE, updatevalue = NA)
plot((testmapw_mask))


C_map_simp <- function(C_map){
    name <- deparse(substitute(C_map))
    gt <- (C_map > -999)
    cell_size <- cellSize(gt, unit = "ha") |> mask(mask = gt)
    carbon_cell <- C_map*cellSize(gt, unit = "ha") # carbon in Mg per cell which is then added up 
    
    area_tot <- sum(values(cell_size), na.rm = T)
    C_mean <- mean(values(C_map), na.rm = T) #mean value of all values of Mg/ha cells
    TotalC_sum <- sum(values(carbon_cell), na.rm =T)
    
    
    return(data.frame("Name" = name, 
                      "Total_area" = area_tot, 
                      "AverageSOC_Mgha" = C_mean, 
                      "total_Carbon_Tg" = TotalC_sum/1e6,
                      stringsAsFactors = T))
}

C_map_simp(testmapw)


C_map_wet_fractions <- function(C_map, WIP){
    name <- deparse(substitute(C_map))
    gt <- (C_map > -999)
    WIP_wetupl <- (WIP >= 0.50)
    WIP_mid <- (WIP >= 0.25 & WIP <= 0.75)
    
    C_map_wet <- mask(C_map, mask = WIP_wetupl, maskvalues = FALSE)
    C_map_upl <- mask(C_map, mask = WIP_wetupl, maskvalues = TRUE)
    C_map_mid <- mask(C_map, mask = WIP_mid, maskvalues = FALSE)
    
    cell_size_all <- cellSize(gt, unit = "ha") |> mask(mask = gt)
    cell_size_wet <- cellSize(gt, unit = "ha") |> mask(mask = C_map_wet)
    cell_size_upl <- cellSize(gt, unit = "ha") |> mask(mask = C_map_upl)
    cell_size_mid <- cellSize(gt, unit = "ha") |> mask(mask = C_map_mid)
    
    carbon_cell_all <- C_map*cell_size_all # carbon in Mg per cell which is then added up 
    carbon_cell_wet <- C_map_wet*cell_size_wet # carbon in Mg per cell which is then added up 
    carbon_cell_upl <- C_map_upl*cell_size_upl # carbon in Mg per cell which is then added up 
    carbon_cell_mid <- C_map_mid*cell_size_mid # carbon in Mg per cell which is then added up 
    
    area_tot <- sum(values(cell_size_all), na.rm = T)
    area_wet <- sum(values(cell_size_wet), na.rm = T)
    area_upl <- sum(values(cell_size_upl), na.rm = T)
    area_mid <- sum(values(cell_size_mid), na.rm = T)
    
    C_mean_all <- mean(values(C_map), na.rm = T) #mean value of all values of Mg/ha cells
    C_mean_wet <- mean(values(C_map_wet), na.rm = T)
    C_mean_upl <- mean(values(C_map_upl), na.rm = T)
    C_mean_mid <- mean(values(C_map_mid), na.rm = T)
    
    TotalC_sum <- sum(values(carbon_cell_all), na.rm =T)
    TotalC_sum_wet <- sum(values(carbon_cell_wet), na.rm =T)
    TotalC_sum_upl <- sum(values(carbon_cell_upl), na.rm =T)
    TotalC_sum_mid <- sum(values(carbon_cell_mid), na.rm =T)
    
    return(data.frame("Name" = c(name, paste0(name, "wet"), paste0(name, "upl"), paste0(name, "mid")), 
                      "Total_area" = c(area_tot, area_wet, area_upl, area_mid),
                      "AverageSOC_Mgha" = c(C_mean_all, C_mean_wet, C_mean_upl, C_mean_mid), 
                      "total_Carbon_Tg" = c(TotalC_sum/1e6, TotalC_sum_wet/1e6, 
                                            TotalC_sum_upl/1e6, TotalC_sum_mid/1e6),
                      stringsAsFactors = T))
}

C_map_wet_fractions(testmap, hoh_stack$WIP)
