#Small modeling script for Teal Carbon Washington
#Anthony J. Stewart
# 2024

library(tidyverse)
library(lme4)
library(lmerTest)
set.seed(11)
setwd('/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')

#### data ####
    # sample_ID is location with repeat measurements at depth 
    # lower_depth is the sampling depth
wa_dat<- read.csv("SOIL CARBON/All_WA/data/dataframes/All_WA_SOC_ExtractPreds_SFcoords.csv")

#Scale predictors 
    # there are other predictors not used here
columns_to_exclude <- c("SOC_stock_spline") # Making sure we don't scale the dependent variable
wa_dat_scale <- wa_dat |> 
    dplyr::select(sample_ID, lower_depth, SOC_stock_spline, 
                  site, CHM, DTM, geomorphons, GEO, WIP,
                  EVI, MNDWI, HLI, lower_depth) |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude),
        ~dplyr::case_when(TRUE ~ scale(.))),
        site = fct_reorder(site, SOC_stock_spline, .fun = "median")) 

hist(wa_dat_scale$SOC_stock_spline)
hist(log10(wa_dat_scale$SOC_stock_spline))


#### model ####
    # most important predictors kept, lowest AIC 
    # CHM = Canopy Height Model
    # DTM = Elevation
    # WIP = Wetland Intrinsic Potential
    # GEO = Geologic Age - categorical
    # site = Watershed study area - categorical (Hoh, Mashel, Colville)
mod <- lmer(log10(SOC_stock_spline) ~ 
                   CHM+WIP*DTM*site+GEO+ 
                   lower_depth + (lower_depth|sample_ID), # allow random slope in sample_ID with depth in model
               data = wa_dat_scale, REML = F, # ML to get AIC 
               control=lmerControl(optimizer="bobyqa")) #optimizer changed to get convergence

#### model diagnostics #### 

plot(mod)
plot(mod, resid(.) ~ fitted(.)|sample_ID, abline = 0) # random effect groupings

qqnorm(wa_dat_agg$SOC, col = as.factor(wa_dat_scale$sample_ID))
qqline(wa_dat_agg$SOC)

lev <- hat(model.matrix(mod)) #leverage 
plot(wa_dat_agg$SOC ~ lev, xlab = "leverage")

cd <- cooks.distance(mod)
plot(lev)
points(cd, col = "blue")

hist(as.vector(unlist(ranef(mod)$sample_ID)))

#### spatial autocorrelation ####
# Create neighbors list using k-nearest neighbors
library(sf)
library(spdep)
#aggregate data by sample ID to remove repeated measures
wa_dat_agg <- wa_dat |> group_by(sample_ID) |> summarise(X = mean(X), 
                                                         Y = mean(Y),
                                                         SOC = sum(SOC_stock_spline))
sf_data <- st_as_sf(wa_dat_agg, coords = c("X", "Y"), crs = 32610)
coords <- st_coordinates(sf_data)
knn_nb <- knearneigh(coords, k = 6)
nb <- knn2nb(knn_nb)

# Create spatial weights
weights <- nb2listw(nb, style = "W")

# semivariogram
inc.lag <- lag.listw(weights, wa_dat_agg$SOC)
plot(inc.lag ~ wa_dat_agg$SOC, pch=16, asp=1)
M1 <- lm(inc.lag ~ wa_dat_agg$SOC)
abline(M1, col="blue")

moran(wa_dat_agg$SOC, weights, length(nb), Szero(weights))

# Calculate Moran's I
moran_test <- moran.test(wa_dat_agg$SOC, weights, alternative="two.sided")
print(moran_test)
# Calculate Moran's I but with MC simulation
MC<- moran.mc(wa_dat_agg$SOC, weights, nsim=999, alternative="two.sided")
MC
plot(MC)

#### DHARMa #### 
    # New package that seems useful to model the residuals of our model for Spatial Autocorrelation
    # Still grouping by sample ID
library(DHARMa)
# simulate residuals from the model 
res <- simulateResiduals(mod)
#recalculate due to the groups in sample location
simulationOutput2 <- recalculateResiduals(res, group = wa_dat$sample_ID)
plot(simulationOutput2)
#aggregating for group coordinates
groupLocations = aggregate(wa_dat[,c("X", "Y")], list(wa_dat$sample_ID), median)

testSpatialAutocorrelation(simulationOutput2, x = groupLocations$X, y = groupLocations$Y)
