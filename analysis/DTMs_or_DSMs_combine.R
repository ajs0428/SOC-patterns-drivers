library(terra)
setwd('/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')

list_rasts <- list.files("Colville/Colville_DSMs/", full.names = TRUE)

for(i in 1:length(list_rasts)){
    r_1 <- rast(list_rasts[[1]]) 
    r <- rast(list_rasts[[i]])
    name_r <- names(r)
    fact = res(r)[1]/res(r_1)[1]
    
    if(res(r)[1]!=res(r_1)[1]){
        r <- disagg(r, fact, filename = paste0("Colville/Colville_DSMs/", 
                                                name_r, "_res.tif"),
                      overwrite = TRUE)
    } else {
        next
    }
}

plot(combined_raster_all)

### Or

for(i in 1:length(list_rasts)){
    r_1 <- rast(list_rasts[[1]]) 
    r <- rast(list_rasts[[i]])
    if(i == 1){
            combined_raster <- r
        } else {
            combined_raster <- terra::mosaic(combined_raster, r)
        }
    }


### ORRRR

f <- lapply(list.files("Colville/Colville_DSMs", pattern = ".tif", full.names = TRUE, include.dirs = FALSE),
            terra::rast)
sprc_f <- terra::sprc(f)

mosaic_f <- mosaic(sprc_f, fun="mean", 
                   filename = "Colville/Colville_WIP_2024/DataExport/Colville_DSMs_mosaic.tif",
                   overwrite = TRUE)
mosaic_f <- rast("Colville/Colville_WIP_2024/DataExport/Colville_DSMs_mosaic.tif") |> project("EPSG:2855")

col_dem <- rast("Colville/Colville_WIP_2024/DataExport/colvilledem_huc2855.tif")
col_poly <- vect("Colville/Colville_WIP_2024/DataExport/ColvilleHUC.shp") |> project("EPSG:2855")

mosaic_f_crp <- terra::mask(mosaic_f, col_poly)
mosaic_f_prj <- project(mosaic_f_crp, col_dem)

canopy_height <- mosaic_f_prj - col_dem
m <- c(-99, 0, 0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
canopy_height_rc1 <- classify(canopy_height, rclmat, include.lowest=TRUE)
writeRaster(canopy_height_rc1, "Colville/Colville_WIP_2024/DataExport/Colville_CanopyHeight_2856.tif")


