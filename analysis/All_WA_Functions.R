# Functions used in All WA SOC Modeling Analysis

#### Libraries ####

library(rgl)
library(terra)
library(lme4)
library(MASS)
library(mgcv)
library(lmerTest)
library(MuMIn)
library(RLRsim)
library(terra)
library(spatialEco)
library(sf)
library(mapview)
library(car)
library(ggplot2)
library(forcats)
library(sjPlot)
library(sjstats)
library(DHARMa)
library(ggeffects)
library(merTools)
library(glmnet)
library(stats)
library(ggcorrplot)
library(RColorBrewer)
library(cowplot)
library(webshot)
library(kableExtra)
library(pdp)
library(vip)
library(formatR)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(ppcor)

#### geospatial ####

terrain_metrics <- function(DEM, path, studyarea){
    stopifnot(file.exists(paste0(path, studyarea, "_HLI", ".tif")) ==FALSE)
    stopifnot(file.exists(paste0(getwd(), "/", path, studyarea, "_geomorphons", ".tif")) ==FALSE)
    
    hli <- spatialEco::hli(DEM, force.hemisphere = "northern")
    writeRaster(hli, names = "HLI", 
                paste0(path, studyarea, "_HLI", ".tif"),
                overwrite = TRUE)
    #geomorphons
    whitebox::wbt_geomorphons(
        dem = sources(DEM),
        output = paste0(getwd(), "/", path, studyarea, "_geomorphons", ".tif"),
        search = 100
    ) 
    
}

stack_resampling <- function(Watershed, WIP, path) {
    stopifnot(file.exists(paste0(getwd(), 
                                 "/", 
                                 path, 
                                 Watershed, 
                                 "_PredictorStack", 
                                 ".tif")) == FALSE)
    
    filelist <- list.files("SOIL CARBON/All_WA/data/Rasters/", 
                           pattern = Watershed, full.names = TRUE,
                           ignore.case = TRUE,
                           include.dirs = FALSE) |> 
        str_subset(pattern = "^(?!.*flood)") |>
        str_subset(pattern = "^(?!.*hillshade)") |>
        str_subset(pattern = "^(?!.*slpasp)") |>
        str_subset(pattern = "^(?!.*fullmodel_v08)")
    f <- lapply(filelist, terra::rast)
    
    rs <- list()
    for (i in 1:length(f)) {
        if(res(f[[i]])[[1]] != res(WIP)[[1]]){
            rs[[i]] <- resample(f[[i]], WIP, method='near', 
                                threads=TRUE)
        } else {
            rs[[i]] <- f[[i]]
        }
    }
    
    s <- rast(rs)
    writeRaster(s, paste0(getwd(), 
                          "/", 
                          path, 
                          Watershed, 
                          "_PredictorStack", 
                          ".tif"), overwrite = TRUE)
    return(s)
}


#### Modeling Analysis ####

# Quick R^2 function
r.sq <- function(y,y.fitted){
    res <- y-y.fitted
    1-sum(res^2)/sum((y-mean(y))^2)
}

fit_vs_var_plot <- function(data, model, x_var, color_var, shape_var){
    if(is.factor(data[x_var][[1]]) == TRUE & is.factor(data[color_var][[1]]) == TRUE ) {
        x_var <- sym(x_var)
        shape_var <- sym(shape_var)
        color_var <- sym(color_var) 
        fitted <- predict(model, newdata = data)
        ggplot(data, aes(y = 10**fitted, x = as.factor(!!x_var))) +
            geom_violin(show.legend = FALSE, scale = "width", linewidth = 0.9)+
            geom_jitter(aes(shape = !!shape_var, , color = !!color_var), 
                        width = 0.2, size = 2, alpha = 0.4) +
            xlab("Study Area\n") + 
            ylab(expression('Fitted Model SOC Stock (g cm'^-2*')')) +
            scale_shape_manual(name = "Study Area", 
                               values = c(16, 17, 18, 19)) +
            scale_colour_viridis_d(name = "Surficial Geology") +
            theme(legend.position = 'right', 
                  legend.key.size = unit(0.6, "cm"),
                  legend.spacing.x = unit(0.8, "cm"),
                  legend.box = "vertical",
                  legend.margin = margin(),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(colour = "grey80"),
                  axis.ticks = element_blank(),
                  text = element_text(size = 9))
    } else if(is.factor(data[x_var][[1]]) == FALSE & is.factor(data[color_var][[1]]) == FALSE) {
        x_var <- sym(x_var)
        shape_var <- sym(shape_var)
        color_var <- sym(color_var) 
        fitted <- predict(model, newdata = data)
        ggplot(data, aes(y = 10**fitted, x = as.factor(!!x_var))) +
            geom_violin(show.legend = FALSE, scale = "width", linewidth = 0.9)+
            geom_jitter(aes(shape = !!shape_var, , color = !!color_var), 
                        width = 0.2, size = 2, alpha = 0.4) +
            xlab("Study Area\n") + 
            ylab(expression('Fitted Model SOC Stock (g cm'^-2*')')) +
            scale_colour_viridis_c(option = "mako", name = deparse(substitute(!!color_var))) +
            theme(legend.position = 'right', 
                  legend.key.size = unit(0.6, "cm"),
                  legend.spacing.x = unit(0.8, "cm"),
                  legend.box = "vertical",
                  legend.margin = margin(),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(colour = "grey80"),
                  axis.ticks = element_blank(),
                  text = element_text(size = 9))
    } else {
        if(x_var == "WIP"){
            label <- "Wetland Intrinsic Potential \n Scaled"
        } else if(x_var == "DTM") {
            label <- "Elevation Scaled\n"
        } else {
            label <- paste0("\n", x_var)
        }
        x_var <- sym(x_var)
        shape_var <- sym(shape_var)
        color_var <- sym(color_var) 
        fitted <- predict(model, newdata = data)
        ggplot(data, aes(y = 10**fitted, x = !!x_var)) +
            geom_point(aes(shape = !!shape_var, color = !!color_var), size = 2, alpha = 0.5)+
            xlab(label) +
            ylab(expression('Fitted Model SOC Stock (g cm'^-2*')')) +
            geom_smooth(aes(y = 10**fitted, x = !!x_var), 
                        method = "lm", color = "#fa3e3e", fill = "#fa3e3e", 
                        linewidth = 0.9, linetype = 5, alpha = 0.3, se = T) + 
            scale_shape_manual(name = "Study Area", 
                               values = c(16, 17, 18, 19)) +
            scale_colour_viridis_d(name = "Surficial Geology") +
            theme(legend.position = 'top', 
                  legend.key.size = unit(0.3, "cm"),
                  legend.spacing.x = unit(0.8, "cm"),
                  legend.box = "vertical",
                  legend.margin = margin(),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(colour = "grey80"),
                  axis.ticks = element_blank(),
                  text = element_text(size = 9))
    }
}

partial_func <- function(vars){
    var_list <- list()
    c_est_list <- list()
    c_pval_list <- list()
    con_var_list <- list()
    p_est_list <- list()
    p_pval_list <- list()
    #site_list <- list()
    
    for(i in 1:length(vars)){
        if(i + 1 <= length(vars)) {
            ct <- cor.test(wa_dat_scale[,vars[[1]]], 
                           wa_dat_scale[,"SOC_stock_spline"])
            pt <- pcor.test(wa_dat_scale[,vars[[1]]], 
                            wa_dat_scale[,"SOC_stock_spline"], 
                            wa_dat_scale[,vars[[i+1]]])
            
            var_list[[i]] <- paste0(vars[[1]], "_NULL")
            c_est_list[[i]] <- ct$estimate[[1]]
            c_pval_list[[i]] <- ct$p.value[[1]]
            con_var_list[[i]] <- vars[[i+1]]
            p_est_list[[i]] <- pt$estimate[[1]]
            p_pval_list[[i]] <- pt$p.value[[1]]
            #site_list[i] <- 
        } else {
            break
        }
    }
    #return(unlist(est_list))
    cdf <- data.frame(var = unlist(var_list), 
                      est = unlist(c_est_list), 
                      pval = unlist(c_pval_list))
    pdf <- data.frame(var = unlist(con_var_list),
                      est = unlist(p_est_list), 
                      pval = unlist(p_pval_list)) 
    df <- rbind(cdf[1,], pdf)
    return(df)
}

#### Visualization #####


ggmaps <- function(data, sm_data, hs){
    rastname <- tools::toTitleCase(str_extract(deparse(substitute(data)), "hoh|mas|col"))
    rasttype <- str_extract(toupper(deparse(substitute(data))), "SOC|WIP")
    
    if(rasttype == "WIP"){
        map <- ggplot() + 
            geom_spatraster(data = data, aes(fill = WIP), maxcell = 2e6) +
            scale_fill_distiller(palette = "YlGnBu",  na.value = NA, direction = 1,
                                 limits = c(0, 1), 
                                 name ='Wetland Intrinsic Potential') +
            ggnewscale::new_scale_fill() + 
            geom_spatraster(data = hs, aes(fill = hillshade), maxcell = 2e6,
                            show.legend = FALSE, alpha = 0.15) +
            scale_fill_distiller(palette = "Greys", na.value = NA)
        sm_map <- ggplot() + geom_spatraster(data = sm_data[[2]], aes(fill = WIP), maxcell = 2e6) +
            scale_fill_distiller(palette = "YlGnBu",  na.value = NA, direction = 1,
                                 limits = c(0, 1)) + 
            theme_void() + theme(legend.position = "none")
        sm_vec <- geom_spatvector(data = sm_data[[1]], fill = NA, colour = "red", linewidth = 0.7)
    } else if(rasttype == "SOC") {
        map <- ggplot() + 
            geom_spatraster(data = data, aes(fill = sum), maxcell = 2e6) +
            scale_fill_viridis_c(option = "inferno",  na.value = NA, 
                                 limits = c(0, 500), name = expression('SOC Stock (Mg ha'^-1*')')) +
            ggnewscale::new_scale_fill() + 
            geom_spatraster(data = hs, aes(fill = hillshade), maxcell = 2e6,
                            show.legend = FALSE, alpha = 0.15) +
            scale_fill_distiller(palette = "Greys", na.value = NA)
        sm_map <- ggplot() + geom_spatraster(data = sm_data[[3]], aes(fill = sum), maxcell = 2e6) +
            scale_fill_viridis_c(option = "inferno",  na.value = NA, 
                                 limits = c(0, 500)) + 
            theme_void() + theme(legend.position = "none")
        sm_vec <- geom_spatvector(data = sm_data[[1]], fill = NA, colour = "cyan", linewidth = 0.7)
    }
    
    legend <- cowplot::get_legend(map + theme(legend.direction = "horizontal",
                                              legend.background = element_blank(),
                                              legend.title.position = "top",
                                              legend.title = element_text(hjust = 0.35), 
                                              legend.text = element_text(size = 8),
                                              legend.key.height = unit(0.5, "cm"),
                                              legend.key.width = unit(1, "cm")))
    
    if (rastname == "Hoh") {
        lbrt <- c("left" = 0.4, "bottom" = 0.15, "right" = 0.75, "top" = 0.5)
        map <- map + 
            coord_sf(expand = FALSE) + 
            annotation_scale(location = "tr", 
                             pad_y = unit(0.15, "cm"), 
                             pad_x = unit(2, "cm"),
                             height = unit(0.1, "cm"))
    } else if (rastname == "Mas"){
        lbrt <- c("left" = -0.15, "bottom" = -0.15, "right" = 0.35, "top" = 0.35)
        map <- map + 
            coord_sf(expand = FALSE) + 
            annotation_scale(location = "tr", 
                             pad_y = unit(0.07, "cm"), 
                             pad_x = unit(2.3, "cm"),
                             height = unit(0.1, "cm"))
    } else if (rastname == "Col"){
        lbrt <- c("left" = 0.6, "bottom" = 0, "right" = 0.9, "top" = 0.3)
        map <- map + 
            coord_sf(expand = FALSE) + 
            annotation_scale(location = "tr", 
                             pad_y = unit(1, "cm"), 
                             pad_x = unit(2, "cm"),
                             height = unit(0.1, "cm"))
    }
    
    mapformat <- map + 
        coord_sf(expand = FALSE) +
        theme(panel.background = element_blank(),
              plot.margin = margin(c(0,0,0,0), "cm"),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank()) 
    
    mapformat_inset <- mapformat + rremove("legend") + sm_vec + 
        inset_element(sm_map, 
                      left = lbrt[["left"]], bottom = lbrt[["bottom"]], 
                      right = lbrt[["right"]], top = lbrt[["top"]], 
                      align_to = "panel")
    
    return(list(mapformat_inset, legend))
}


smallmapcrop <- function(wip, soc){
    extent <- ext(wip) 
    
    if (tools::toTitleCase(str_extract(deparse(substitute(wip)), "hoh|mas|col")) == "Hoh"){
        midx <- extent$xmax[[1]] - ((extent$xmax[[1]] - extent$xmin[[1]])/2)
        midy <- extent$ymax[[1]] - ((extent$ymax[[1]] - extent$ymin[[1]])/2) +5000
    } else if (tools::toTitleCase(str_extract(deparse(substitute(wip)), "hoh|mas|col")) == "Col"){
        midx <- extent$xmax[[1]] - ((extent$xmax[[1]] - extent$xmin[[1]])/2)
        midy <- extent$ymax[[1]] - ((extent$ymax[[1]] - extent$ymin[[1]])/2) -4700
    } else {
        midx <- extent$xmax[[1]] - ((extent$xmax[[1]] - extent$xmin[[1]])/2)
        midy <- extent$ymax[[1]] - ((extent$ymax[[1]] - extent$ymin[[1]])/2) 
    }
    
    df <- data.frame(x = midx, y = midy)
    pt <- vect(df, geom = c("x", "y"), crs = crs(wip))
    buffpt <- buffer(pt, 1000) 
    
    cropwip <- crop(wip, buffpt, mask = TRUE)
    cropsoc <- crop(soc, buffpt, mask = TRUE)
    return(list(buffpt, cropwip, cropsoc))
}

predict_plot <- function(model, data) {
    
    predict_graph <- ggplot(data, 
                            aes(y = 10**(fitted(model)), x = (SOC_stock_spline))) +
        geom_jitter(color='black', 
                    aes(fill = (data$WIP*100), 
                        shape = as.factor(site)),
                    size = 4, stroke = 0.9, alpha = 0.8) +
        scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"),
                             name = "Wetland \nProbability %", n.breaks = 5, limits = c(0, 100)) +
        scale_shape_manual(name = "Study Area", 
                           values = c(21, 22, 23, 24)) +
        geom_smooth(aes(y = 10**(fitted(model)), x = (SOC_stock_spline)), 
                    method = "lm", color = "#fa3e3e", fill = "#fa3e3e", 
                    linewidth = 0.9, linetype = 5, alpha = 0.3, se = T) +
        xlab(expression('Sampled SOC Stock (g cm'^-2*')')) + 
        ylab(expression('Predicted SOC Stock (g cm'^-2*')')) + 
        geom_abline(intercept = 0, slope = 1, linewidth = 0.9, linetype = "dashed") +
        annotate("text", 
                 label = paste("R^{2} == ", 
                               signif(r.sq((data$SOC_stock_spline),
                                           10**fitted(model)), 3)), 
                 x = 0.5, y = 3.5, size = 4, parse = T) + 
        annotate("text", label = "Model Fit", 
                 x = 1.5, y = 3.5, size = 4) +
        annotate("segment", color = "#fa3e3e",
                 x = 1.3, y = 3.3, linewidth = 0.9, linetype = 5, 
                 xend = 1.75) + 
        annotate("text", label = "1:1", 
                 x = 2.5, y = 3.5, size = 4) +
        annotate("segment", color = "black",
                 x = 2.3, y = 3.3, linewidth = 0.9, linetype = "dashed",
                 xend = 2.75) +
        xlim(0, 4) +
        ylim(0, 4)  +
        theme(legend.position = 'right', 
              legend.key.size = unit(0.6, "cm"),
              legend.spacing.x = unit(1.2, "cm"),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "grey80"),
              axis.ticks = element_blank(),
              text = element_text(size = 12)) +
        guides(guide_legend(byrow = TRUE))
}
