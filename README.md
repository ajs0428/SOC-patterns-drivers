Code for analysis, figures, tables, writing of the manuscript: \
Patterns and drivers of newly revealed wetland and wet soil carbon stocks

This manuscript is under development and all results are preliminary

Data files needed for reproducing. Other data may be included in scripts\
but these are the necessary files
- Tabulated
	- PNW Soil Carbon Stocks (pnw_stocks_data.csv)
	- PNW Soil Carbon Concentration (pnw_lab_data.csv) 
- Vector 
	- Lithology (LITH)
	- Landform (geomorphons)
- Raster
	- Canopy Height (CHM)
	- Elevation (DTM)
	- Enhanced Vegetation Index (EVI)
	- Heat Load Index (HLI)
	- Mean Annual Precipitation (MAP)
	- Mean Annual Temperature (MAT) 
	- Potential Evapotranspiration (PET)
	- Normalized Difference Vegetation Index (NDVI)
	- Wetland Intrinsic Potential (WIP)


The majority of the R code is contained in Quarto (.qmd) files \
which is similar to Rmarkdown and can be run in RStudio

R workflow for reproducing:

- The first steps are to take raw data and process to create \
splined observations. Then extract geospatial data for modeling 

1) PNW_SOC_spline_processing.qmd
2) PNW_geospatial_processing.qmd
3) PNW_Training_Testing_Data_Split_Scale.qmd

- Once the complete dataframes are created, modeling can be done \
through the following but in no particular order:
	- PNW_SOC_Linear_Modeling_Analysis.qmd
	- PNW_SOC_RF_Modeling.qmd
	- PNW_SOC_Edaphic_Linear_Modeling_Analysis.qmd
	- PNW_SOC_Edaphic_RF_Modeling.qmd
	- PNW_SOC_Linear_Modeling_Analysis_noWIP.qmd
	- PNW_SOC_RF_Modeling_noWIP.qmd	

	- .qmd scripts with "Edaphic" in the name refer to modeling \
SOC % with additional predictors including soil pH and texture. \
Otherwise the scripts model stocks 

 	- Additional .qmd documents with "noWIP" are similar \
modeling scripts but without the WIP parameter

- Finally mapped SOC stock predictions are produced with 
	- PNW_SOC_Prediction.qmd 
	- QGIS for stylistic rendering
	
Figures and tables for the manuscript are created in
    - PNW_Figures_Tables.qmd
There are many figures and tables which require derived tables \
or models from the earlier scripts. However these are also \
in a zip folder
- required_figurestables.zip


Authors: 
Anthony J Stewart, Meghan Halabisky, David V. D’Amore, \
Diogo Spinola, Chad Babcock, L. Monika Moskal, David Butman
