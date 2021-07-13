# Panel Spatial regression
# Daniel Silva, 2021

# Load packages and data (pseudo Columbus_crime from Anselin, adapted to Belem-PA)
  library("spdep")
  #library("rgdal") # in order to import vectors (i.e., .shp)
  library("spatialreg")
  library("readr")
  library("dplyr")
  library("splm")
  library("plm")

  df_panel = read_csv("data/panel_regrowth_Cerrado.csv")
    df_panel = subset(df_panel, drop_4balance==0)
    coord = read_csv("data/coord_panelSR_Cerrado.csv") #separated file bc  we're holding the weights constant
    coord = subset(coord, drop_4balance==0)

## 1. Initial data analysis ####
  # 1.1. Spatial weights matrix and neighbors, based on contiguity or distance
    # From coordinates in a dataframe
    xy = SpatialPoints(coord[,3:2]) #from long to lat
    nbs = make.sym.nb(knn2nb(knearneigh(xy,k=6))) # symmetric weight matrix needed for panel
      plot(nbs, xy)

    lws <- nb2listw(nbs, zero.policy = T) # using symmetric weight matrix
  
  # 1.2. simple panel OLS and FE regression for reference
    ols = plm(regrowth_ha ~ lpw +embargos_def +veg_pers, data = df_panel, index = c("geocode", "year"), model="pooling")
      summary(ols)

# 2. Panel spatial regression (2001-2016) ####
# functions 'splm', 'spreml', and 'spgm'; the last run IV for panel (see the options method and instruments in the function)

  # Fixed effect
    sfe = spgm(regrowth_ha ~ lpw +embargos_def +veg_pers, data = df_panel, index = c("geocode", "year"),
               listw = lws, model = "within", lag = FALSE, spatial.error = TRUE)
      summary.splm(sfe)
  
  # Random effect
    sre = spreml(regrowth_ha ~ lpw +embargos_def +veg_pers, data = df_panel, index = c("geocode", "year"),
                 w = lws, errors = "sem2re")
      summary.splm(sre)
  
  # Marginal effect (demands a splm model)
    impac1 <- impacts(sfe, listw = lws, time = 16)
      summary(impac1, zstats=TRUE, short=TRUE)
      
  # Hausman test
    sphtest(regrowth_ha ~ lpw +embargos_def +veg_pers, data = df_panel,
                   listw = lws, spatial.model="error")
   
    sphtest(sre, sfe) #compare the previous models (e.g., FE and RE)
    