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

## 1. Initial data analysis ####
  # 1.1. Spatial weights matrix and neighbors, based on contiguity or distance
    # From coordinates in a dataframe
    xy = SpatialPoints(df_panel[,14:13]) #from long to lat
    nb = knn2nb(knearneigh(xy))
    nbs <- make.sym.nb(nb) #symetric weights
      plot(nb, xy)
    
    # From polygons (i.e., .shp), use 'poly2nb()'
    
    w = nb2listw(nb, style="W", zero.policy = T)
    lws <- nb2listw(nbs) #symetric
  
  # 1.2. simple panel OLS and FE regression for reference
    ols = plm(regrowth ~ lpw, data = df_panel, index = c("geocode", "year"), model="pooling")
      summary(ols)
    
    fe = plm(regrowth ~ lpw, data = df_panel, index = c("geocode", "year"), model="within")
      summary(fe)
  
    # Lagrange multiplier test for spatial lag and spatial error dependencies
    lm.LMtests(ols, w, test=c("LMlag", "LMerr")) #if not signif for spatial, dump it

  # 1.3. Moran's I test
    moran.test(df$lpw2016, w) #Consider to use 'zero.policy=T'
    moran.plot(df$lpw2016, w)

    #local = localmoran(x = df$lpw2016, listw = w) # Ii: local moran statistic; E.Ii: expectation of local moran statistic; Var.Ii: variance of local moran statistic; Z.Ii: standard deviate of local moran statistic

# 2. Panel spatial regression (2001-2016) ####
  #df_panel = pdata.frame(df_panel, index = c("geocode", "year"))
  
  # Random effect
    sre = spml(regrowth ~ lpw, data = df_panel, index = c("geocode", "year"),
               listw = lws, model = "random", lag = FALSE, spatial.error = "b")
      summary(sre)
  
  # Fixed effect
    sfe <- spml(regrowth ~ lpw +embargos_def, data = df_panel, listw = w,
              model="within", spatial.error="b", Hess = FALSE)
      summary(sfe)
  
  # Marginal effect
    impac1 <- impacts(sre, listw = w, time = 16)
      summary(impac1, zstats=TRUE, short=TRUE)
