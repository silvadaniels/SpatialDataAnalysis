# Spatial Econometrics in R
# Daniel Silva, 2021

# Load packages and data (pseudo Columbus_crime from Anselin, adapted to Belem-PA)
  library("spdep")
  library("rgdal")
  library("spatialreg")
  
  df = readOGR(dsn="/Users/dan/Downloads/Geostatist", layer = "FTP_Soy-munAmaz2013-2017")
  # last columns are: ton_soy, invest, custeio, labor, soy_ha
  
## Initial data analysis ####
  # Spatial weights matrix and neighbors, based on contiguity or distance
    # From coordinates in a dataframe
      #xy = SpatialPoints(df[,8:9])
      #nb = knearneigh(xy, k=1, longlat = T) # command for points; check d range in summary
      #  nb = knn2nb(nb)
    
    # From polygons (i.e., .shp)
      nb = poly2nb(df, row.names = df$CD_GEOCODI, queen=T) # used in polygons. For contiguity; and you can include a snap distance (500km) or remove distance
    
    w = nb2listw(nb, style="W", zero.policy = T)
    
  # OLS regression for reference
    ols = lm(Veg_PastAg ~ lpw, data = df)
      summary(ols)
      
  # Lagrange multiplier test for spatial lag and spatial error dependencies
    lm.LMtests(ols, w, test=c("LMlag", "LMerr"))

## Spatial analysis/regression and Moran's I ####
# run for contiguity and distance weight matrix (just change the 'w'), then compare the results
  # Moran's I test
    moran.test(df$lpw, w) #Consider to use 'zero.policy=T'
    moran.plot(df$lpw, w)
  
  # Spatial lag model
    slm = lagsarlm(Veg_PastAg ~ lpw, data = df, w)
      summary(slm)
  
  # Spatial error model
    sem = errorsarlm(Veg_PastAg ~ lpw, data = df, w)
      summary(sem)
      
  # Spatial lag of X
    slx = lmSLX(Veg_PastAg ~ lpw, data = df, w)
      summary(slx)
  
  # Other models ####
    # SDEM Spatial Durbin Error Model (add lag X to SEM): y=XB+WxT+u, u=LWu+e
      reg5=errorsarlm(Veg_PastAg ~ lpw, data = df, w, etype="emixed")
        summary(reg5)
        
    # SDM Spatial Durbin Model (add lag X to SAR): y=pWy+XB+WXT+e 
      reg6=lagsarlm(Veg_PastAg ~ lpw, data = df, w, type="mixed")
      
    # Manski All-inclusive Model: y=pWy+XB+WXT+u, u=LWu+e (not recommended)
      reg7=sacsarlm(Veg_PastAg ~ lpw, data = df, w, type="sacmixed") #listw2 allows for a different weights matrix for the error structure if desired
      
    # SARAR a.k.a. Kelejian-Prucha, Cliff-Ord, or SAC If all T=0,y=pWy+XB+u, u=LWu+e
      reg8=sacsarlm(Veg_PastAg ~ lpw, data = df, w, type="sac") #listw2 allows for a different weights matrix for the error structure if desired

# Marginal effects (test aforementioned models) ####
  impacts(slm,listw=w) # Calculates de direct (changes on i affect i), indirect (changes on j affect i), and global effect (i+j)
    summary(impacts(slm,listw=w,R=500),zstats=TRUE) #Add zstats,pvals; R=500 not needed for SLX

# Hausman test
  Hausman.test(reg4)
  
# Also, check the estimates of equilibrium effect (section 6a)