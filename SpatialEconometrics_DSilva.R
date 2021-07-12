# Spatial Econometrics
# Daniel Silva, 2021

# Load packages and data (pseudo Columbus_crime from Anselin, adapted to Belem-PA)
  library("spdep")
  #library("rgdal") # in order to import vectors (i.e., .shp)
  library("spatialreg")
  library("readr")
  library("dplyr")

  df = read_csv("data/regrowth-cs_Cerrado.csv")

## 1. Initial data analysis (cross-section 2001-2016) ####
  # 1.1. Spatial weights matrix and neighbors, based on contiguity or distance
    # From coordinates in a dataframe
      xy = SpatialPoints(df[,15:14]) #from long to lat
      nb = knn2nb(knearneigh(xy, k=6))
        plot(nb, xy) #visualize to double check
    
    # From polygons (i.e., .shp), use 'poly2nb()'
  
    w = nb2listw(nb, style="W", zero.policy = T) # zero.policy to consider 'island' of points
    
  # 1.2. OLS regression for reference
    ols = lm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df)
      summary(ols)
      
  # 1.3. Lagrange multiplier test for spatial lag and spatial error dependencies
    lm.LMtests(ols, w, test=c("LMlag", "LMerr"))

## 2. Spatial analysis/regression and Moran's I (cross-section 2001-2016) ####
# run for contiguity and distance weight matrix (just change the 'w'), then compare the results
  # Moran's I test
    moran.test(df$lpw2016, w) #Consider to use 'zero.policy=T'
    moran.plot(df$lpw2016, w)
    
    #local = localmoran(x = df$lpw2016, listw = w) # Ii: local moran statistic; E.Ii: expectation of local moran statistic; Var.Ii: variance of local moran statistic; Z.Ii: standard deviate of local moran statistic
  
  # Spatial lag model
    slm = lagsarlm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w)
      summary(slm)
  
  # Spatial error model
    sem = errorsarlm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w)
      summary(sem)
      
  # Spatial lag of X
    slx = lmSLX(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w)
      summary(slx)
  
  # 2.1. Other models ####
    # SDEM Spatial Durbin Error Model (add lag X to SEM): y=XB+WxT+u, u=LWu+e
      reg5=errorsarlm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w, etype="emixed")
        summary(reg5)
        
    # SDM Spatial Durbin Model (add lag X to SAR): y=pWy+XB+WXT+e 
      reg6=lagsarlm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w, type="mixed")
      
    # Manski All-inclusive Model: y=pWy+XB+WXT+u, u=LWu+e (not recommended)
      reg7=sacsarlm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w, type="sacmixed") #listw2 allows for a different weights matrix for the error structure if desired
      
    # SARAR a.k.a. Kelejian-Prucha, Cliff-Ord, or SAC If all T=0,y=pWy+XB+u, u=LWu+e
      reg8=sacsarlm(log(regrowth_2001_2016) ~ log(lpw2016) +embargos_def, data = df, w, type="sac") #listw2 allows for a different weights matrix for the error structure if desired

# 3. Marginal effects (test aforementioned models) ####
  impacts(slx,listw=w) # Calculates de direct (changes on i affect i), indirect (changes on j affect i), and global effect (i+j)
    summary(impacts(slx,listw=w,R=500),zstats=TRUE) #Add zstats,pvals; R=500 not needed for SLX

# Hausman test
  
# Also, check the estimates of equilibrium effect (section 6a)
