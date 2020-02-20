Soil data
=========

We start by importing the results of a soil texture analysis from a
table. The table contains a row for each layer with the following
columns:

-   Name
-   lab number
-   thickness in \[cm\]
-   texture in \[%\]
-   organic matter in \[%\]
-   comment

``` r
d <- read.csv("../data/soil_profile.csv", sep=";")
```

    ##           name labno thickness     clay    csilt    fsilt      csand     msand
    ## 1 Masotsheni 1  FS51        50 36.98730 13.00652 19.10333 11.5700224 2.4227825
    ## 2 Masotsheni 2  FS52       220 26.92769 20.03921 21.29166 21.3302412 4.6215523
    ## 3 Masotsheni 3  FS53        70 56.47385 10.34635 20.90826  0.8390795 0.6293097
    ## 4 Masotsheni 4  FS54       250 24.01935 15.30436 14.02900 35.5801914 5.4818318
    ##      fsand     vfsand som           comment
    ## 1 3.213895 13.6961376 4.8           Topsoil
    ## 2 5.383347  0.4062903 3.2         Paleosols
    ## 3 3.094106  7.7090433 1.2 Tansitional layer
    ## 4 5.016393  0.5688693 0.8        Base layer

Here are the results plotted in a USDA soil texture triangle
![](soil_parameter_preparation_files/figure-markdown_github/texture%20triangle-1.png)

Rill erodibility
================

Is calculated after Alberts et al. (1995) using the proportion of very
fine sand *V**F**S* \[%\] and the content of organic matter *O**M* \[%\]
with the equation

``` r
erodibility <- function(vfsand, om){
  0.00197+0.00030*vfsand+0.038633*exp(-1.84*om)
}
d$Kr <- erodibility(d$vfsand, d$som)
```

    ## [1] 0.006084481 0.002198998 0.008529263 0.011005624

Critical velocity
=================

Mean particle size diameter
---------------------------

We calculate the mean particle size diameter $Dg% in \[m\] using the
equation by Renard et al. (1997) with the primary fraction of the grain
size classes *f*<sub>*i*</sub> in \[%\] and the arithmetic mean of the
borders of the grain size classes *m*<sub>*i*</sub> in \[mm\].

``` r
mean_particle_diameter <- function(clay=0, silt=0, sand=0, csilt=0, fsilt=0, csand=0, msand=0, fsand=0, vfsand=0, method="geometric", classification="USDA"){
  if(classification=="USDA"){
    gsize <- c(mean(c(0, 0.002)),    # Clay
               mean(c(0.002, 0.05)), # Silt
               mean(c(0.05, 2)),   # Sand
               mean(c(0.05, 0.02)),  # Coarse silt
               mean(c(0.02, 0.002)), # Fine silt
               mean(c(2, 0.5)),      # Coarse sand
               mean(c(0.5, 0.25)),   # Medium sand
               mean(c(0.25, 0.1)),   # Fine sand
               mean(c(0.1, 0.05)))   # Very fine sand    
  }
  if(method=="geometric"){
    mean_diameter <- exp(0.01 * (
        clay * log(gsize[1]) +
        silt * log(gsize[2]) +
        sand * log(gsize[3]) +
        csilt * log(gsize[4]) +
        fsilt * log(gsize[5]) +
        csand * log(gsize[6]) +  
        msand * log(gsize[7]) +  
        fsand * log(gsize[8]) +  
        vfsand * log(gsize[9])
    ))
  }
  return(mean_diameter/1000) # from mm to m
}

d$Dg <- mean_particle_diameter(d$clay, d$csilt, d$fsilt, d$csand, d$msand, d$fsand, d$vfsand, method="geometric", classification="USDA")
```

Note, that the Mean particle diameter is converted from \[mm\] to \[m\]
for further processing

    ## [1] 2.601007e-05 3.014672e-05 1.229285e-05 2.602513e-05

Shear stress
------------

There are multiple ways to estimate Shear stress *τ*<sub>*c*</sub> from
soil texture. However, you’ll notice the large differences between the
results of the different equations. We’ll continue the method from
Alberts et al. (1995).

### after Smerdon (1961)

``` r
shear_stress <- function(clay){
  shear_stress <- 3.11*10^(0.0182*clay) # in Pa (Pascal)
}

d$tau <- shear_stress(d$clay)
```

    ## [1] 14.653094  9.612678 33.157705  8.509667

### after Alberts et al. (1995)

``` r
shear_stress <- function(clay, vfsand){
  shear_stress <- 2.67 + 0.065 * clay - 0.058 * vfsand
}

d$tau <- shear_stress(d$clay, d$vfsand)
```

    ## [1] 4.279799 4.396735 5.893676 4.198263

Critical velocity after Bogomolov & Mikhailov (1988)
----------------------------------------------------

Finally, we’ve collected all parameters to calculate the critical
velocity of erosion initiation *V*<sub>*c**r**i**t*</sub> \[m/s\] using
the equation from Bogomolov & Mikhailov (1988): The equation includes
some predefined parameter values, obtained from Märker (2001) and
Sidrochuk (2019, pers. comm.), comprising the colloidal suspended load
factor *m*<sub>1</sub>, the turbulence parameter *n*<sub>1</sub>, the
density of water *p*, the gravity constant *g* and the soil mechanical
coefficent *K*<sub>0</sub>.

``` r
v_crit <- function(m1=1,     # 1=clear water,...,1.4=colloidal suspended load > 0.1 kg/m3
                   n1=4,     # Turbulence: 4
                   ps,       # Density of soil kg/m3
                   p=1000,   # Density of water: ~1000 kg/m3
                   g=9.81,   # Newton/kg
                   d,        # mean particle size diameter (m)
                   K0=1,     # Soil mechanical coefficent: 0.5 (Maerker 2002) or 1 (Sidorchuk 2019)
                   Cnf       # Soil shear strength (Pa)
){
  # Critical velocity of bed erosion initiation (after 1pers. comm. Sidorchuck 2019)
  v_crit <- 1.25 * sqrt(
    (2*m1/1.54*p*n1)*
      ((ps-p)*g*d + 1.25*Cnf*K0)
  )
  return(v_crit)
}

d$v_crit <- v_crit(ps=2500, d=d$Dg, Cnf=d$tau)
```

Critical slope
==============

The critical slope is simply the slope measured in the field or obtained
from a high-res DEM coverted to radians.

``` r
s_crit <- function(slope, unit="degrees"){
  if(unit=="degrees"){
    s_crit <- slope*pi/180
  }else if(unit=="percent"){
    s_crit <- atan(slope/100)
  }
}
# Since all layers have the same slope in this example, we can assign the same value. Otherwise get the value from the table.
d$s_crit <- s_crit(57, "degrees")
```

    ## [1] 0.9948377 0.9948377 0.9948377 0.9948377

Export
======

Now we can export all relevant variables of each layer to an individual
PARAM.txt parameter file.

``` r
layer_export <- function(years, lnumber=1, v_crit, s_crit, Kr, Qdat="S_150.txt"){
  for (i in seq(1,lnumber)){
    filename <- paste("../output/PARAM", toString(i), ".txt", sep="")
    Layer <- file(filename, open="w")
    lines <- c("0\t0", 
               i,              # Number of layer
               v_crit[i],      # Critical Velocity
               "0",            # Index of texture layer (unused)
               s_crit[i],      # Stable slope
               Kr[i],          # Erodibility coefficent
               "S_150.txt",
               years,
               seq(1, years, 1))
    writeLines(lines, con=Layer)
    close(Layer)
    rm(i, filename, Layer, lines)
}}


# Apply
layer_export(years=59, lnumber=nrow(d), v_crit=d$v_crit, s_crit=d$s_crit, Kr=d$Kr)
```

References
==========

Alberts E. E., Nearing M. A., Weltz M. A., Risse L. M., Pierson F. B.,
Zhang X. C., Laflen J. M., and Simanton J. R. 1995, Soil component.
USDA-Water Erosion Prediction Project: Hillslope Profile and Watershed
Model Documentation. Chapter 7, pp. 7.1–7.45. Rep. No. 10. USDA-ARS
National Soil Erosion Research Laboratory, West Lafayette, IN.
[Open](./resources/Alberts_1995.pdf)

Bogomolov, A. I. and K. A. Mikhaylov (1972). Hydravlica. Moscow,
Stroyizdat.

Maerker, M. (2001). Regionale Erosionsmodellierung unter Verwendung des
Konzepts der Erosion Response Units (ERU) am Beispiel zweier
Flusseinzugsgebiete im südlichen Afrika. Chemisch-Geowissenschaftlichen
Fakultät. Jena, Friedrich-Schiller-Universität Jena. Dr. rer. nat.: 212.

Smerdon, E. T., & Beasley, R. P. (1961). Critical tractive forces in
cohesive soils. Agricultural Engineering, 42(1), 26-29.
