KwaThunzi runoff time-series simulation
================

## Data import and preparation

Import the runoff dataset 1960-2018 of the Umkomazi gauging station (ID:
U1H005) available at the [Department of Water
Affairs](http://www.dwa.gov.za/Hydrology/). The runoff Q is in unit
\[m³/s\].

``` r
runoff <- read_table2("../data/U1H005_daily.txt", 
                      col_types = cols(DATE = col_integer()), 
                      skip = 10)
head(runoff, 3)
```

    ## # A tibble: 3 x 3
    ##       DATE     Q  QUAL
    ##      <int> <dbl> <dbl>
    ## 1 19600901  2.61     1
    ## 2 19600902  2.58     1
    ## 3 19600903  2.56     1

Clean the data and set column formats. The quality code indicating the
validity of each observation is available [here]() Create seperate
columns for ‘year’ and ‘doy’(day-of-year) for each observation.

``` r
runoff[is.na(runoff$QUAL),2] <- NA # remove Quality measures that shifted in the Q column
runoff <- ts.format(runoff, format="%Y%m%d", cols=c(1,2)) # set date format
runoff$year <- as.integer(substr(runoff$Date, 1,4)) # derive year
runoff$doy <- as.numeric(strftime(runoff$Date, format = "%j"))  # derive day-of-year
runoff <- na.omit(runoff) # delete empty fields
head(runoff, 3)
```

    ## # A tibble: 3 x 5
    ##   Date                    Q  QUAL  year   doy
    ##   <dttm>              <dbl> <dbl> <int> <dbl>
    ## 1 1960-09-01 00:00:00  2.61     1  1960   245
    ## 2 1960-09-02 00:00:00  2.58     1  1960   246
    ## 3 1960-09-03 00:00:00  2.56     1  1960   247

![](runoff_endless_experiment_files/figure-gfm/plot%20time%20series-1.png)<!-- -->

## Runoff transformation and aggregation

Transform discharge by base 10 logarithm.

``` r
runoff$Qlog <- log10(runoff$Q)  # calculate Log base 10 of runoff
runoff <- runoff[!(runoff$Qlog=='-Inf'),] # delete -Inf
head(runoff,3)
```

    ## # A tibble: 3 x 6
    ##   Date                    Q  QUAL  year   doy  Qlog
    ##   <dttm>              <dbl> <dbl> <int> <dbl> <dbl>
    ## 1 1960-09-01 00:00:00  2.61     1  1960   245 0.416
    ## 2 1960-09-02 00:00:00  2.58     1  1960   246 0.412
    ## 3 1960-09-03 00:00:00  2.56     1  1960   247 0.408

![](runoff_endless_experiment_files/figure-gfm/log%20tranformation%20plot-1.png)<!-- -->

Calculate daily mean and standard deviation for each day-of-year.

``` r
#Create empty data frame
runoff_aggregates <- data.frame()

#aggregate mean runoff per day
runoff_aggregates <- aggregate(runoff, by=list(runoff$doy), FUN=mean)[,c('doy','Qlog')]

#aggregate mean runoff per day
runoff_aggregates$sd <- aggregate(runoff, list(runoff$doy), FUN=sd)[,c('Qlog')]

names(runoff_aggregates) <- c('doy','mean','sd')
```

## Runoff simulation

Simulate expected runoff of each day of any timeseries. Create
randomized values, based on the observed annual discharge distribution,
which was derived from the time series 1960-2018. The random values are
estimated from mean and standard deviation of each daily runoff.

``` r
set.seed(2020)
number_of_years <- 500 # years to be simulated
runoff_pred <- data.frame()

for(i in runoff_aggregates$doy){
  x <- rnorm(number_of_years, mean=runoff_aggregates$mean[i], sd=runoff_aggregates$sd[i])
  pred <- data.frame(cbind(seq(1,500,1), # year
                         rep(i,500),   # doy
                         x))           # predicted values
  runoff_pred <- rbind(runoff_pred, pred)
}
names(runoff_pred) <- c('year','doy','Qlog')
```

    ## Warning: Removed 5 rows containing non-finite values (stat_bin2d).

![](runoff_endless_experiment_files/figure-gfm/simulation%20plot-1.png)<!-- -->
Backtransformation of the logarithmic runoffs to regular Q \[m³/s\]

``` r
runoff_pred$Q <- 10^runoff_pred$Qlog
```

## Calculation of area specific runoff and export to file

Calculate area specific runoff (q) after Baumgartner and Liebscher
(1996) and Casper and Bormann (2016) using  with Q = runoff at the
gauging station and A = catchment area

``` r
catchment_area <- 1744 * 1000000 # area in [m²]
runoff_pred$q <- runoff_pred$Q/catchment_area
```

Order rows by year and day-of-year

``` r
runoff_pred <- runoff_pred[order(runoff_pred$year, runoff_pred$doy),]
head(runoff_pred, 6)
```

    ##      year doy      Qlog         Q            q
    ## 1       1   1 1.4499444 28.180219 1.615838e-08
    ## 501     1   2 1.8259061 66.973978 3.840251e-08
    ## 1001    1   3 1.3288005 21.320654 1.222515e-08
    ## 1501    1   4 1.9601689 91.236559 5.231454e-08
    ## 2001    1   5 0.8988826  7.922871 4.542931e-09
    ## 2501    1   6 1.5452669 35.096750 2.012428e-08

Export data to a structured file with the two columns ‘year’ and ‘q’,
which is ready to use with the gully model

``` r
out <- runoff_pred[,c("year", "q")] 
out$q <- format(out$q, scientific = FALSE) # change to decimals
write.table(out, "../output/S_150_lts.txt", sep="\t", row.names=F, col.names=F, quote=F)
```

## References

BAUMGARTNER, A. & LIEBSCHER, H. J. 1996. Lehrbuch der Hydrologie, Band
1: Allgemeine Hydologie, Quantitative Hydrologie, Stuttgart, Germany,
Schweizerbarth. CASPER, M. & BORMANN, H. 2016. Abfluss im
Gewässersystem. In: FÖHRER, N., BORMANN, H., MIEGEL, K., CASPER, M.,
BRONSTERT, A., SCHUMANN, A. & WEILER, M. (eds.) Hydrologie. Bern,
Switzerland: Haupt.
