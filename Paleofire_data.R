#install.packages("paleofire",repo="http://cran.r-project.org")
library(paleofire)

### select charcoal series between 30° and 90° latitude
### and -110° and -41° longitude, and include only those with 
### at least one geochronological (14C or 210Pb dating, tephra 
### layer, etc.) control point each 2500 year.

ID <- pfSiteSel(lat>30 & lat<90, long>-110 & long<(-41), date_int<=2500, num_version<400)
length(ID$id_site)

sumID <- summary(ID)
ID <- pfSiteSel(id_site %in% ID$id_site & num_samp>=20)
length(ID$id_site)

TR1 <- pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"))

####################
### Method 1
############
### Data-binning procedure used to calculate a composite curve.
### Binning sequence must be supplied for 'pfComposite' function to calculate
### 'mean' charcoal value per series within each binning interval
### Mean is then calculated for all series
### Below 500 yr bin width is used for span of 11,000 yrs

### CI for 'pfComposite' function calculated by bootstrap re-sampling of
### binned charcoal series and calculation of the mean for each bin n times,
### using the argument 'nboot' (by default => nboot=1000).

COMP1 <- pfComposite(TR2, binning=TRUE, bins=seq(from=0,to=11000, by=500))


####################
### Method 2
############
### Two-stage smoothing method:
### (i) "pre-bins" individual charcoal series using non-overlapping
### bins in order to ensure that records with high sample resolution do
### not have a disproportionate influence on the composite record, and
### that data will be not interpolated for records with a lower resolution
### (ii) smooth the pre-binned series using a locally-weighted scatter plot
### smoother "LOWESS" with a pre-defined constant bandwidth (half-width)
### given in the years (hw argument). Note: 20 yr non-overlapping bins and
### LOWESS smoother with 500 yr window half-width ==> 1000-yr smoothing window.
### 'tarAge' argument defines the center of each bin (binhw being the pre-binning
### bin half width) and ages for LOWESS estim. takes place (i.e. -50 to 12,000)

### 'pfCompositeLF' function ==> CIs calculated by bootstrap re-sampling
### the pre-binned series (as opposed to individual samples) and then applying
### the LOWESS curve fitting. This is repeated n times and is defined using 'nboot'
### CIs are estimated based on the distrib. of bootstrap replicates.

COMP2 <- pfCompositeLF(TR2, tarAge=seq(-50,12000,20), binhw=10, hw=500, nboot=100)


### Basic plot of data

par(mfrow=c(2,1))
plot(COMP1,conf=c(0.025,0.975),main="(a)")
plot(COMP2,conf=c(0.05,0.95),main="(b)")


### because the values in charcoal series are auto-correlated, testing the
### significance of their temporal variations is not straight-forward.
### Circular function ==> "moving" or "circular" block bootstrap, which
### consists of splitting each charcoal series into (n-b+1) overlapping blocks
### of data (n = sample size, b = block size). Blocks are randomly selected
### (with replacement) to reconstruct new individual charcoal series that are
### then used to estimate the CIs around the charcoal series composite mean.
### This procedure may dampen the long-term trends in data, particularly if
### individual sites record opposing trends
### If observed trend does not exceed the CIs then there the composite values
### at that time are not greater than by chance if records are not at all
### synchronized. However, a composite curve not exceeding the CIs does not 
### exclude the occurrence of important trends in the composite series.
### This procedure is used for testing significance of local minima or maxima
### in the composite time series.

### Circular block bootstrap method commonly used in Superposed Epoch Analysis
### a compositing technique that aims to find the response of a variable to
### one or multiple particular events (e.g. insect outbreaks, fire) using
### dendrochronological series or paleoecological proxy series

### Note: b=NULL argument indicates that the block size is automatically
### calculated for each series using the formulation above.

circboot <- pfCircular(COMP1, b=NULL, nboot=100, conf=c(0.005,0.025,0.975,0.995))
plot(circboot)



