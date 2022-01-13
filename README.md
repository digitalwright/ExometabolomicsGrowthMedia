# ExometabolomicsGrowthMedia

<h2><b>ExometMedia_MS</b>: code for DIMS data analysis for  "Balancing trade-offs imposed by growth media and mass spectrometry for bacterial exometabolomics"</h2>

<hr>

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)

<hr> 

## Description

ExometMedia_MS is processing script designed to read in tables of m/z features (resulting from spectral stitching via [DIMSpy](https://github.com/computational-metabolomics/dimspy/blob/master/README.rst) and output metrics about that well/plate (e.g. number of unique m/z features, number of poor quality wells). ExometMedia_MS is intended to replicate the results found in the paper "Balancing trade-offs imposed by growth media and mass spectrometry for bacterial exometabolomics"; alterations will be required for general use.

## Installation

This script has the following minimum requirements:

1. R 3.6 or greater
2. In order to plot beeswarm plots, the package ['beeswarm'](https://cran.r-project.org/web/packages/beeswarm/index.html) is required. 

Input files used for the paper "Balancing trade-offs imposed by growth media and mass spectrometry for bacterial exometabolomics" can be found on on [GitHub](https://link). 


## Usage

This code can be run through terminal by setting the working directory, then running the R script:

```{bash}
#cd /Users/aed68/Desktop/ExometMedia_MS
#Rscript ExometMedia_MS.R
```


If you choose to run this in RStudio, you can run the script using source():
```{r}
#source("ExometMedia_MS.R")
```

Or alternatively, you can open the script in RStudio, set the working directory in line 13 of the script, and run it.
<hr>
The following global variables are defined in the first part of the script; the descriptions of each can be found below.

<b>max_mzG</b>: set this to the maximum m/z feature you want in your analysis. Default is 1400 m/z.

<b>min_SNR</b>: this sets the inclusion criteria for a given m/z feature in a plate. For example, if min_SNR = 10 (default), a feature must occur with an SNR >= 10 at least once in a plate to be included in the analysis.

<b>thresh</b>: this sets the threshold for inclusion as a background feature in blank wells. For example, if thresh = 0.5 (default), a feature must occur in at least half of the blank wells in order to be counted as a media feature, which is later subtracted as background.

<b>ppm</b>: Set the ppm tolerance window (+/- ppm) for searching compounds in NPAtlas. The default tolerance is set as +/- 7 ppm.

<hr>
Additionally, the script will print out several tables of data; descriptions of each can be found below.

<b>inj_summary</b> is a table summarizing the number of injection types.

<b>z_summary</b> is a table that summarizes the number of common features produced by each Streptomyces isolate.
A feature is included if it is detected in at least 50% of the isolate wells, after backgrdound subtraction.
Blank wells should have 0 features after background subtraction.

<b>samp_feats_means</b> is a table that summarizes the mean number of features that fall between the indicated m/z ranges.
This data is visualized in Figure 4E.

<b>bad</b> and <b>bad_frac</b> are tables of injection data generated from section 1D inj_summary.
These tables show the number of total injections (excluding QC injections), how many have empty windows of data,
and how many have windows of poor quality data (poor win SNR <1%).
This is used to generate Figure 4A.
