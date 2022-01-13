#####
## TITLE: 
## Balancing trade-offs imposed by growth media and mass spectrometry for bacterial exometabolomics 
##
## AUTHORS: 
## Ann E. Donnelly, Nithya Narayanan, Caroline M. C. Birer-Williams, Travis J. DeWolfe, 
## Rosalie K. Chu, Christopher R. Anderton, Erik S. Wright
##
## Script to process MS data and generate related figures.
#####

#set repo, if this is not already done
#setwd("~/Desktop/ExometabolomicsGrowthMedia")

#The package 'beeswarm' is required to generate figures
library(beeswarm)

#### Setting variables ####

#features must occur at least once in a plate with SNR greater than this value to be 
#included in analysis (Section 2B)
min_SNR <- 10

#for background subtraction: a feature must be present in at least this fraction of blank wells 
#to be counted (Section 5)
thresh <- 0.50   

#tolerance for NPAtlas search
ppm <- 7  #mass +/- this ppm

#truncating data - set max m/z here (for NPAtlas)
max_mzG <- 1400

#### A few other things ####
#pull in master sample organization on plates 2-8
s <- read.csv("./Samples.csv",
              stringsAsFactors=FALSE)

#name plates by media type
SAMPLES <- c(AED2="G-N-P-T-NA", AED3="S-N-P-T-NA", AED4="CAA-A", 
             AED5="AA-P-NA", AED6="ISP1-A", AED7="ISP2-A", AED8="NB-A")

#coloring for plots
#Sequential colors for 16 media types, by classification (Fig 2)
#If you want to change colors, do this here. 
#This is compatible with opt_analysis code for well darkness.
col16s <- c('#bdc9e1','#67a9cf','#02818a','#66c2a4','#238b45','#00441b',            #complex
            '#fec44f','#ec7014','#993404',                                          #defined + A
            '#B7A4A8','#fa9fb5','#BBA9BE', '#df65b0','#ce1256','#980043','#4d004b') #defined + NA

#reordering plates to match the rest of the OPT paper
re_ord <- c(1,2,4,3,5,6,7)
#these are the colors used for reordered plates analyzed in this MS analysis
plate_col <- c(14,15,16,4,2,3,1)




st <- Sys.time()

##### 1 Load in plate data #####
load("./dataR/pm2.RData")
load("./dataR/pm3.RData")
load("./dataR/pm4.RData")
load("./dataR/pm5.RData")
load("./dataR/pm6.RData")
load("./dataR/pm7.Rdata")
load("./dataR/pm8.Rdata")

load("./dataR/features.RData")

plate <- list(pm2,pm3,pm4,pm5,pm6,pm7,pm8)


#known issue with false start on plate 2, remove errant QC and A02 2 rows in the pm so only 105 rows not 107
plate[[1]] <- plate[[1]][c(-11, -82), ]



##### 1B Global alignment to QC exact masses #####
intense_peaks <- list()
intense_features <- list()
for (i in 1:length(plate)) {
  g <- grep("QC", rownames(plate[[i]]))
  intense_peaks[[i]] <-
    apply(plate[[i]][g, ], 1, function(x)
      which(x > max(x) * .2))
  intense_features[[i]] <-
    features[[i]][unique(unlist(intense_peaks[[i]]))]
}

tune_mix <- c(118.086255, 922.009798, 1221.990637)

tune_mix_ppm <- matrix(NA, nrow = length(plate), ncol = 3)
for (i in 1:length(plate)) {
  a1 <-
    intersect(which(intense_features[[i]] > 118),
              which(intense_features[[i]] < 119))
  a2 <-
    intersect(which(intense_features[[i]] > 922),
              which(intense_features[[i]] < 923))
  a3 <-
    intersect(which(intense_features[[i]] > 1221.9),
              which(intense_features[[i]] < 1222))
  tune_mix_ppm[i, 1] <-
    (intense_features[[i]][a1[1]] - tune_mix[1]) / (tune_mix[1]) * 1e6
  tune_mix_ppm[i, 2] <-
    (intense_features[[i]][a2[1]] - tune_mix[2]) / (tune_mix[2]) * 1e6
  tune_mix_ppm[i, 3] <-
    (intense_features[[i]][a3[1]] - tune_mix[3]) / (tune_mix[3]) * 1e6
}
row.names(tune_mix_ppm) <- SAMPLES
colnames(tune_mix_ppm) <- tune_mix

#mean ppm shift
tmp <- apply(tune_mix_ppm, 1, mean)

#plate correction by mean of 3 tune mix peaks
for (i in 1:length(plate)) {
  colnames(plate[[i]]) <- features[[i]] + (features[[i]] * (tmp[i] * 10 ^
                                                              -6))
  features[[i]] <- features[[i]] + (features[[i]] * (tmp[i] * 10 ^ -6))
}


##### 1C PEG removal #####
print("Cleaning up spectra...")
pBar <- txtProgressBar(style = 3)
# work with original plates prior to background subtraction or SNR choice
p44_plate <- rep(list(list()), length(plate))
p44_feats <- rep(list(list()), length(plate))

for (i in 1:length(plate)) {
  #create a list of features for each well of a plate
  p1 <-
    apply(plate[[i]], 1, function(x)
      which(x != 0)) #this keeps the original indexing so you can pull the original feats back
  p2 <-
    lapply(p1, function(x)
      as.numeric(names(x))) #these are just m/z
  
  p44 <- list()
  for (j in seq_along(p2)) {
    # go through rows in plate
    p3 <-
      round(outer(p2[[j]], p2[[j]], FUN = "-"), 3) # subtract all the m/z features
    p3[which(p3 < 44.01)] <- 0 #zero in on the cadence at 44 m/z
    p3[which(p3 > 44.03)] <- 0
    
    #find all the consecutive starts and stops, by row
    p44 <- which(p3 != 0, arr.ind = TRUE)
    
    #intensity matrix - get rid of outliers (next peak must be at least 50% of previous)
    p44_int <- matrix(plate[[i]][j, p1[[j]][p44[, 1:2]]], ncol = 2)
    if (nrow(p44_int) != 1) {
      #need at least 2 rows to apply
      pp <-
        which(apply(p44_int, 1, function(x)
          min(x) > max(x) / 2) == TRUE) #adjacent peak must be >50% intensity
      p44 <- p44[pp, ]
    }
    
    #find all the features that occur both as a start and a stop
    p44_feats <- intersect(p44[, 1], p44[, 2])
    
    #find first and last features of a distribution
    p44_other <-
      p44[-match(p44_feats, p44)] #remove all the features that have been called already
    starts <-
      p44_other[which(round(p2[[j]][p44_other] + 88.052, 2) %in% round(p2[[j]][p44], 2) == TRUE)]
    stops <-
      p44_other[which(round(p2[[j]][p44_other] - 88.052, 2) %in% round(p2[[j]][p44], 2) == TRUE)]
    
    p44_feats <- unique(c(p44_feats, starts, stops))
    
    p44_plate[[i]][[j]] <-
      p1[[j]][p44_feats] #original feature numbers
  }
  
  setTxtProgressBar(pBar, i / length(plate))
}
close(pBar)

#unique features by plate
p_poly <- list()
for (i in 1:length(plate)) {
  p_poly[[i]] <-
    unique(unlist(p44_plate[[i]])) #features[[i]][p_poly[[i]]] if you want mw
}

#remove polymer from plates - subtract any polymer feature BY EACH WELL (not global across plate)
for (i in 1:length(plate)) {
  for (j in 1:nrow(plate[[i]])) {
    plate[[i]][j, p44_plate[[i]][[j]]] <- 0
  }
}


##### 1D ID failed injections and poor quality data #####
d <- list()
win <-
  c(70 + 170 * 0:(ceiling((max_mzG - 70) / 170))) #this will change depending on what max_mzG is set to

for (i in 1:length(plate)) {
  missingwindows <- vector()
  min_mz <- vector()
  max_mz <- vector()
  mw_start <- vector()
  mw_stop <- vector()
  winSNR <- matrix(NA, nrow = nrow(plate[[i]]), ncol = length(win) - 1)
  
  for (j in 1:nrow(plate[[i]])) {
    w <- features[[i]][which(plate[[i]][j, ] != 0)]
    
    if (length(w) < 2) {
      next
    } else{
      missingwindows[j] <- max(diff(w)) #maxiumum distance between 2 peaks
      mw_start[j] <- w[which.max(diff(w))] #start of each mw
      mw_stop[j] <- w[which.max(diff(w)) + 1] #stop of each mw
      min_mz[j] <- min(w) #overall start
      max_mz[j] <- max(w) #overall stop
      
      for (k in 1:(length(win) - 1)) {
        r <-
          range(which(features[[i]] >= win[k] &
                        features[[i]] < win[k + 1])) #start and stop of features
        winSNR[j, k] <- sum(plate[[i]][j, r[1]:r[2]]) #sum of the SNR
      }
    }
  }
  
  d[[i]] <-
    data.frame(rownames(plate[[i]]), missingwindows) #empty windows
  d[[i]]$min <- min_mz #first feature
  d[[i]]$max <- max_mz #last feature
  d[[i]]$mw_start <- mw_start #start of missing window
  d[[i]]$mw_stop <- mw_stop #stop of missing window
  
  d[[i]]$poorSNR_1percent <-
    apply(winSNR, 1, function(x)
      any(x < apply(winSNR, 2, median) * 0.01)) #median by window
}


#look through all of the wells in each plates for gaps that are >170 (stiched windows) or missed stops/starts
d2 <- list()
d3 <-
  list() #which wells have windows with poor SNR (sum <1% of mean of any window, by media)
mw <- 170
for (i in seq_along(d)) {
  u <- unique(c(
    which(d[[i]]$missingwindows > mw),
    which(d[[i]]$min > (win[1] + 170 - 1)),
    which(d[[i]]$max < win[length(win)] - 170 + 1)
  )) #the way the windows are set up
  d2[[i]] <- d[[i]][u, ]
  d3[[i]] <- d[[i]][which(d[[i]]$poorSNR_1percent == TRUE), ]
}

####correction!!
d4 <- list()
for (i in seq_along(d)) {
  u <- unique(c(
    which(d[[i]]$missingwindows > 170),
    which(d[[i]]$min > (win[2] - 1)),
    which(d[[i]]$max < win[length(win) - 1] - 1)
  )) #the way the windows are set up
  d4[[i]] <- d[[i]][u, ]
}

#a list of all of the bad samples for each plate
bad_samples <- lapply(lapply(d2, rownames), as.numeric)
poorSNR_samples_1percent <- lapply(lapply(d3, rownames), as.numeric)
poorSNR_samples_noempty <- list() #1percent only
for (i in 1:length(plate)) {
  poorSNR_samples_noempty[[i]] <-
    poorSNR_samples_1percent[[i]][which(poorSNR_samples_1percent[[i]] %in% bad_samples[[i]] == FALSE)]
}


#looking at samples/blank only
bad_samples_sb <- list()
poorSNR_samples_1percent_sb <- list()
poorSNR_samples_noempty_sb <- list()
for (i in 1:length(plate)) {
  g <- grep("QC", rownames(plate[[i]]))
  if (length(which(bad_samples[[i]] %in% g == TRUE)) == 0) {
    bad_samples_sb[[i]] <- bad_samples[[i]]
  } else{
    bad_samples_sb[[i]] <-
      bad_samples[[i]][-which(bad_samples[[i]] %in% g == TRUE)]
  }
  if (length(which(poorSNR_samples_1percent[[i]] %in% g == TRUE)) == 0) {
    poorSNR_samples_1percent_sb[[i]] <- poorSNR_samples_1percent[[i]]
  } else{
    poorSNR_samples_1percent_sb[[i]] <-
      poorSNR_samples_1percent[[i]][-which(bad_samples[[i]] %in% g == TRUE)]
  }
  if (length(which(poorSNR_samples_noempty[[i]] %in% g == TRUE)) == 0) {
    poorSNR_samples_noempty_sb[[i]] <- poorSNR_samples_noempty[[i]]
  } else{
    poorSNR_samples_noempty_sb[[i]] <-
      poorSNR_samples_noempty[[i]][-which(bad_samples[[i]] %in% g == TRUE)]
  }
}

#summarizing feature counts, by well
feat_bs <- rep(list(vector()), length(plate)) #empty windows
feat_bspsnr <-
  rep(list(vector()), length(plate)) #empty windows and low snr 1%
feat_psnr <- rep(list(vector()), length(plate)) #low snr 1% only
for (i in 1:length(plate)) {
  m = 1
  n = 1
  o = 1
  for (j in bad_samples[[i]]) {
    feat_bs[[i]][m] <- length(which(plate[[i]][j, ] != 0)) # features
    m = m + 1
  }
  for (k in poorSNR_samples_1percent[[i]]) {
    feat_bspsnr[[i]][n] <- length(which(plate[[i]][k, ] != 0)) # features
    n = n + 1
  }
  
  rrr <-
    poorSNR_samples_1percent[[i]][which(poorSNR_samples_1percent[[i]] %in% bad_samples_sb[[i]] == FALSE)]
  for (l in rrr) {
    feat_psnr[[i]][o] <- length(which(plate[[i]][l, ] != 0)) # features
    o = o + 1
  }
}

####only exclude empty window runs that have < 500 total features
feats500 <- lapply(feat_bs, function(x)
  which(x > 500))
feat_psnr500 <- lapply(feat_psnr, function(x)
  which(x > 500))

bad_samples500 <- list()
poorSNR_samples_noempty500 <- list()
for (i in 1:length(plate)) {
  bad_samples500[[i]] <- bad_samples[[i]][-feats500[[i]]]
  poorSNR_samples_noempty500[[i]] <-
    poorSNR_samples_noempty[[i]][-feat_psnr500[[i]]]
}

bad_samples500sb <- list()
poorSNR_samples_noempty500sb <- list()
for (i in 1:length(plate)) {
  bad_samples500sb[[i]] <- bad_samples_sb[[i]][-feats500[[i]]]
  poorSNR_samples_noempty500sb[[i]] <-
    poorSNR_samples_noempty_sb[[i]][-feat_psnr500[[i]]]
}

#which injections to remove
bad_remove <- list()
for (i in 1:length(plate)) {
  bad_remove[[i]] <-
    c(bad_samples500[[i]], poorSNR_samples_noempty500[[i]])
}

bad_remove[which(lapply(bad_remove, length) == 0)] <- NA


##### 2 Remove poor quality injections from analysis #####

#injection summary table (build part before QC, and fill in rest later)
inj_summary <- matrix(NA, ncol = 10, nrow = length(plate))
inj_summary[, 1] <- sapply(plate, nrow)
inj_summary[, 2] <-
  sapply(plate, function(x)
    length(grep("QC", rownames(x))))
inj_summary[, 3] <-
  sapply(plate, function(x)
    length(grep("blank", rownames(x))))
inj_summary[, 4] <-
  inj_summary[, 1] - (inj_summary[, 2] + inj_summary[, 3])
inj_summary[, 5] <- inj_summary[, 1] - inj_summary[, 2]
colnames(inj_summary) <-
  c(
    "injections",
    "QC",
    "blank",
    "samples",
    "s+b",
    "clean injections",
    "clean QC",
    "clean blank",
    "clean samples",
    "clean s+b"
  )
row.names(inj_summary) <- SAMPLES

####adjust plate 6 and 8 - two duplicate sample injections.
#the bad runs get removed in bad_samples
inj_summary[c(5, 7), 4:5] <- inj_summary[c(5, 7), 4:5] - 1

#get rid of bad runs on each plate according to bad_samples row list - 
#this is missed windows >70 m/z (including starts/stops)
for (i in 1:length(plate)) {
  if (is.na(bad_remove[[i]][1]) == TRUE) {
    plate[[i]] <- plate[[i]]
  } else{
    plate[[i]] <- plate[[i]][-bad_remove[[i]], ]
  }
}








#finish the summary table
inj_summary[, 6] <- sapply(plate, nrow)
inj_summary[, 7] <-
  sapply(plate, function(x)
    length(grep("QC", rownames(x))))
inj_summary[, 8] <-
  sapply(plate, function(x)
    length(grep("blank", rownames(x))))
inj_summary[, 9] <-
  inj_summary[, 6] - (inj_summary[, 7] + inj_summary[, 8])
inj_summary[, 10] <- inj_summary[, 6] - inj_summary[, 7]

inj_summary <- as.data.frame(inj_summary)

##### 2A removing known expt issues #####
#specific data points with issues upon experimental setup
plate[[4]] <-
  plate[[4]][-grep("blank_D01", row.names(plate[[4]])), ] #something grew in this blank well, remove
plate[[6]] <-
  plate[[6]][-grep("SC_G04", row.names(plate[[6]])), ] #this agar pad fell out (S. coelicolor)


print("inj_summary")
print(inj_summary)
print("inj_summary is a table summarizing the number of injection types.")

##### 2B setting low SNR, removing features that don't meet this threshold #####
#remove features that don't occur at least once in a plate with SNR > min_SNR
for (i in 1:length(plate)) {
  w <- which(apply(plate[[i]], 2, max) < min_SNR)
  plate[[i]] <- plate[[i]][, -w]
  features[[i]] <- features[[i]][-w]
}


##### 3 deisotoping #####
iso1 <- list()
for (i in seq_along(plate)) {
  ppp <- outer(features[[i]], features[[i]], FUN = "-")
  ppp[which(ppp > 1.00341)] <- 0
  ppp[which(ppp < 1.00331)] <- 0
  
  iso1[[i]] <- which(ppp != 0, arr.ind = TRUE)[, 1]
  
  ppp <- 0
}

#get rid of 13C +1 isotopes
for (i in seq_along(plate)) {
  plate[[i]][, iso1[[i]]] <- 0 #set all 13C isotopes to 0
}



##### 4 BG/QC summary #####
#creating common features/all features lists for each plate in QC and blanks wells

#for each plate, create a summary features vector as 1 and 0 (presence/absence)
#if ANY of the samples records a given feature, count as a 1
QC_features_all <- features
BG_features_all <- features

for (i in 1:length(plate)) {
  QC <- grep("QC", rownames(plate[[i]]))
  
  for (k in seq_along(features[[i]])) {
    if (any(plate[[i]][QC, k] > 0)) {
      QC_features_all[[i]][k] <- 1
    } else{
      QC_features_all[[i]][k] <- 0
    }
  }
  
  BG <- grep("blank", rownames(plate[[i]]))
  
  for (k in seq_along(features[[i]])) {
    if (any(plate[[i]][BG, k] > 0)) {
      BG_features_all[[i]][k] <- 1
    } else{
      BG_features_all[[i]][k] <- 0
    }
  }
}




##### 5 BG/QC subtraction #####
#do this on the plate level, features level
#initialize lists identical to previous lists, so we can subtract features out accordingly without messing up original lists
plate_QCBG_all <- plate
features_QCBG_all <- features

#the intersect of the QC and BG peaks doesn't matter (e.g. don't need to take unique only)
for (i in 1:length(plate)) {
  plate_QCBG_all[[i]] <-
    plate[[i]][, c(-which(QC_features_all[[i]] == 1), -which(BG_features_all[[i]] == 1))]
  features_QCBG_all[[i]] <-
    features[[i]][c(-which(QC_features_all[[i]] == 1), -which(BG_features_all[[i]] == 1))]
}


#fraction of features from media
inj_ba <-
  list() #number of features in each well that were coming from the media
inj_ba_before <- list()
inj_ba_after <- list()
inj_ba_media <- list()
for (i in 1:length(plate)) {
  g <- grep("QC|blank", rownames(plate_QCBG_all[[i]]))
  
  inj_ba_before[[i]] <-
    apply(plate[[i]][-g, ], 1, function(x)
      length(which(x > 0))) #features before sub
  inj_ba_after[[i]] <-
    apply(plate_QCBG_all[[i]][-g, ], 1, function(x)
      length(which(x > 0))) #after media sub
  inj_ba_media[[i]] <-
    inj_ba_before[[i]] - inj_ba_after[[i]] #features due to media
  
  inj_ba[[i]] <- 1 - (inj_ba_after[[i]] / inj_ba_before[[i]])
}



##### 6 Media features summary #####
#look at common features in at least 50% of injections
samps <- s[25:48, 2]
z <- rep(list(rep(list(NA), 24)), length(plate)) #raw feature counts
z_feats <-
  rep(list(rep(list(NA), 24)), length(plate)) #common features vector
z_feats_count <-
  rep(list(rep(list(NA), 24)), length(plate)) #common features count

#chose which subtracted plate/features you want to use
feat_set <- features_QCBG_all
plate_set <- plate_QCBG_all


for (i in 1:length(plate)) {
  for (j in seq_along(samps)) {
    zz <-
      grep(samps[j], row.names(plate[[i]])) #call rows of a given sample j
    if (length(zz) == 0) {
      #if there are no rows, go to the next j
      print(paste(i, j))
      next
    } else{
      #find feature counts in individual injections
      for (k in seq_along(zz)) {
        z[[i]][[j]][k] <- length(which(plate_set[[i]][zz[k], ] > 0))
      }
      
      
      #find common features
      for (m in seq_along(feat_set[[i]])) {
        if ((length(which(plate_set[[i]][zz, m] > 0)) / length(zz)) >= thresh) {
          z_feats[[i]][[j]][m] <- 1
        } else{
          z_feats[[i]][[j]][m] <- 0
        }
      }
    }
  }
  
  z_feats_count[[i]] <-
    lapply(z_feats[[i]], function(x)
      length(which(x == 1)))
  
}

#summary of common feature counts
z_summary <- matrix(NA, ncol = 24, nrow = length(plate))
for (i in 1:length(plate)) {
  z_summary[i, ] <- unlist(z_feats_count[[i]])
}

row.names(z_summary) <- SAMPLES
colnames(z_summary) <- samps

print("z_summary")
print(z_summary)
print("z_summary is a table that summarizes the number of common features produced by each Streptomyces isolate.")
print("A feature is included if it is detected in at least 50% of the isolate wells, after backgrdound subtraction.")
print("Blank wells should have 0 features after background subtraction.")



#summarizing by media type

#making a list of all QC/BG subtracted injections, by media, excluding background injections
z_media <- list()
for (i in 1:length(plate)) {
  z_media[[i]] <- unlist(z[[i]][-13])
}

#create a list of common features for each plate, by sample
ss <- list()
for (i in 1:length(plate)) {
  ss[[i]] <- lapply(z_feats[[i]], function(x)
    which(x == 1))
}



#summary of common features
samp_feats <- matrix(NA, ncol = length(plate), nrow = 24)
colnames(samp_feats) <- SAMPLES
row.names(samp_feats) <- s[25:48, 2]
for (i in 1:length(plate)) {
  for (j in 1:24) {
    samp_feats[j, i] <- length(features_QCBG_all[[i]][ss[[i]][[j]]])
  }
}




#any other feature m/z thresholds
samp_feats2 <- matrix(NA, ncol = length(plate), nrow = 24)
colnames(samp_feats2) <- SAMPLES
row.names(samp_feats2) <- s[25:48, 2]

greater_than <- c(70, 350, 700, 1050)
less_than <- c(350, 700, 1050, 1400)
samp_feats_ranges <- list()

for (a in seq_along(greater_than)) {
  for (i in 1:length(plate)) {
    for (j in 1:24) {
      samp_feats2[j, i] <-
        length(intersect(
          which(as.numeric(features_QCBG_all[[i]][ss[[i]][[j]]]) >= greater_than[a]),
          which(as.numeric(features_QCBG_all[[i]][ss[[i]][[j]]]) < less_than[a])
        ))
    }
  }
  #save a selected bunch of ranges
  samp_feats_ranges[[a]] <- samp_feats2
}

#calculating mean features in common feature vectors, that fall within a range of masses
samp_feats_means <- list()
for (i in seq_along(greater_than)) {
  samp_feats_means[[i]] <- colMeans(samp_feats_ranges[[i]][-13, ])
}
samp_feats_means <- do.call(rbind, samp_feats_means)
row.names(samp_feats_means) <- paste(greater_than, "to", less_than)

print("samp_feats_means")
print(samp_feats_means)
print("samp_feats_means is a table that summarizes the mean number of features that fall between the indicated m/z ranges.")
print("This data is visualized in Figure 4E.")

##### 7 Adduct correction #####

pAd_plate <- list()

for (i in 1:length(plate)) {
  #create a list of features for each well of a plate
  p1 <-
    apply(plate_QCBG_all[[i]], 1, function(x)
      which(x != 0)) #this keeps the original indexing so you can pull the original feats back
  p2 <-
    lapply(p1, function(x)
      as.numeric(names(x))) #these are just m/z
  
  pAd <- list()
  for (j in seq_along(p2)) {
    # go through rows in plate
    p3 <-
      round(outer(p2[[j]], p2[[j]], FUN = "-"), 3) # subtract all the m/z features
    p3[which(p3 < 15)] <- 0
    p3[which(p3 > 40)] <- 0
    
    #find all adduct pairs
    pg1 <- which((p3 > 15.97 & p3 < 15.98) == TRUE, arr.ind = TRUE)
    pg2 <- which((p3 > 21.98 & p3 < 21.99) == TRUE, arr.ind = TRUE)
    pg3 <- which((p3 > 38.95 & p3 < 38.96) == TRUE, arr.ind = TRUE)
    
    #how many features should be removed by collapsing adducts?
    col1 <- length(unique(c(pg2[, 2], pg3[, 2]))) #number of unique H+ pairs
    
    #number of unique that don't occur as an H+ pair
    col2 <- nrow(pg1) - #number of K-Na pairs
      length(intersect(pg1[, 1], pg3[, 1])) - #K occurs as an H+ pair
      length(intersect(pg1[, 2], pg2[, 1])) #Na occurs as an H+ pair
    
    pAd[[j]] <- col1 + col2 #adduct count that should be removed from total feature #
    
  }
  pAd_plate[[i]] <- pAd
}


#create corrected feature counts
pAd_corrected <- rep(list(list()), length(plate))


for (i in 1:length(plate)) {
  for (j in 1:nrow(plate_QCBG_all[[i]])) {
    if (pAd_plate[[i]][[j]] == 0) {
      pAd_corrected[[i]][[j]] <- NA
    } else{
      pAd_corrected[[i]][[j]] <-
        length(which(plate_QCBG_all[[i]][j, ] != 0)) - #how many features in the well
        pAd_plate[[i]][[j]] #adduct collapsing - subtract this number of features
    }
  }
}




##### 8 NPAtlas #####
print("Searching Natural Product Atlas...")
pBar <- txtProgressBar(style = 3)
##### Load NPAtlas
load("./NPAtlas.RData")

w <- which(x$genus == "Streptomyces")
x <- x[w, ]

x2 <-
  data.frame(x$compound_m_plus_h, x$compound_m_plus_na) # combine adducts
x2$x.compound_m_plus_k <-
  x$compound_accurate_mass + 38.9632 #add a column for k adducts
x2 <- apply(x2, 2, unique) #unique features only
x2 <- x2[-which(x2[, 3] > max_mzG), ]

#x2 <- nonstrepto

##### find NPAtlas hits in each well of each plate
npa_rows5 <- rep(list(rep(list(NA), 1)), length(plate))

pd <- 1 / ppm * 1e6

for (i in 1:length(plate)) {
  a <-
    lapply(apply(plate_QCBG_all[[i]], 1, function(x)
      which(x != 0)),
      function(x)
        features_QCBG_all[[i]][x]) #find all features in each well of a plate
  
  
  mm1 <- matrix(NA, ncol = 105, nrow = nrow(x2))
  mm2 <- matrix(NA, ncol = 105, nrow = nrow(x2))
  mm3 <- matrix(NA, ncol = 105, nrow = nrow(x2))
  for (k in 1:nrow(x2)) {
    #for each of the sets of adducts in NPAtlas (rows)
    ah <- x2[k, 1] #H+
    ana <- x2[k, 2] #Na+
    ak <- x2[k, 3] #K+
    
    #find which plate features a match and NPA feature k
    m1 <-
      which(sapply(a, function(p)
        any(ah - (ah / pd) <= p & ah + (ah / pd) >= p)) == TRUE)
    m2 <-
      which(sapply(a, function(p)
        any(ana - (ana / pd) <= p & ana + (ana / pd) >= p)) == TRUE)
    m3 <-
      which(sapply(a, function(p)
        any(ak - (ak / pd) <= p & ak + (ak / pd) >= p)) == TRUE)
    
    #take the trues, find the NPA features
    mm1[k, m1] <- TRUE
    mm2[k, m2] <- TRUE
    mm3[k, m3] <- TRUE
  }
  mm11 <- apply(mm1, 2, function(x) which(x == TRUE))
  mm22 <- apply(mm2, 2, function(x) which(x == TRUE))
  mm33 <- apply(mm3, 2, function(x) which(x == TRUE))
  
  mmm <- list()
  for (j in 1:nrow(plate_QCBG_all[[i]])) {
    mmm[[j]] <-
      table(c(mm11[[j]], mm22[[j]], mm33[[j]])) #table what NPA features occur with multipe adducts in a well
  }
  
  npa_rows5[[i]] <-
    lapply(mmm, function(x)
      as.numeric(names(which(x >= 2)))) #at least 2 features
  
  setTxtProgressBar(pBar, i / length(plate))
  
}

close(pBar)

#save this as "real" to compare with next section
npa_rows5_real <- npa_rows5


##### fraction of NPAtlas covered in a plate
npa_frac5 <- list()
for (i in 1:length(plate)) {
  npa_frac5[[i]] <- length(unique(unlist(npa_rows5[[i]]))) / nrow(x2)
}


##### total feature count per well
#samples only
npa_feats5 <- list()
for (i in 1:length(plate)) {
  g <- grep("QC|blank", rownames(plate_QCBG_all[[i]]))
  npa_feats5[[i]] <- unlist(lapply(npa_rows5[[i]], length))[-g]

}


#NPAtlas features as a fraction of total feature count by well
#by well

npa_feats5_perwell <- list()
for (i in 1:length(plate)) {
  g <- grep("QC|blank", rownames(plate_QCBG_all[[i]]))
  a <- unlist(pAd_corrected[[i]][-g]) #used the feature count adjusted for adducts
  npa_feats5_perwell[[i]] <- npa_feats5[[i]] / a
}




#### 8A NPAtlas selectivity ####
load("./NPAtlas.RData") #reload to start with original unedited NPA, this is x
x2 <- data.frame(x$compound_m_plus_h, x$compound_m_plus_na) # combine adducts
x2$x.compound_m_plus_k <- x$compound_accurate_mass + 38.9632 #add a column for k adducts

#divide into two sets - streptos, and non-streptos
strepto <- x2[which(x$genus == "Streptomyces"),] #4830
nonstrepto <- x2[-which(x$genus == "Streptomyces"),] #24176

#remove the compounds that are over 1400 m/z
strepto <- strepto[-which(strepto[, 3] > max_mzG), ] #4693
nonstrepto <- nonstrepto[-which(nonstrepto[, 3] > max_mzG), ] #23665

#keep only unique masses in each set
strepto <- apply(strepto, 2, unique) #3400
nonstrepto <- apply(nonstrepto, 2, unique) #9801

#identify masses that are in both sets
s1 <- which(strepto[,1] %in% nonstrepto[,1] == TRUE)
n1 <- which(nonstrepto[,1] %in% strepto[,1] == TRUE)

#remove masses
#strepto <- strepto[-s1,] #2324
nonstrepto <- nonstrepto[-n1,] #8725

x2 <- nonstrepto

print("Still determining selectivity...")
pBar <- txtProgressBar(style = 3)

##### find NPAtlas hits in each well of each plate
npa_rows5 <- rep(list(rep(list(NA), 1)), length(plate))

pd <- 1 / ppm * 1e6

for (i in 1:length(plate)) {
  a <-
    lapply(apply(plate_QCBG_all[[i]], 1, function(x)
      which(x != 0)),
      function(x)
        features_QCBG_all[[i]][x]) #find all features in each well of a plate
  
  
  mm1 <- matrix(NA, ncol = 105, nrow = nrow(x2))
  mm2 <- matrix(NA, ncol = 105, nrow = nrow(x2))
  mm3 <- matrix(NA, ncol = 105, nrow = nrow(x2))
  for (k in 1:nrow(x2)) {
    #for each of the sets of adducts in NPAtlas (rows)
    ah <- x2[k, 1] #H+
    ana <- x2[k, 2] #Na+
    ak <- x2[k, 3] #K+
    
    #find which plate features a match and NPA feature k
    m1 <-
      which(sapply(a, function(p)
        any(ah - (ah / pd) <= p & ah + (ah / pd) >= p)) == TRUE)
    m2 <-
      which(sapply(a, function(p)
        any(ana - (ana / pd) <= p & ana + (ana / pd) >= p)) == TRUE)
    m3 <-
      which(sapply(a, function(p)
        any(ak - (ak / pd) <= p & ak + (ak / pd) >= p)) == TRUE)
    
    #take the trues, find the NPA features
    mm1[k, m1] <- TRUE
    mm2[k, m2] <- TRUE
    mm3[k, m3] <- TRUE
  }
  mm11 <- apply(mm1, 2, function(x) which(x == TRUE))
  mm22 <- apply(mm2, 2, function(x) which(x == TRUE))
  mm33 <- apply(mm3, 2, function(x) which(x == TRUE))
  
  mmm <- list()
  for (j in 1:nrow(plate_QCBG_all[[i]])) {
    mmm[[j]] <-
      table(c(mm11[[j]], mm22[[j]], mm33[[j]])) #table what NPA features occur with multipe adducts in a well
  }
  
  npa_rows5[[i]] <-
    lapply(mmm, function(x)
      as.numeric(names(which(x >= 2)))) #at least 2 features
  
  setTxtProgressBar(pBar, i / length(plate))
  
}

close(pBar)

npa_rows5_fake <- npa_rows5





strepto_feats <- lapply(lapply(npa_rows5_real, unlist), unique)
nonstrepto_feats <- lapply(lapply(npa_rows5_fake, unlist), unique)

strepto_table <- matrix(0, nrow=3400, ncol=7)
for(j in 1:length(plate)){
strepto_table[strepto_feats[[j]],j] <- 1
}
  

nonstrepto_table <- matrix(0, nrow=8725, ncol=7)
for(j in 1:length(plate)){
  nonstrepto_table[nonstrepto_feats[[j]],j] <- 1
}
  

npa_hits <- data.frame(colSums(strepto_table)[re_ord], colSums(nonstrepto_table)[re_ord])
colnames(npa_hits) <- c("strepto", "nonstrepto")
row.names(npa_hits) <- SAMPLES

r <- npa_hits

#bootstrap replication
tot <- c(3400, 8725) # total compounds
reps <- 1e7 # bootstrap replicates
quants <- c(Lower=0.025, Median=0.5, Upper=0.975)

stats <- array(NA_real_,
               c(nrow(r), length(quants), 3),
               dimnames=list(rownames(r),
                             names(quants),
                             c("TPRs", "FPRs", "Log-odds")))
for (i in seq_len(nrow(r))) {
  #print(rownames(r)[i])
  
  TPRs <- rpois(reps, r[i, 1])/tot[1]
  FPRs <- rpois(reps, r[i, 2])/tot[2]
  logodds <- log(TPRs/FPRs)
  #print(sum(logodds <= 0)/reps*nrow(r)) # p-adj
  
  stats[i,, "TPRs"] <- quantile(TPRs, quants)
  stats[i,, "FPRs"] <- quantile(FPRs, quants)
  stats[i,, "Log-odds"] <- quantile(logodds, quants)
}





##### FIGURES #####

##### Fig 4A  #####
#bad spec summary, samples/blank only
bad <- data.frame(inj_summary[,5],
                  unlist(lapply(bad_samples_sb,length)), 
                  unlist(lapply(poorSNR_samples_1percent_sb, length)), #these can inclue empty windows
                  unlist(lapply(poorSNR_samples_noempty_sb, length)))
colnames(bad) <- c("injections", "empty windows", "poor win SNR <1%", "p1only")
row.names(bad) <- SAMPLES

bad <- as.matrix(bad)

#what fraction of 95 injections
bad_frac <- data.frame(1-bad[,1]/95,
                       bad[,2]/95,
                       bad[,4]/95) #pull only the poor SNR specs that are not included in empty windows
colnames(bad_frac) <- c("injections","empty windows", "poor win SNR <1%")
row.names(bad_frac) <- SAMPLES

print("bad")
print(bad)
print("bad_frac")
print(bad_frac)
print("bad and bad_frac are tables of injection data generated from section 1D inj_summary")
print("These tables show the number of total injections (excluding QC injections), how many have empty windows of data,")
print("and how many have windows of poor quality data (poor win SNR <1%).")
print("This is used to generate Figure 4A.")


#out of 95 injections (excludes QC runs)
dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 4, 2, 1), xpd=FALSE, lwd=8)
barplot(colSums(t(bad_frac[re_ord,1:3])), ylim=c(0,1), col=NA, border=col16s[plate_col], axisnames=FALSE)

barplot(t(bad_frac[re_ord,1:3]), ylim=c(0,1), border=NA, col=c("black","gray50", "gray80"), add=TRUE, axisnames=FALSE)
par(lwd=2)
barplot(height=(bad_frac[,1]+sapply(bad_samples500sb, length)/95+sapply(poorSNR_samples_noempty500sb, length)/95)[re_ord], col="#ffffff80", add=TRUE, 
        density=rep(15,7), angle=45, names="", border="#ffffff80")

mtext("Fraction of compromised injections", side=2, padj=-4, font=1, cex=1)
legend("topleft", legend=c("Failed injection", "Empty windows", "Low SNR windows", "Irrecoverable injections"), 
       col= c("black", "gray50", "gray70", "black"), pt.cex=1.5, 
       cex=1, pch=15, box.lty=0, y.intersp=1.2, bg=NULL)
#irrecoverable injections legend box needs to be changed by hand


##### Fig 4B #####
#all of the features in every well, nothing subtracted
dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 4, 2, 1), xpd=TRUE, lwd=1)
beeswarm(inj_ba_before[re_ord], method="swarm",corral="wrap", corralWidth = 0.7, 
         pch=21, col=col16s[plate_col], bg=paste0(col16s[plate_col],"70"),
         spacing=0.7, cex=0.7, side=0, labels="", ylim=c(0,2000))
bxplot(inj_ba_before[re_ord], probs = 0.5, lwd=4, add = TRUE, col="gray70")
bxplot(inj_ba_before[re_ord], probs = 0.5, lwd=2, add = TRUE, col="black")
mtext("Raw feature count", side=2, padj=-4, font=1, cex=1)
axis(side=1, at=1:length(plate), labels=rep("",length(plate)))



##### Fig 4C #####
####fraction of features in each well coming from media
dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 4, 2, 1), xpd=TRUE, lwd=1)
beeswarm(inj_ba[re_ord], method="swarm",corral="wrap", corralWidth = 0.7, 
         pch=21, col=col16s[plate_col], bg=paste0(col16s[plate_col],"70"),
         spacing=0.7, cex=0.7, side=0, labels="", ylim=c(0,1))
bxplot(inj_ba[re_ord], probs = 0.5, lwd=4, add = TRUE, col="gray70")
bxplot(inj_ba[re_ord], probs = 0.5, lwd=2, add = TRUE, col="black")
mtext("Fraction of features from media", side=2, padj=-4, font=1, cex=1)
axis(side=1, at=1:length(plate), labels=rep("",length(plate)))


##### Fig 4D #####
#corrected for multiple adducts per unique compound
dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 4, 2, 1), xpd=TRUE, lwd=1)
beeswarm(lapply(pAd_corrected[re_ord], unlist), method="swarm",corral="wrap", corralWidth = 0.7,
         pch=21, col=col16s[plate_col], bg=paste0(col16s[plate_col],"70"),
         spacing=0.7, cex=0.7, side=0, labels="",main="", ylim=c(0,1000))
bxplot(lapply(pAd_corrected[re_ord], unlist), add=TRUE, probs=0.5, lwd=5, col="gray90")
bxplot(lapply(pAd_corrected[re_ord], unlist), add=TRUE, probs=0.5, lwd=3)
mtext("Number of unique features per well", side=2, padj=-4, font=1, cex=1)
axis(side=1, at=1:length(plate), labels=rep("",length(plate)))


##### Fig 4E #####
#barplot of mean features in common feature vectors, that fall within a range of masses
grays <- c(paste0("gray", c(10,40,60,90)))

dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mar = c(2, 4, 2, 1), xpd=FALSE, lwd=8)
barplot(colSums(samp_feats_means[,re_ord]), col=NA, border=col16s[plate_col], ylim=c(0,1200), xaxt="n")
mtext("Mean m/z features in ranges", side=2, padj=-4, font=1, cex=1)
par(lwd=1)
barplot(samp_feats_means[,re_ord],  col=grays, border=NA,  add=TRUE, xaxt="n")
legend("topright",legend=c("< 350","350-700","700-1050","> 1050" ), pch=15, pt.cex=1.5, 
       box.lty = 0, cex=1, y.intersp=1.1, col=grays)


#### Fig 5A ####
#barplot of 7ppm  unique features by medium
plotVar <- function(var) {
  xloc <- barplot(stats[, "Median", var]*100,
                  col=col16s[plate_col],
                  ylab="",
                  #ylim=range(0, stats[,, var]*100),
                  ylim=c(0,5),
                  names="")
  arrows(xloc,
         stats[, "Lower", var]*100,
         xloc,
         stats[, "Upper", var]*100,
         code=3,
         angle=90,
         length=0.03)
  
  stats[, "Median", var]*100
}


dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 5, 2, 1), xpd=FALSE, lwd=1)
plotVar("TPRs") # Figure 5A
abline(h=0)

mtext(expression(
  atop(paste("Percent of ", italic("Streptomyces"), "-associated"), 
       paste("NPA compounds detected"))), 
  side=2, padj=-1, font=1, cex=1) 





##### Fig 5B #####
dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 5, 2, 1), xpd=TRUE, lwd=1)
beeswarm(npa_feats5[re_ord[1:7]], method="swarm",corral="wrap", corralWidth = 0.7,
         pch=21, col=col16s[plate_col], bg=paste0(col16s[plate_col],"70"), ylim=c(0,30),
         spacing=0.7, cex=0.7, side=0, labels="",main="")
bxplot(npa_feats5[re_ord], add=TRUE, probs=0.5, lwd=5, col="gray90")
bxplot(npa_feats5[re_ord], add=TRUE, probs=0.5, lwd=3)
#mtext("NPAtlas features detected per well", side=2, padj=-4, font=1, cex=1)
axis(side=1, at=1:length(plate), labels=rep("",(length(plate))))
mtext(expression(
  atop(paste(italic("Streptomyces"), "-associated NPA"), 
       paste("compounds detected per well"))), 
  side=2, padj=-1, font=1, cex=1) 


#### Fig 5C ####
dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 5, 2, 1), xpd=TRUE, lwd=1)
beeswarm(npa_feats5_perwell[re_ord], method="swarm",corral="wrap", corralWidth = 0.7,
         pch=21, col=col16s[plate_col], bg=paste0(col16s[plate_col],"70"), ylim=c(0,0.10),
         spacing=0.7, cex=0.7, side=0, labels="",main="")
bxplot(npa_feats5_perwell[re_ord], add=TRUE, probs=0.5, lwd=5, col="gray90")
bxplot(npa_feats5_perwell[re_ord], add=TRUE, probs=0.5, lwd=3)
#mtext("Fraction of NPAtlas features per well", side=2, padj=-4, font=1, cex=1)
axis(side=1, at=1:(length(plate)), labels=rep("",length(plate)))
mtext(expression(
  atop(paste("Fraction of ", italic("Streptomyces"), "-associated"), 
       paste("NPA compounds detected per well"))), 
  side=2, padj=-1, font=1, cex=1) 


#### Fig 5D ####


plotVar2 <- function(var) {
  xloc <- barplot(stats[, "Median", var],
                  col=col16s[plate_col],
                  ylab="",
                  #ylim=range(0, stats[,, var]),
                  ylim=c(-0.6,1),
                  names="")
  arrows(xloc,
         stats[, "Lower", var],
         xloc,
         stats[, "Upper", var],
         code=3,
         angle=90,
         length=0.03)
  
  stats[, "Median", var]
}


dev.new(width = 5, height = 4,  noRStudioGD = TRUE)
par(mfrow=c(1,1), mar = c(2, 5, 2, 1), xpd=FALSE, lwd=1)
plotVar2("Log-odds") # Figure 5D
abline(h=0)


mtext(expression(
  atop(paste("Selectivity for ", italic("Streptomyces"), "-produced"), 
       paste("compounds vs. other natural products (log-odds)"))), 
  side=2, padj=-1, font=1, cex=1) 


et <- Sys.time() 
et-st
