###############################################################################
# calcRatio                                                                   #
###############################################################################
calcRatio <- function( # ARGUMENTS
                       object,
                       file             = character(),
                       featurelist,
                       pairs,
                       plot             = TRUE,
                       plotSpectrum     = FALSE,
                       dir              = character(),
                       filter           = FALSE,
                       filter_tol       = 0.15,
                       mz_width         = 0.05,
                       target_intensity = 2/3,
                       rt_tolerance     = 100000,
                       mode             = 1
                      )
{
	# object      : xcmsSet
	# (after going through steps, group(), retcor(), group(), and fillPeaks())
	# file        : featurelist file [.csv format]
	# featurelist : features in pairs form
	# note        : if there is no featurelist, then calculate features based on only pairs
	# mode        1) only MixtoMix
	#             2) MixtoMix and Light Labeled sample
	#             3) MixtoMix and Heavy Labeled sample
	#             4) MixtoMix, Light and Heavy Labeled sample
	#             5) MixtoMix, Light and Heavy Labeled sampleA, Light and Heavy Labeled sampleA
	
	
	# requires xcms
	require(xcms) || stop("Couldn't load xcms")
	
	
	# get class information
	classinfo <- classInfo(object = object, mode = mode)
	
		
	# wrong input
	if(missing(featurelist) && missing(file) && missing(pairs))
	{ cat('Please check your input\n') }
	
	if( missing(mode) )
	{ cat('mode is missing?\n') }
	
		
	# mode 1
	else if(missing(featurelist) && missing(file))
	{
		cat('Mode 1: MixtoMix\nGenerating featurelist...\n')
		
		# calculate ratio
		calculateRatio(
		                object           = object,
		                features         = pairs,
		                plot             = plot,
		                plotSpectrum     = plotSpectrum,
		                dir              = dir,
		                filter           = filter,
		                filter_tol       = filter_tol,
		                mz_width         = mz_width,
		                target_intensity = target_intensity,
		                classinfo        = classinfo,
		                mode = mode
		               )
	}
	
	
	# mode 2
	else if(missing(pairs) && missing(file))
	{
		cat('Mode 1: MixtoMix and Light Labeled sample\nCalculating Ratio...\n')
				
		
		# calculate ratio
		calculateRatio(object=object, features=featurelist, plot=plot, plotSpectrum=plotSpectrum, dir=dir,
		               filter=filter, filter_tol=filter_tol, mz_width=mz_width, target_intensity=target_intensity)
	}
	else if(missing(file))
	{
		cat('Calculating ratio based on featurelist and pairs\n')
		
		
		# match featurelist and pairs
		matched <- union_features(featurelist=featurelist, XtoYfeatures=pairs, tol_mz_ppm=50, tol_rt=50)
		
		
		# calculate ratio
		calculateRatio(object=object, features=matched$features, plot=plot, plotSpectrum=plotSpectrum, dir=dir,
		               filter=F, filter_tol=filter_tol, mz_width=mz_width, target_intensity=target_intensity)
	}
	else if(missing(pairs))
	{
		cat('Calculating ratio based on featurelist using', file, '\n')
		
		
		# match featurelist and pairs
		matched <- union_features(featurelistFile=file, XtoYfeatures=pairs, tol_mz_ppm=50, tol_rt=50)
		
		
		# calculate ratio
		calculateRatio(object=object, features=matched$features, plot=plot, plotSpectrum=plotSpectrum, dir=dir,
		               filter=F, filter_tol=filter_tol, mz_width=mz_width, target_intensity=target_intensity)
		
	} # end_if
} # end_function



###############################################################################
# classInfo                                                                   #
###############################################################################
classInfo <- function( # ARGUMENTS
                       object,
                       mode = numeric()
                      )
{
	# object : xcmsSet() - should be grouped
	
	
	# input ok?
	if( missing(object) )
	{ stop('Input(object) missing at classInfo') }
	if( missing(mode) )
	{ stop('Input(mode) missing at classInfo') }
  
  
	# initialize
	re         <- NULL
	samples    <- sampnames(object)
	n          <- length(samples)
	classlabel <- sampclass(object)
	classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]
	mode1      <- ( mode == 1 && length(levels(sampclass(object))) == 1 )
	mode2      <- ( mode == 2 && length(levels(sampclass(object))) == 2 )
	mode3      <- ( mode == 3 && length(levels(sampclass(object))) == 2 )
	mode4      <- ( mode == 4 && length(levels(sampclass(object))) == 3 )
	mode5      <- ( mode == 5 && length(levels(sampclass(object))) == 5 )
	
	
	if( mode1 || mode2 || mode3 || mode4 || mode5)
	{
		class1             <- levels(sampclass(object))[1]
		class1at           <- which(classlabel %in% class1)
		numSampInClass1    <- length(class1at)
		c1end              <- class1at[length(class1at)]
		re$class1          <- class1
		re$class1at        <- class1at
		re$numSampInClass1 <- numSampInClass1
		re$c1end           <- c1end
		
	
		if( mode2 || mode3 || mode4 || mode5 )
		{
			class2             <- levels(sampclass(object))[2]
			class2at           <- which(classlabel %in% class2)
			numSampInClass2    <- length(class2at)
			c2end              <- class2at[length(class2at)]
			re$class2          <- class2
			re$class2at        <- class2at
			re$numSampInClass2 <- numSampInClass2
			re$c2end           <- c2end
		}
	
		
		if( mode4 || mode5 )
		{
			class3             <- levels(sampclass(object))[3]
			class3at           <- which(classlabel %in% class3)
			numSampInClass3    <- length(class3at)
			c3end              <- class3at[length(class3at)]
			re$class3          <- class3
			re$class3at        <- class3at
			re$numSampInClass3 <- numSampInClass3
			re$c3end           <- c3end
		}
		
		
		if( mode5 )
		{
			class4             <- levels(sampclass(object))[4]
			class4at           <- which(classlabel %in% class4)
			numSampInClass4    <- length(class4at)
			c4end              <- class4at[length(class4at)]
			re$class4          <- class4
			re$class4at        <- class4at
			re$numSampInClass4 <- numSampInClass4
			re$c4end           <- c4end
			
			class5             <- levels(sampclass(object))[5]
			class5at           <- which(classlabel %in% class5)
			numSampInClass5    <- length(class5at)
			c5end              <- class5at[length(class5at)]
			re$class5          <- class5
			re$class5at        <- class5at
			re$numSampInClass5 <- numSampInClass5
			re$c5end           <- c5end
		}
		
		
		re$classlabel      <- classlabel
	}
	else
	{
		stop('\nCheck with your mode and classes
	        mode 1 - needs only one class   (MixtoMix)
	        mode 2 - needs only two class   (MixtoMix and one sample)
	        mode 3 - needs only two class   (MixtoMix and one sample)
	        mode 4 - needs only three class (MixtoMix and two samples)
	        mode 5 - needs only five  class (MixtoMix and five samples)\n')
	}
		
	invisible(re)
} # end_function



###############################################################################
#                                                                             #
###############################################################################
calculateRatio <- function( # ARGUMENTS
                            object,
                            features,
                            plot             = TRUE,
                            plotSpectrum     = boolean(),
                            dir              = character(),
                            filter           = boolean(),
                            filter_tol       = 0.15,
                            mz_width         = 0.05,
                            target_intensity = 2/3,
                            rt_tolerance     = 100000,
                            classinfo,
                            mode
                           )
{
	# object      : xcmsSet
	# (after going through steps, group(), retcor(), group(), and fillPeaks())

	
	# INPUT OK?
	if( missing(object) )
	{ stop('Input(object) missing at calculateRatio') }
	if( missing(features) )
	{ stop('Input(features) missing at calculateRatio') }
	if( missing(classinfo) )
	{ stop('Input(classinfo) missing at calculateRatio') }
	if( missing(mode) )
	{ stop('Input(mode) missing at calculateRatio') }


	groupmat <- groups(object)
	if ( length(groupmat) == 0 )
	{ stop("No group information found") }
	

	# get mzrange AND rtrange using features
	mzrange          <- cbind(features$mz-mz_width, features$mz+mz_width)
	rtrange          <- cbind(features$rt-100,  features$rt+100)
	rtrange_accurate <- cbind(features$rtmin-5,  features$rtmax+5)
	

	# get Extracted Ion Chromatograms
	EICs1 <- getEIC(object, mzrange, rtrange,          rt="corrected")
	EICs2 <- getEIC(object, mzrange, rtrange_accurate, rt="corrected")
	

	# get object information
	groupidx   <- groupnames(EICs2)
	sampleidx  <- sampnames (EICs2)

	
	if (is.numeric(groupidx))
	{ groupidx <- EICs2@groupnames[groupidx] }
	
	grpidx <- match(groupidx, EICs2@groupnames)

	if (is.numeric(sampleidx))
	{ sampleidx <- names(EICs2@eic)[sampleidx] }

	sampidx <- match(sampleidx, names(EICs2@eic))
	numSamples  <- length(sampleidx)
	numFeatures <- dim(EICs2@rtrange)[1]


	# get new retention time range using target --------------------------------#
	# Description : 1. Take the max intensity
	#               2. Get target_intensity of the nearest distance from center
	#                  (default = 2/3) 
	#               3. Duplicate
	count        <- 1
	max_samp_val <- 0
	max_samp     <- 1
	rtrange_c    <- NULL
	new_rtrange  <- NULL
	min_val      <- Inf
	while( count < numFeatures )
	{
		for( samp in 1:numSamples )
		{
			maxoEIC1 <- data.frame(EICs2@eic[[samp]][count])
			maxoEIC2 <- data.frame(EICs2@eic[[samp]][count+1])
			
			maxo1 <- max(maxoEIC1$intensity)
			maxo2 <- max(maxoEIC2$intensity)
			if( maxo1 >= maxo2 )
			{ max <- maxo1 }
			else
			{ max <- maxo2 }
			if( max_samp < max )
			{
				max_samp_val <- max
				max_samp <- samp
			}
		}
		
		
		pairedEIC1    <- data.frame(EICs1@eic[[max_samp]][count])
		pairedEIC2    <- data.frame(EICs1@eic[[max_samp]][count+1])
		maxoEIC1      <- data.frame(EICs2@eic[[max_samp]][count])
		maxoEIC2      <- data.frame(EICs2@eic[[max_samp]][count+1])
		
		
		maxo1         <- max(maxoEIC1$intensity)
		maxo2         <- max(maxoEIC2$intensity)
		maxo1at       <- which.min(abs(pairedEIC1$intensity-maxo1))
		maxo2at       <- which.min(abs(pairedEIC2$intensity-maxo2))
		maxo1scantime <- pairedEIC1$rt[maxo1at]
		maxo2scantime <- pairedEIC2$rt[maxo2at]
		ratio         <- NULL
		stdev         <- NULL
		target_caught_at <- 0; t1 <- 0; tz1 <- 0; tz2 <- 0;
		delta_t <- 0; t_left <- 0; t_right <- 0;
		
		
		if( abs(maxo1scantime-maxo2scantime) <= rt_tolerance )
		{
			if( maxo1 > maxo2 )
			{
				intensity_half <- maxo1 * target_intensity
				t1 <- maxo1scantime
				for(c in maxo1at:1)
				{
					if( (pairedEIC1$intensity[c]-intensity_half) < 0 )
					{
						target_caught_at <- c
						break()
					}
				}
				if( target_caught_at <= 0 )
				{ tz1 <- 0 }
				else
				{ tz1 <- pairedEIC1$rt[target_caught_at] }
				for( c in maxo1at:length(pairedEIC1$intensity) )
				{
					if( (pairedEIC1$intensity[c]-intensity_half) < 0 )
					{
						target_caught_at <- c
						break()
					}
				}
				if( target_caught_at == 0 )
				{ tz2 <- t1 - 5 }
				else
				{ tz2 <- pairedEIC1$rt[target_caught_at] }
				if(tz2 >= tz1)
				{
					delta_t <- abs(t1 - tz1)
				}
				else
				{
					delta_t <- abs(t1 - tz2)
				}
			}
			else
			{
				intensity_half <- maxo2 * target_intensity
				t1 <- maxo2scantime
				for( c in maxo2at:1 )
				{
					if( (pairedEIC2$intensity[c]-intensity_half) < 0 )
					{
						target_caught_at <- c
						break()
					}
				}
				if( target_caught_at <= 0 )
				{ tz1 <- 0 }
				else
				{ tz1 <- pairedEIC2$rt[target_caught_at] }
				for(c in maxo2at:length(pairedEIC2$intensity))
				{
					if( (pairedEIC2$intensity[c]-intensity_half) < 0 )
					{
						target_caught_at <- c
						break()
					}
				}
				if( target_caught_at <= 0 )
				{ tz2 <- t1 - 5 }
				else
				{ tz2 <- pairedEIC2$rt[target_caught_at] }
				if(tz2 >= tz1)
				{
					delta_t <- abs(t1 - tz1)
				}
				else
				{
					delta_t <- abs(t1 - tz2)
				}
			} # end_if( maxo1 > maxo2 )
		} # end_if( abs(maxo1scantime-maxo2scantime) <= rt_tolerance )
		
		
		# new time range
		t_left      <- t1 - delta_t
		t_right     <- t1 + delta_t
		rtrange_c   <- cbind(t_left, t_right)
		new_rtrange <- rbind(new_rtrange, rtrange_c) # for lightEIC
		new_rtrange <- rbind(new_rtrange, rtrange_c) # for heavyEIC
		count       <- count + 2
	} # end_while
	
	# get new retention time ended ---------------------------------------------#
	
	
	# update features with new retention time range
	features$rtmin <- new_rtrange[,1]
	features$rtmax <- new_rtrange[,2]
	
	
	# with new retention time range
	#EICs1       <- getEIC(object, mzrange, rtrange,     rt="corrected")
	EICs3       <- getEIC(object, mzrange, new_rtrange, rt="corrected")
	count       <- 1
	ratio_list  <- NULL
	while( count < numFeatures )
	{
		intensities <- NULL
		rts         <- NULL
		ratio_c     <- NULL
		for( samp in 1:numSamples )
		{
			pairedEIC1 <- data.frame(EICs1@eic[[samp]][count])
			pairedEIC2 <- data.frame(EICs1@eic[[samp]][count+1])
			ratioEIC1  <- data.frame(EICs3@eic[[samp]][count])
			ratioEIC2  <- data.frame(EICs3@eic[[samp]][count+1])
			pwid       <- (t_right-t_left)/length(ratioEIC1$intensity)
			
			# calculate area
			into1  <- pwid*sum(ratioEIC1$intensity)
			into2  <- pwid*sum(ratioEIC2$intensity)
			
			# calculate ratio
			ratios <- ratioEIC1$intensity / ratioEIC2$intensity
			ratio  <- mean(ratios, na.rm=T)
			ratio_c <- cbind(ratio_c, ratio)
		} # end_for

		
		if( count == 1 )
		{
			ratio_list <- ratio_c
			ratio_list <- rbind(ratio_list, 1/ratio_c)
		}
		else
		{
			ratio_list <- rbind(ratio_list, ratio_c)
			ratio_list <- rbind(ratio_list, 1/ratio_c)
		}
		
		count <- count + 2			
	} # end_while


	ratio_list <- getRatio(
	                        mode      = mode,
	                        classinfo = classinfo,
	                        ratiolist = ratio_list,
	                        stdevlist = ratio_list
	                       )
	if( mode == 4 )
	{ stdevlist <- ratio_list$stdevlist }
	
	ratio_list <- ratio_list$ratiolist


	# final list
	featurelist <- features
	features    <- cbind(features, ratio_list)

	
	if( filter )
	{
		features <- features[features$feature_ratio < 1+filter_tol,]
		features <- features[features$feature_ratio > 1-filter_tol,]
	}

	if( plot )
	{
		final_ratio <- plotEICs(
		          object,
		          features,
		          dir          = dir,
		          mz_width     = mz_width,
		          classinfo    = classinfo,
		          plotSpectrum = plotSpectrum,
		          stdevlist    = stdevlist,
		          mode         = mode
		         )
	} # end_if( plot )


	# modify column names
	if( mode == 1 )
	{
		columnnames          <- c(sampnames(EICs3), classinfo$class1)
		colnames(ratio_list) <- paste(columnnames,  '.ratio', sep='')
	}
	else if( mode == 2 || mode == 3 )
	{
		columnnames          <- c(sampnames(EICs3), classinfo$class1, classinfo$class2)
		colnames(ratio_list) <- paste(columnnames,  '.ratio', sep='')
	}
	else if( mode == 4 )
	{
		columnnames          <- c(sampnames(EICs3), classinfo$class1, classinfo$class2, classinfo$class3)
		colnames(ratio_list) <- paste(columnnames,  '.ratio', sep='')
		colnames(stdevlist)  <- paste(columnnames,  '.stdev', sep='')
	}
	else if( mode == 5 )
	{
		columnnames          <- c(sampnames(EICs3), classinfo$class1, classinfo$class2, classinfo$class3, classinfo$class4, classinfo$class5, 'SampL/ML', 'SampL/ML')
		colnames(ratio_list) <- paste(columnnames,  '.ratio', sep='')
	}
	else # if( !(mode == 1 or 2 or 3 or 4 or 5) )
	{
		stop('\nCheck with your mode and classes
	        mode 1 - needs only one class   (MixtoMix)
	        mode 2 - needs only two class   (MixtoMix and one sample)
	        mode 3 - needs only two class   (MixtoMix and one sample)
	        mode 4 - needs only three class (MixtoMix and two samples)
	        mode 5 - needs only five  class (MixtoMix and five samples)\n')
	}


	# final list
	featurelist <- cbind(featurelist, ratio_list)
	
	if( mode == 4 )
	{
		featurelist <- cbind(featurelist, stdevlist[dim(stdevlist)[2]])
		list( ratio = ratio_list, features = featurelist, stdev = stdevlist, final = final_ratio$final )
	}
	else
	{
		list( ratio = ratio_list, features = featurelist )
	}
}



###############################################################################
#                                                                             #
###############################################################################
getRatio <- function( # ARGUMENTS
                      mode,
                      classinfo,
                      ratiolist,
                      stdevlist
                     )
{
	# INPUT OK?
	if( missing(mode) )
	{ stop('Input(mode) missing at getRatio') }
	if( missing(classinfo) )
	{ stop('Input(classinfo) missing at getRatio') }
	if( missing(ratiolist) )
	{ stop('Input(ratiolist) missing at getRatio') }
	if( missing(stdevlist) )                   
	{ stop('Input(stdevlist) missing at getRatio') }
  
	# mode1
	if( mode == 1 )
	{
		numRatios   <- dim(ratiolist)[1]
		count       <- 1
		ML_MH_ratio <- NULL
		while( count <= numRatios )
		{ # should calculate only mixtoMix ratios (mean values)
			if( count == 1 )
			{
				ML_MH_ratio <- mean(c(ratiolist[count, classinfo$class1at]), na.rm=T)
			}
			else
			{
				ML_MH_ratio <- rbind(ML_MH_ratio, mean(c(ratiolist[count, classinfo$class1at]), na.rm=T))
			}
			count <- count + 1
		} # end_while
		ratiolist <- cbind(ratiolist, ML_MH_ratio)
	}
	
	
	# mode2, 3, 4, and 5: Light-Labeled (and Heavy-Labeled respectively)
	else if( mode == 2 || mode == 3 || mode == 4 )
	{
		numRatios       <- dim(ratiolist)[1]
		count           <- 1
		ML_MH_ratio     <- NULL
		Samp_Mix_ratio1 <- NULL
		Mix_Samp_ratio1 <- NULL
		Samp_Mix_ratio2 <- NULL
		Mix_Samp_ratio2 <- NULL
		
		ML_MH_sd        <- NULL
		Samp_Mix_sd1    <- NULL
		Mix_Samp_sd1    <- NULL
		Samp_Mix_sd2    <- NULL
		Mix_Samp_sd2    <- NULL

		finalRatio      <- NULL
		
		while( count <= numRatios )
		{
			if( count == 1 )
			{
				ML_MH_ratio     <- mean(c(ratiolist[count, classinfo$class1at]), na.rm=T)
				Samp_Mix_ratio1 <- mean(c(ratiolist[count, classinfo$class2at]), na.rm=T)
				
				if( !is.na(mean(c(stdevlist[count, classinfo$class2at]))) && !is.na(mean(c(stdevlist[count, classinfo$class1at]))) )
				{
					ML_MH_sd        <- sd(c(stdevlist[count, classinfo$class1at]), na.rm=T)
					Samp_Mix_sd1    <- sd(c(stdevlist[count, classinfo$class2at]), na.rm=T)
				}
				else
				{
					ML_MH_sd     <- NA
					Samp_Mix_sd1 <- NA
				}

				if( mode == 4 )
				{
					Mix_Samp_ratio1 <- mean(c(ratiolist[count, classinfo$class3at]), na.rm=T)
					if( !is.na(mean(c(stdevlist[count, classinfo$class3at]))) )
					{ Mix_Samp_sd1  <- sd(c(stdevlist[count, classinfo$class3at]), na.rm=T) }
					else
					{ Mix_Samp_sd1 <- NA }
				}
			}
			else
			{
				ML_MH_ratio     <- rbind(ML_MH_ratio,     mean(c(ratiolist[count, classinfo$class1at]), na.rm=T))
				Samp_Mix_ratio1 <- rbind(Samp_Mix_ratio1, mean(c(ratiolist[count, classinfo$class2at]), na.rm=T))
				
				if( !is.na(mean(c(stdevlist[count, classinfo$class1at]))) && !is.na(mean(c(stdevlist[count, classinfo$class2at]))) )
				{
					ML_MH_sd      <- rbind(ML_MH_sd,        sd(c(stdevlist[count, classinfo$class1at]), na.rm=T))
					Samp_Mix_sd1  <- rbind(Samp_Mix_sd1,    sd(c(stdevlist[count, classinfo$class2at]), na.rm=T))
				}
				else
				{
					ML_MH_sd      <- rbind(ML_MH_sd,     NA)
					Samp_Mix_sd1  <- rbind(Samp_Mix_sd1, NA)
				}
				
				if( mode == 4 )
				{
					Mix_Samp_ratio1 <- rbind(Mix_Samp_ratio1, mean(c(ratiolist[count, classinfo$class3at]), na.rm=T))
					if( !is.na(mean(c(stdevlist[count, classinfo$class3at]))) )
					{ Mix_Samp_sd1    <- rbind(Mix_Samp_sd1, sd(c(stdevlist[count, classinfo$class3at]), na.rm=T)) }
					else
					{	Mix_Samp_sd1      <- rbind(Mix_Samp_sd1, NA) }
				}
			}
			count <- count + 1
		} # end_while
		
		if( mode == 2 || mode == 3 )
		{ ratiolist <- cbind(ratiolist, ML_MH_ratio, Samp_Mix_ratio1) }
		else if( mode == 4 )
		{
			ratiolist <- cbind(ratiolist, ML_MH_ratio, Samp_Mix_ratio1, Mix_Samp_ratio1)
			stdevlist <- cbind(stdevlist, ML_MH_sd,    Samp_Mix_sd1,    Mix_Samp_sd1)
		}
		else if( mode == 5 )
		{ ratiolist <- cbind(ratiolist, ML_MH_ratio, Samp_Mix_ratio1, Mix_Samp_ratio1, Samp_Mix_ratio2, Mix_Samp_ratio2) }
	}
	
	if( mode == 4 )
	{ list( stdevlist = stdevlist, ratiolist = ratiolist ) }
	else
	{ list( ratiolist = ratiolist ) }
}



###############################################################################
# plotEICs                                                                    #
###############################################################################
plotEICs <- function( # ARGUMENTS
                      object,
                      features,
                      plotSpectrum = boolean(),
                      dir          = character(),
                      mz_width     = mz_width,
                      classinfo,
                      mode,
                      stdevlist,
                      ...
                     )
{
	# object   : xcmsSet ( retention corrected )
	# features : pairs information
	
	
	# INPUT OK?
	if( missing(object) )
	{ stop('Input(object) missing at plotEICs') }
	if( missing(features) )
	{ stop('Input(features) missing at plotEICs') }
	if( missing(plotSpectrum) )
	{ stop('Input(plotSpectrum) missing at plotEICs') }
	if( missing(dir) )
	{ stop('Input(dir) missing at plotEICs') }
	if( missing(mz_width) )
	{ stop('Input(mz_width) missing at plotEICs') }
	if( missing(classinfo) )
	{ stop('Input(classinfo) missing at plotEICs') }
	if( missing(mode) )
	{ stop('Input(mode) missing at plotEICs') }
	if( missing(stdevlist) )                   
	{ stop('Input(stdevlist) missing at plotEICs') }
	
	
	# create directory for EIC files
	if (length(dir))
	{
		eicdir <- paste(dir, "_eic", sep="")
		dir.create(eicdir)
	}
	
	
	# get mzrange AND rtrange using features
	mzrange <- cbind(features$mz-mz_width, features$mz+mz_width)
	rtrange <- cbind(features$rt-100,  features$rt+100)
	
	
	# get Extracted Ion Chromatograms
	plotEIC <- getEIC(object, mzrange=mzrange, rtrange=rtrange)
	
	
	# get object information
	groupidx    <- groupnames(plotEIC)
	sampleidx   <- sampnames (plotEIC)
	
	
	if (is.numeric(groupidx))
	{ groupidx  <- plotEIC@groupnames[groupidx] }
	grpidx  <- match(groupidx, plotEIC@groupnames)
	if (is.numeric(sampleidx))
	{ sampleidx <- names(plotEIC@eic)[sampleidx] }
	sampidx <- match(sampleidx, names(plotEIC@eic))


	numSamples  <- length(sampleidx)
	numFeatures <- dim(plotEIC@rtrange)[1]
	count       <- 1
	while(count < numFeatures)
	{
		min_val  <- Inf
		temp_min <- Inf
		for(samp in 1:numSamples)
		{
			pairedEIC1 <- data.frame(plotEIC@eic[[samp]][count])
			pairedEIC2 <- data.frame(plotEIC@eic[[samp]][count+1])
			
			temp_min <- min(dim(pairedEIC1)[1], dim(pairedEIC2)[1])
			if( temp_min < min_val )
			{ min_val <- temp_min }
			
			if( samp == 1 )
			{
				intensities <- cbind(                         pairedEIC1[1:min_val,]$intensity)
				intensities <- cbind(intensities[1:min_val,], pairedEIC2[1:min_val,]$intensity)
				rts         <- cbind(                         pairedEIC1[1:min_val,]$rt)
				rts         <- cbind(rts        [1:min_val,], pairedEIC2[1:min_val,]$rt)
			}
			else
			{
				intensities <- cbind(intensities[1:min_val,], pairedEIC1[1:min_val,]$intensity)
				intensities <- cbind(intensities[1:min_val,], pairedEIC2[1:min_val,]$intensity)
				rts         <- cbind(rts        [1:min_val,], pairedEIC1[1:min_val,]$rt)
				rts         <- cbind(rts        [1:min_val,], pairedEIC2[1:min_val,]$rt)
			}
		} # end_for
		
		
		if( mode == 1 )
		{
			plotEIC(
			         y       = intensities,
			         x       = rts,
			         eicdir  = eicdir,
			         peaknum = features$pairID[count],
			         mass    = features$mass[count],
			         mzrange = c(features$mz[count], features$mz[count+1]),
			         rtrange = c(features$rtmin[count], features$rtmax[count]),
			         ratio   = features[count, length(features)],
			         classinfo = classinfo,
			         mode    = mode
			        )
		}
		else if( mode == 2 )
		{
			plotEIC(
			         y         = intensities,
			         x         = rts,
			         eicdir    = eicdir,
			         peaknum   = features$pairID[count],
			         mass      = features$mass[count],
			         mzrange   = c(features$mz[count], features$mz[count+1]),
			         rtrange   = c(features$rtmin[count], features$rtmax[count]),
			         ratio1    = features[count,  length(features)],    # = L/H
			         ratio2    = features[count, (length(features)-1)], # = L/H
			         classinfo = classinfo,
			         mode      = mode
			        )
		}
		else if( mode == 3 )
		{
			plotEIC(
			         y         = intensities,
			         x         = rts,
			         eicdir    = eicdir,
			         peaknum   = features$pairID[count],
			         mass      = features$mass[count],
			         mzrange   = c(features$mz[count], features$mz[count+1]),
			         rtrange   = c(features$rtmin[count], features$rtmax[count]),
			         ratio1    = features[count  ,  length(features)],    # = H/L
			         ratio2    = features[count+1, (length(features)-1)], # = L/H
			         classinfo = classinfo,
			         mode      = mode
			        )
		}
		else if( mode == 4 )
		{
			SampH_div_MixL   <- features [count+1, length(features)]
			SampL_div_MixH   <- features [count,  (length(features)-1)]
			MixL_div_MixH    <- features [count,  (length(features)-2)]
			
			SampH_div_MixL_e <- stdevlist[count+1, dim(stdevlist)[2]]
			SampL_div_MixH_e <- stdevlist[count,  (dim(stdevlist)[2]-1)]
			MixL_div_MixH_e  <- stdevlist[count,  (dim(stdevlist)[2]-2)]
			
			SampL_div_MixL   <- SampL_div_MixH/MixL_div_MixH
			
			SampL_div_MixL_e <- sqrt(SampL_div_MixH_e^2 + MixL_div_MixH_e^2)
			
			Samp_div_Mix     <- mean(c(SampH_div_MixL, SampL_div_MixL), na.rm = T)
			Samp_div_Mix_e   <- SampH_div_MixL_e + SampL_div_MixL_e
			
			temp <- cbind(
			               features$pairID[count],
			               features$mass[count],
			               SampH_div_MixL,
			               SampH_div_MixL_e,
			               SampL_div_MixH,
			               SampL_div_MixH_e,
			               MixL_div_MixH,
			               MixL_div_MixH_e,
			               SampL_div_MixL,
			               SampL_div_MixL_e,
			               Samp_div_Mix,
			               Samp_div_Mix_e
			              )

			if( mode == 4 )
			{
				if( count == 1 )
				{ final_result_list <- temp }
				else
				{ final_result_list <- rbind(final_result_list, temp) }
			}
			
			plotEIC(
			         y            = intensities,
			         x            = rts,
			         eicdir       = eicdir,
			         peaknum      = features$pairID[count],
			         mass         = features$mass[count],
			         mzrange      = c(features$mz[count], features$mz[count+1]),
			         rtrange      = c(features$rtmin[count], features$rtmax[count]),
			         ratio1       = features [count+1, length(features)],     # = H/L
			         ratio2       = features [count,  (length(features)-1)],  # = L/H
			         ratio3       = features [count,  (length(features)-2)],  # = L/H
			         stdev1       = stdevlist[count+1, dim(stdevlist)[2]],    # = H/L
			         stdev2       = stdevlist[count,  (dim(stdevlist)[2]-1)], # = L/H
			         stdev3       = stdevlist[count,  (dim(stdevlist)[2]-2)], # = L/H
			         finalratio   = Samp_div_Mix,
			         finalratio_e = Samp_div_Mix_e,
			         classinfo    = classinfo,
			         mode         = mode
			        )
		}
		
		if( plotSpectrum )
		{
			files <- cdfpaths( object )
			xraw  <- xcmsRaw( files[1] )
			xraw@scantime <- object@rt$corrected[[1]]
			plotMSpec(
			           xraw,
			           mass    = features$mass[count],
			           rtrange = c(features$rtmin[count], features$rtmax[count]),
			           eicdir  = eicdir,
			           peaknum = features$pairID[count]
			          )
		}
		count <- count + 2
	} # end_while


	if( mode == 4 )
	{	
		columnnames <- c(
		                  'peaknum',
		                  'mass',
		                  'SampH/MixL',
		                  'SampH/MixL_err',
		                  'SampL/MixH',
		                  'SampL/MixH_err',
		                  'MixL/MixH',
		                  'MixL/MixH_err',
		                  'SampL/MixL',
		                  'SampL/MixL_err',
		                  'Samp/Mix',
		                  'Samp/Mix_err'
		                 )
		colnames(final_result_list) <- columnnames


		list( final = final_result_list )
	}
}



###############################################################################
#                                                                             #
###############################################################################
plotEIC <- function( # ARGUMENTS
                     y,
                     x,
                     eicdir       = character(),
                     peaknum      = numeric(),
                     mzrange,
                     rtrange,
                     class1       = classinfo$class1,
                     class2       = classinfo$class2,
                     ratio        = numeric(),
                     ratio1       = numeric(),
                     ratio2       = numeric(),
                     ratio3       = numeric(),
                     stdev1       = numeric(),
                     stdev2       = numeric(),
                     stdev3       = numeric(),
                     finalratio   = numeric(),
                     finalratio_e = numeric(),
                     mass         = numeric(),
                     classinfo,
                     mode,
                     ...
                    )
{
	# y : matrix
	# x : matrix
	
	
	if (length(eicdir))
	{
		eicpicfile <- paste(peaknum, "_eic.png", sep="")
		eicpdf     <- paste(peaknum, "_eic.pdf", sep="")
		if( capabilities("png") )
			png(file.path(eicdir, eicpicfile), width = 800, height = 600)
		else
			pdf(file.path(eicdir, eicpdf), width = 800/72,
			height = 600/72, onefile = FALSE)
	}
	
	legendtext1 <- paste('Light -', round(mzrange[1], 4))
	legendtext2 <- paste('Heavy -', round(mzrange[2], 4))
	linewidth   <- c(1, 2)
	
	if( mode == 1 )
	{
		linetype         <- NULL
		linetype[1:((classinfo$numSampInClass1)*2)] <- 1
		matplot(y=y, x=x, type='l', lty=linetype, lwd=linewidth, col=c('red', 'blue'), xlab = "Retention Time", ylab = "Intensity")
		legend ("topright", c(" Light ", " Heavy "), pch='RB', col=c('red', 'blue'),)
		legend ("topleft",  c(legendtext1, legendtext2), lty=c(1, 1), lwd=linewidth)
	}
	else if( mode == 2 || mode == 3)
	{
		n1               <- classinfo$numSampInClass1
		n2               <- classinfo$numSampInClass2
		color1           <- NULL
		color2           <- NULL
		color1[1:(2*n1)] <- c('red',   'blue') # red   and blue
		color2[1:(2*n2)] <- c('green', 'gray') # green and gray
		colorAll         <- c(color1, color2)
		c1L              <- paste(c(classinfo$class1), '- Light')
		c1H              <- paste(c(classinfo$class1), '- Heavy')
		c2L              <- paste(c(classinfo$class2), '- Light')
		c2H              <- paste(c(classinfo$class2), '- Heavy')

		matplot(y=y, x=x, type='l', lty='11', col=colorAll, xlab = "Retention Time", ylab = "Intensity")
		legend ("topright", c(c1L, c1H, c2L, c2H), pch='RBGR', col=c(2, 4, 'green', 'gray'))
	}
	else if( mode == 4 )
	{
		n1               <- classinfo$numSampInClass1
		n2               <- classinfo$numSampInClass2
		n3               <- classinfo$numSampInClass3
		totn             <- n1 + n2 + n3
		color1           <- NULL
		color2           <- NULL
		color3           <- NULL
		color1[1:(2*n1)] <- 'gray'#c('red',   'blue')
		color2[1:(2*n2)] <- 'blue'#c('green', 'gray')
		color3[1:(2*n2)] <- 'red'#c('pink',  'orange')
		linetype         <- NULL
		linetype[1:totn] <- 1
		colorAll         <- c(color1, color2, color3)
		c1L              <- paste(c(classinfo$class1), '- Light')
		c1H              <- paste(c(classinfo$class1), '- Heavy')
		c2L              <- paste(c(classinfo$class2), '- Light')
		c2H              <- paste(c(classinfo$class2), '- Heavy')
		c3L              <- paste(c(classinfo$class3), '- Light')
		c3H              <- paste(c(classinfo$class3), '- Heavy')

		matplot(y=y[,(n1*2+1):dim(y)[2]], x=x[,(n1*2+1):dim(x)[2]], type='l', lwd=linewidth, lty=linetype, col=c(color2, color3), xlab = "Retention Time", ylab = "Intensity")
		legend ("topright", c(classinfo$class2, classinfo$class3), pch='RB', col=c(2, 4))
		legend ("topleft",  c(legendtext1, legendtext2), lty=c(1, 1), lwd=linewidth)

		#matplot(y=y, x=x, type='l', lwd=linewidth, lty=linetype, col=colorAll, xlab = "Retention Time", ylab = "Intensity")
		#legend ("topright", c(c1L, c1H, c2L, c2H, c3L, c3H), pch='RBGRPO', col=c(2, 4, 'green', 'gray', 'pink', 'orange'))
		#legend ("topleft",  c(legendtext1, legendtext2), lty=c(1, 1), lwd=linewidth)
	}

		
	if( !missing(rtrange) )
	{
		abline(v=rtrange[1], col='black', lty=3)
		abline(v=rtrange[2], col='black', lty=3)
	}
	
	
	if( mode == 1 )
	{
		mtext(paste('Ratio L/H =', round(ratio, 2), ' '), cex = 1)
	}
	else if( mode == 2 || mode == 3 )
	{
		mtext(paste('Light -', round(mzrange[1], 3), '   ',
	              'Heavy -', round(mzrange[2], 3), '\n',
	              'Ratio Samp/Mix', '=', round(ratio1, 2), ' ',
	              'Ratio MixtoMix', '=', round(ratio2, 2), ' '), cex = 1)
	}
	else if( mode == 4 )
	{
		mtext(paste('SampH/MixL', '=', round(ratio1, 2),     '+/-', round(stdev1, 2), ', ',
	              'SampL/MixH', '=', round(ratio2, 2),     '+/-', round(stdev2, 2), ', ',
	              'MixL/MixH',  '=', round(ratio3, 2),     '+/-', round(stdev3, 2), '  ',
	              'Samp/Mix',   '=', round(finalratio, 2), '+/-', round(finalratio_e, 2)), cex = 1.2,)
	}
	
	              
	if( !missing(mass) )
	{ title(main=paste('Peak #', peaknum, ' - Neutral Mass [', round(mass, 3), ']'), font.main = 2, cex.main = 2) }
	else
	{ title(main=paste('Peak #', peaknum)) }
	
	if (length(eicdir)) { dev.off() }
}



###############################################################################
#                                                                             #
###############################################################################
plotMSpec <- function( # ARGUMENTS
                       object,
                       mass,
                       rtrange,
                       mzsize,
                       eicdir   = character(),
                       peaknum  = numeric()
                      )
{
	# object : xcmsRaw
	if (length(eicdir))
	{
		eicpicfile <- paste(peaknum, "_massSpec.png", sep="")
		eicpdf <- paste(peaknum, "_massSpec.pdf", sep="")
		
		if( capabilities("png") )
		{ png(file.path(eicdir, eicpicfile), width = 800, height = 600) }
		else
		{ pdf(file.path(eicdir, eicpdf), width = 800/72, height = 600/72, onefile = FALSE) }
	}
	
	mzs <- c(mass-mzsize, mass+mzsize)
	
	plotSpec(object, mzrange=mzs, , timerange=rtrange)
	
	if (length(mass) && length(peaknum))
	{ mtext(paste('Peak #', peaknum, ' - Neutral Mass [', round(mass, 4), ']\n',
                   'Light -', round(mass[1], 4),
                '   Heavy -', round(mass[2], 4)), cex = 0.8) }
   	
	if (length(eicdir))
	{ dev.off() }
}



###############################################################################
#                                                                             #
###############################################################################
setGeneric("plotSpec", function(object, ...) standardGeneric("plotSpec"))

setMethod("plotSpec", "xcmsRaw", function(object, ident = FALSE, 
                                          vline = numeric(0), real_mz, ...) {

    sel <- profRange(object, ...)

    title = paste("Averaged Mass Spectrum: ", sel$timelab, " (", 
                  sel$scanlab, ")",  sep = "")

    points <- cbind(profMz(object)[sel$massidx], 
                    rowMeans(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
    
		light_mz_at <- which.min(abs(points[,1]-real_mz[1]))
		heavy_mz_at <- which.min(abs(points[,1]-real_mz[2]))
		ylim        <- max(points[light_mz_at,2], points[heavy_mz_at,2])
		plot(points, ylim=c(0, ylim*1.1), type="h", main = title, xlab="m/z", ylab="Intensity")
    
    if (length(vline))
        abline(v = vline, col = "red")
    
    if (ident)
        return(identify(points, labels = round(points[,1], 6)))
		
		num_points <- dim(points)[1]
		if( num_points > 6 )
		{
			labels <- points[(light_mz_at-1):(light_mz_at+1),]
			labels <- rbind(labels, points[(heavy_mz_at-1):(heavy_mz_at+1),])
			text(labels[,1], labels[,2], labels[,1], adj=c(.5, -.5), col='blue3')
		}
    
    invisible(points)
})
