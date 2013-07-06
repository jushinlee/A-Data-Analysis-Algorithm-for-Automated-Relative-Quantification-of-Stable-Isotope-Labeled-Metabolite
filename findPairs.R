findPairs <- function(object, peaklist, mz_shift=2.0067, mz_tol=0.0005, rt_tol=10, light_label=28.0312, heavy_label=30.0379)
{
	# object : xcmsSet
	# peaklist from the xcmsSet (or xcmsRaw through findPeaks() function)
	# default value of mz_shift = 2.0067 (dalton)
	# default value of mz_tol   = 0.0005 (dalton)
	# default value of rt_tol   = 10     (seconds)
	# light_label               = 28.0312 (dalton)
	# heavy_label               = 30.0379 (dalton)
	
	# make data frame consistent
	if (missing(object) && missing(peaklist))
	{ cat('missing peaklist!') }
	else if(missing(object))
	{
		cnames   <- colnames(peaklist)
		peaklist <- data.frame(peaklist)
		colnames(peaklist) <- cnames
	}
	else
	{
		peaklist <- peaktable(object)
	}
	
	# get number of peaks
	numPeaks <- dim(peaklist)[1]
	templist <- NULL

	# allocate annotation column
	Labels <- NA
	peaklist <- cbind(peaklist, Labels)
	
	# initialization
	end    <- numPeaks - 1
	pairID <- 1
	temp1  <- NULL
	temp2  <- NULL
	for(i in 1:end)
	{
		start <- i+1
		for(j in start:numPeaks)
		{
			mz_shift_max <- mz_shift + mz_tol
			actual_mz_tol <- abs(abs(peaklist$mz[i] - peaklist$mz[j]) - mz_shift)
			actual_rt_tol <- abs(peaklist$rt[i] - peaklist$rt[j])

			# if distance between peaks are further than the given mz_shift, then break the loop
			if(actual_mz_tol > mz_shift_max) break

			# store paired peaks
			if((actual_mz_tol <= mz_tol) && (actual_rt_tol <= rt_tol))
			{
				mass <- peaklist[i,]$mz - light_label
				temp1 <- cbind(pairID, mass, peaklist[i,])
				temp2 <- cbind(pairID, mass, peaklist[j,])
				templist <- rbind(templist, temp1, temp2)
				annotation <- paste('[', pairID, '] [M+Light]  ([M]=', round(mass, 6), ')', sep='')
				
				# annotate here

				pairID <- pairID + 1 # increment pairID
			}
		} # end_for
	}   # end_for
	
	# remove labels
	 templist <- templist[-which(colnames(templist)=="Labels")]
	
	list(pairs=templist, peaklist=peaklist)
}

peaktable<-function(xs)
{
	if (nrow(xs@groups) > 0)
	{
		groupmat <- groups(xs)
		ts <- data.frame(cbind(groupmat, groupval(xs, "medret", "mz"), groupval(xs, "medret", "into"), groupval(xs, "medret", "maxo"), groupval(xs, "medret", "rtmin"), groupval(xs, "medret", "rtmax")), row.names = NULL)
		cnames <- colnames(ts)
		if (cnames[1] == "mzmed")
			cnames[1] <- "mz"
		else stop("Peak information ?!?")
		if (cnames[4] == "rtmed") 
			cnames[4] <- "rt"
		else stop("Peak information ?!?")
			colnames(ts) <- cnames
	}
	else if (length(xs@sampnames) == 1)
	{
		cnames   <- colnames(xs@peaks)
		ts <- data.frame(xs@peaks)
		colnames(ts) <- cnames
	}
	else stop("First argument must be a xcmsSet with group information or contain only one sample.")
	ts
}
