#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//LL & TR 15/4/2019
//Input is the ROI mask already generated from Suite2P and the registered movie.
//Does interpolation to make movie square.  Corrects background and calculates DF/F0.  
//Leaves behind pop wave with subscript "_DFF0" and the squared up movie ("_sqr") and corresponding ROI mask ("sqrROI").
//Approach is to measure local bkg and implement Svoboda factor.  This prevents -ve signals and is sensible.
//The local bkg is calculated form mask that approximates all ROIs as circles.  
//Measures average radius of ROIs and excludes pixels that are within 1.5*radius (i.e those that are too close) 
//but includes those that are > 1.5 r and < 2.25 r, as log as they are not on a neighbouring ROI exclusion zone.

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Function DFF0(ROIMask, registeredMovie)

	wave ROIMask, registeredMovie

	wavestats/q ROIMask
	variable nROIs = (V_min*(-1))		//nROIs is the number of ROIs - edit - TR - 15/04/19 ROI index begins at 0 in Suite2P: +1 if using suite2P
	variable nx = dimsize(ROIMask, 0)
	variable ny = dimsize(ROIMask,1)
	variable nFrames = dimsize(registeredMovie,2)
	variable deltat = dimdelta(registeredMovie, 2)
	variable exclusionradius						//Distance from center of each ROI to be excluded. Currently in pixels
	variable maxradius								//Max distance from center of each ROI to be included. Currently in pixels
	variable SvobodaFactor = 0.5					//Otherwise, bkg seems to be overestimated. 
	variable scalefactor

	Make/O/N=(nFrames, nROIs) ROI_pop=0, BKG_pop=0
	Make/O/N=(nFrames) ROI
	SetScale /P x,0,deltat, ROI_pop, BKG_pop, ROI
	Make/O/N=(nROIs) ROInPixels, SurroundnPixels
	Make/O/N=(nROIs, 2)  ROIcentres=0

	String inputwavename = NameofWave(registeredMovie)
	String DFFOwavename = inputwavename + "_DFF0"					//This will be the output wave
	String Squarewavename = inputwavename + "_sqr"	
	String ROIMaskwavename = inputwavename + "_sqrROI"

	variable i, j, k, l, n, ROInum
	variable m = 0									//Pixel value of ROI number x will be -x
	variable F0 									//This will be the fluorescence at rest
	variable, dx, dy


	//Rescale movie to be square
	//Make a temp movie which is downscaled appropriately
	Make/O/N=(nx, ny) tempFrameMovie
	if(nx>ny)
		Make/O/N=(ny, ny, nFrames) tempMovie
		Make/O/N=(ny, ny) tempMask, tempFrame, destFrame
		scalefactor=ny/nx
	else
		Make/O/N=(nx, nx, nFrames) tempMovie
		Make/O/N=(nx, nx) tempMask, tempFrame, destFrame
		scalefactor=nx/ny
	endif

	for(i=0; i<nFrames; i+=1)
		tempFrameMovie[][]=registeredMovie[p][q][i]
		if(nx>ny)
			ImageInterpolate/F={scalefactor, 1}/DEST=destFrame Bilinear tempFrameMovie
		else
			ImageInterpolate/F={1, scalefactor}/DEST=destFrame Bilinear tempFrameMovie
		endif
		tempMovie[][][i]= destFrame[p][q]
	endfor

	Duplicate/O tempMovie regMovie

	//Rescale ROI mask to be square
	if(nx>ny)
		ImageInterpolate/F={scalefactor, 1}/DEST=destFrame Bilinear ROIMask
	else
		ImageInterpolate/F={1, scalefactor}/DEST=destFrame Bilinear ROIMask
	endif
	tempMask[][]= destFrame[p][q]
	
	Duplicate/O tempMask ROIwave, surround_ROI
	surround_ROI = 1
	KillWaves/Z/F tempMask

	// insert TR - 15/04/19 - reset nx and ny variables to be the dimensions of the squared waves

	nx = dimsize(ROIWave, 0)
	ny = dimsize(ROIWave,1)

	//Find the center of mass of each ROI
	ROInPixels=0
	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			ROInum = ROIwave[i][j] *-1
			if(ROInum >= 0)
				ROInPixels[ROInum]+=1	
				ROIcentres[ROInum][0]+=i
				ROIcentres[ROInum][1]+=j
			endif
		endfor
	endfor

	for (i=0;i<nROIs;i+=1)	
		ROIcentres[i][]/=ROInPixels[i]	
	endfor

	//Calculate the concentric circles over each ROI for the exclusion mask
	variable aveROIradius = sqrt(mean(ROInPixels)/Pi)			//mean radius in pixels
	exclusionradius = aveROIradius*1.25
	maxradius	 = exclusionradius*1.5

	//Make the exclusion mask
	Duplicate/O ROIwave ExclusionMask
	ExclusionMask = 1

	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			for (k=0;k<nROIs;k+=1)
				dx = abs(i-ROIcentres[k][0])
				dy = abs(j-ROIcentres[k][1])
				if( (sqrt(dx^2 + dy^2) < exclusionradius))
					ExclusionMask[i][j]=0
				endif
			endfor
		endfor
	endfor

	
	// edit TR - 15/04/19 - loop frames
	//Make a version of the movie in which exclusion zones are blanked to zero
	Duplicate/O regMovie ExcluderegMovie

	for(i=0;i<(dimsize(ExcluderegMovie,2));i+=1)
		duplicate/o/r=[][][i] ExcluderegMovie tempExc
		tempExc*=ExclusionMask
		ExcluderegMovie[][][i]=tempExc[p][q]
	endfor

	//Now calculate signal in each ROI
	ROInPixels[]=0
	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			ROInum = ROIwave[i][j]
			if(ROInum <= 0)
				ROInPixels[ROInum*-1]+=1
				MatrixOP/O ROI = beam(regMovie,i,j)
				ROI_pop[][(ROInum*-1)] += ROI[p]    
			endif
		endfor
	endfor

	for(i=0; i<nROIs; i+=1)
		ROI_pop[][i]= ROI_pop[p][i]/ROInPixels[i]
	endfor

	//Measure the surround
	SurroundnPixels[]=0
	for (i=0; i< nx; i+=1)
		for (j=0; j< ny; j+=1)
			if(ExclusionMask[i][j]==1)
				for (k=0;k<nROIs;k+=1)	
					dx = abs(i-ROIcentres[k][0])
					dy = abs(j-ROIcentres[k][1])
					if((sqrt(dx^2 + dy^2) < maxradius))
						SurroundnPixels[k]+=1
						MatrixOP/O ROI = beam(regMovie,i,j)   
						BKG_pop[][k] += ROI[p] 						 
					endif
				endfor
			endif
		endfor
	endfor

	BKG_pop*=SvobodaFactor
	for(i=0; i<nROIs; i+=1)
		BKG_pop[][i]= BKG_pop[p][i]/SurroundnPixels[i]
	endfor	

	//Now background correct (i.e subtract background)
	Duplicate/O ROI_pop BKGCorr_pop DFF0_pop
	BKGCorr_pop = ROI_pop - BKG_pop

	//Now find F0 and calculate DF/F0
	// edit TR - 15/04/19 - now computes f0 as the mode rather than the average
	for(i=0; i<nROIs; i+=1)
		duplicate/free/o/r=[][i] BKGCorr_pop ROI
		
		Make/o/n=100 Histo
		Histogram /b=1 ROI, Histo	// Histogram it
		wavestats /q/m=1 histo				// compute the wavestats... includes V_maxLoc -> the x location of the largest bin
		//if(V_MaxLoc>0)
		F0 = V_MaxLoc
		//	print F0
		//else
		//	F0 = 1
		//endif
		DFF0_pop[][i]=(ROI[p] -F0)/ F0
	
	endfor
	
		// If you are imaging synapses, some of the ROIs are lost in the downsampling step... this happens when the ROI is a single pixel in any dimension.
		// These traces come out as a series of NaNs. Since these are spurious ROIs as it is, we shall simply zero them here. Since correcting for axial 
		// motion is neccessary with synapses, these will be excluded later in the correctqa() procedure
	
	
	
	
	
	Duplicate/O DFF0_pop $DFFOwavename
	Duplicate/O regMovie $squarewavename
	Duplicate/O ROIwave $ROIMaskwavename

	KillWaves/Z/F SurroundnPixels, ROI, ROInPixels, ROIWave, regMovie, BKGCorr_pop,ROI_pop, BKG_pop
	KillWaves/Z/F ROInPixels, ROIcentres, ROIcentres, surround_ROI, DFF0_pop, ROIwave

End
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

