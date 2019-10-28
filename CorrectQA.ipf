#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "Sarfia"
#include "Advanced ROI tools"


// Thomas Ryan 27/02/19 - 

// Summary: This procedure uses existing waves representing the estimated axial displacement over time 
// (shiftsName), the zProfile of each ROI arranged in rows (zProfilesName) and the experimental 
// activity traces from the same ROIs. For synapses, each profile is normalised to the interpolated 
// fluorescence value expected at the imaging plane. These profiles are then fit with a Moffat function;
// a variation of a t-distribution which is well suited to describe how synapses recorded in-vivo are 
// convolved with the PSF of the microscope in the axial dimension. These fits are then used as a 
// look-up table, to generate a vector describing how we predict changes to fluorescence may occur due 
// to axial displacements over the course of an experiment. This vector is then used to correct the 
// fluorescent activity traces of the ROI to more closely reflect the true intensity in our intended plane. 

// Corrected traces are passed to a new 3D wave; CorrectedQA_DivMoff, with a new vector 
// CorrectedROI_Id containing the original ROI indexes for back referencing. Another wave, 
// CorrectedQA_SubMoff is similar, but computed by estimating a new baseline for the ROI and subtracting 
// this from the raw trace. In this case, the original baseline fluorescence is added back to the time 
// series following subtraction. Fitted curves and coefficients are stored in zFitCurves and zFitCoefs
// respectively.

// The first flag within the code below allows you to select between performing the above procedure using a 
// fitted Moffat function, or using a more empirical fitted spline. This may be more appropriate for ROIs
// which have arbitrary size and shape in the z-dimension. Set to: 

// -> 1 for Moffat fitting
// ->2 for spline fitting
// ->3 for both.

// In the case of empirical correction through a spline the equivalent results
// waves listed above will be generated with the suffix "Emp"

// Several exclusion criteria are included so that only synapses we are confident of correcting are passed. 
// Those that are excluded are still "corrected" and can be found in the corrected waves prefixed with
// "All_", the "CorrectedQA" waves only contain the passed and corrected waves. Original ROI index of passed 
// ROIs is stored in "ROI_Id_Orig"

function CorrectQA()

	string shiftsName=wavelist("zestimate","","")
	string zProfilesName=wavelist("CroppedRef_Ch1_reg_DFF0","","")
	string QAName=wavelist("CroppedCh1_DFF0","","")
	
	wave shifts=$shiftsName, zProfiles=$zProfilesName, QA=$QAName
	duplicate/o/free $shiftsName shifts
	duplicate/o/free $zProfilesName zProfiles
	duplicate/o/free $QAName QA
	//////// Choose Moffat (1) vs. spline (2) or both (3)\\\\\\\\\
	
	variable flag=1
	
	//////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	/////////////// Set distance between z-slices \\\\\\\\\\\\\\\\
	
	variable/g zstep=0.5
	
	/////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
		
	// Find some useful variables to be used as global variables to pass to different procedures
	variable/g numRois=dimsize(zProfiles,1)
	variable/g numSlices=dimsize(zProfiles,0)
	variable/g recLength=dimsize(shifts,0)
	
	
	variable i, j, R
	
	// Compute the range of z-positions the zProfiles cover
	variable/g startZ=(((floor(numslices/2))*zstep)*-1)+zstep
	make/o/N=(numSlices-2) xwave = startZ+(p*zstep)

	// Compute the starting position as the average of the estimated position of the first 5 time points
	variable/g startPos=mean(shifts,0,4)
	
	// Smooth the estimate of shifts (similar to lopass filter) 
	smooth 5,shifts
	
	if(flag==1) // If we are correcting using a Moffat only...

		CoQAMoff(shifts,zProfiles,QA,xwave,startZ,numSlices,numRois,reclength,startPos)

	else
		if(flag==2) // If we are correcting using Spline fits...

			CoQAEmp(shifts,zProfiles,QA,xwave,startZ,numSlices,numRois,reclength,startPos)

		else
			if(flag==3) // If we are correcting using both a Moffat and spline...

				CoQAMoff(shifts,zProfiles,QA,xwave,startZ,numSlices,numRois,reclength,startPos)
				CoQAEmp(shifts,zProfiles,QA,xwave,startZ,numSlices,numRois,reclength,startPos)
				
			else
			
				print "Please select a valid flag -> 1 for Moffat fitting, 2 for empirical spline fitting, 3 for both"
				
			endif
		endif
	endif
end
	
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Thomas Ryan 26/02/19

// Summary - Fits a spline to the profiles of each ROI and uses this as a look up table to predict what
// fluorescent changes we expect from axial displacement. Correction is performed using both divisive and 
// subtractive  methods. These are stored in waves CorrectedQA_EmpDiv & CorrectedQA_EmpSub respectively. 
// Plotted spline fits are stored in zFitCurves_Emp

// Inputs:
// shifts			-	wave		-	1D vector describing the shifts in the axial plane over time in microns
// zProfiles		-	wave		-	2D wave containing the raw axial profiles (columns) of each ROI (rows)
// QA				-	wave 		-	2D wave containing the fluorescent activity over time (columns) of each ROI (rows)
//	xwave			-	wave		-	1D vector containing the x-scaling of profiles
// startZ			-	variable	-	the axial position of the first slice of the zstack (changes depending on whether stack was taken from the bottom or top of the volume)
//	numSlices		-	variable	-	the number of slices in the reference stack
// numRois		-	variable	-	the number of ROIs
//	recLength		-	variable	-	the length of the time series
//	startPos		-	variable	-	the estimated starting axial position of the time series in microns

//	Generates:
// 

function CoQAEmp(shifts,zProfiles,QA,xwave,startZ,numSlices,numRois,reclength,startPos)

	wave shifts,zProfiles,QA,xwave
	variable numSlices, numRois, recLength,startZ,startPos
	variable i,j,startF

	// Create some objects to accept the results
	make/o /n=(dimsize(QA,0),dimsize(QA,1)) DivTrace_emp, All_CorrectedQA_DivEmp, SubTrace_emp, CorrectedQA_SubPri, All_CorrectedQA_SubEmp
	make/o/n=(numROIs) ROI_Log_Keep_Emp=1
	// keep these in the same time scaling as the original traces
	copyscaling(QA,DivTrace_emp)
	copyscaling(QA,SubTrace_emp)

	for(i=0;i<numRois;i+=1) // cycle through ROIs
			
		duplicate/free/o/r=[1,(numSlices-2)][i] zProfiles tempZ 	// grab this ROI's zProfile (not first or last slice)
		interpolate2/T=2 /Y=SS_tempZ tempZ								// interpolate the zProfile with a cubic spline
		setscale/i x,xwave[0],xwave[inf], tempZ,zProfiles,SS_tempZ	// reset scale as we now have 200 points
		
		// if we're on the first profile then make the waves to store zprofile curves and the new xwave 
		// for interpolated data
		
		if(i==0) 
			make/o/n=(dimsize(SS_tempZ,0),numRois) zFitCurves_Emp
			variable delt=dimdelta(SS_tempZ,0)
			make/o/n=(200) xwaveInterp=startZ+(p*delt)
		endif
		
		duplicate/free/o SS_tempZ tempZ_norm
		variable wmin=wavemin(tempZ_norm) 	// find the minimum in the profile and subtract
		tempZ_norm-=wmin
		variable normaliser=interp(startPos,xwaveInterp,tempZ_norm)	// find the intesity at the starting plane and divide
		tempZ_norm/=normaliser
			
		zFitCurves_Emp[][i]=tempZ_norm[p]				// store the fitted splines normalised and minimum subtracted
		
		// In case of variations in intensity between the zstack and the activiy trace,(e.g. because of 
		// filtering and averaging etc), scale the fit to the starting fluorescence.
		
		duplicate/free/o/r=[][i] QA QAtemp
		redimension/n=(-1) QAtemp
		setscale/p x,0,1,QAtemp
		startF=mean(QAtemp,0,14) // find the starting fluorescence for this ROI (baseline) based n first 1.5 sec)
		
		duplicate/o tempZ_norm, tempZ_Scaled
		tempZ_Scaled*=startF		// scale the stack to the activity trace
		setscale/i x,xwave[0],xwave[inf], tempZ_Scaled, tempZ_Norm // rescale: may be redundant 

		// Now cycle through every time frame and compute the change in fluorescence expected due to shifts
		// Done on the normalised profiles, this will give us the divisor trace
		// Done on the scaled profiles, this will give us the subtractor trace
		
		for(j=0;j<reclength;j+=1) 
		
			DivTrace_emp[j][i]=interp((shifts[j]),xwaveInterp,tempZ_norm) // build divisor trace
			SubTrace_emp[j][i]=interp((shifts[j]),xwaveInterp,tempZ_Scaled) // and subtractor trace
			
		endfor
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////Exclusion Module //////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////
		
		variable UpperThreshWidth=15
		variable LowerThreshWidth=7
		
		// find levels at a half the maximum value of the synapse (FWHM)
		
		variable HalfMax=wavemax(tempZ_norm) // vaguely similar to 2*SD
		HalfMax/=2
		
		FindLevels/q/B=5 /D=LevelRisew /edge=1 tempZ_Norm,HalfMax // find rising level at halfmax
		
		variable FifthMax=wavemax(tempZ_norm) 
		FifthMax/=5

		// If we don't find a rising level, or more than one, then exclude...

		if(V_Flag==2)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif

		if(V_levelsFound!=1)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif

		FindLevels/q/B=5 /D=LevelFallw /edge=2 tempZ_Norm,HalfMax // find falling level
		
		// If we don't find a falling level, or more than one, then exclude...

		if(V_Flag==2)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif

		if(V_levelsFound!=1)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif	
		
		variable LevelFall=LevelFallw[0]
		variable LevelRise=LevelRisew[0]
		variable LevelWidth=LevelFall-LevelRise
		killwaves LevelFallw,LevelRisew
		// if this width is wider than threshWidth (typical synapse) then exclude...

		if(LevelWidth>UpperThreshWidth)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif

		// if this width is narrower than threshWidth (typical synapse) then exclude...
		if(LevelWidth<LowerThreshWidth)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif
		
		// If shifts cause intensity to fall by more than 90% (division factor of < 0.1 ) exclude as this is close to noise
		duplicate/o/free/r=[][i] Divtrace_emp temp
		variable maxDrop=wavemin(temp)
		
		if(maxDrop<0.1)
			ROI_Log_Keep_Emp[i]=0
			continue
		endif
		
		
		// if the shifts do not fall within the levelRange, defined by width at the "fifthMax", exclude... (this does not work well with the empirical computation of the correction factor as we have no model to eliminate noise entirely)
		// This removes synapses whose intensity falls to below 20% of the maximum
		
		//variable HighShift=wavemax(shifts)
		//variable LowShift=wavemin(shifts)

		//if(LowShift<LevelRise)
		//	ROI_Log_Keep_Emp[i]=0
		//	continue
		//endif

		//if(Highshift>LevelFall)
		//	ROI_Log_Keep_Emp[i]=0
		//	continue
		//endif
		
	endfor
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	setscale/i x,xwave[0],xwave[inf], ZFitCurves_Emp // rescale: may be redundant 

	// Now use the divisor and subtractive traces to correct the raw QA wave
		
	MatrixOP/o All_CorrectedQA_DivEmp=QA/DivTrace_emp				// divide by the divisor trace
	Matrixop/o All_CorrectedQA_SubPri=QA-SubTrace_emp				// subtract by the subtractor trace
	
	// If there is no axial movement, then the difference between QA and SubTrace will = 0
	// Add baseline fluorescence for this ROI back
		 
	for(i=0;i<numRois;i+=1)
		duplicate/free/o/r=[][i] All_CorrectedQA_SubPri tempSub
		redimension/n=(-1)  tempSub
		
		duplicate/free/o/r=[][i] QA QAtemp
		redimension/n=(-1) QAtemp
		setscale/p x,0,1,QAtemp
		startF=mean(QAtemp,0,14) // find the starting fluorescence for this ROI (baseline)
		tempSub+=startF
		All_CorrectedQA_SubEmp[][i]=tempSub[p]
		
	endfor
	copyscaling(QA,All_CorrectedQA_DivEmp)					// rescale corrected_Div wave
	copyscaling(QA,All_CorrectedQA_SubEmp)					// rescale corrected_Sub wave


	////////////////////////////////// Exclusion ///////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Split the ROIs we can correct confidently

	variable numCorrected=sum(ROI_Log_Keep_Emp)
	make/o/n=(reclength,numCorrected) CorrectedQA_DivEmp, CorrectedQA_SubEmp
	make/o/n=(numCorrected) ROI_ID_Orig_Emp
	variable count=0

	for(i=0;i<numRois;i+=1)

		if(ROI_Log_Keep_Emp[i]==1)

			duplicate/o/r=[][i] All_CorrectedQA_DivEmp CorrectedQA_DivEmpTemp
			redimension/n=(-1) CorrectedQA_DivEmpTemp

			duplicate/o/r=[][i] All_CorrectedQA_SubEmp CorrectedQA_SubEmpTemp
			redimension/n=(-1) CorrectedQA_SubEmpTemp

			CorrectedQA_DivEmp[][count]=CorrectedQA_DivEmpTemp[p]
			CorrectedQA_SubEmp[][count]=CorrectedQA_SubEmpTemp[p]
			ROI_ID_Orig_Emp[count]=i
			count+=1

		endif
	endfor
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
end

// Thomas Ryan 28/01/19 - Thanks to Benjamin James for help with the curve evaluation step

// Summary - Fits a Moffat function to the zProfile of each ROI, normalised to the intensity at the starting plane
// Reads off this fit at x=shifts for each time point to create a divisor trace
// Corrects the activity traces by dividing each trace by this divisor trace 
 
// Inputs:
// shifts			-	wave		-	1D vector describing the shifts in the axial plane over time in microns
// zProfiles		-	wave		-	2D wave containing the raw axial profiles (columns) of each ROI (rows)
// QA				-	wave 		-	2D wave containing the fluorescent activity over time (columns) of each ROI (rows)
//	xwave			-	wave		-	1D vector containing the x-scaling of profiles
// startZ			-	variable	-	the axial position of the first slice of the zstack (changes depending on whether stack was taken from the bottom or top of the volume)
//	numSlices		-	variable	-	the number of slices in the reference stack
// numRois		-	variable	-	the number of ROIs
//	recLength		-	variable	-	the length of the time series
//	startPos		-	variable	-	the estimated starting axial position of the time series in microns

//	Generates:
// 


function CoQAMoff(shifts,zProfiles,QA,xwave,startZ,numSlices,numRois,reclength,startPos)

	wave shifts, zProfiles,QA,xwave
	variable startZ,numSlices,numRois,reclength,startPos

	variable i,j,R

	// Make some waves to store results
	make/o/n=(numSlices-2,numRois) NormalisedzProfiles=0
	make/o/n=(numRois) err,chi2Mat, chi2Mat_FirstFit,zPosition, ROI_Log_Keep_Moff=1
	
	string message
	for(i=0;i<numRois;i+=1) // cycle through ROIs

		// grab the profile (not the first and last slice)
		duplicate/o/r=[1,(numSlices-2)][i] zProfiles tempZ
		redimension/s/n=(-1) tempZ	// make it 1D (because it is)
		setscale/i x,xwave[0],xwave[inf], tempZ // set the x scale to match the scale of shifts
	
		// Make an xwave corresponding to the fitted curve (200pnts). This is used for interpolation of the starting fluorescence
		interpolate2/T=2 /Y=smoothed tempZ	
		variable delta=dimdelta(smoothed,0)
		make/o/n=(200) xwaveinterp = startZ+(p*delta)
		killwaves smoothed
		
		// fit the Moffat Function to profiles
		duplicate/o tempZ w
		
		// We require an initial guess at the coefficients. find the peak and its position here
		variable peak =wavemax(w)
		findlevel/q w, peak
		variable r0 = v_levelx // find x-index that corresponds to the peak
	
		// initial guesses at coefficients in order: peak amplitude.....position of peak 
		Make/D/N=5/O W_coef
		W_coef[0] = {peak,1,1,r0,0}						
	
		// fit the Moffat function 
		try
			FuncFit/q radialMoffat W_coef w /d ; AbortOnRTE
		catch
			message="ROI# " + num2str(i) + " cannot be fit with a Moffat... skipping"
			print message
			err=GetRTError(1)
			continue
		endtry
		wave fit_w
		
		// save chi2 value for this fit
		chi2Mat_FirstFit[i]=V_chisq
				
		// duplicate and normalise to the fluorescence of the Moffat fit at imaging plane 
	
		variable v=W_coef[4]
	
		// If the baseline is below 0, then we do not want to subtract a negative number. Set this to the minimum value in the trace if so
		// This usually happens if you don't have enough of the profile to fit properly (i.e. the edges are still falling at the limits of your zstack)
		
		if(v<0)
	
			v=wavemin(tempZ)
			message="WARNING: ROI#" + num2str(i) + " has a baseline below zero. Consider a larger zstack"
			print message
	
		endif	
	
		// subtract the baseline of the first Moffat fit -> better than using the minimum unless it's below zero
		duplicate/o tempZ w_sub
		w_sub-=v	
	
		// Normalise to the fluorescence at the imaging plane
		duplicate/o w_sub w_norm
		variable normalizer=interp(startPos,xwaveInterp,fit_w)
		normalizer-=v
		w_norm/=normalizer
	
		// save the coefficients of the first fit
		duplicate/o W_coef W_CoefOrig
	
		// rebuild W_Coef to reflect the previous fit but with adjusted values to account for subtraction and normalization
		peak=wavemax(w_norm)
		findlevel/q w_norm, peak
		r0 = v_levelx
	
		// Fit No. 2 after normalising and subtracting the baseline (W_Coef[4])
		W_coef[0] = {peak,1,1,r0,0}
		try
			FuncFit/q/h="00001" radialMoffat W_coef w_norm /d 	; AbortOnRTE	// perform this fit with the baseline fixed at 0
		catch
			message="ROI# " + num2str(i) + " cannot be fit with a Moffat... skipping"
			print message
			err=GetRTError(1)
			continue
		endtry 		// perform this fit with the baseline fixed at 0
	
		wave fit_w_norm
		chi2Mat[i]=V_chisq
		
		/////////////////////////////////////////////////////////////
				
		
		NormalisedzProfiles[][i]=w_norm[p]			// save the normalised profile (baseline subtracted)
		zPosition[i]=W_coef[3]							// save the peak positions in z for each ROI
		 
		
		if(i==0) // if we're on the first profile then make the zFit waves to accept results
		
			make/o/n=(dimsize(fit_w,0),numRois) zFitCurves, zFitNormCurves
			make/o/n=(5,numRois) zFitCoeffs, zFitNormCoeffs
			
		endif
		
		zFitCurves[][i]=fit_w[p]
		zFitCoeffs[][i]=W_CoefOrig[p] 
		zFitNormCoeffs[][i]=W_Coef[p]
		zFitNormCurves[][i]=fit_w_norm[p]
		setscale/i x,xwave[0],xwave[inf],NormalisedzProfiles, zFitCurves, zFitNormCurves, fit_w, fit_w_norm
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////Exclusion Module //////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		// First check the fit using chi2 values
		variable fitMax=wavemax(fit_w_norm)
		variable chi2thresh=.6 // .6 estimated by eye from the distribution of chi2 values in a typical experiment. This distribution is usually gaussian, with a mean of c.0.15, and SD of 0.1
		variable crit= V_chisq/fitMax	// chi square is sensitive to the scale of the fits, which are not neccessarily the same. Here I normalise by the maximum value in the fit. Since the fit's 
		// baseline is fixed at zero, this represents the full range of the values in the fit. Not sure if this is "right", but it seems to have the desired effect. 
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		if(crit>chi2thresh) 
		
			ROI_Log_Keep_Moff[i]=0
		
			continue
		endif
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		
		variable UpperThreshWidth=10
		variable LowerThreshWidth=4
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		
		// find levels at a half the maximum value of the synapse
		
		variable HalfMax=wavemax(fit_w_norm) 
		HalfMax/=2
		
		FindLevels/q /D=LevelRisew /edge=1 fit_w_norm,HalfMax // find rising level
		

		// If we don't find a rising level, or more than one, then exclude...

		if(V_Flag==2)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif

		if(V_levelsFound!=1)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif

		FindLevels/q /D=LevelFallw /edge=2 fit_w_norm,HalfMax // find falling level
		
		// If we don't find a falling level, or more than one, then exclude...

		if(V_Flag==2)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif

		if(V_levelsFound!=1)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif	
		
		variable LevelFall=LevelFallw[0]
		variable LevelRise=LevelRisew[0]
		variable LevelWidth=LevelFall-LevelRise

		// if this width is wider than threshWidth (typical synapse) then exclude...

		if(LevelWidth>UpperThreshWidth)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif
		
		// if this width is narrower than threshwidth (typical synapse) then exclude...

		if(LevelWidth<LowerThreshWidth)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif
		
		// Use the deriviative to find the range we can correct for. This rests on the assumption that when the derivative approaches zero, we are within noise (have reached the edge of the synapse) and therefore cannot glean useful informaion as to its expected fluorescence)

		differentiate fit_w_norm /D=tempZ_DIV
		killwaves fit_w, fit_w_norm
		// set thresholds of the slope of the Moffat fit to find the bounds of the synapse in x to match with the shifts. If the shifts fall within this range, we can correct
		variable lowerBound=.15
		variable upperBound=-.15
		
		variable LevelUpperW=xwave[inf],LevelLowerW=xwave[0]
		
		// Find the first rise
		FindLevel/q /edge=1 tempZ_DIV,lowerBound 
		
		// if we don't find a rising level then set it to -9.5um, since this means that we do not return to noise with negative deflections, provided we remain inside the zstack. The derivative typically peaks at around 0.3. So we will check that 
		if(V_Flag!=0)
			LevelLowerW=V_levelX
		endif

		// Find the second rising level				
		FindLevel/q/edge=1 tempZ_DIV,upperBound

		// Again, if we don't find a rising level then set it to +9.5um, since this means that we do not return to noise with positive deflections, provided we remain inside the zstack. The derivative typically peaks at around 0.3. So we will check that 
		if(V_Flag!=0)
			LevelUpperW=V_levelX
		endif

		// if the shifts do not fall within the levelRange, exclude...
		//variable HighShift=wavemax(shifts)
		//variable LowShift=wavemin(shifts)

		//if(LowShift<LevelLowerW)
		//	ROI_Log_Keep_Moff[i]=0
		//	continue
		//endif

		//if(Highshift>LevelUpperW)
		//	ROI_Log_Keep_Moff[i]=0
		//	continue
		//endif
		
	endfor
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	// Compute divisor trace
	valFromCoeffs(zFitNormCoeffs,shifts)  
	wave DivTrace
	
	
	
	// check that division factor never falls below 0.1
	for(i=0;i<numRois;i+=1) // cycle through ROIs	
		duplicate/o/free/r=[][i] Divtrace temp
		variable maxDrop=wavemin(temp)
		
		if(maxDrop<0.1)
			ROI_Log_Keep_Moff[i]=0
			continue
		endif 
	endfor
	
	
	// create subtrace from DivTrace ; scaled to starting fluorescence
	//duplicate/o DivTrace subtrace
	//for(i=0;i<numRois;i+=1)
	//	duplicate/o/r=[][i] subtrace SubTracetemp
	//	redimension/n=(-1) SubTracetemp
	//	
	//	duplicate/free/o/r=[][i] QA QAtemp
	//	redimension/n=(-1) QAtemp
	//	setscale/p x,0,1,QAtemp
	//	variable startF=mean(QAtemp,0,14) // find the starting fluorescence for this ROI (baseline)
		
	//	SubTracetemp*=startF
	//	subtrace[][i]=SubTracetemp[p]
	//endfor
	
	// Correct through subtraction
	//MatrixOP/o All_CorrectedQA_SubMoff=QA-SubTrace
	
	// If no axial motion has occured then subtraction will bring the trace to zero. Add starting fluorescence back here
	//for(i=0;i<numRois;i+=1)
	
	//	duplicate/free/o/r=[][i] All_CorrectedQA_SubMoff tempSub
	//	redimension/n=(-1)  tempSub
		
	//	duplicate/free/o/r=[][i] QA QAtemp
	//	redimension/n=(-1) QAtemp
	//	setscale/p x,0,1,QAtemp
	//	startF=mean(QAtemp,0,14) // find the starting fluorescence for this ROI (baseline)
		
	//	tempSub+=startF
	//	All_CorrectedQA_SubMoff[][i]=tempSub[p]
	//endfor

	// Divide the raw trace by the divisor trace to correct
	MatrixOP/o All_CorrectedQA_DivMoff=QA/DivTrace
	
	// Rescale both corrected QA waves
	copyscaling(QA,All_CorrectedQA_DivMoff)
	//copyscaling(QA,All_CorrectedQA_SubMoff)
	
	////////////////////////////////// Exclusion ////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Final check of the time series to ensure there are no NaNs (generated by the DFF0 procedure)//
	// Remove these from the corrected QA, along with any we cannot correct confidently//////////////

	for(i=0;i<numRois;i+=1)

		duplicate/o/free/r=[][i] All_CorrectedQA_DivMoff tempCorr

		if (numtype(tempCorr[0]) == 2) // is this trace NaNs?
			ROI_Log_Keep_Moff[i]=0
		endif
		
	endfor

	variable numCorrected=sum(ROI_Log_Keep_Moff)
	make/o/n=(reclength,numCorrected) CorrectedQA_DivMoff, CorrectedQA_SubMoff, Pre_CorrectedQA
	make/o/n=(numCorrected) ROI_ID_Orig_Moff
	make/o/n=((dimsize(NormalisedzProfiles,0)),numCorrected) zProfilesCorrected=0
	make/o/n=((dimsize(zFitNormCurves,0)),numCorrected) zMoffatFitCorrected=0
	variable count=0

	for(i=0;i<numRois;i+=1)

		if(ROI_Log_Keep_Moff[i]==1)

			duplicate/free/o/r=[][i] All_CorrectedQA_DivMoff CorrectedQA_DivMoffTemp
			redimension/n=(-1) CorrectedQA_DivMoffTemp

		//	duplicate/free/o/r=[][i] All_CorrectedQA_SubMoff CorrectedQA_SubMoffTemp
		//	redimension/n=(-1) CorrectedQA_SubMoffTemp
			
			duplicate/o/free/r=[][i] QA temp
			redimension/n=(-1) temp
			
			duplicate/o/free/r=[][i] zFitNormCurves tempCorrZFit
			redimension/n=(-1) tempCorrZFit
			
			duplicate/o/free/r=[][i] NormalisedzProfiles tempCorrZ
			redimension/n=(-1) tempCorrZ
			
			zProfilesCorrected[][count] = tempCorrZ[p]
			zMoffatFitCorrected[][count]	= tempCorrZFit[p]	
			Pre_CorrectedQA[][count] = temp[p]
			CorrectedQA_DivMoff[][count]=CorrectedQA_DivMoffTemp[p]
		//	CorrectedQA_SubMoff[][count]=CorrectedQA_SubMoffTemp[p]
			ROI_ID_Orig_Moff[count]=i
			count+=1

		endif
	endfor
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
end

// Thomas Ryan 2019 Fitting wizard GUI and drop-down menu functionality for Moffat
// edit TR - 07/02/19 - added 5th coefficient for baseline. change made in both the fitMoffat and the valsfromCoeff functions
Function radialMoffat(w,r) : FitFunc
	Wave w
	Variable r

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(r) = B + A * (1+(R-r0)^2/alpha^2)^-beta
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ r
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = alpha
	//CurveFitDialog/ w[2] = beta
	//CurveFitDialog/ w[3] = r0
	//CurveFitDialog/ w[4] = B (baseline)

	return w[4]+(w[0] * (1+(R-w[3])^2/w[1]^2)^-w[2])
End

Window Graph0() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(343.5,90.5,738,299) wSmooth,blah,fit_wSmooth
EndMacro

// Thomas Ryan 28/01/19
// Fits a Moffat function to the input wave, computing two distinct fits to the unnormalised and the
// normalised profiles for subtractive and divisive methods of correction respectively. Stores chi
// squared values for exclusion later

function fitMoffat(w)

	wave w
	
	
	// We require an initial guess at the coefficients. find the peak and its position here
	variable peak =wavemax(w)
	findlevel w, peak
	variable r0 = v_levelx // find x-index that corresponds to the peak
	
	// initial guesses at coefficients in order: peak amplitude.....position of peak 
	Make/D/N=5/O W_coef
	W_coef[0] = {peak,1,1,r0,0}						
	
	// fit the Moffat function 
	FuncFit radialMoffat W_coef w /d
	make/o/n=1 chi=V_chisq
	wave fit_w
end

// Thomas Ryan 28/01/19
////////// Short function to find Y at a given X using the equation //////////

function valFromCoeffs(w,xW)
	wave w, xW
	
	variable nR = dimsize(w,1) 
	variable nS = dimsize(xW,0)
	
	make/o/n=(nS,nR) DivTrace=0
	
	DivTrace =  (w[0][q] * (1+(xW[p]-w[3][q])^2/w[1][q]^2)^-w[2][q])+w[4][q] // This plugs the x position into the Moffat equation using the coefficients found during fitting
	
	
end


//############# Utility Functions Thomas Ryan 2019###############\\

// Utility function to create a locomotion trace from a rotary encoder recorded through a scanimage channel 
// Can be scaled to convert Piezo monitor recorded through a scanimage channel with different input voltages 
function LocTrace(WName)

	string WName

	wave w=$WName

	variable dimX=dimsize(w,0),dimY=dimsize(w,1),dimZ=dimsize(w,2),i
	make/o/free/n=(dimZ) VTrace
	for(i=0;i<dimZ;i+=1)

		duplicate/free/o/r=[0][0][i] w temp

		wavestats/q temp

		VTrace[i]=V_avg
	endfor

	variable baseline=wavemin(VTrace)//Vtrace[0]//mean(Vtrace,0,14)
	MatrixOP/o/free VoltTrace=(VTrace-baseline)
	duplicate/o VoltTrace LocoTrace
	smooth 3,locotrace

	// as voltage drops, objective moves down
	//PiTrace*=-0.647668 // Input range -20V, 20V
	//PiTrace*=-0.198432 // Input range -10V, 10V 
end

// Utility function to compute deltaF over F - 
//	WARNING - this is bare bones: operates on the mode as baseline, no correction for bleach or neuropil at this stage

function DFF0Trace(WName)

	string WName
	string outName=WName+"_DFF0"

	duplicate/o $WName ROI
	duplicate/o $WName DFF0w

	Make/o/n=100 Histo
	Histogram /b=1 ROI, Histo	// Histogram it
	wavestats /q/m=1 histo				// compute the wavestats... includes V_maxLoc -> the x location of the largest bin
	
	variable F0 = V_MaxLoc
	
	DFF0w[]=(ROI[p] -F0)/ F0
		
	duplicate/o DFF0w $outName
		
end
		
		

