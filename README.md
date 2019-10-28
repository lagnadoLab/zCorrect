# zCorrect

Collection of procedures for the estimation of z-displacement and correction of population activity for resulting intensity changes

Experimental Acquisition and Data Handling Note;

A reference stack is required for estimation of z-displacements and modelling of intensity as a function of z position. We recommend taking at least 10 volumes 
covering 10um above and below the imaging plane at 0.5um intervals. Both anatomical and reporter chanels should be present in both the time series and the reference
volume images. Volume and time series data should all be in this same folder. The procedure assumes there is a 1:1 ratio of recordings and z-stack series. 
The file list is sorted and time series and zstacks are assumed to sort into two sequential lists (not interleaved). Easiest way to ensure this works is to have the 
filename begin with "z..." if it is a z-stack, and be suffixed with the experiment number corresponding to the time series.

E.G.

TimeSeries_1.tiff
Zstack_1.tiff
TimeSeries_2.tiff
Zstack_2.tiff

becomes:

TimeSeries_1.tiff
TimeSeries_2.tiff
Zstack_1.tiff
Zstack_2.tiff


If another naming convention is applied, please ensure that files are sorted in this way before running the procedure.

Once acquired, follow the steps below to generate and correct population activity.


GENERAL USE WITH BATCH SCRIPT

Procedure file: zCorrectMasterBatchLite.ipf

	All image processing and displacement estimation is performed in the function "zCorrectMaster_Batch". Input parameters are as follows:

		i. numVols	=	[the number of volumes in your reference stack image series]
		ii. ImPlane	=	the estimated position of your imaging plane for initial realignment. N.B. This need not be perfect. 
		iii. nFolders	=	OPTIONAL: the number of folders you wish to process. To enter this parameter use the syntax "nFolders=n"


	Input in procedure itself:	

		i. path1 = OPTIONAL: a string variable specifying the full path to the folder containing your raw data. If not provided, you will be asked to navigate to a folder.

	Additional paths may be added for the processing of multiple folders of data in sequence. If desired, then add further string variables as demonstrated in the provided 
	procedure. These will be called path2, path3 etc.



	Example:

	If we have one reference stack image series consisting of 10 volumes each of 41 slices, and one time series measured in the central position in the folder ...\\MyDataSession\, run:

	/////////////////////////////

	zCorrectMaster_Batch(10,20,nFolders=1)

	/////////////////////////////

	Here, 10 is the number of volumes, 20 is the slice number we expect the time series to most closely match (central plane), and 1 is the number of folders.

	Additionally you may specify...

	string path2 =  "...\\MySecondDataFolder\" 
	string path3 = "...\\MyThirdDataFolder\"

	...within the procedure if desired. Only [nFolders] will be processed



By default these functions will utilise SARFIA procedures already loaded into IGOR. SARFIA is freely available from the IgorExchange and is described in Dorostkar et al. 2010.
N.B. Once SARFIA has been downloaded, search for the LoadScanImage.ipf and replace it with the procedure file in this repository.

The procedure will call a series of functions to perform the following processes:

1. Loading and channel splitting of time series recording
2. Loading and channel splitting of the reference stack image series
3. Construction of average reference stack by realigning each volume and taking an average projection of the result
4. Registration of the time series images to itself
5. Registration of the reference stack to the realigned time seriesy
6. Cropping of time series and reerence stack images to remove zero-padded boundaries
7. Displacement estimation

These procedures generate these waves:

- CroppedCh1 = registered and cropped raw time series data for Ch1
- CroppedCh2 = same for Ch2

- CroppedRef_Ch1_reg = registereed and cropped reference stack for Ch1
- CroppedRef_Ch2_reg = same for Ch2

- zestimate = the estimated axial displacement over time


Additionally, the following waves and variables are saved

-ZShiftY/X = the largest shift in the y or x dimensions for all registered volumes
-lowestY/X = the largest negative shift from all registrations procedures... used to crop images
-highestY/X = the largest positive shift  

- CroppedRefCh1/2 = cropped reference stack before realignment for Ch1 and 2
- CorrMap = the 2D correlelogram for the zestimation procedures

- The original movies are also retained in their split channels


The experiment is then named as per the time series recording and saved in the same folder as an IGOR experiment file beore removing all waves and variables and moving on to the next image pair


Reference stack functions:

CreateRefStack_3dReg = will create an average reference stack from a series of volumes as described. Volume series are split into individual volumes and registered in 3D with the central (chronologically) volume. 

Inputs: 	string variables:	zStackCh1Name	=	string variable, the name of Ch1 of the volume series
					zStackCh2Name	=	string variable, the name of Ch2 of the volume series 
		numeric variable:	numVols		=	numeric variable informing the number of volumes in the series

Outputs:	waves:			VolumeCh1	=	Ch1 reference realigned stack
					VolumeCh2	=	Ch2 reference realigned stack

		numeric variables:	zShiftX		= 	largest shift in x across all realignment steps 
					zShifty		= 	largest shift in y across all realignment steps

CreateRefStack_PlaneToTS = will create an average reference stack from a series of volumes as described. Planes of each volumes are first realigned to each other, then each plane is realigned to the adjacent plane closest to the central plane of its own volume

Inputs and outputs are the same as above, but outputs are suffixed with "_PlReg"

N.B. the scripts also generates maps of the SD of the reference volumes for each channel which are used later. These are not currently stored but this code can be easily modified to store this information.

Other functions:

RegisterAll		=	automatically registers the time series images, and subsequently the reference stack is shifted to match the realigned time series (faster than realigning the entire time series (1000s of frames) to the reference stack)

Inputs:		string variables: 	Ch1MovieName		=	Ch1 time series wave name	 
					Ch2MovieName		=	Ch2 time series wave name
					RefStack_Ch1		=	Ch1 reference stack wave name	
					RefStack_Ch2		=	Ch2 reference stack wave name
					VolCh1SD		=	Standard deviation volume produced in CreateRefStack scripts
					Ch1PosSDMovieName	=	Average stacks plus S.D. produced in CreateRefStack scripts
					Ch1NegSDMovieName	=	Average stacks minus S.D. produced in CreateRefStack scripts

		numerical variable:	ImPlane		=	numerical variable estimating the initial imaging plane (usually the central plane)

					The following are the largest positive and negative shifts in x and y dimensions over all registration proceudres

					highestX		
					highestY
					lowestX
					lowestY

Inputs:			Realigned versions of all ImageInputs suffixed with "_reg"


CropScalarShifts	=	following registration, this script takes all time series and reference stacks and crops them by the single largest shift in each dimension from the realignment steps applied to reference stack images and time series. 

Inputs:		string variables:	RegCh1Movie	=	time series Ch1 wave name
					RegCh2Movie	=	time series Ch2 wave name
					RegCh1Ref	=	reference stack Ch1 wave name
					RegCh2Ref	=	reference stack Ch2 wave name
					RegPosSDRef	=	reference stack plus S.D.
					RegNegSDRef	=	reference stack minus S.D.
					RegVolSDRef	=	reference stack S.D. projection

		wave variables:		shiftX		=	Vector of x transforms from time series	
					shiftY		=	Vector of y transforms from time series	

		numeric variables:	ZshiftX		=	largest shift during registration of reference stack in x
					ZshiftY		=	largest shift during registration of reference stack in y

		The following are the largest positive and negative shifts in x and y dimensions over all registration proceudres

					highestX
					highestY	
					lowestX
					lowestY

Z_estimate		=	automatically estimates the axial location of each time frame of a time series from within a reference stack through a cross correlation procedure as described. 

Inputs:		wave variables:		MovieStackIn	=	time series after realignment and pre-processing
					RefStackIn	=	reference stack after realignment and pre-processing

Outputs:	wave variable:		zestimate	=	vector variable describing the estimated axial offset from the central plane for each time frame. 

This function operates much faster on downsampled images, options for which are provided within the script. The parameters for this downsampling step are defined in the first few lines of the procedure  
LL - Note the variables that you should check: zstep, downsamplez (= 1), downsamplexy (= 2).  
For beads do not downsample.  For vessels, downsampling xy by factor 2 has no effect on Zestimate.  
Even downsampling by 4 works (and saves a lot of time!).  The downsample variables should be 1 or a multiple of 2. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CorrectQA = 	This procedure uses existing waves representing the estimated axial displacement over time generated by "z_estimate", and the zProfile of each ROI and the experimental 
 		activity traces from the same ROIs. Each profile is normalised to the interpolated fluorescence value expected at the imaging plane. These profiles are then fit with
		a Moffat function;a variation of a t-distribution which is well suited to describe how synapses recorded in-vivo are convolved with the PSF of the microscope in the 
		axial dimension. These fits are then used as a look-up table, to generate a vector describing how we predict changes to fluorescence may occur due to axial displacements 
		over the course of an experiment. This vector is then used to correct the fluorescent activity traces of the ROI to more closely reflect the true intensity in our intended plane. 

Inputs: 	No direct inputs; instead the script searches the existing waves in the first few lines.

		shiftsName	=	wavelist("zestimate","","")			//	default name for the axial shifts over time generated using the "z_estimate" procedure.
		zProfilesName	=	wavelist("CroppedRef_Ch1_DFF0","","")		//	activity matrix generated by user's procedures following ROI extraction (e.g. SARFIA, 
		QAName		=	wavelist("CroppedCh1_DFF0_QA","","")

Outputs:	All_CorrectedQA_DivMoff	=	All Corrected traces
 		CorrectedQA_DivMoff	=	Corrected traces that have passed exclusion criteria
		PreCorrectedQA		=	Original, uncorrected traces for those that passed exclusion criteria
 		CorrectedROI_Id		=	containing the original ROI indexes for back referencing
		zFitCurves		=	profile curve fits (moffat or gaussian
		zFitCoefs		=	the coefficients for all profile fits

Another wave, CorrectedQA_SubMoff is similar to CorrectedQA_DivMoff, but computed by estimating a new baseline for the ROI and subtracting this from the raw trace. In this case, the original baseline fluorescence is added back to the time	series following subtraction. 


The first flag within the code below allows you to select between performing the above procedure using a fitted Moffat function, or using a more empirical fitted spline. This may be more appropriate for ROIs which have arbitrary size and shape in the z-dimension. Set to: 

		-> 1 for Moffat fitting
		-> 2 for spline fitting
		-> 3 for both.

In the case of empirical correction through a spline the results waves listed above will be generated with the suffix "Emp"

Several exclusion criteria are included so that only synapses we are confident of correcting are passed. Those that are excluded are still "corrected" and can be found in the corrected waves prefixed with "All_", the "CorrectedQA" waves only contain the passed and corrected 
waves. Original ROI index of passed ROIs is stored in "ROI_Id_Orig"


CoQAMoff


CoQAEmp

Other Functions:

radialMoffat	=	Fitting wizard GUI and drop-down menu functionality for Moffat	
fitMoffat	=	Fits a Moffat function to the input wave, computing two distinct fits to the unnormalised and the normalised profiles for subtractive and divisive methods of correction respectively. Stores chi squared values for exclusion later
valFromCoeffs	=	Function to find Y at a given X using the Moffat equation
LocTrace	=	Utility function to create a locomotion trace from a rotary encoder recorded through a scanimage channel. Can be scaled to convert Piezo monitor recorded through a scanimage channel with different input voltages
DFF0Trace	=	Utility function to compute simple deltaF over F. Should only be used for basic observation and with great care: computes baseline flourescence as the mode of the trace and does NOT correct for bleach or neuropil. 
