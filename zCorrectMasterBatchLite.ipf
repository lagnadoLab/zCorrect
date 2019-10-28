#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "Sarfia"
#include "LoadScanImage"
//###################################################################################################################
// Thomas Ryan	-	2019

// zCorrectMaster_Batch:	This batch script will cycle through a number of folders containing raw data from ScanImage and automatically create a registered and cropped average reference stack and time series. 
// 								It then uses these to estimate displacements in z, and saves the experiment using the name of each time series file. 

// INPUTS 	: 			numVols 				= 		the number of volumes collected for construction of the average reference stack
//							ImPlane 				= 		the approximate frame number within the acquired volume that the time series most closely matches; usually the central plane i.e. floor(number of slices / 2)
//							nFolders				= 		OPTIONAL parameter specifying number of folders defined as string variables called path1...n

// OUTPUTS	:			CroppedCh1 			= 		registered and cropped raw time series data for Ch1
//							CroppedCh2 			= 		same for Ch2
//							CroppedRef_Ch1_reg 	= 		registereed and cropped reference stack for Ch1
//							CroppedRef_Ch2_reg 	= 		same for Ch2
//							zestimate 			= 		the estimated axial displacement over time


//	Additionally, the following waves and variables are saved

//							ZShiftY/X 			= 		the largest shift in the y or x dimensions for all registered volumes
//							lowestY/X 			= 		the largest negative shift from all registrations procedures... used to crop images
//							highestY/X 			= 		the largest positive shift  
//							CroppedRefCh1/2 		= 		cropped reference stack before realignment for Ch1 and 2
//							CorrMap 				= 		the 2D correlelogram for the zestimation procedures

//	The original movies are also retained in their split channels

// N.B. The method by which the reference stack is generated depends on which procedure you choose to use.
// Please see lines 126 & 127 and choose 

function zCorrectMaster_Batch(numVols,ImPlane,[nFolders])

	variable numvols, ImPlane, nFolders
	
	string path1="\MyFirstDataFolder\""
	// String variable with full path to directory containing raw ScanImage files. Can add more paths here called path2, path3 etc, but ensure you specifiy how many paths there are in the input parameter -  nFolders
	//	string path2=
			
	string S_path
	variable FolderNum
	 	
	if(paramisdefault(nFolders))										//	Check if multiple folders
		nFolders=1
	endif
	
	variable chkPath=strlen(path1)
	
	if(chkPath==0)
		getfilefolderinfo/D													//Open a dialogue choose the folder to process if no path is specified
		path1=S_path
	endif
		
	FOR(FolderNum=0;FolderNum<nFolders;FolderNum+=1)				// cycle through the number of folders
		switch(FolderNum)
			case 0:
				S_path=path1
				break
			case 1:				
				//			S_path=path2
				break
			case 2:
				//			S_path=path3
				break
			case 3:
				//			S_path=path4
				break
		endswitch
	
		newpath/O pathstr, S_path			
		string filelist= indexedfile(pathstr,-1,".Tif")	
		filelist=SortList(filelist, ";", 16)				  	  	//Create a list of all the .tif files in the folder. -1 parameter addresses all files.
		Variable Length=itemsinlist(filelist)
		string pathstr=s_path
		print pathstr
		string RecFileName, zstackFileName
	
		variable  ZnChannels,Rnchannels,rL,zL,k,i
	
		/// cycle through experimental files and zstacks... assumes that there is a 1:1 ratio of recordings and zstacks

		for(i=0;i<(Length/2);i+=1)
	
			RecFileName= StringFromList(i, filelist)
			zstackFileName=StringFromList((i+(Length/2)), filelist)
			rL=strlen(RecFileName)
			zL=strlen(zstackFileName)
			string message="Loading Recording " + num2str(i) + " of " + num2str(Length/2)
			print message
			AutoLoadScanImage(pathstr, RecFileName)

 
			string waveN= RecFileName[0,(rL-5)]
			string c1=waveN+"_ch1"
			string c2=waveN+"_ch2"
			string c3=waveN+"_ch3"
			string c4=waveN+"_ch4"
			duplicate/FREE/o/r=[][][0] $waveN firstrecpic
			redimension/n=(-1,-1) firstrecpic     
			string zwaveN= zstackFileName[0,(zL-5)]
			string zc1=zwaveN+"_ch1"
			string zc2=zwaveN+"_ch2"
			string zc3=zwaveN+"_ch3"
			RnChannels =  nChannelsFromHeader(firstrecpic)
		
			print "Splitting channels of recording"

			SplitChannels($waveN,RnChannels) 
			killwaves $waveN    		     		    	      // Split the channels in n channels and name it automatically  'blabla_ch1'and 'blabla_ch2'
		
			print "Loading zStack for this experiment"
			AutoLoadScanImage(pathstr, zstackFileName)
			duplicate/FREE/o/r=[][][0] $zwaveN firstzpic  
			redimension/n=(-1,-1) firstzpic
			ZnChannels = nChannelsFromHeader(firstzpic)
			print "Splitting the zstack channels"
			SplitChannels($zwaveN,ZnChannels) 	
			killwaves $zwaveN firstzpic
		
			variable/g ZshiftX,ZshiftY,highestX,highestY,lowestX,lowestY
			Print "Creating a reference stack"
			
			// Replace one function with the other depending on which method of reference stack creation you wish to use
			//CreateRefStack_3dReg(zc1,zc2,numVols)
			CreateRefStack_PlaneToTS(zc1,zc2,numVols)
		
			//outputs: shiftXReg, shiftYReg////////////////////////
			string Ch1Ref="VolumeCh1_PlReg"
			string Ch2Ref="VolumeCh2_PlReg"
			string VolCh1SD="VolumeStDevCh1_Plreg" // N.B.  remove "_PlReg" suffix if using CreateRefStack_3DReg
			string Ch1PosSDRef="VolumePosSDCh1"
			string Ch1NegSDRef="VolumeNegSDCh1"		
			///////////////////////////////////////////////////////
			
			print "Registering all Images"
			
			RegisterAll(ImPlane,c1, c2, Ch1Ref, Ch2Ref,VolCh1SD,Ch1PosSDRef,Ch1NegSDRef)
			//outputs: ////////////////////////////////////////////
			wave Ch1MovieOut,Ch2MovieOut, shiftX,shiftY
			string RegCh1Movie = c1 + "_reg"
			string RegCh2Movie = c2 + "_reg"
			string RegCh1Ref = Ch1Ref + "_reg"
			string RegCh2Ref = Ch2Ref + "_reg"
			string RegPosSDRef = Ch1PosSDRef + "_reg"
			string RegNegSDRef = Ch1NegSDRef + "_reg"
			string RegVolSDRef = VolCh1SD + "_reg"
			///////////////////////////////////////////////////////
			
			killwaves Ch1MovieOut,Ch2MovieOut
			
			CropScalarShifts(RegCh1Movie,RegCh2Movie,RegCh1Ref,RegCh2Ref,RegPosSDRef,RegNegSDRef,RegVolSDRef,shiftX,shiftY,ZshiftX,ZshiftY,highestX,highestY,lowestX,lowestY)
			// outputs: ///////////////////////////////////////////
			wave CroppedCh2,CroppedRefCh2,CroppedCh1, CroppedRefCh1
			///////////////////////////////////////////////////////
			killwaves $RegCh1Movie, $RegCh2Movie,$RegCh1Ref,$RegCh2Ref
			Print "Estimating z-displacements"
		
			Z_estimate(CroppedCh2, CroppedRefCh2)
			// outputs: ///////////////////////////////////////////
			wave Zestimate
			
			wave M_stdvImage, W_coef, W_sigma
			killwaves M_stdvImage,W_coef, W_sigma

			// realign zstack to actual rather than presumed starting position so that ROI maps match
		
			//first find the start position
			duplicate/o/r=[0,5] zEstimate zEstToMean
			variable startPos=mean(zEstToMean)  
		
			// then find the slice this most closely corresponds to
			variable sliceN=(round(startPos/0.5))+ImPlane
		
			//now register the zstack to the time series using this plane to determine the shifts
			// make the testwave from the reference stack
		
			duplicate/FREE/o/r=[][][sliceN] CroppedRefCh2 tempTestWave
			redimension/s/n=(-1,-1) tempTestWave
			MatrixOP/FREE/o tempTestWave_NaNBusted = ReplaceNaNs(tempTestWave, 0)
		
			// make a reference wave from the average of the first 5 frames of the recording
			make/o/FREE /n=(dimsize(CroppedCh2,0),dimsize(CroppedCh2,1),5) reftoA 
						
			for(k=0;k<5;k+=1)
				duplicate/FREE/o/r=[][][k] CroppedCh2, tt
				redimension/s /n=(-1,-1) tt
				reftoA[][][k]=tt[p][q]
			endfor
			
			ImageTransform averageImage reftoA
			wave m_AveImage
			duplicate/FREE/o m_AveImage Ch2refm
			killwaves m_AveImage
			redimension/s Ch2Refm
			MatrixOP/FREE/o Ch2Ref_NaNBusted = ReplaceNaNs(Ch2Refm, 0)
			
			// now register
			imageregistration /q /rot={0,0,0} /stck /pstk /csnr=0 /refm=0 /tstm=0 testwave=tempTestWave_NaNBusted, refwave=Ch2Ref_NaNBusted 
			wave m_regParams
			duplicate/o m_regParams pwaveRef
			killwaves m_regparams

			//now apply the shifts to every slice of the zstack in both channels
			make/free/o/n=(dimsize(CroppedRefCh2,0),dimsize(CroppedRefCh2,1),dimsize(CroppedRefCh2,2)), RefOutCh1_Collected, RefOutCh2_Collected
		
			for (k=0;k<dimsize(CroppedRefCh2,2);k+=1)
				duplicate/FREE/o /r=[][][k] CroppedRefCh2 StackReg_Ch2
				duplicate/FREE/o /r=[][][k] CroppedRefCh1 StackReg_Ch1
				redimension /s /n=(-1,-1) StackReg_Ch1, StackReg_Ch2

				imageregistration /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=StackReg_Ch2 // apply shifts to rest of the z-stack in Ch2 (assumes rigid shift through whole volume)
				wave M_RegOut
			
				duplicate/FREE/o M_RegOut RefOutCh2
				killwaves M_regout

				imageregistration /stck /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=StackReg_Ch1 // apply shifts to rest of the z-stack in Ch1 (assumes rigid shift through whole volume) 
				wave M_RegOut
			
				duplicate/FREE/o M_RegOut RefOutCh1
				MatrixOP/FREE/o RefOutCh1_NaNBusted = ReplaceNaNs(RefOutCh1, 0)
				MatrixOP/FREE/o RefOutCh2_NaNBusted = ReplaceNaNs(RefOutCh2, 0)	

			
				RefOutCh1_Collected[][][k]=RefOutCh1_NaNBusted[p][q]
				RefOutCh2_Collected[][][k]=RefOutCh2_NaNBusted[p][q]
			endfor
			
			string RefStackOutName_Ch2="CroppedRef_Ch2_reg",RefStackOutName_Ch1="CroppedRef_Ch1_reg"		
			duplicate/o RefOutCh1_Collected $RefStackOutName_Ch1
			duplicate/o RefOutCh2_Collected $RefStackOutName_Ch2
			copyscaling(CroppedRefCh1, $RefStackOutName_Ch1)
			copyscaling(CroppedRefCh2, $RefStackOutName_Ch2)
						
			// tidy waves
			wave M_regMaskOut,xwave,absshiftY,absshiftX,shiftXY,VolumeStDevCh2_Plreg,VolumeStDevCh1_Plreg, tempsubA, tempSubB,VolumeStDevCh2,VolumeStDevCh1,M_StdvImage,CroppedSDRefCh1,CroppedNegSD,CroppedPosSD,VolumeStDevCh1_Plreg_reg,VolumeNegSDCh1_reg,VolumePosSDCh1_reg,shiftx,shifty,VolumeNegSDCh2,VolumeNegSDCh1,VolumePosSDCh1,VolumeCh2_Plreg,VolumeCh1_Plreg,VolumeCh1,VolumeCh2
			killwaves M_regOut,M_regMaskOut,pwaveref,M_stdvImage,zesttomean,xwave,absshiftY,absshiftX,shiftXY,VolumeStDevCh2_Plreg,VolumeStDevCh1_Plreg, tempsubA, tempSubB,VolumeStDevCh2,VolumeStDevCh1,M_StdvImage,CroppedSDRefCh1,CroppedNegSD,CroppedPosSD,VolumeStDevCh1_Plreg_reg,VolumeNegSDCh1_reg,VolumePosSDCh1_reg,shiftx,shifty,VolumeNegSDCh2,VolumeNegSDCh1,VolumePosSDCh1,VolumeCh2_Plreg,VolumeCh1_Plreg,VolumeCh1,VolumeCh2
			
			// now save the experiment		
			string saveName=WaveN+".pxp"
			message="Saving experiment as "+saveName
			print message
			saveExperiment/P=pathstr as saveName
			print "done"
			print "Killing all waves before we move on to the next file" 
			killAllWaves()
			KillVariables/A/Z
			
		endfor
	ENDFOR
end
//###################################################################################################################

//###################################################################################################################
// Thomas Ryan - 2019
//	CreateRefStack_3dReg = will create an average reference stack from a series of volumes as described. 
// Volume series are split into individual volumes and registered using a volumetric sub-pixel affine 
// transformation algorithm using the central (chronologically) volume as a template 

// Inputs: 	zStackCh1Name		=		string variable, the name of Ch1 of the volume series
//				zStackCh2Name		=		string variable, the name of Ch2 of the volume series 
//				numVols				=		numeric variable informing the number of volumes in the series

// Outputs:	VolumeCh1				=		Ch1 reference realigned stack
//				VolumeCh2				=		Ch2 reference realigned stack
//				zShiftX				= 		largest shift in x across all realignment steps 
//				zShifty				= 		largest shift in y across all realignment steps

function CreateRefStack_3DReg(zStackCh1Name,zStackCh2Name,numVols)

	string zStackCh1Name, zStackCh2Name
	variable numVols
	variable filter=0

	wave Ch1W=$zStackCh1Name
	wave Ch2W=$zStackCh2Name
	
	variable wCh1Min=wavemin(Ch1W)
	variable wCh2Min=wavemin(Ch2W)	
	
	// offset image values to remove negative values if they exist - takes an extremely long time!!!
	//if(wCh1Min<0)
	//Ch1W+=(wavemin(Ch1W)+1) 
	//endif
	//if(wCh2Min<0)
	//Ch2W+=(wavemin(Ch2W)+1) 
	//endif
		
	variable dimX=dimsize(Ch2W,0)
	variable dimY=dimsize(Ch2W,1)
	variable dimZConcat=dimsize(Ch2W,2)
	variable dimZVol=(dimZConcat/numVols)

	string VolCh1Name="Ch1Volume_"
	string VolCh2Name="Ch2Volume_"

	string PlaneCh1Name="Ch1Plane_"
	string PlaneCh2Name="Ch2Plane_"

	string PlaneCh2WaveName,PlaneCh2AvWaveName,VolCh2WaveName,VolCh2RegWaveName
	string PlaneCh1WaveName,PlaneCh1AvWaveName,VolCh1WaveName,VolCh1RegWaveName
	
	make/o/n=(dimX,dimY,DimZVol) VolumeCh1=0	
	make/o/n=(dimX,dimY,DimZVol) VolumeCh2=0 // make matrices to accept the final volumes

	make/o/n=(dimX,dimY,numVols) PlaneTempCh1=0 // make matrices to accept each plane later
	make/o/n=(dimX,dimY,numVols) PlaneTempCh2=0

	variable i,j

	for(i=0;i<numVols;i+=1) // cycle through the number of volumes to preallocate each volume
		if(i==(floor(numVols/2))) // then this is the middle volume. We'll use this one to register the others
		
			string midCh2VolWaveName=VolCh2Name + num2str(i)
			VolCh1WaveName=VolCh1Name + num2str(i)
			
			duplicate/o/r=[][][(i*dimZVol),(((i*dimZVol)+dimZVol)-1)] Ch2W $midCh2VolWaveName
			duplicate/o/r=[][][(i*dimZVol),(((i*dimZVol)+dimZVol)-1)] Ch1W $VolCh1WaveName
			
			setscale/p z,0,1, $midCh2VolWaveName
			redimension/s $midCh2VolWaveName
		else							
			VolCh1WaveName=VolCh1Name + num2str(i)
			VolCh2WaveName=VolCh2Name + num2str(i)

			duplicate/o/r=[][][(i*dimZVol),(((i*dimZVol)+dimZVol)-1)] Ch1W $VolCh1WaveName
			duplicate/o/r=[][][(i*dimZVol),(((i*dimZVol)+dimZVol)-1)] Ch2W $VolCh2WaveName

			setscale/p z,0,1, $VolCh1WaveName
			setscale/p z,0,1, $VolCh2WaveName
						
			redimension/s $VolCh1WaveName
			redimension/s $VolCh2WaveName
		endif
	endfor

	// now we've made the split waves we can register them one by one to the middle volume
	// while we do this we can build the waves for each plane we will average later
	
	for(i=0;i<numVols;i+=1) // cycle through volumes
		
		VolCh1WaveName=VolCh1Name + num2str(i)
		VolCh2WaveName=VolCh2Name + num2str(i)
		
		VolCh1RegWaveName=VolCh1WaveName+"_reg"
		VolCh2RegWaveName=VolCh2WaveName+"_reg"
		
		if(i==(numVols/2))
		
			duplicate/o $VolCh2WaveName $VolCh2RegWaveName
			duplicate/o $VolCh1WaveName $VolCh1RegWaveName
			
		else
		
			Imageregistration/q /stck /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=$VolCh2WaveName, refWave=$midCh2VolWaveName // register volume with the middle volume in Ch2
			wave M_RegOut, M_RegParams
			MatrixOP/o/free temp = ReplaceNaNs(M_Regout, 0)	
			copyscaling($VolCh2WaveName, temp)
			duplicate/o temp $VolCh2RegWaveName
			Imageregistration/q /stck /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=M_RegParams testwave=$VolCh1WaveName // apply transformation to Ch1

			wave M_RegOut
			MatrixOP/o/free temp = ReplaceNaNs(M_Regout, 0)	//replace NaN's with 0
			copyscaling($VolCh1WaveName, temp)
			if(filter!=0)
				MatrixFilter/o/N=(filter) gauss temp
			endif
			duplicate/o temp $VolCh1RegWaveName
			
		endif
		
		
		for(j=0;j<dimZVol;j+=1) // cycle through the planes
		
			PlaneCh1WaveName=PlaneCh1Name+num2str(j)
			PlaneCh2WaveName=PlaneCh2Name+num2str(j)
			
			if(i==0)
			
				duplicate/o PlaneTempCh1 $PlaneCh1WaveName	// if this is the first volume we need to make the waves for each plane
				duplicate/o PlaneTempCh2 $PlaneCh2WaveName
				
			endif
			
			duplicate/o $PlaneCh1WaveName tempCh1Plane			
			duplicate/o $PlaneCh2WaveName tempCh2Plane
			
			duplicate/o/r=[][][j] $VolCh1RegWaveName tempCh1
			duplicate/o/r=[][][j] $VolCh2RegWaveName tempCh2
			
			setscale/p z,0,1, tempCh1
			setscale/p z,0,1, tempCh2
						
			tempCh1Plane[][][i] = tempCh1[p][q]
			tempCh2Plane[][][i] = tempCh2[p][q]
			
			duplicate/o tempCh1Plane $PlaneCh1WaveName
			duplicate/o tempCh2Plane $PlaneCh2WaveName
			
		endfor
	endfor
	
	/// now filter and average project each plane
	for(j=0;j<DimZVol;j+=1)
	
		PlaneCh1WaveName=PlaneCh1Name+num2str(j)
		PlaneCh2WaveName=PlaneCh2Name+num2str(j)
		
		ImageTransform averageImage $PlaneCh1WaveName
		wave M_AveImage
		duplicate/o M_AveImage, tempCh1
		if(filter!=0)
			MatrixFilter/o/N=(filter) gauss tempCh1
		endif
		VolumeCh1[][][j] = tempCh1[p][q]
		
		ImageTransform averageImage $PlaneCh2WaveName
		wave M_AveImage
		duplicate/o M_AveImage, tempCh2
		if(filter!=0)
			MatrixFilter/o/N=(filter) gauss tempCh2
		endif
		VolumeCh2[][][j] = tempCh2[p][q]
		
		killwaves $PlaneCh1WaveName, $PlaneCh2WaveName, M_AveImage
	endfor
	
	// clean up after all this
	killwaves M_RegOut,PlaneTempCh1, tempCh1Plane, tempCh1, M_RegOut, PlaneTempCh2, tempCh2Plane, tempCh2

	for(i=0;i<numVols;i+=1)
	
		VolCh1WaveName=VolCh1Name + num2str(i)
		VolCh2WaveName=VolCh2Name + num2str(i)
		
		VolCh1RegWaveName=VolCh1WaveName+"_reg"
		VolCh2RegWaveName=VolCh2WaveName+"_reg"
			
		killwaves $VolCh1WaveName,$VolCh1RegWaveName, $VolCh2WaveName,$VolCh2RegWaveName
		
	endfor
	duplicate/FREE/o/r=[0][]	M_regParams tempshiftX
	duplicate/FREE/o/r=[1][]	M_regParams tempshiftY	
	setscale/p x,0,1, tempshiftX, tempshiftY 
	setscale/p y,0,1, tempshiftX, tempshiftY
	variable/g lowestX=0,lowestY=0,highestX=0,highestY=0
	lowestX=wavemin(tempshiftX)
	lowestY=wavemin(tempshiftY)
	highestX=wavemax(tempshiftX)
	highestY=wavemin(tempshiftY)
end
//###################################################################################################################

//###################################################################################################################
// Thomas Ryan - 01/02/19

// Creates a 3D average volume to be used as a reference for cross correlation and z-motion correction
// Channel 2 is assumed to be the structural marker and all alignments are performed in this channel
// Each z-stack should be a concatenated stack of >=3 volumes. 
// z-stacks are split into their individual planes, filtered with a gaussian, and realigned to the central (chronological) image.
// Following this, each plane is then aligned to the adjacent plane closest to the central plane.
// all shifts in X and Y dimensions are stored, and may be used to search for shearing artefacts, but not in this procedure.
// These shifts are also used later in "CropShifts" procedures, used to crop the zero padded edges introduced by imageregistration.
// Alignments are passed to Channel 1, and a reference stack of each channel contrusted from the average projection of each plane.

// Check filter variable to change the size of the gaussian filter

function CreateRefStack_PlaneToTS(zStackCh1Name,zStackCh2Name,numVols)

	string zStackCh1Name, zStackCh2Name
	variable numVols

	variable filter=0 // set the SD of the gaussian filter applied to each plane before averaging. Set to zero for no filtering
	variable regCh=2


	variable templowestX,templowestY,temphighestX,temphighestY
	variable/g highestX,highestY,lowestX,lowestY
		
	wave Ch1W=$zStackCh1Name
	wave Ch2W=$zStackCh2Name
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// OPTIONAL MODULE TO OFFSET NEGATIVE VALUES IN YOUR IMAGES IF THEY EXIST. TAKES AN EXTREMELY LONG TIME!!
	//variable wCh1Min=wavemin(Ch1W)
	//variable wCh2Min=wavemin(Ch2W)	
	//if(wCh1Min<0)
	//Ch1W+=(wavemin(Ch1W)+1) 
	//endif
	//if(wCh2Min<0)
	//Ch2W+=(wavemin(Ch2W)+1) 
	//endif
	///////////////////////////////////////////////////////////////////////////////////////////////
	
		
	variable dimX=dimsize(Ch2W,0)
	variable dimY=dimsize(Ch2W,1)
	variable dimZConcat=dimsize(Ch2W,2)
	variable numPlanes=(dimZConcat/numVols)

	string shiftXName
	string shiftYName
	
	make/o/n=(dimX,dimY,numPlanes) VolumeCh1=0	, VolumeStDevCh1=0
	make/o/n=(dimX,dimY,numPlanes) VolumeCh2=0, VolumeStDevCh2=0 // make matrices to accept the final volumes

	make/FREE/o/n=(dimX,dimY,numVols) PlaneTempCh1=0 // make matrices to accept each plane
	make/FREE/o/n=(dimX,dimY,numVols) PlaneTempCh2=0

	variable i,j,k, imcount, centralPlane
	if(regCh==2)
		for(i=0;i<numPlanes;i+=1) // cycle through the planes
			imcount=0
			for(j=0;j<dimZConcat;j+=numPlanes) // cycle through the frames of each plane
		
				if(imcount==(floor(numVols/2))) // use middle time frame to align rest of the plane
			
					duplicate/FREE/o/r=[][][(j+i)] Ch2W tempRefCh2
					redimension/s/n=(-1,-1) tempRefCh2
			
				endif
			
				duplicate/FREE/o/r=[][][(j+i)] Ch1W tempCh1
				duplicate/FREE/o/r=[][][(j+i)] Ch2W tempCh2
	
				redimension/s/n=(-1,-1) tempCh1,tempCh2
				setscale/p x,0,1, tempCh1,tempCh2
			
		
				PlaneTempCh1[][][imcount]=tempCh1[p][q]
				PlaneTempCh2[][][imcount]=tempCh2[p][q]
			
				imcount+=1
			
			endfor
			
			// Filter (if applicable) 
			if(filter!=0)
				MatrixFilter/o/N=(filter) gauss PlaneTempCh2
			endif
		
			// Register Ch2 of this plane 			
			Imageregistration/q /stck /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=PlaneTempCh2, refWave=tempRefCh2 // register volume with the middle volume in Ch2
			wave M_RegOut, M_RegParams
			MatrixOP/o/free temp = ReplaceNaNs(M_Regout, 0)	
			copyscaling(PlaneTempCh2, temp)
			killwaves M_RegOut		
		
			// Average Project this plane and insert into final Volume
			Imagetransform averageImage temp
			wave M_AveImage, M_StdvImage
			duplicate/o/free M_AveImage tempAv
			duplicate/o/free M_StdvImage tempRtStD
						
			// construct the Reference volume and standard deviation volumes
			VolumeStDevCh2[][][i]=tempRtStD[p][q]
			VolumeCh2[][][i]=tempAv[p][q]		
			killwaves tempRtStD

			// apply transformations to Ch1 for this plane
			Imageregistration/q /stck /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=M_RegParams testwave=PlaneTempCh1 // apply transformation to Ch1
			wave M_RegOut
			MatrixOP/o/free temp = ReplaceNaNs(M_Regout, 0)	//replace NaN's with -50000 to make them easy to spot when i miss them cropping
			copyscaling(PlaneTempCh1, temp)
		
		
			// Filter (if applicable) then average Project this plane and insert into final Volume
			if(filter!=0)
				MatrixFilter/o/N=(filter) gauss temp 
			endif
		
		
			Imagetransform averageImage temp
			wave M_AveImage, M_StdvImage
			duplicate/o M_AveImage tempAv,tempSubA, tempSubB
			duplicate/o/free M_StdvImage, tempRtStD
			
			// construct the Reference volume and standard deviation volumes
			VolumeStDevCh1[][][i]=tempRtStD[p][q]
			VolumeCh1[][][i]=tempAv[p][q]		
			killwaves tempRtStD
								
			// take the shifts and keep the highest ones		
			duplicate/FREE/o/r=[0][]	M_regParams tempshiftX
			duplicate/FREE/o/r=[1][]	M_regParams tempshiftY	
			setscale/p x,0,1, tempshiftX, tempshiftY 
			setscale/p y,0,1, tempshiftX, tempshiftY
				
			templowestX=wavemin(tempshiftX)
			templowestY=wavemin(tempshiftY)
			temphighestX=wavemax(tempshiftX)
			temphighestY=wavemin(tempshiftY)
		
			lowestX=lowestX<templowestX ? lowestX : templowestX
			lowestY=lowestY<templowestY ? lowestY : templowestY
			highestX=highestX>temphighestX ? highestX : temphighestX
			highestY=highestY>temphighestY ? highestY : temphighestY
		
			// Housekeeping
			wave M_AveImage,M_RegMaskOut,M_RegOut, M_StdvImage, M_RegParams
			killwaves M_AveImage, M_AveImage,M_RegMaskOut,M_RegOut, M_StdvImage, M_RegParams
		endfor
	
		// Because we align the time series to itself (ignoring z movement at this stage), the alignment of the z-stack and the time series are not neccessarily like for like. 
		// specifically, translations in z can also result in translations in x and y. These are 'corrected' by the imageregistration steps in the time series, but not the reference stack. 
		// This can result in a reference stack that suggests a fall/rise in fluorescence much steeper than the equivalent z-translation actually causes in the time series. 
		// To account for this, we register each plane of the reference stack (Ch2 or 1 depending on regCh parameter above) to it's adjacent plane, starting from the central plane and working our way out. 
		centralPlane=floor(numPlanes/2)
		make/o/n=(dimsize(VolumeCh2,0),dimsize(VolumeCh2,1),dimsize(VolumeCh2,2)), VolumeCh2_Plreg, VolumeCh1_Plreg, VolumeStDevCh2_Plreg, VolumeStDevCh1_Plreg

	
		for(i=0;i<floor(numPlanes/2);i+=1) 
			//grab image i+1 before and after the central plane 
			duplicate/free/o/r=[][][centralPlane+i+1] VolumeCh2 tempCh2_rhs
			duplicate/free/o/r=[][][centralPlane-i-1] VolumeCh2 tempCh2_lhs
			duplicate/free/o/r=[][][centralPlane+i+1] VolumeCh1 tempCh1_rhs
			duplicate/free/o/r=[][][centralPlane-i-1] VolumeCh1 tempCh1_lhs
			duplicate/free/o/r=[][][centralPlane+i+1] VolumeStDevCh2 tempStDevCh2_rhs
			duplicate/free/o/r=[][][centralPlane-i-1] VolumeStDevCh2 tempStDevCh2_lhs
			duplicate/free/o/r=[][][centralPlane+i+1] VolumeStDevCh1 tempStDevCh1_rhs
			duplicate/free/o/r=[][][centralPlane-i-1] VolumeStDevCh1 tempStDevCh1_lhs
	
			//grab reference images for both rhs and lhs. If i==0 then the central plane is corect for both
			duplicate/free/o/r=[][][centralPlane+i] VolumeCh2 tempCh2Ref_rhs
			duplicate/free/o/r=[][][centralPlane-i] VolumeCh2 tempCh2Ref_lhs
	
			redimension/s /n=(-1,-1) tempCh2_rhs,tempCh2_lhs,tempCh2Ref_rhs,tempCh2Ref_lhs,tempCh1_rhs,tempCh1_lhs,tempStDevCh2_rhs, tempStDevCh2_lhs, tempStDevCh1_rhs, tempStDevCh1_lhs
			// rhs first
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=tempCh2_rhs refwave=tempCh2Ref_rhs
			wave W_regParams, M_RegOut
			duplicate/o W_RegParams rhs_pwave	
			VolumeCh2_Plreg[][][centralPlane+1+i]=M_Regout[p][q]
			killwaves W_regParams, M_RegOut
	
			// apply shifts to Ch1 and StDev Images
			// Ch1
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=rhs_pwave testwave=tempCh1_rhs 
			wave M_RegOut
			VolumeCh1_Plreg[][][centralPlane+1+i]=M_Regout[p][q]
	
			// stdev Ch1
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=rhs_pwave testwave=tempStDevCh1_rhs 
			wave M_RegOut
			VolumeStDevCh1_Plreg[][][centralPlane+1+i]=M_Regout[p][q]
	
			// stdev Ch2
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=rhs_pwave testwave=tempStDevCh2_rhs 
			wave M_RegOut
			VolumeStDevCh2_Plreg[][][centralPlane+1+i]=M_Regout[p][q]
			
			// now lhs
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=tempCh2_lhs refwave=tempCh2Ref_lhs
			wave W_regParams, M_RegOut
			duplicate/o W_RegParams lhs_pwave	
			VolumeCh2_Plreg[][][centralPlane-1-i]=M_Regout[p][q]
			killwaves W_regParams, M_RegOut
	
			// apply shifts to Ch1 and StDev Images
			// Ch1
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=lhs_pwave testwave=tempCh1_lhs 
			wave M_RegOut
			VolumeCh1_Plreg[][][centralPlane-1-i]=M_Regout[p][q]
	
			// stdev Ch1
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=lhs_pwave testwave=tempStDevCh1_lhs 
			wave M_RegOut
			VolumeStDevCh1_Plreg[][][centralPlane-1-i]=M_Regout[p][q]
	
			// stdev Ch2
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=lhs_pwave testwave=tempStDevCh2_lhs 
			wave M_RegOut
			VolumeStDevCh2_Plreg[][][centralPlane-1-i]=M_Regout[p][q]
	
		endfor
		
		// fill in the central (unshifted reference) frame
		VolumeCh1_Plreg[][][centralPlane]=VolumeCh1[p][q][centralPlane]
		VolumeCh2_Plreg[][][centralPlane]=VolumeCh2[p][q][centralPlane]
		VolumeStDevCh1_Plreg[][][centralPlane]=VolumeStDevCh1[p][q][centralPlane]
		VolumeStDevCh2_Plreg[][][centralPlane]=VolumeStDevCh2[p][q][centralPlane]

	elseif(regCh==1)	// if regCh==1, then do the registration exactly as above but using channel 1 instead of channel 2 to register
		centralPlane=floor(numPlanes/2)
		for(i=0;i<numPlanes;i+=1) // cycle through the planes
			imcount=0
			for(j=0;j<dimZConcat;j+=numPlanes) // cycle through the frames of each plane
		
				if(imcount==(floor(numVols/2))) // use middle time frame to align rest of the plane
			
					duplicate/FREE/o/r=[][][(j+i)] Ch1W tempRefCh1
					redimension/s/n=(-1,-1) tempRefCh1
			
				endif
			
				duplicate/FREE/o/r=[][][(j+i)] Ch1W tempCh1
				duplicate/FREE/o/r=[][][(j+i)] Ch2W tempCh2
	
				redimension/s/n=(-1,-1) tempCh1,tempCh2
				setscale/p x,0,1, tempCh1,tempCh2
					
				PlaneTempCh1[][][imcount]=tempCh1[p][q]
				PlaneTempCh2[][][imcount]=tempCh2[p][q]
			
				imcount+=1
			
			endfor
			
			// Register Ch1 of this plane 			
			Imageregistration/q /stck /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=PlaneTempCh1, refWave=tempRefCh1 // register volume with the middle volume in Ch2
			wave M_RegOut, M_RegParams
			MatrixOP/o/free temp = ReplaceNaNs(M_Regout, 0)	
			copyscaling(PlaneTempCh1, temp)
			killwaves M_RegOut		
				
			// Filter (if applicable) and average Project this plane and insert into final Volume
			if(filter!=0)
				MatrixFilter/o/N=(filter) gauss temp
			endif
		
			if(imcount>1)
				Imagetransform averageImage temp
				wave M_AveImage
				duplicate/o M_AveImage temp
			endif
			VolumeCh1[][][i]=temp[p][q]		

			// apply transformations to Ch2 for this plane
			Imageregistration/q /stck /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=M_RegParams testwave=PlaneTempCh2 // apply transformation to Ch1
			wave M_RegOut
			MatrixOP/o/free temp = ReplaceNaNs(M_Regout, 0)	//replace NaN's with -50000 to make them easy to spot when i miss them cropping
			copyscaling(PlaneTempCh2, temp)
				
			// Filter (if applicable) then average Project this plane and insert into final Volume
			if(filter!=0)
				MatrixFilter/o/N=(filter) gauss temp 
			endif
		
			if(imcount>1)
				Imagetransform averageImage temp
				wave M_AveImage
				duplicate/o M_AveImage temp
			endif
			VolumeCh2[][][i]=temp[p][q]	
					
			// take the shifts and keep the highest ones		
			duplicate/FREE/o/r=[0][]	M_regParams tempshiftX
			duplicate/FREE/o/r=[1][]	M_regParams tempshiftY	
			setscale/p x,0,1, tempshiftX, tempshiftY 
			setscale/p y,0,1, tempshiftX, tempshiftY
				
			templowestX=wavemin(tempshiftX)
			templowestY=wavemin(tempshiftY)
			temphighestX=wavemax(tempshiftX)
			temphighestY=wavemin(tempshiftY)
		
			lowestX=lowestX<templowestX ? lowestX : templowestX
			lowestY=lowestY<templowestY ? lowestY : templowestY
			highestX=highestX>temphighestX ? highestX : temphighestX
			highestY=highestY>temphighestY ? highestY : temphighestY
		
			// Housekeeping
			wave M_AveImage,M_RegMaskOut,M_RegOut, M_StdvImage, M_RegParams
			killwaves M_AveImage,M_RegMaskOut,M_RegOut, M_StdvImage, M_RegParams
		endfor
		
		duplicate/o VolumeCh1 VolumeCh1_Plreg
		duplicate/o VolumeCh2 VolumeCh2_Plreg
			
		// Because we align the time series to itself (ignoring z movement at this stage), the alignment of the z-stack and the time series are not neccessarily like for like. 
		// specifically, translations in z can also result in translations in x and y. These are 'corrected' by the imageregistration steps in the time series, but not the reference stack. 
		// This can result in a reference stack that suggests a fall/rise in fluorescence much steeper than the equivalent z-translation actually causes in the time series. 
		// To account for this, we register each plane of the reference stack (Ch2 or 1 depending on regCh parameter above) to it's adjacent plane, starting from the central plane and working our way out. 

		for(i=0;i<floor(numPlanes/2);i+=1)
			//grab image i+1 before and after the central plane
			duplicate/free/o/r=[][][centralPlane+i+1] VolumeCh2_Plreg tempCh2_rhs
			duplicate/free/o/r=[][][centralPlane-i-1] VolumeCh2_Plreg tempCh2_lhs
			duplicate/free/o/r=[][][centralPlane+i+1] VolumeCh1_Plreg tempCh1_rhs
			duplicate/free/o/r=[][][centralPlane-i-1] VolumeCh1_Plreg tempCh1_lhs
	
			//grab reference images for both rhs and lhs. If i==0 then the central plane is corect for both
			duplicate/free/o/r=[][][centralPlane+i] VolumeCh1_Plreg tempCh1Ref_rhs
			duplicate/free/o/r=[][][centralPlane-i] VolumeCh1_Plreg tempCh1Ref_lhs

			redimension/s /n=(-1,-1) tempCh2_rhs,tempCh2_lhs,tempCh1Ref_rhs,tempCh1Ref_lhs,tempCh1_rhs,tempCh1_lhs
			// rhs first
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=tempCh1_rhs refwave=tempCh1Ref_rhs
			wave W_regParams, M_RegOut
			duplicate/o W_RegParams rhs_pwave	
			VolumeCh1_Plreg[][][centralPlane+1+i]=M_Regout[p][q]
			killwaves W_regParams, M_RegOut
			// apply shifts to Ch1
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=rhs_pwave testwave=tempCh2_rhs 
			wave M_RegOut
			VolumeCh2_Plreg[][][centralPlane+1+i]=M_Regout[p][q]
		
			// now lhs
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} testwave=tempCh1_lhs refwave=tempCh1Ref_lhs
			wave W_regParams, M_RegOut
			duplicate/o W_RegParams lhs_pwave	
			VolumeCh1_Plreg[][][centralPlane-1-i]=M_Regout[p][q]
			killwaves W_regParams, M_RegOut
			// apply shifts to Ch2
			Imageregistration/q /csnr=0 /refm=0 /tstm=0 /rot={0,0,0} /user=lhs_pwave testwave=tempCh2_lhs 
			wave M_RegOut
			VolumeCh2_Plreg[][][centralPlane-1-i]=M_Regout[p][q]
	
		endfor
	else
		print "Choose channel 1 or 2 to register to"
	endif

	killwaves lhs_pwave, rhs_pwave, M_RegOut, M_RegMaskOut
	
	// Create volumes of reference stack + and - SD of each pixel.
	
	MatrixOP/o VolumePosSDCh1=VolumeCh1_Plreg+VolumeStDevCh1
	MatrixOP/o VolumeNegSDCh1=VolumeCh1_Plreg-VolumeStDevCh1
	MatrixOP/o VolumePosSDCh2=VolumeCh2_Plreg+VolumeStDevCh2
	MatrixOP/o VolumeNegSDCh2=VolumeCh2_Plreg-VolumeStDevCh2
	
end
//###################################################################################################################

//###################################################################################################################
// Thomas Ryan - 2019
// RegisterAll	:	automatically registers the time series images, and subsequently the reference stack is shifted to match the realigned time series (faster than realigning the entire time series (1000s of frames) to the reference stack)

// ImageInputs:		Ch1MovieName			=	string variable, Ch1 time series wave name	 
//						Ch2MovieName			=	string variable, Ch2 time series wave name
//						RefStack_Ch1			=	string variable, Ch1 reference stack wave name	
//						RefStack_Ch2			=	string variable, Ch2 reference stack wave name
//						VolCh1SD				=	string variable, Standard deviation volume produced in CreateRefStack scripts
//						Ch1PosSDMovieName	=	string variable, Average stacks plus S.D. produced in CreateRefStack scripts
//						Ch1NegSDMovieName	=	string variable, Average stacks minus S.D. produced in CreateRefStack scripts


//	VariableInputs:	ImPlane				=	numerical variable, estimate of the initial imaging plane (usually the central plane)

//	The following inputs are the largest positive and negative shifts in x and y dimensions over all registration proceudres

//						highestX		
//						highestY
//						lowestX
//						lowestY

// 	OUTPUTS:			Realigned versions of all ImageInputs suffixed with "_reg"
//						shiftX			=	Vector of x transforms from time series	
//						shiftY			=	Vector of y transforms from time series
//						ZshiftX		=	largest shift during registration of reference stack in x
//						ZshiftY		=	largest shift during registration of reference stack in y

function RegisterAll(ImPlane,Ch1MovieName, Ch2MovieName, RefStack_Ch1, RefStack_Ch2, VolCh1SD,Ch1PosSDMovieName, Ch1NegSDMovieName)
	
	variable ImPlane
	string Ch1MovieName, Ch2MovieName, RefStack_Ch1, RefStack_Ch2,VolCh1SD,Ch1PosSDMovieName, Ch1NegSDMovieName
	
	variable i, j, dimX, dimY, dimZ, nFramesStack
	
	//////// Change this to alter the SD of the gaussian filter applied to images for registration ////////////
	variable filter=0
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	dimX=dimsize($Ch1MovieName,0)
	dimY=dimsize($Ch1MovieName,1)
	dimZ=dimsize($Ch1MovieName,2)
	nFramesStack=dimsize($RefStack_Ch1,2)
	
	// Register the Anatomical Channel to Self ////////////
	Make/o/n=(dimX,dimY,dimZ), Ch1MovieOut, Ch2MovieOut //create empty matrices to store output
	make/FREE/o/n=(dimX,dimY,dimZ), Ch2tempReg, Ch2tempFiltReg //create empy matrices for the filtered and toRegister versions of Ch2
	duplicate/free/o $Ch1MovieName, Ch1tempReg //create toRegister version of Ch1
	duplicate/free/o $Ch1PosSDMovieName, Ch1PosSD
	duplicate/free/o $Ch1NegSDMovieName, Ch1NegSD
	duplicate/free/o $VolCh1SD, VolCh1SDw
	
	redimension/s Ch1tempReg, VolCh1SDw, Ch1PosSD, Ch1NegSD, Ch2tempReg, Ch2tempFiltReg //redimension to single float for imageRegistration
	
	// create filtered image to perform initial registration on
	for(i=0;i<dimZ;i+=1) 
		duplicate/FREE/o/r=[][][i] $Ch2MovieName tempst //pull each slice of Ch2 movie one by one
		redimension/s /n=(-1,-1) tempst
		Ch2tempReg[][][i]=tempst[p][q] //use unadulterated image to populate the toRegister Ch2
		MatrixFilter/n=(filter) gauss,tempst // filter frame with a 2D gaussian kernel with SD of "filter" pixels
		Ch2tempFiltReg[][][i]=tempst[p][q] // populate filtered image with filtered frame
	endfor
	killwaves tempst
			
	// Currently registers to the first frame, uncomment below lines to average first 5 frames of the movie to use as a reference 			
		
	//make/o/FREE /n=(dimX,dimY,5) reftoA 
						
	i=0//	for(i=0;i<5;i+=1)								
	duplicate/FREE/o/r=[][][i] Ch2tempFiltReg, Ch2ref//tt
	//		redimension/s /n=(-1,-1) tt
	//		reftoA[][][i]=tt[p][q]
	//	endfor
	//	ImageTransform averageImage reftoA
	//	wave m_AveImage
	//	duplicate/FREE/o m_AveImage Ch2ref
	redimension/s/n=(-1,-1) Ch2Ref
	//	killwaves m_AveImage
		
	//register the filtered movie with the average of first 5 frames in ch2 (or simply the first frame)
		
	imageregistration /q /rot={0,0,0} /stck /pstk /csnr=0 /refm=0 /tstm=0 testwave=Ch2tempfiltReg, refwave=Ch2ref 
	wave m_regParams 
	
	make/o /n=(dimZ) shiftXY, shiftX, shiftY
	duplicate/o/r=[0][]	m_regParams shiftX, absshiftX
	duplicate/o/r=[1][]	m_regParams shiftY, absshiftY
	absshiftX=abs(shiftX)
	absshiftY=abs(shiftY)
	shiftXY[]=absshiftX[p]+absshiftY[p]
	
	
	duplicate/FREE/o m_regParams pwave //grab shifts
	killwaves m_regParams, Ch2tempfiltReg
	
	// apply these shifts to the unfiltered version of Ch2
	imageregistration /q /rot={0,0,0} /stck /pstk /csnr=0 /refm=0 /tstm=0 /user=pwave testwave=Ch2tempReg
	wave m_regout
	duplicate/o m_regout Ch2MovieOut // copy registration result for Ch2

	//apply the same shifts to the unfiltered Ch1
	imageregistration /q /rot={0,0,0} /stck /pstk /csnr=0 /refm=0 /tstm=0 /user=pwave testwave=Ch1tempReg //apply previous registration results to the other channel
	duplicate/o m_regout Ch1MovieOut // copy registration result for Ch1
	
	killwaves m_regout, pwave
	
	//////////// While we have these waves handy, register the zstack to the registered Ch2 ///////////

	make/free/o/n=(dimX,dimY) tempStack
	duplicate/FREE/o/r=[][][ImPlane] $RefStack_Ch2 tempStack // take central plane of volume to register against
	redimension/s/n=(-1,-1) tempStack
	MatrixOP/FREE/o stack_NaNBusted = ReplaceNaNs(tempStack, 0)
	MatrixFilter/n=(filter) gauss,stack_NaNBusted

	//register slice of the zstack (corresponding to provided ImPlane) with the same referece image as the time series
	imageregistration/q /rot={0,0,0} /csnr=0 /refm=0 /tstm=0 testwave=stack_NaNBusted, refwave=Ch2ref 
	wave W_regParams
	duplicate/FREE/o W_regParams pwaveRef
	killwaves W_regParams
	
	variable/g ZshiftX, ZshiftY
	ZshiftX=pwaveRef[0][0]
	ZshiftY=pwaveRef[1][0]
		
	//apply these shifts to both channels of the original zstacks and the SD Maps
	make/free/o/n=(dimsize($RefStack_Ch2,0),dimsize($RefStack_Ch2,1),dimsize($RefStack_Ch2,2)), RefOutCh1_Collected, RefOutCh2_Collected, RefOutPosSDCh1_Collected,RefOutNegSDCh1_Collected, SDOutCh1_Collected // Create matrices to accept each slice (Probably not needed)
	
	for (i=0;i<dimsize($RefStack_Ch2,2);i+=1)
		duplicate/FREE/o /r=[][][i] $RefStack_Ch2 StackReg_Ch2
		duplicate/FREE/o /r=[][][i] $RefStack_Ch1 StackReg_Ch1
		duplicate/FREE/o /r=[][][i] Ch1PosSD StackRegPosSD_Ch1
		duplicate/FREE/o /r=[][][i] Ch1NegSD StackRegNegSD_Ch1
		duplicate/FREE/o /r=[][][i] VolCh1SDw StackSD_Ch1
		redimension /s /n=(-1,-1) StackReg_Ch1, StackReg_Ch2,StackRegPosSD_Ch1, StackRegNegSD_Ch1, StackSD_Ch1
		
		MatrixOP/FREE/o stackCh1_NaNBusted = ReplaceNaNs(StackReg_Ch1, 0)
		MatrixOP/FREE/o stackCh2_NaNBusted = ReplaceNaNs(StackReg_Ch2, 0)
		MatrixOP/FREE/o stackPosSDCh1_NaNBusted = ReplaceNaNs(StackRegPosSD_Ch1, 0)
		MatrixOP/FREE/o stackNegSDCh1_NaNBusted = ReplaceNaNs(StackRegNegSD_Ch1, 0)
		MatrixOP/FREE/o StackSD_Ch1_NaNBusted = ReplaceNaNs(StackSD_Ch1, 0)
		
		imageregistration /stck /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=stackCh2_NaNBusted // apply shifts to rest of the z-stack in Ch2 (assumes rigid shift through whole volume)
		wave M_RegOut
		string RefStackOutName_Ch2=RefStack_Ch2+"_reg"
		duplicate/FREE/o M_RegOut RefOutCh2

		imageregistration /stck /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=stackCh1_NaNBusted // apply shifts to rest of the z-stack in Ch1 (assumes rigid shift through whole volume) 
		string RefStackOutName_Ch1=RefStack_Ch1+"_reg"		
		duplicate/FREE/o M_RegOut RefOutCh1
		
		imageregistration /stck /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=StackSD_Ch1_NaNBusted // apply shifts to rest of the z-stack in SD (assumes rigid shift through whole volume) 
		string SDStackOutName_Ch1=VolCh1SD+"_reg"		
		duplicate/FREE/o M_RegOut SDOutCh1
		
		imageregistration /stck /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=stackPosSDCh1_NaNBusted // apply shifts to rest of the z-stack in PosSD (assumes rigid shift through whole volume) 
		string PosSDOutName_Ch1=Ch1PosSDMovieName+"_reg"		
		duplicate/FREE/o M_RegOut RefOutPosSDCh1
		
		imageregistration /stck /rot={0,0,0} /q /csnr=0 /refm=0 /tstm=0 /user=pwaveRef testwave=stackNegSDCh1_NaNBusted // apply shifts to rest of the z-stack in PosSD (assumes rigid shift through whole volume) 
		string NegSDOutName_Ch1=Ch1NegSDMovieName+"_reg"		
		duplicate/FREE/o M_RegOut RefOutNegSDCh1
		
		MatrixOP/FREE/o RefOutCh1_NaNBusted = ReplaceNaNs(RefOutCh1, 0)
		MatrixOP/FREE/o RefOutCh2_NaNBusted = ReplaceNaNs(RefOutCh2, 0)	
		MatrixOP/o/free RefOutPosSDCh1_NaNBusted = ReplaceNaNs(RefOutPosSDCh1,0)
		MatrixOP/o/free RefOutNegSDCh1_NaNBusted = ReplaceNaNs(RefOutNegSDCh1,0)
		MatrixOP/FREE/o SDOutCh1_NaNBusted = ReplaceNaNs(SDOutCh1, 0)

		RefOutCh1_Collected[][][i]=RefOutCh1_NaNBusted[p][q]
		RefOutCh2_Collected[][][i]=RefOutCh2_NaNBusted[p][q]
		RefOutPosSDCh1_Collected[][][i]=RefOutPosSDCh1_NaNBusted[p][q]
		RefOutNegSDCh1_Collected[][][i]=RefOutNegSDCh1_NaNBusted[p][q]
		SDOutCh1_Collected[][][i]=SDOutCh1_NaNBusted[p][q]
	endfor

	duplicate/o RefOutCh1_Collected $RefStackOutName_Ch1
	duplicate/o RefOutCh2_Collected $RefStackOutName_Ch2
	duplicate/o RefOutPosSDCh1_Collected $PosSDOutName_Ch1
	duplicate/o RefOutNegSDCh1_Collected $NegSDOutName_Ch1
	duplicate/o SDOutCh1_Collected $SDStackOutName_Ch1
	
	///////////////// Zap some NaNs, scale and name waves appropriately //////////////////////
	string Ch2MovieOutName = Ch2MovieName + "_reg"
	string Ch1MovieOutName = Ch1MovieName + "_reg"
	MatrixOP/FREE/o Ch1w_NaNBusted = ReplaceNaNs(Ch1MovieOut, 0)	
	MatrixOP/FREE/o Ch2w_NaNBusted = ReplaceNaNs(Ch2MovieOut, 0)	
	copyscaling(Ch1tempReg, Ch1w_NaNBusted)
	copyscaling(Ch2tempReg, Ch2w_NaNBusted)
	duplicate /o Ch1w_NaNBusted, $Ch1MovieOutName
	duplicate /o Ch2w_NaNBusted, $Ch2MovieOutName

	//Housekeeping
	wave M_RegMaskOut
	killwaves M_RegOut, M_RegMaskOut, Ch2tempReg, Ch1tempReg
end
//#######################################################################################################


// Thomas Ryan	-	2019
//	CropScalarShifts	:	following registration, this script takes all time series and reference stacks and crops them by the single largest shift in each dimension from the realignment steps applied to reference stack images and time series. 
// ImageInputs are generated by the RegisterAll procedure

//	ImageInputs:		RegCh1Movie	=	time series Ch1 wave name
//						RegCh2Movie	=	time series Ch2 wave name
//						RegCh1Ref		=	reference stack Ch1 wave name
//						RegCh2Ref		=	reference stack Ch2 wave name
//						RegPosSDRef	=	reference stack plus S.D.
//						RegNegSDRef	=	reference stack minus S.D.
//						RegVolSDRef	=	reference stack S.D. projection

//	VariableInputs:	shiftX			=	Vector of x transforms from time series	
//						shiftY			=	Vector of y transforms from time series	
//						ZshiftX		=	largest shift during registration of reference stack in x
//						ZshiftY		=	largest shift during registration of reference stack in y

//		The following inputs are the largest positive and negative shifts in x and y dimensions over all registration proceudres

//						highestX
//						highestY	
//						lowestX
//						lowestY


function CropScalarShifts(RegCh1Movie,RegCh2Movie,RegCh1Ref,RegCh2Ref,RegPosSDRef,RegNegSDRef,RegVolSDRef,shiftX,shiftY,ZshiftX,ZshiftY,highestX,highestY,lowestX,lowestY)

	string RegCh1Movie, RegCh2Movie, RegCh1Ref, RegCh2Ref,RegPosSDRef,RegNegSDRef,RegVolSDRef
	wave shiftX, shiftY
	variable highestX,highestY,lowestX,lowestY, ZshiftX, ZshiftY
	
	// find biggest shifts in movie
	variable highestMovX, lowestMovX, highestMovY, lowestMovY
	highestMovX=wavemax(shiftX)
	lowestMovX=wavemin(shiftX)
	highestMovY=wavemax(shiftY)
	lowestMovY=wavemin(shiftY)

	highestX=highestX>highestMovX ? highestX : highestMovX
	lowestX=lowestX<lowestMovX ? lowestX : lowestMovX
	highestY=highestY>highestMovY ? highestY : highestMovY
	lowestY=lowestY<lowestMovY ? lowestY : lowestMovY

	variable remRowsStart, remRowsEnd, remColsStart, remColsEnd

	remRowsStart=highestX>ZshiftX ? highestX : ZshiftX
	remRowsEnd=lowestX<ZshiftX ? lowestX : ZshiftX
	remColsStart=highestY>ZshiftY ? highestY : ZshiftY
	remColsEnd=lowestY<ZshiftY ? lowestY : ZshiftY

	remRowsStart=ceil(remRowsStart)
	remRowsEnd=floor(remRowsEnd)
	remColsStart=ceil(remColsStart)
	remColsEnd=floor(remColsEnd)
	
	remRowsEnd*=-1
	remColsEnd*=-1
	
	variable dimRows=dimsize($RegCh1Movie,0)
	variable dimCols=dimsize($RegCh1Movie,1)

	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegCh1Movie, CroppedCh1
	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegCh2Movie, CroppedCh2
	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegCh1Ref, CroppedRefCh1
	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegCh2Ref, CroppedRefCh2
	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegPosSDRef, CroppedPosSD
	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegNegSDRef, CroppedNegSD
	duplicate/o/rmd=[remRowsStart+1,(dimRows-(remRowsEnd+1))][remColsStart+1,(dimCols-(remColsEnd+1))][] $RegVolSDRef, CroppedSDRefCh1
	setscale/p x,0,1, CroppedCh1, CroppedCh2, CroppedRefCh1, CroppedRefCh2, CroppedPosSD, CroppedNegSD, CroppedSDRefCh1
	setscale/p y,0,1, CroppedCh1, CroppedCh2, CroppedRefCh1, CroppedRefCh2, CroppedPosSD, CroppedNegSD, CroppedSDRefCh1

end



//  Leon Lagnado, 25/1/2019


//Input is the movie and the reference stack (both after x-y alignment and with the reference stack 
//appropriately averaged and pre-processed).

//Output is a wave called  Zestimate which shows estimated z displacement per point.

//Note the variables that you should check: zstep, downsamplez (= 1), downsamplexy (= 2).  
//For beads do not downsample.  For vessels, downsampling xy by factor 2 has no effect on Zestimate.  
//Even downsampling by 4 works (and saves a lot of time!).  The downsample variables should be 1 or a multiple of 2. 

//	Check for offsets // DONE TR 27/01/19
// TR- various edits/bugs 01/19
// TR speed per frame - 28/01/19 - after optimizing
// - Sample		Downsamplexy	Dimensions		Laptop/Comp		Result 
// - Beads		1					256x100		 	Laptop			0.82s
// - Vessel		2					512x100			Laptop			0.69s
// TR 12/03/19 - added option to downsample in x and y dimensions seperately (images are typically 512 x 100)
// TR- 01/06/19modified to generate a correlation map over time. 
			
function Z_estimate(MovieStackIn, RefStackIn)

	wave MovieStackIn, RefStackIn
	
	// DOWNSAMPLING PARAMETERS
	variable deltat=dimdelta(MovieStackIn,0)
	variable i, j, k
	variable MovieFrame_sd
	variable zstep =0.5					//The distance between adjacent planes in the reference stack in microns
	variable downsamplez = 1
	variable downsamplex = 2
	variable downsampley = 1
	
	//Some downsampling can be introduced here to speed up cross-correlation, but be careful!
	//Down-sampling does not work well with images that have high spatial frequencies (such as beads)
	//DOES work well with images of vessels.  Downsampling x4 seems to work almost as well as straight 400x100 images.
	Variable zdim = DimSize(RefStackIn, 2)
	
	Duplicate/FREE/O MovieStackIn, MovieStack
	Resample/DIM=0/DOWN=(downsamplex) MovieStack
	Resample/DIM=1/DOWN=(downsampley) MovieStack
	Resample/DIM=2/DOWN=(downsamplez) MovieStack
	Duplicate/FREE/O/r=[][][1,(zdim-2)] RefStackIn, RefStack // edit - TR 29/01/19 - Remove first and last z-slice because of piFoc shearing artefact
	Resample/DIM=0/DOWN=(downsamplex) RefStack
	Resample/DIM=1/DOWN=(downsampley) RefStack
	Resample/DIM=2/DOWN=(downsamplez) RefStack
		
	zdim = DimSize(RefStack, 2)
	Variable nFrames = DimSize(MovieStack, 2)
	variable xdim=DimSize(RefStack, 0)
	Variable ydim = DimSize(RefStack, 1)
	
	variable startZ=(((floor(zdim/2))*zstep)*-1)
	make/o/N=(zdim) xwave = startZ+(p*zstep)
	Duplicate/FREE/O RefStack, RefStack_meansub

	Make/FREE /O/N=(xdim, ydim) MovieFrame, RefFrame
	Make/FREE /O/N=(zdim) Ref_avg, Ref_sd
	Make /O/N=(zdim) Corr_values
	Make /O/N=(nFrames) Zestimate
	Make /O/N=(200) fit_Corr_values
	Make/FREE /O/N=5 W_coef

	//for timing while getting code right and fast

	//  Make array containing sd for each frame of ref stack and make the mean subtracted ref stack
	for(i=0; i<zdim; i+=1)
		duplicate/FREE/o/r=[][][i] RefStack RefFrame
		redimension/n=(-1,-1) RefFrame
		//RefFrame[][]=RefStack[p][q][i] // duplicate is twice as fast than wave assignment: this for loop Wave Assignment: 0.085s, duplicating: 0.046s
		WaveStats/Q RefFrame
		RefStack_meansub[][][i]=RefFrame[p][q]-V_avg
		Ref_sd[i] = V_sdev
	endfor

	//To compensate for variations in intensity (e.g due to bleaching or changes in laser power) 
	//after cross-correlating two images, divide by the product of the sd of each image
	////https://en.wikipedia.org/wiki/Cross-correlation#Zero-normalized_cross-correlation_(ZNCC)
	make/o/n=(zdim,nFrames) CorrMap

	for(j=0; j<nFrames; j+=1)
		// edited TR 28/01/19///////////////////////////////
		duplicate/FREE/o/r=[][][j] MovieStack MovieFrame
		redimension/n=(-1,-1) MovieFrame // edit TR 28/01/19: duplicating and redimension twice as fast as wave assignment
		//MovieFrame[][]=MovieStack[p][q][j]
		WaveStats/Q MovieFrame			//using WaveSats because it is multithreaded
		MovieFrame-=V_avg					//subtract the mean
		MovieFrame_sd = V_sdev			//find sd of that movie frame
		for(i=0; i<zdim; i+=1)
			duplicate/FREE/o/r=[][][i] RefStack_meansub RefFrame
			redimension/n=(-1,-1) RefFrame
			MatrixOP/FREE/O MovieFrame_corr=Correlate(RefFrame,MovieFrame,4)
			variable dividor=Ref_sd[i]*MovieFrame_sd
			MovieFrame_corr/=dividor
			Corr_values[i]=wavemax(MovieFrame_corr)
		endfor
		/////////////////////////////////////////////////////
		setscale/i x,xwave[0],xwave[inf], Corr_values
		CorrMap[][j]=corr_values[p]
		setscale/i x,xwave[0],xwave[inf], CorrMap
		setscale/p y,0,deltat, CorrMap
		variable wmin=wavemin(Corr_values)
		Corr_values-=wmin // just for fitting
		CurveFit/Q gauss Corr_values /D
		WaveStats/Q fit_Corr_values
		Zestimate[j]=V_maxloc*-1
		doupdate
		
	endfor
	
	
	killwaves Corr_values,fit_Corr_values
	

end
//#################################################################################################

//########### Utility functions #################\\
			
// Utility function that pools waves of a given name from multiple experiment files and concatenates them into a target wave for bulk analysis
function PoolData()

	string pathstr="[pathToFolder]",subDataFolder="[subFolderNAme]",objList="waveName",targetStr="targetWaveName"	
			
	//	prompt pathstr, "type the full path of the folder containing your experiments"
	//	prompt subDataFolder, "Type the name of the subdatafolder (if any) your data is in"
	//	prompt objList, "Type the name of the wave you want to concatenate"
	//	prompt targetStr, "Type the name of the wave you want to concatenate into"
	//	doprompt "File selection", pathstr, subDataFolder, objList, targetStr
			
			
	newpath/O path, pathstr			
	string filelist= indexedfile(path,-1,".pxp")	
			
	Variable numFiles=itemsinlist(filelist)
		
	variable i,j
		
	for(i=0;i<numFiles;i+=1)	// cycle through files in the folder
		
		string expName= StringFromList(i, filelist)
		string message = "Loading from Experiment " + expName + " ( " + num2str(i+1) + "/" +num2str(numFiles) + " in this folder"
		print message
		
		string	expPath = pathstr + "\\" + expName
		
		// load objects in the object string list		
		LoadData/S=subDataFolder/J=objList   /Q expPath
		duplicate/free/o $objList objW
		if(i==0)
			duplicate/free objW targetW
			killwaves $objList, objW
		else
			concatenate/np/Kill {objW}, targetW
			killwaves $objList
		endif
	endfor
	duplicate/o targetW $targetstr
end

// Utility functions to attempt loading waves from files with given prefix (wildcard workaround)	
function LoadMatchesCh3 (prefixstr, sourceFile)
	string prefixstr, sourceFile
	variable i
	string resultstr = ""
    
	i = 0
	do
		if (tryload (prefixstr + num2str (i) + "_Ch3", sourceFile))
			resultstr = AddListItem (prefixstr + num2str (i), resultstr, ";", inf)
		else
			break
		endif
		i += 1
	while (1)
end

function LoadMatchesCh4 (prefixstr, sourceFile)
	string prefixstr, sourceFile
	variable i
	string resultstr = ""
    
	i = 0
	do
		if (tryload (prefixstr + num2str (i) + "_Ch4", sourceFile))
			resultstr = AddListItem (prefixstr + num2str (i), resultstr, ";", inf)
		else
			break
		endif
		i += 1
	while (1)
end

// Utiity function to search experiment files for a wave with a full or partial stringmatch
function TryLoad (WName, sourceFile)
	string WName, sourceFile
    
	if (strlen (sourceFile) == 0 || strlen (WName) == 0)
		return 1
	endif
	if (exists (WName))
		return 0
	endif
	loaddata /q/j=WName sourcefile
	return (!exists (WName))
end


// Utility function to kill all waves in the databrowser
function killAllWaves()

	string list=wavelist("*",";","")
	Variable numList = ItemsInList(list),j
			
	for (j=0;j<numList;j+=1)
		KillWaves/Z $(StringFromList(j, list))
	endfor
end

function killAllWavesList(List)
	string List
	Variable numList = ItemsInList(list),j
			
	for (j=0;j<numList;j+=1)
		KillWaves/Z $(StringFromList(j, list))
	endfor
	
end


// Utility Macro to load specific waves from a given experiment, and save under a different name. Can be used to perform computations and other functions using loaded waves
macro load&sav()

	string savename="[saveFilePath]"
	string expPath="[expFilePath]"
	string objlist="waveName1;waveName2;waveNameN"




	string namedWave1="nameOfWave1"
	string namedWave2="nameOfWave2"

	objlist = AddListItem(namedWave1,objList,";")
	objlist = AddListItem(namedWave2,objList,";")

	LoadData/J=objList   /Q expPath

	saveexperiment as saveName

endmacro