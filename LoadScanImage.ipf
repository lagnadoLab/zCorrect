#pragma rtGlobals=1		// Use modern global access method.

constant Z_factor = 4, ImageLength = 665.6

//Z_factor is the factor by which distances in the z dimension have to be multiplied in order to return actual distances in µm. 
//This is important only for the correct scaling of image stacks, not for single images or movies.
//ImageLength is the side length in µm of an image taken at zoom 1. This will, among other factors, depend on the objective used.

//For images saved with a different program, only basic functionality (i.e. loading of image stacks or file sequences) is available.

//02/02/2010 MMD: added version support for msPerLineFromHeader
//26/04/2010 MMD: added AverageFrames(imStack, n) and AutoAverageFrames(imStack)
//30/6/2018 LL:modified several functions to work with ScanImage 5 using Igor 8.  In particular, the new big tiff format for Scanimage files.
//13/02/2019 TR:added new function AutoLoadScanImageSingle to load images as single precision float for opimal memory usage in specific applications (e.g. z-correct)
//////////////////////////////////////////////////


Function /t LoadScanImage()

String ImgWaveName, FirstWave
string header, s_info = "No header info available\r"
Variable PointPos

ImageLoad/Q/O/C=-1/LR3D/BIGT=1/LTMD

Wave/t tagWave = root:Tag0:T_Tags
String s1=tagWave[12]		
Variable len=strlen(s1)
String s2=s1[18,len]			// to skip "imageDescription:"

if (v_flag == 0)
	return "-1"
endif

	//header = s_info
	header = s2
	PointPos = strsearch(S_Filename, ".tif", 0)
	ImgWaveName = S_FileName[0,PointPos-1]
	ImgWaveName = ReplaceString("-", ImgWaveName, "_")
	
	PointPos = strsearch(S_Wavenames, ";", 0)
	FirstWave =S_Wavenames[0,PointPos-1]
	
if (waveexists($ImgWaveName))
	killwaves /z $ImgWaveName
endif


duplicate /o $FirstWave, $ImgWaveName
Killwaves /z $FirstWave


Note $ImgWaveName, s2
Note $ImgWaveName, "file.path="+s_path
Note $ImgWaveName, "file.name="+s_filename

redimension /s $ImgWaveName		//convert to single precision floating point
									//Comment: large resolution files may exceed the system's memory when converted to double-precision FP

Return ImgWaveName
End

//////////////////////////////////////////////////////

Function /wave AutoLoadScanImage(pathstr, filenamestr)			//Does not open a file dialogue, takes path and filename as input parameters
	String pathstr, filenamestr

	String ImgWaveName, FirstWave
	string header, s_info = "No header info available\r"
	Variable PointPos
	
	NewPath/o/q path, pathstr

	ImageLoad/Q/O/C=-1/LR3D/BIGT=1/LTMD/P=path filenamestr
	
	Wave/t tagWave = root:Tag0:T_Tags
	String s1=tagWave[12]		
	Variable len=strlen(s1)
	String s2=s1[18,len]			// to skip "imageDescription:"
	
	if (v_flag == 0)
		Abort
	endif
	
		header = s_info
		PointPos = strsearch(S_Filename, ".tif", 0)
		ImgWaveName = S_FileName[0,PointPos-1]
		ImgWaveName = ReplaceString("-", ImgWaveName, "_")
		
		PointPos = strsearch(S_Wavenames, ";", 0)
		FirstWave =S_Wavenames[0,PointPos-1]
		
	if (waveexists($ImgWaveName))
		killwaves /z $ImgWaveName
	endif
	
	
	duplicate /o $FirstWave, $ImgWaveName
	Killwaves /z $FirstWave
	
	redimension /d $ImgWaveName		//convert to double precision floating point
	
	Note $ImgWaveName, s2
	Note $ImgWaveName, "file.path="+s_path
	Note $ImgWaveName, "file.name="+s_filename
	
	Wave ReturnWv =  $ImgWaveName
	Return ReturnWv
End

///////////////////////////// Added TR 13/02/19 for memory optimisation /////////////////////////////
// Conversion to double precision massively expands the memory required to load images but is unneccessary for many applications

Function /wave AutoLoadScanImageSingle(pathstr, filenamestr)			//Does not open a file dialogue, takes path and filename as input parameters
	String pathstr, filenamestr

	String ImgWaveName, FirstWave
	string header, s_info = "No header info available\r"
	Variable PointPos
	
	NewPath/o/q path, pathstr

	ImageLoad/Q/O/C=-1/LR3D/BIGT=1/LTMD/P=path filenamestr
	
	Wave/t tagWave = root:Tag0:T_Tags
	String s1=tagWave[12]		
	Variable len=strlen(s1)
	String s2=s1[18,len]			// to skip "imageDescription:"
	
	if (v_flag == 0)
		Abort
	endif
	
		header = s_info
		PointPos = strsearch(S_Filename, ".tif", 0)
		ImgWaveName = S_FileName[0,PointPos-1]
		ImgWaveName = ReplaceString("-", ImgWaveName, "_")
		
		PointPos = strsearch(S_Wavenames, ";", 0)
		FirstWave =S_Wavenames[0,PointPos-1]
		
	if (waveexists($ImgWaveName))
		killwaves /z $ImgWaveName
	endif
	
	
	duplicate /o $FirstWave, $ImgWaveName
	Killwaves /z $FirstWave
	
	redimension/s /d $ImgWaveName		//convert to single precision floating point
	
	Note $ImgWaveName, s2
	Note $ImgWaveName, "file.path="+s_path
	Note $ImgWaveName, "file.name="+s_filename
	
	Wave ReturnWv =  $ImgWaveName
	Return ReturnWv
End

//////////////////////////////////////////////////////

Function ZoomFromHeader(PicWave)
wave  Picwave
variable ZoomFactor
string header = note(PicWave)

ZoomFactor = str2num(StringByKey("SI.hRoiManager.scanZoomFactor",TrimString(ReplaceString("\n",header,";"),1)," = "))

Return ZoomFactor
End

//////////////////////////////////////////////////////

Function FramesFromHeader(PicWave)
wave PicWave
variable Frames
string header = note(PicWave)

//Frames = NumberByKey("state.acq.numberOfFrames", Header, "=","\r")
Frames = str2num(StringByKey("SI.hStackManager.framesPerSlice",TrimString(ReplaceString("\n",header,";"),1)," = "))

return Frames
End

//////////////////////////////////////////////////////

Function XRelFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

String Pos=StringByKey("SI.hMotors.motorPosition",TrimString(ReplaceString("\n",header,";"),1)," = ")
Variable len=strlen(Pos) 
String xposstring = Pos[1, len]
Variable xrel = str2num(xposstring)
return xrel

End
//////////////////////////////////////////////////////

Function YRelFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

String Pos=StringByKey("SI.hMotors.motorPosition",TrimString(ReplaceString("\n",header,";"),1)," = ")
Variable len=strlen(Pos) 
String xposstring = Pos[1, len]
Variable xrel = str2num(xposstring)
String yposstring = Pos[(1+strlen(num2str(xrel))), len]
Variable yrel = str2num(yposstring)

return yrel

End


//////////////////////////////////////////////////////
Function /t FilePathFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

return StringByKey("file.path", Header, "=","\r")

End

//////////////////////////////////////////////////////
Function /t FilenameFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

return StringByKey("file.name", Header, "=","\r")

End

//////////////////////////////////////////////////////

Function ZRelFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

String Pos=StringByKey("SI.hMotors.motorPosition",TrimString(ReplaceString("\n",header,";"),1)," = ")
Variable len=strlen(Pos) 
String xposstring = Pos[1, len]
Variable xrel = str2num(xposstring)
String yposstring = Pos[(1+strlen(num2str(xrel))), len]
Variable yrel = str2num(yposstring)
String zposstring = Pos[(2+strlen(num2str(xrel))+strlen(num2str(yrel))), len]
Variable zrel = str2num(zposstring)

return zrel

End
//////////////////////////////////////////////////////

Function msPerLineFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

String sPerLineStr=StringByKey("SI.hRoiManager.linePeriod",TrimString(ReplaceString("\n",header,";"),1)," = ")
Variable msPerLine
print sPerLineStr
msPerLine = str2num(sPerLineStr)*1000

return msPerLine

End

//////////////////////////////////////////////////////

Function sPerLineFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

String sPerLineStr= StringByKey("SI.hRoiManager.linePeriod",TrimString(ReplaceString("\n",header,";"),1)," = ")
Variable sPerLine=str2num(sPerLineStr)

return sPerLine

End

//////////////////////////////////////////////////////

Function /t ExpDateFromHeader(PicWave)
wave PicWave
string header = note(PicWave)
string ExpDate, month, day, year, ExpTime
variable pointer

ExpDate = StringByKey("state.internal.triggerTimeString", Header, "=","\r")
if (stringmatch(expdate, ""))
	return ""
endif

ExpDate = ReplaceString("'", ExpDate, "")

if(stringmatch(ExpDate[1], "/"))
	Month = "0"+ExpDate[0]
	pointer = 1
Elseif(stringmatch(ExpDate[2], "/"))
	Month = ExpDate[0,1]
	pointer = 2
else
	DoAlert 0, "Something's wrong with the date..."
	return ""
endif

if(stringmatch(ExpDate[pointer+2], "/"))
	Day = "0"+ExpDate[pointer+1]
	pointer += 1
Elseif(stringmatch(ExpDate[pointer+3], "/"))
	Day = ExpDate[pointer+1,pointer+2]
	pointer += 2
else
	DoAlert 0, "Something's wrong with the day..."
	return ""
endif

year = ExpDate[pointer+2, pointer+5]

ExpTime=ExpDate[pointer+7, pointer+14]

string result = Day+"/"+Month+"/"+Year+" "+ExpTime

Return Result
End

//////////////////////////////////////////////////////

Function ZSlicesFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

Variable ZSlices= str2num(StringByKey("SI.hStackManager.numSlices", TrimString(ReplaceString("\n",header,";"),1)," = "))

Return ZSlices

End

//////////////////////////////////////////////////////

Function ZStepsizeFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

Variable ZStepSize= str2num(StringByKey("SI.hStackManager.stackZStepSize", TrimString(ReplaceString("\n",header,";"),1)," = "))

return ZStepSize

End

//////////////////////////////////////////////////////

Function ScanImageVersionFromHeader(PicWave)
wave PicWave
string header = note(PicWave)

Return NumberByKey("state.software.version", Header, "=","\r")

End

//////////////////////////////////////////////////////

function ApplyHeaderInfo(Wave3D)
	wave wave3D
	variable x_res,y_res,z_res, slices, stepsize
	
	variable Zoom=ZoomFromHeader(Wave3D)
	
	if(NumType(Zoom) == 2)			//break, if no ScanImage header info available
		return -1
	endif
	
	variable timeperline=sPerLineFromHeader(Wave3D), nChannels
	
	x_res = ImageLength/ zoom / dimsize(wave3d,0)*1e-6	//in m
	y_res = ImageLength / zoom / dimsize(wave3d,1)*1e-6	//in m
	z_res = timeperline *  dimsize(wave3d,1)	//in s
	slices = ZSlicesFromHeader(wave3d)
	stepsize= ZStepsizeFromHeader(wave3d) * z_factor*1e-6 //in m
	
	setscale /P x, 0, x_res,"m",wave3d
	setscale /P y, 0, y_res,"m",wave3d
	
	if(slices > 1)
		setscale /P z, 0, stepsize,"m",wave3d
	else
		setscale /P z, 0, z_res,"s",wave3d
	endif

end

////////////////////////////////////////////////

function /t LoadMovie()		//loads a ScanImage movie from multiple files

String ImgWaveName, FirstWave, FileName, FNTrunc, FNNum, FNPath, StrNum
	string  header, s_info = "No header info available\r"
	Variable PointPos, startnum, cont = 1, newframes,currentframes, frames
	
	ImageLoad /Q /O /C=-1
	
	if (v_flag == 0)
		return "-1"
	endif
	
		header = s_info
		PointPos = strsearch(S_Filename, ".tif", 0)
		ImgWaveName = S_FileName[0,PointPos-1]
		FileName = S_FileName[0,PointPos-1]
		ImgWaveName = ReplaceString("-", ImgWaveName, "_")
		PointPos = strsearch(S_Wavenames, ";", 0)
		FirstWave =S_Wavenames[0,PointPos-1]
		FNPath = S_Path
	
		StrNum=FileName[strlen(FileName)-3,strlen(FileName)-1]	//last 3 characters in FileName
		StartNum=Str2Num(StrNum)
		if(NumType(StartNum==2))
			DoAlert 0, "The last three characters of the filename are not numbers."
			return "-1"
		endif
	
		FNTrunc = FileName[0,strlen(FileName)-4]
		
	if (waveexists(:tw0))
		killwaves /z tw0
	endif
		
	duplicate /o/free $FirstWave, TempImage
	killwaves /z  $FirstWave
	
	
	variable refnum
	
	string fn2
	
	Do
		startnum +=1
		currentframes=dimsize(TempImage,2)
		
		if (startnum < 10)
			FNNum = "00"+Num2Str(startnum)
		elseif (startnum < 100)
			FNNum = "0"+Num2Str(startnum)
		else
			FNNum = Num2Str(startnum)
		endif
			
		FileName = 	FnPath+FNTrunc+FNNum+".tif"
		fn2 = 	FnPath+FNTrunc+FNNum
		open /z=1 /r  refnum as FileName
		
		if (v_flag==0)
			ImageLoad /c=-1 /n=tw /o /q FileName
			
			
			PointPos = strsearch(S_Wavenames, ";", 0)
			FirstWave =S_Wavenames[0,PointPos-1]
			Wave sec = $firstwave
			
			
			frames=dimsize(sec,2)
			if(frames==0)
				frames = 1
			endif
			
			newframes=currentframes+frames
			
			redimension/n=(-1,-1,newframes) TempImage
			
			TempImage[][][currentframes,newframes-1] = sec[p][q][r-currentframes]
			
									
			close refnum
			killwaves $firstwave
		else
			cont = 0
		endif
		
		
	While(cont)
	
	
	if (waveexists($ImgWaveName))
		killwaves /z $ImgWaveName
	endif
	
	
	duplicate /o TempImage, $ImgWaveName
	Killwaves /z $FirstWave, TempImage
	
	redimension /s $ImgWaveName		//convert to single precision floating point
	
	
	Note $ImgWaveName, header
	Note $ImgWaveName, "file.path="+s_path
	Note $ImgWaveName, "file.name="+s_filename
	
	Print "Saved as "+ ImgWaveName
	
	Return ImgWaveName
End


//////////////////////////////////////////////////////
// LoadFrap is a modification of LoadScanImage to load FRAP experiments
// saved in iVision (script?) on one specific microscope. It can be used as
// a template for similar procedures.

Function /t LoadFrap()

String ImgWaveName, FirstWave
Variable PointPos

ImageLoad /Q /O /C=-1

if (v_flag == 0)
	return "-1"
endif

	//removing the extension from the filename

	PointPos = strsearch(S_Filename, ".tif", 0)
	ImgWaveName = S_FileName[0,PointPos-1]
	ImgWaveName = ReplaceString("-", ImgWaveName, "_")
	
	PointPos = strsearch(S_Wavenames, ";", 0)
	FirstWave =S_Wavenames[0,PointPos-1]
	
if (waveexists($ImgWaveName))
	killwaves /z $ImgWaveName
endif


duplicate /o $FirstWave, $ImgWaveName
Killwaves /z $FirstWave

redimension /d $ImgWaveName		//convert to double precision floating point

setscale /p x,0,63e-6/512,"m" $ImgWaveName	//63µm is the side length of an image
setscale /p y,0,63e-6/512,"m" $ImgWaveName

variable zdim = 1

prompt  zdim, "Time between frames (s)"
doprompt "Enter variables", zdim
if(v_flag)
	setscale /p z,0,1,"Frame", $ImgWaveName
else
	setscale /p z,0,zdim,"s", $ImgWaveName
endif

Note $ImgWaveName, "file.path="+s_path
Note $ImgWaveName, "file.name="+s_filename

Return ImgWaveName

end


//////////////////////////////////////////////////////

Function nChannelsFromHeader(PicWave)
	wave PicWave
	string header = note(PicWave)
	
	Variable nChannels=0
	
	String nChannelsStr = StringByKey("SI.hChannels.channelsActive", TrimString(ReplaceString("\n",header,":"),1)," = ", ":")
	//nChannels = NumberByKey("state.acq.savingChannel1", Header, "=","\r") + NumberByKey("state.acq.savingChannel2", Header, "=","\r") + NumberByKey("state.acq.savingChannel3", Header, "=","\r")

	if(strlen(nChannelsStr)==3)
		nChannels=1
	elseif(strlen(nChannelsStr)==5)
		nChannels=2
	elseif(strlen(nChannelsStr)==7)	
		nChannels=3
	elseif(strlen(nChannelsStr)==9)
		nChannels=4
	endif
	
	return nChannels

End

//////////////////////////////////////////////////////

Function SplitChannels(PicWave,nChannels)
	wave PicWave
	variable nChannels
	
	variable nFrames = DimSize(PicWave,2), FramesPerChannel, Rest, ii
	String wvName
	FramesPerChannel=nFrames/nChannels
	Rest=FramesPerChannel-trunc(FramesPerChannel)
	
	if(Rest)
		Print "WARNING: inequal number of frames per channel."
		FramesPerChannel=trunc(FramesPerChannel)
	endif
	
	For(ii=0;ii<nChannels;ii+=1)
		wvName=NameOfWave(PicWave)+"_Ch"+Num2Str(ii+1)
		Duplicate /o PicWave $wvName
		wave w=$wvName
		Redimension/n=(-1,-1,FramesPerChannel) w
		MultiThread w=PicWave[p][q][r*nChannels+ii]	
	EndFor	
End

//////////////////////////////////////////////////////

Function/wave InvertImage(image)			//inverts an image
	Wave image
	
	Variable imax = WaveMax(image)
	String ResName = NameOfWave(image)+"_inv"
	
	Duplicate/o/free image inv_image
	
	FastOP inv_image = -1*image+(imax)
	
	Duplicate/o inv_image $ResName
	
	Return $ResName
End


//////////////////////////////////////////////////////

Function/wave AverageFrames(imStack, n)		//averages every n frames in a stack
	wave imStack
	variable n
	
	variable nFrames, ii, newn, start, stop
	string newname = NameOfWave(imStack)+"_av"+num2str(n)
	
	nFrames=DimSize(imStack,2)
	
	newn = round(nFrames/n)
	
	duplicate/o/free imStack, av_stack
	redimension/n=(-1,-1,newn) av_stack
	
	For(ii=0;ii<newn;ii+=1)
	
		start=ii*n
		stop=start+n-1
		
		duplicate/o/free /r=[0,*][0,*][start,stop] imStack, tobeAv
		
		MatrixOP /o/free avImage=sumBeams(tobeAv)/n
	
		av_stack[][][ii]=avImage[p][q]
	
	
	EndFor
	
	duplicate/o av_stack, $newName
	return $newName
End

//////////////////////////////////////////////////////

Function /wave AutoAverageFrames(imStack)
	wave imStack

	string header = note(imStack)
	variable n = NumberByKey("state.acq.numberOfFrames", Header, "=","\r")		//get n from Header

	if(numtype(n))
		Abort "unablo to retrieve state.acq.numberOfFrames from wave note."
	else
		wave result= AverageFrames(imStack, n)
	endif

	return result
End
	
	