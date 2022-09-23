print("\n> TIF_2_Segmentable");

setupFile=getArgument();
print("setupFile:"+setupFile);

/********************************************************************/
//					PARSE PARAMETERS
/********************************************************************/

//Parse user-defined Parameters
expeLabel=getParam(setupFile, "expeLabel");
drugLabel=getParam(setupFile, "drugLabel");
list_pos=split(getParam(setupFile, "list_pos"),',');
list_frames=getParam(setupFile, "list_frames");
list_channel=split(getParam(setupFile, "list_channel"),',');
str_channel=split(getParam(setupFile, "str_channel"),',');
segmentable_channel=getParam(setupFile, "segmentable_channel");

numPos=list_pos.length;
if(indexOf(list_frames,"-")>0){
    from_frame=substring(list_frames, 0, indexOf(list_frames,"-"));
    to_frame=substring(list_frames, indexOf(list_frames,"-")+1);
    list_frames=""+from_frame;
    for(i=from_frame;i<=to_frame; i++){
	  if(i>from_frame)
        list_frames=list_frames+","+i;
    }
}
frame_names=split(list_frames,",");
list_frames=split(list_frames,",");
//processes onyl DIC
	//str_channel=newArray(str_channel[0]);
///////////
niters=25;
//niters=1;
print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/

pathDATA=getParam(setupFile, "pathDATA");
dirTIF=pathDATA+getParam(setupFile, "dirTIF");  //Input
dirSEGMENTABLE=pathDATA+getParam(setupFile, "dirSEGMENTABLE");  //Output
exec("mkdir "+dirSEGMENTABLE);  //Create output directory (level1)

for(i=0; i<list_pos.length; i++){
	exec("mkdir "+dirSEGMENTABLE+list_pos[i]);  //Create output directory (level2)
	for(j=0; j<str_channel.length; j++){
		exec("mkdir "+dirSEGMENTABLE+list_pos[i]+"/"+str_channel[j]);  //Create output directory (level3)
		}
}

/********************************************************************/
//							MACRO
/********************************************************************/
psf=pathDATA+"PSF_Real-16bit.tif";
setBatchMode(true);
for(p=0; p<numPos; p++){
	for(c=0; c<str_channel.length; c++){
	
		//Loads stack
		
        dirInput=dirTIF+""+list_pos[p]+"/"+str_channel[c]+"/";
		print("Input: "+dirInput);
		run("Image Sequence...", "open="+dirInput+" sort");
        
        //Export substack as separate files
        dirOutput=dirSEGMENTABLE+""+list_pos[p]+"/"+str_channel[c]+"/";
        
		for (s=1; s<=frame_names.length; s++){ 		
  			run("Make Substack...", " slices="+s); 
  			index=IJ.pad(frame_names[s-1],3);     
  			fileName=expeLabel+"_"+list_pos[p]+"_"+str_channel[c]+"_"+index+".tif";
  			
  			decon(dirInput,dirOutput,fileName,psf,niters);
  			open(dirOutput+fileName);
			run("Clear Results");
			run("Measure");
			mode=getResult("Mean",0);
			
			if(mode==0){
				close(fileName);
				decon(dirInput,dirOutput,fileName,psf,1);
				open(dirOutput+fileName);	
				}
			selectWindow("Results");
			run("Close");
			//close("Results");
			selectWindow(fileName);
			
			run("Enhance Contrast...", "saturated=0 normalize equalize");
			if(segmentable_channel=="DIC"){
				run("Gaussian Blur...", "sigma=1");
				run("Subtract Background...", "rolling=100  sliding");
				}
			//run("Gaussian Blur...", "sigma=1");
			//run("Invert");
			run("8-bit");
			//run("Subtract Background...", "rolling=50 light sliding");
			if(segmentable_channel=="GFP"||segmentable_channel=="DsRed"){
				run("Subtract Background...", "rolling=100  sliding");
				}
			//run("Subtract Background...", "rolling=100  sliding");
			saveAs("Tiff", dirOutput+fileName);
			close();
  			
  			
			//print(dirOutput+fileName);
			
  			close(); 	
 		} //end of slices
 		close();
	} //end of channels

	if(segmentable_channel=="DsRed+GFP"){
		
		seg_channels=split(segmentable_channel,"+");
		dirInput=dirSEGMENTABLE+""+list_pos[p]+"/";
		dirOutput=dirSEGMENTABLE+""+list_pos[p]+"/"+segmentable_channel+"/";
		exec("mkdir "+dirOutput);  //Create output directory (level1)
		
		for(f=0;f<list_frames.length;f++){
				 	
			index=IJ.pad(list_frames[f],3);
			imgs=newArray(seg_channels.length);
			for(c=0;c<seg_channels.length;c++){
				fileName=expeLabel+"_"+list_pos[p]+"_"+seg_channels[c]+"_"+index+".tif";
				img=dirInput+seg_channels[c]+"/"+fileName;
				//Array.concat(imgs,img);
				imgs[c]=fileName;
				open(img);
			}
			
			fileName=expeLabel+"_"+list_pos[p]+"_"+segmentable_channel+"_"+index+".tif";
			imageCalculator("Min create", imgs[0],imgs[1]);	
			//run("Subtract Background...", "rolling=50 light sliding");
			run("Subtract Background...", "rolling=10 sliding");   //this one is the good one AND NO gb
			saveAs("Tiff", dirOutput+fileName);
			close(imgs[0]);
			close(imgs[1]);
			close(fileName);
			print(dirOutput+fileName);
			}
		
		
		
		}


	
} //end of pos
setBatchMode(false);



/********************************************************************/
//						uJ_functions
/********************************************************************/
function getParam(setupFile, varName){
	filestring=File.openAsString(setupFile);
	rows=split(filestring, "\n"); 
	ret="";
	for(i=0; i<rows.length; i++){ 
		columns=split(rows[i],"="); 
		if(columns[0]==varName)
			ret=columns[1];
	} 
	return ret;
}

function decon(input,output,filename,pathToPsf,iters) {
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// MRNSD
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//print("Entering decon with:",input,"\n",output,"\n",filename,"\n",pathToPsf,"\n",subname);
	pathToBlurredImage = input + filename;
	filename=split(filename,".");
	filename=filename[0];
	pathToDeblurredImage = output + filename;
	preconditioner = "FFT"; //available options: FFT, NONE
	preconditionerTol = "-1"; //if -1, then GCV is used to compute preconditionerTol
	boundary = "REFLEXIVE"; //available options: REFLEXIVE, PERIODIC, ZERO
	resizing = "AUTO"; // available options: AUTO, MINIMAL, NEXT_POWER_OF_TWO
	output = "SAME_AS_SOURCE"; // available options: SAME_AS_SOURCE, BYTE, SHORT, FLOAT  
	precision = "SINGLE"; //available options: SINGLE, DOUBLE
	stoppingTol = "-1"; //if -1, then stoppingTol is computed automatically
	threshold = "-1"; //if -1, then disabled
	logConvergence = "false";
	//maxIters = "25";
	maxIters = iters;
	nOfThreads = "8";
	showIter = "false";

	//print("maxIters: "+maxIters);
	
	call("edu.emory.mathcs.restoretools.iterative.ParallelIterativeDeconvolution2D.deconvolveMRNSD", pathToBlurredImage, pathToPsf, pathToDeblurredImage, preconditioner, preconditionerTol, boundary, resizing, output, precision, stoppingTol, threshold, logConvergence, maxIters, nOfThreads, showIter);
}






exit;