print("\n> TIF_2_Composite");

args=split(getArgument(),",");
Array.print(args);
setupFile=args[0];
print("setupFile:"+setupFile);
currentPos=args[1]; 
print("currentPos:"+currentPos);
currentChannel=args[2]; 
print("currentChannel:"+currentChannel);

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
list_frames=split(list_frames,",");
print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/

pathDATA=getParam(setupFile, "pathDATA");
dirTIF=pathDATA+getParam(setupFile, "dirTIF"); //Input
dirROIS=pathDATA+getParam(setupFile, "dirROIS"); //Input
dirCOMPOSITE=pathCOMPOSITE+getParam(setupFile, "dirCOMPOSITE"); //Output

File.makeDirectory(dirCOMPOSITE);  //Create output directory (level1)


/********************************************************************/
//							MACRO
/********************************************************************/

setBatchMode(true); 

if(currentPos!=""){
	str_rois=getROIFiles(dirROIS+currentPos+"/", "");
}else{
	str_rois=getROIFiles(dirROIS, "");
}


list_rois=split(str_rois,",");
print("Loading "+list_rois.length+" Images from "+dirROIS);


for(i=1; i<list_rois.length; i++){
	run("Clear Results");

	thisPos=substring(list_rois[i],indexOf(list_rois[i], "_xy")+1,indexOf(list_rois[i], "_xy")+5);
	thisFrame=substring(list_rois[i],indexOf(list_rois[i], ".zip")-3,indexOf(list_rois[i], ".zip"));	
	print(i+": Analyzing "+" pos="+thisPos+", frame="+thisFrame);

	//Open images to extract fluorescence intensity
	fileNameGFP=expeLabel+"_"+thisPos+"_GFP_"+thisFrame+".tif";
	open(dirTIF+""+thisPos+"/GFP/"+fileNameGFP);
	rename("GFP");
	
	fileNameDsRed=expeLabel+"_"+thisPos+"_DsRed_"+thisFrame+".tif";
	open(dirTIF+""+thisPos+"/DsRed/"+fileNameDsRed);
	rename("DsRed");
	
	//Loads ROI Manager
	roiManager("reset");
	thisChannel="DsRed+GFP";  //TMP
	roiManager("Open",dirROIS+""+thisPos+"/"+thisChannel+"/"+list_rois[i]);
	
	//Now exports data
	File.makeDirectory(dirDATA+thisPos);  //Create output directory (level2)
	fileNameData=expeLabel+"_"+thisPos+"_"+thisFrame+".txt";
	numRois = roiManager("count");
	print("   Exporting "+numRois+" cells from pos:"+thisPos+" frame:"+thisFrame+"");
	run("Clear Results"); 
	for (j=0; j<numRois; j++) {
	
		selectWindow("GFP");
		roiManager("select", j);
		label=call("ij.plugin.frame.RoiManager.getName", j);
  
		getStatistics(area, mean);
		setResult("id",j, ""+j);
		//setResult("imageDIC", j, fileNameDIC);
		setResult("imageGFP", j, fileNameGFP);
		setResult("imageDsRed", j, fileNameDsRed);
		setResult("pos",j, thisPos);
		setResult("frame",j, thisFrame);
		setResult("roi_label",j, label);
		setResult("GFP", j, mean);

		selectWindow("DsRed");
		roiManager("select", j);
		label=call("ij.plugin.frame.RoiManager.getName", j);

		getStatistics(area, mean);
		setResult("DsRed", j, mean);	
  
		updateResults();
		saveAs("Results", dirDATA+""+thisPos+"/"+fileNameData);
		//print("   Exporting "+numRois+" cells from pos:"+thisPos+" frame:"+thisFrame+"");
		
	}  //for rois
	run("Clear Results"); 
	run("Close All");
	roiManager("reset");
	
	
}

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

function getROIFiles(dirROI, str_rois) {
     list = getFileList(dirROI);
     for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/"))
           str_rois=str_rois+getROIFiles(""+dirROI+list[i], str_rois);
        else{
        	if (endsWith(list[i], ".zip")){
  				str_rois = str_rois+","+list[i];
  				print(">"+list[i]);
  			}
        }
     }
     
     return str_rois;
  }
print("Exit");
exit;