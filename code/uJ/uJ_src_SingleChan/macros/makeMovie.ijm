print("\n> displayLayers");

args=split(getArgument(),",");
Array.print(args);
setupFile=args[0];
print("setupFile:"+setupFile);

currentPos=args[1]; 
print("currentPos:"+currentPos);

list_channels=split(args[2],"+");
print("channels: "+args[2]);

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
dirFIGURES=pathDATA+getParam(setupFile, "dirFIGURES"); //Input
dirMOVIES=pathDATA+getParam(setupFile, "dirMOVIES"); //Output

File.makeDirectory(dirMOVIES);  //Create output directory (level1)

/********************************************************************/
//							MACRO
/********************************************************************/

setBatchMode(true); 
	
print("*"+list_channels.length);
for(c=0; c<list_channels.length; c++){	
	currentChannel=list_channels[c];
	File.makeDirectory(dirMOVIES+"/"+currentChannel);  //Create output directory (level3)

	print(c+": makeMovie: "+" pos="+currentPos+"");

	
	dirStack="";
	
	dirStack=dirTIF+currentPos+"/"+currentChannel+"/";	
	prefileNameStack=dirStack+expeLabel+"_"+currentChannel+"_"+currentPos+"_";
	
	
	for(t=0; t<list_frames.length; t++){
		thisFrame=IJ.pad(list_frames[t],3);
		
		fileNameFrame=prefileNameStack+""+thisFrame+".tif";
			print("> "+fileNameFrame);

		existFrame=File.exists(fileNameFrame);
		if((existFrame)==1){
		
			if(t==0){
				open(fileNameFrame);
				rename(currentChannel);
			}else{
				open(fileNameFrame);
				run("Copy");
				close();

				selectWindow(currentChannel);
				run("Add Slice");
				run("Paste");
			}
		}

	}	
	
	if(currentChannel=="GFP")
		run("Merge Channels...", "c2=GFP");
	if(currentChannel=="DsRed")
		run("Merge Channels...", "c1=DsRed");
			
	fileNameMovie=expeLabel+"_"+currentPos+"_"+currentChannel+".avi";
	run("AVI... ", "compression=JPEG frame=7 save="+dirMOVIES+""+currentChannel+"/"+fileNameMovie);
	print("Saving "+dirMOVIES+""+currentChannel+"/"+fileNameMovie);	

	close("*");
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
  				//print(">"+list[i]);
  			}
        }
     }
     
     return str_rois;
  }
print("Exit");
exit;