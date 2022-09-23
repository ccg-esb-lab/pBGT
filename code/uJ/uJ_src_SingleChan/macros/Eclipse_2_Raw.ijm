print("> Eclipse_2_Raw");

setupFile=getArgument();
print("   setupFile:"+setupFile);

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

frame_names=split(list_frames,",");

//list_frames=split(list_frames,",");
//print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/

pathDATA=getParam(setupFile, "pathDATA");
dirECLIPSE=pathDATA+getParam(setupFile, "dirECLIPSE");  //Input
dirRAW=pathDATA+getParam(setupFile, "dirRAW");  //Output
exec("mkdir "+dirRAW);  //Create output directory (level1)

for(i=0; i<list_pos.length; i++){
	exec("mkdir "+dirRAW+list_pos[i]);  //Create output directory (level2)
	for(j=0; j<str_channel.length; j++){
		exec("mkdir "+dirRAW+list_pos[i]+"/"+str_channel[j]);  //Create output directory (level3)
	}
}

/********************************************************************/
//							MACRO
/********************************************************************/


setBatchMode(true);
for(p=0; p<numPos; p++){
	for(c=0; c<list_channel.length; c++){
		
		//Loads stack
		str_filter=list_pos[p]+""+list_channel[c];  
        run("Image Sequence...", "open="+dirECLIPSE+" file="+str_filter+" sort");
        
        //Export substack as separate files
        run("Make Substack...", "  slices="+list_frames);
        dirName=dirRAW+""+list_pos[p]+"/"+str_channel[c]+"/";
        print("   Saving "+dirName);
		for (s=1; s<=nSlices; s++){   	
  			run("Make Substack...", " slices="+s); 
  			index=IJ.pad(frame_names[s-1],3);                  
  			fileName=expeLabel+"_"+list_pos[p]+"_"+str_channel[c]+"_"+index+".tif"; 
			saveAs("Tiff", dirName+fileName); 
  			close(); 	
 		} 
 		
 		//Close all windows
 		while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      	} 
	}
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

exit;