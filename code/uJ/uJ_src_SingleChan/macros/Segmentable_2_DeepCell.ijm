print("\n> Segmentable_2_DeepCell");

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
list_frames=split(list_frames,",");

//str_channel=newArray(str_channel[0]); //processes onyl DIC
if(segmentable_channel!=0){
str_channel=newArray(""+segmentable_channel); //processes the segmentable channel
}
print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/

pathDATA=getParam(setupFile, "pathDATA");

dirSEGMENTABLE=pathDATA+getParam(setupFile, "dirSEGMENTABLE");  //Input
dirDEEPCELL=pathDATA+getParam(setupFile, "dirDEEPCELL");  //Output
exec("mkdir "+dirDEEPCELL);  //Create output directory (level1)
dirDEEPCELL_seg=dirDEEPCELL+"segmentable/";
exec("mkdir "+dirDEEPCELL_seg);  //Create output directory (level1)

for(i=0; i<list_pos.length; i++){
	exec("mkdir "+dirDEEPCELL_seg+list_pos[i]);  //Create output directory (level2)
	for(j=0; j<str_channel.length; j++){
		exec("mkdir "+dirDEEPCELL_seg+list_pos[i]+"/"+str_channel[j]);  //Create output directory (level3)
		}
}




/********************************************************************/
//							MACRO
/********************************************************************/

setBatchMode(true);
for(p=0; p<numPos; p++){
	for(c=0; c<str_channel.length; c++){
	
		//Loads stack
		
        dirInput=dirSEGMENTABLE+""+list_pos[p]+"/"+str_channel[c]+"/";
		print(dirInput);
		run("Image Sequence...", "open="+dirInput+" sort");
        
        //Export substack as separate files
        dirOutput=dirDEEPCELL_seg+""+list_pos[p]+"/"+str_channel[c]+"/";
        
		for (s=1; s<=nSlices; s++){ 
		
  			run("Make Substack...", " slices="+s); 
  			ss=IJ.pad(s-1,3);
  			fileName=ss+"_"+expeLabel+"_"+list_pos[p]+"_"+str_channel[c]+"_aligned_"+ss+".tif";
  			
  			saveAs("Tiff", dirOutput+fileName);
			close();
  			
  			
			print(dirOutput+fileName);
		
  		
 			} 
 		close();
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