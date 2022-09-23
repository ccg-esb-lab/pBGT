print("\n> TIF_2_Montage");

args=split(getArgument(),",");
Array.print(args);
setupFile=args[0];
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
pathDATA=getParam(setupFile, "pathDATA");
signalFile=pathDATA+getParam(setupFile, "signal_file");

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

end_frame=parseInt(to_frame);
print(signalFile);

print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/


dirTIF=pathDATA+getParam(setupFile, "dirTIF"); //Input

dirMontage=pathDATA+getParam(setupFile, "dirMontage"); //Output
l=lastIndexOf(dirMontage, "/");
dirMontage=substring(dirMontage,0,l);
dirMontage=dirMontage+"_single/";
print(dirMontage);




File.makeDirectory(dirMontage);  //Create output directory (level1)



/********************************************************************/
//							MACRO
/********************************************************************/

setBatchMode(true); 
//setBatchMode(false); 

	//Read the Signal file 

	fileStr=File.openAsString(signalFile); 
	fileArray=split(fileStr, "\n"); 
	
subset=newArray();
subset=Array.concat(subset,1);
flag=1;


	for(i=1;i<fileArray.length;i++){
		line=split(fileArray[i],"\t");
		frame=parseInt(line[0]);
		if(frame<=end_frame){
		
			red=255*parseFloat(line[2]);
			green=255*parseFloat(line[3]);
			blue=255*parseFloat(line[1]);

			if((green==255*0.5) && (red==255*0.5) && (flag==1)){
				
				subset=Array.concat(subset,i);
				flag=0;
				}
			if((green==255) || (red==255)){
				subset=Array.concat(subset,i);
				}
			if(frame==24){
				print(green,red,255*0.5);
				}
		
			}
		}
Array.print(subset);



	
for(c=0;c<str_channel.length;c++){
	channel=str_channel[c];
	flag=0;
	for(i=0; i<numPos; i++){
	//for(i=0; i<2; i++){
		thisPos=list_pos[i];
		thisDir=dirMontage+thisPos+"/";
		//File.makeDirectory(thisDir);
		
		//Open fluorescence images and edit them according to signal
		fileName=expeLabel+"_"+thisPos+"_"+channel+"_001.tif";
		run("Image Sequence...", "open="+dirTIF+thisPos+"/"+channel+"/"+fileName+" sort");
		
			
		
		
		for(j=0;j<subset.length;j++){
			s=IJ.pad(subset[j],3);
			selectWindow(channel);
			run("Make Substack...", " slices="+s); 
			run("Canvas Size...", "width=640 height=542 position=Bottom-Center");
		/*
			line=split(fileArray[s],"\t");
			red=255*line[2];
			green=255*line[3];
			blue=255*line[1];
			setForegroundColor(red, green, blue);
			fillRect(0, 0, 640, 30);
			*/
			frame_name=expeLabel+"_"+thisPos+"_"+channel+"_Frame_"+s;
			
			rename(frame_name);
		//	saveAs("Tiff", thisDir+frame_name+".tif");
			
			}
		
		selectWindow(channel);
		close();
	
		run("Images to Stack", "name="+channel+" title=[Frame] use");
	/*	if(flag==0){
			selectWindow(channel);
			run("Duplicate...", "duplicate");
			selectWindow(channel+"-1");
			run("Canvas Size...", "width=640 height=30 position=Top-Center");
			run("Make Montage...", "columns="+nSlices+" rows=1 scale=1 ");
			rename("xy00");
			saveAs("Tiff", dirMontage+expeLabel+"_Montage_"+channel+"_"+"xy00"+".tif");
			flag=1;
			selectWindow(channel+"-1");
			close();
			}*/
		
		selectWindow(channel);
		run("Canvas Size...", "width=640 height=512 position=Bottom-Center");
		//run("Make Montage...", "columns="+nSlices+" rows=1 scale=1 label");
		run("Make Montage...", "columns="+nSlices+" rows=1 scale=1 ");
		rename(thisPos);
		saveAs("Tiff", dirMontage+expeLabel+"_Montage_"+thisPos+"_"+channel+".tif");
		selectWindow(channel);
		//rename(expeLabel+"_Montage_"+thisPos);
		//saveAs("Tiff", dirMontage+expeLabel+"_Montage_"+thisPos+".tif");
		close();
	
		
	}

	run("Images to Stack", "name=StackXY title=xy use");
	
	
	run("Make Montage...", "columns=1 rows="+nSlices+" scale=1 ");
	rename(expeLabel+"_Montage");
	saveAs("Tiff", dirMontage+expeLabel+"_Montage_ALL"+"_"+channel+".tif");
	selectWindow(expeLabel+"_Montage_ALL"+"_"+channel+".tif");
	close();

	selectWindow("StackXY");
	close();

	//Close all windows
 		while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      	} 

	
}



setBatchMode(false); 

//Close all windows
 		while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      	} 


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



print("Exit");
exit;






