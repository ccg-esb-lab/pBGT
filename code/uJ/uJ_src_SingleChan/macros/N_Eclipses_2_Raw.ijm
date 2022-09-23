print("> N-Eclipses_2_Eclipse");

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
signal_file=getParam(setupFile, "signal_file");
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
//exec("mkdir "+dirRAW);  //Create output directory (level1)
//exec("mkdir "+dirECLIPSE);  //Create output directory (level1)
File.makeDirectory(dirRAW);
//File.makeDirectory(dirECLIPSE);

/********************************************************************/
//							MACRO
/********************************************************************/

setBatchMode(true);
signal_file=pathDATA+""+signal_file;


print(signal_file);

signal=File.openAsString(signal_file);
signal=split(signal,"\n");


dirs_list=newArray();
frames_plus_list=newArray();
frame_plus=0;
	last_dir="";

for(i=1;i<signal.length;i++){
	
	line=split(signal[i],"\t");
	frame=parseInt(line[0]);
	dir=line[4];
	
	if(last_dir!=dir){
		frame_plus=frame;
		dirs_list=Array.concat(dirs_list,dir);
		frames_plus_list=Array.concat(frames_plus_list,frame_plus);

			
			}
last_dir=dir;	
	}
Array.print(dirs_list);
Array.print(frames_plus_list);




for(posi=0;posi<list_pos.length;posi++){
//for(posi=0;posi<2;posi++){
	pos=list_pos[posi];
	this_dir=dirRAW+pos;
	File.makeDirectory(this_dir);
	
	for(ci=0;ci<list_channel.length;ci++){
	//for(ci=0;ci<1;ci++){
		channel=list_channel[ci];
		channel_name=str_channel[ci];
		this_dir=dirRAW+pos+"/"+channel_name+"/";
		File.makeDirectory(this_dir);
		
		frame_plus=0;
		last_dir="";
		for(di=0;di<dirs_list.length;di++){
			dir=dirs_list[di];
			frame_plus=frames_plus_list[di];
								
			
			n_start=frame_plus;
			print(this_dir,"--->",pos,channel,n_start);	
			n_star=IJ.pad(n_start,3);
			img_prefix=expeLabel+"_"+pos+"_"+channel_name+"_";
			save_name=this_dir+img_prefix+n_start+".tif";
			
		

			if(dir=="x"){
				next_plus=frames_plus_list[di+1];
				
				for(x=frame_plus;x<next_plus;x++){
				
				newImage("Untitled", "8-bit noise", 640, 512, 1);
				run("16-bit");
				}
			run("Images to Stack", "name=Stack title=[] use");	
			//print(save_name);
			}else{
				run("Image Sequence...", "open="+pathDATA+dir+" file="+pos+channel+" sort");
				}
			
			run("Image Sequence... ", "format=TIFF name="+img_prefix+" start="+n_start+" digits=3  save="+save_name);
			close();
				
		
		
			
		
		
		}
	
	}
	
	
		
	
	
}



//Close all windows
 		while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      	} 

exit;






















/*
 
 newImage("Untitled", "8-bit black", 640, 512, 1);
run("Tiff...");
saveAs("Tiff", "/home/charly/Lab/Projects/uJ/uJ_data/HT-Sines/20180422_Charly_HT13_Sine_8hrs_rep2/data_eclipse_002/Untitled.tif");

 


 
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


*/
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