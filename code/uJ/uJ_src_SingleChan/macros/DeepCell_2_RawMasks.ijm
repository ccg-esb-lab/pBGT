print("\n> DeepCell_2_RawMasks");

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

if(segmentable_channel!=0){
str_channel=newArray(""+segmentable_channel); //processes the segmentable channel
}


list_frames=split(list_frames,",");
print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/

pathDATA=getParam(setupFile, "pathDATA");
dirDEEPCELL=pathDATA+getParam(setupFile, "dirDEEPCELL");  //Input
dirDEEPCELL_masks=dirDEEPCELL+"Masks/";  //Input
dirRAWMASKS=pathDATA+getParam(setupFile, "dirRAWMASKS"); //Output

//processes onyl DIC
	str_channel=newArray(str_channel[0]);
///////////


exec("mkdir "+dirRAWMASKS);  //Create output directory (level1)
for(i=0; i<list_pos.length; i++){
	exec("mkdir "+dirRAWMASKS+list_pos[i]);  //Create output directory (level2)
	for(j=0; j<str_channel.length; j++){
		exec("mkdir "+dirRAWMASKS+list_pos[i]+"/"+str_channel[j]);  //Create output directory (level3)
		}
}

/********************************************************************/
//							MACRO
/********************************************************************/

//threshold_Down=10;
threshold_Down=50;
n_range=10;


for(p=0; p<numPos; p++){
	for(c=0; c<str_channel.length; c++){
		dirInput=dirDEEPCELL_masks+""+list_pos[p]+"/"+str_channel[c]+"/";
		dirOutput=dirRAWMASKS+""+list_pos[p]+"/"+str_channel[c]+"/";
		print(dirInput);
		
		run("Image Sequence...", "open="+dirInput+" file=masks_ sort");
		rename(list_pos[p]+"-"+str_channel[c]);
		selectWindow(list_pos[p]+"-"+str_channel[c]);
		a=getImageInfo();
		a=split(a,"\n");
		//Array.print(a);
		w=split(a[8]," ");
		h=split(a[9]," ");
		d=split(a[10]," ");
		w=parseInt(w[1]);
		h=parseInt(h[1]);
		d=parseInt(d[1]);


		print(w,h,d);
		
		selectWindow(list_pos[p]+"-"+str_channel[c]);
		
		
		run("Orthogonal Views");
		selectWindow("XZ "+h/2);
		run("Copy");
		run("Internal Clipboard");
		selectWindow("Clipboard");
		run("8-bit");
		rename(list_pos[p]+"_X");
		//saveAs("Tiff", dirTemp+list_pos[p]+"_X.tif");

		selectWindow("YZ "+w/2);
		run("Copy");
		run("Internal Clipboard");
		run("Rotate 90 Degrees Left");
		run("8-bit");
		rename(list_pos[p]+"_Y");
		//saveAs("Tiff", dirTemp+list_pos[p]+"_Y.tif");
		
		selectWindow("XZ "+h/2);
		close();
		
	

		selectWindow(list_pos[p]+"_X");
		makeRectangle(1,1,w-2,d-2);

				
		run("Plot Profile");
		xpoints=newArray();
		ypoints=newArray();
		xvalue_up=0;
		xvalue_down=0;
		Plot.getValues(xpoints, ypoints);
		//Array.print(xpoints);
		//Array.print(ypoints);
		
		Array.getStatistics(ypoints, min, max, mean, stdDev);
		threshold_Down=mean;
		//threshold_Down=mean-mean/2;
		print(threshold_Down);
		flagUp=1;
		flagDown=1;
		start=floor(ypoints.length/2);
		
		for(i=0;i<start-n_range;i++){
			
			i_up=start+i;
			i_down=start-i;
			//y_up=ypoints[i_up];
			y_up=getArrayAverage(ypoints,i_up,i_up+n_range);
			//y_down=ypoints[i_down];
			y_down=getArrayAverage(ypoints,i_down-n_range,i_down);
			//print(i,start,i_up,i_down,y_up,y_down);
			if(y_up<threshold_Down && flagUp ){
				print("--->"+i,start,i_up,i_down,y_up,y_down);
				xvalue_down=xpoints[i_up];
				flagUp=0;	
				}	
			if(y_down<threshold_Down && flagDown ){
				print("<---"+i,start,i_up,i_down,y_up,y_down);
				xvalue_up=xpoints[i_down];
				flagDown=0;	
				}	
		
			}
	
	
		
		close(list_pos[p]+"_X");
		close("Plot of "+list_pos[p]+"_X");
		//close("Clipboard");

		selectWindow(list_pos[p]+"_Y");
 		makeRectangle(1,1,h-2,d-2);
		run("Plot Profile");
	
		xpoints=newArray();
		ypoints=newArray();
		yvalue_up=0;
		yvalue_down=0;
		Plot.getValues(xpoints, ypoints);
		flagUp=1;
		flagDown=1;
		start=floor(ypoints.length/2);
		
		Array.getStatistics(ypoints, min, max, mean, stdDev);
		threshold_Down=mean;
		//threshold_Down=mean-mean/2;
		
		
		for(i=0;i<start-n_range;i++){
			i_up=start+i;
			i_down=start-i;
			//y_up=ypoints[i_up];
			//y_down=ypoints[i_down];
			
			y_up=getArrayAverage(ypoints,i_up,i_up+n_range);
			y_down=getArrayAverage(ypoints,i_down-n_range,i_down);
			
						
			if(y_up<threshold_Down && flagUp){
				yvalue_down=xpoints[i_up];
				flagUp=0;	
				}	
			if(y_down<threshold_Down && flagDown){
				yvalue_up=xpoints[i_down];
				flagDown=0;	
				}	
		
			}
	
		close(list_pos[p]+"_Y");
		close("Plot of "+list_pos[p]+"_Y");
		//close("Clipboard");
		rx=xvalue_up;
		ry=yvalue_up;
		rw=xvalue_down-xvalue_up;
		rh=yvalue_down-yvalue_up;
			
		//print("Values: "+xvalue_up+" "+yvalue_up+" "+xvalue_down+" "+yvalue_down);
		if(rx==0 || ry==0 || rw<350 || rh<2000){
			print("Rectangle Not Found");
			
			waitForUser("Help!","We clould not find the trap,\nplease set it yourself with a rectangle\nthen click OK");
			setTool("rectangle");
			getSelectionBounds(rx, ry, rw, rh);
			run("Select None");
			//exit;
			}
		//print(rxx,ryy,rww,rhh);
		print(rx,ry,rw,rh);		
			
		print("Rectangle: "+xvalue_up+" "+yvalue_up+" "+xvalue_down-xvalue_up+" "+yvalue_down-yvalue_up);
			
			selectWindow(list_pos[p]+"-"+str_channel[c]);
			for (s=1; s<=nSlices; s++){ 
			
  				run("Make Substack...", " slices="+s); 
  				index=IJ.pad(list_frames[s-1],3); 
  				fileName=expeLabel+"_"+list_pos[p]+"_"+str_channel[c]+"_"+index+".tif";
  				
				makeRectangle(rx,ry,rw,rh);
				run("Clear Outside");
				run("8-bit");
				run("Canvas Size...", "width=640 height=512 position=Center");
  				saveAs("Tiff", dirOutput+fileName);
				close();
  				print(dirOutput+fileName);
		
 				} 
 		//close(list_pos[p]+"-"+str_channel[c]);
 		//Close all windows
 			while(nImages>0){ 
          		print("Closing..."+nImages);
          		selectImage(nImages); 
 			    close(); 
 			}
 		
 		
		}
		
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


function getArrayAverage(array,from,to){
	avg=0;
	n=to-from;
	for(i=from;i<=to;i++){
		avg=avg+array[i];
		}
	//avg=avg/n;


	return avg;
}





exit;