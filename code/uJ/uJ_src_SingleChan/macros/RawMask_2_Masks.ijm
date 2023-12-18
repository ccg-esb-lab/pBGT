/*
Requires ImageJ > 1.49
Last Update: 17/05/2017 (RPM)

Use:
<Space>: Re-color 
<Shift> <Space>: Toggle View
<Left-click>: If two cells selected, combine them. Else extend cell.
<Right-click>: Split cells
<Shift> <Left-click>: insert cell
<Shift> <Right-click>: delete cell.
<arrow up>: zoom out
<arrow down>: zoom in
<arrow left>: save/prev frame
<arrow right>: save/next frame
<s>: save
<q>: save/exit

NOTES:
- In RoiManager/Options: Uncheck Associate Show All ROIs with Slices
- Keyboard Listener is an ImageJ Plugin. Copy uJ_keyListener.class in Plugins folder


*/

print("\n> RawMask_2_Masks");

args=split(getArgument(),",");
Array.print(args);

setupFile=args[0];
/*
currentPos=args[1]; 
currentChannel=args[2]; //TMP 
doAll=args[3];
*/
print(setupFile);
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

/********************************************************************/
//					Setup Run
/********************************************************************/

	doAll=1;
	currentPos="All";
	list_frames_reindexed=list_frames;
	
	labMode=newArray("All","By trap");
  
	Dialog.create("Procedure mode");
	Dialog.addRadioButtonGroup("Mode", labMode, 1, 2, "By trap");
	Dialog.show();

	pMode=Dialog.getRadioButton;
	print("Mode: "+pMode);
  
	if(pMode!="All"){
		//doAll=0;	
		Dialog.create("Select trap");
		columns = 4;
		rows =-floor(-list_pos.length/columns );
		labels_index = newArray("natural","by layers");
		Dialog.addRadioButtonGroup("Reindex",labels_index,1,2,"natural");
		labels_index = newArray("Automatic","Manual");
		Dialog.addRadioButtonGroup("Correction",labels_index,1,2,"Manual");
		Dialog.addRadioButtonGroup("Trap",list_pos,rows,columns,list_pos[0]);
		Dialog.show();
		reindex=Dialog.getRadioButton;
		doAll=Dialog.getRadioButton;
		currentPos=Dialog.getRadioButton;
		list_pos=newArray(currentPos);
		
		if(reindex=="by layers"){
			list_frames_reindexed=makeNewIndexes(list_frames);	
			}
		
		if(doAll=="Manual"){
			doAll=0;
			}else{
				doAll=1;
				}
		
		
		}





	

if(segmentable_channel!=0){
str_channel=newArray(""+segmentable_channel); //processes the segmentable channel
currentChannel=str_channel[0];
}
//currentChannel=str_channel[0];
print("setupFile:"+setupFile);
print("currentPos:"+currentPos);
print("currentChannel:"+currentChannel);

print("Loading: "+list_pos.length+" positions/"+list_channel.length+" channels/"+list_frames.length+" frames");

print("List Frames:");
Array.print(list_frames);
print("List Frames Reindexed:");
Array.print(list_frames_reindexed);


/********************************************************************/
//					CREATE DIRECTORY STRUCTURE
/********************************************************************/
dirROOT=getParam(setupFile, "dirROOT");
pathDATA=getParam(setupFile, "pathDATA");
dirRAWMASKS=pathDATA+getParam(setupFile, "dirRAWMASKS"); //Input
dirSEGMENTABLE=pathDATA+getParam(setupFile, "dirSEGMENTABLE"); //Input
dirMASKS=pathDATA+getParam(setupFile, "dirMASKS"); //Output
dirROIS=pathDATA+getParam(setupFile, "dirROIS"); //Output





File.makeDirectory(dirROIS);
File.makeDirectory(dirMASKS);
for(i=0; i<list_pos.length; i++){
	File.makeDirectory(dirROIS+list_pos[i]);
	File.makeDirectory(dirMASKS+list_pos[i]);
	
	File.makeDirectory(dirROIS+list_pos[i]+"/"+currentChannel);
	File.makeDirectory(dirMASKS+list_pos[i]+"/"+currentChannel);
	
	
}

/********************************************************************/
//							MACRO
/********************************************************************/

lutPath=dirROOT+"uJ_src_pBGT/luts/glasbey.lut";
roiColors = loadColorsfromLut(lutPath);

r_cell=2.0;
//max_area_cell=5000;
//min_area_cell=50;

max_area_cell=10000;
min_area_cell=30;


num_thetas=60;
r_cut=15;



shift=1; ctrl=2; alt=8; leftButton=16; rightButton=4; insideROI = 32;
x2=-1; y2=-1; z2=-1; flags2=-1;

/************************************************/
// Loads image and setups environment


for(pos=0;pos<list_pos.length;pos++){	
/*
if(doAll==1){
	setBatchMode(true);
}
*/

	currentPos=list_pos[pos];
	currentFrame=0;
	currentFrame_reindexed=list_frames_reindexed[0];
	
	numFrames=list_frames.length;
	print("numFrames="+numFrames);
	
	while (currentFrame>-1) { 
		
		if(currentFrame==0){
			currentFrame=1;
		}else{
			uJCloseFrame();
		}
		
		uJLoadFrame(currentFrame_reindexed, currentPos, currentChannel, dirSEGMENTABLE, dirRAWMASKS, dirMASKS);
		thisFrame=IJ.pad(currentFrame_reindexed,3);
		selectWindow(expeLabel+"_"+currentPos+"_"+thisFrame+".tif");
		run("Original Scale");
		scy=screenHeight;
		scx=screenWidth;
		setLocation(scx/4, scy/4);
		//setLocation((scx/2), (scy/2));
		run("In [+]");run("In [+]");run("In [+]");
		
		//run("Set... ", "zoom=200 x=320 y=256");
		run("Set... ", "zoom=100");
		/************************************************/
		// Initial mask segmentation
		/************************************************/
		
		setSlice(1);  //Go to mask
		n=detect_all();
		print("  ... ("+n+" cells)");
		

		
		reColorRois(2);
		setSlice(2);  //Go to segmentable
		roiManager("Show All");
	
		/************************************************/
		// Correct Mask Segmentation
		/************************************************/
	
		selected_rois = newArray("0");
		prev_x=-1; prev_y=-1;
		prev_slice=2;
	
		overlay_on=false;
		left_clicked=false;
		correctingFrame=1;
	
		setOption("DisablePopupMenu", true);
		while (correctingFrame!=0) { 
		
			getCursorLoc(x, y, z, flags);
		
			//******* Plugin Key_Listener ********
			// Communicates with plugin via drawing pixels
			// (1,1): Save
			// (0,0): Quit
			// (1,0): Prev slice
			// (0,1): Next slice
			if(getPixel(1, 1)==255){  //SAVE
				setPixel(1, 1, 0);
				
				uJSaveFrame(currentFrame_reindexed, currentPos, currentChannel, dirMASKS, dirROIS);
				
				print("[save]");
			}
			if(getPixel(0, 0)==255){  //QUIT
				
				uJSaveFrame(currentFrame_reindexed, currentPos, currentChannel, dirMASKS, dirROIS);
			
				setPixel(0, 0, 0);
				correctingFrame=0; //Exit correctingFrame loop
				currentFrame=-1;   //Exit correctingScript loop
				
				print("[quit]");
			}
			if(getPixel(0, 1)==255){  //NEXT FRAME
				setPixel(0, 1, 0);
				print("**"+currentFrame+"*"+numFrames);
				if(currentFrame<numFrames && currentFrame>-1){
					correctingFrame=0;
					if(currentFrame==0)
						currentFrame=1;
					currentFrame=currentFrame+1;
					
					uJSaveFrame(currentFrame_reindexed, currentPos, currentChannel, dirMASKS, dirROIS);
					currentFrame_reindexed=list_frames_reindexed[currentFrame-1];
					print("[next] "+currentFrame_reindexed+"/"+numFrames);
				}else{  //End of movie
					setPixel(0, 0, 255); 
					uJSaveFrame(currentFrame_reindexed, currentPos, currentChannel, dirMASKS, dirROIS);
					
				}
				
				
			}
			if(getPixel(1, 0)==255){ //PREVIOUS FRAME
				setPixel(1, 0, 0);
				if(currentFrame>1){
					correctingFrame=0;
					currentFrame=currentFrame-1;
				}
				
				uJSaveFrame(currentFrame_reindexed, currentPos, currentChannel, dirMASKS, dirROIS);
				currentFrame_reindexed=list_frames_reindexed[currentFrame-1];
				print("[prev] "+currentFrame_reindexed+"/"+numFrames);
			}
			
			/******* <space>: Toggle view ********/
			if(isKeyDown("space") && isKeyDown("shift")){
					s = s + "<space>";
				  
					//Toggle overlay
					if(overlay_on==true){
						roiManager("Show None");
						overlay_on=false;
					}else{
						roiManager("Show All");
						overlay_on=true;
					}
				
					if(this_slice==1){
						reColorRois(2);
						setSlice(2);  //Go to dic
						roiManager("Show None");
						prev_slice=2;
					}else{
						reColorRois(1);
						setSlice(1);  //Go to mask
						prev_slice=1;
					}
					setKeyDown("None");
					roiManager("deselect");	
					//setKeyDown("Shift");
			}
			
			/******* <space>: Toggle overlay ********/
			if(isKeyDown("space") && isKeyDown("shift")==0){
					reColorRois(this_slice);
					setSlice(this_slice);
					roiManager("Show All");
					setKeyDown("None");
					roiManager("deselect");	
			}
			
			if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
			
				this_slice=getSliceNumber();
			
				s = " ";
				if (flags&shift!=0)
					 s = s + "<shift>";
				else
					selected_cell=-1;
				if (flags&ctrl!=0)
					s = s + "<ctrl>";
				if (flags&insideROI!=0)	s = s + "<inside>";
			
			
				/******* <left>: ADD/JOIN CELLS ********/
				if (flags&leftButton!=0 && flags&ctrl==0){
					s = s + "<left>";
				
					if(left_clicked==false){
						left_clicked=true;
						//print(s+": click");
						selected_rois = newArray();
					}
				
					this_roi=getRoiIndex(x,y,z);
					if(this_roi>-1){
						iarray=array_index(selected_rois, this_roi);
						if(iarray==-1){
							selected_rois = Array.concat(selected_rois,this_roi); 
						}
					}
				
					//print("(x"+x+","+y+"): "+this_roi);
					setForegroundColor(255,255,255);
					if(prev_x>-1 && prev_y>-1){
						setSlice(1);
						setLineWidth(2*r_cell);
						drawLine(prev_x, prev_y, x, y);
					}
					prev_x=x;
					prev_y=y;
					setSlice(prev_slice);	
					run("Select None");
				}else{
					if(left_clicked==true){
						//print(" <left>: release");
						left_clicked=false;
						prev_x=-1;
						prev_y=-1;
					
						if(lengthOf(selected_rois)>0){
							print("Merging Cells: ");
							Array.print(selected_rois); 
							roiManager("select", selected_rois);
							roiManager("delete");
						}else{  //insert new cell
							//print("  Inserting Cell: "+roi_label);
							setSlice(1);
							setForegroundColor(255,255,255);
							fillOval(x-r_cell, y-r_cell, 2*r_cell, 2*r_cell);
						}
					
						roi_label="_n"+n+"_x"+x+" y"+y;
						roi_label="roi_f"+currentFrame_reindexed+"_n"+n+""+"_x"+x+"_y"+y;
						icell=detectCell(x,y);
						if(icell>-1){
							print("  Updating Cell: "+roi_label);
							roiManager("update");
							lastRoi = roiManager("count")-1;
					
							roiManager("select", lastRoi);
							roiManager("rename", roi_label);
					
							n=n+1;
							
							
						}
						reColorRois(prev_slice);
							setSlice(prev_slice);	
							roiManager("Show All");
							roiManager("deselect");			
						
					}
				}
			
			
				/***************/
				if (flags&rightButton!=0){
					s = s + "<right>";
				
					/******* <shift><right>: DELETE CELL ********/
					if(flags&shift!=0){
						iRoi=getRoiIndex(x,y,z);
						if(iRoi>-1){
							roiManager("select", iRoi);
							roiManager("Set Fill Color", "black");
							setForegroundColor(0,0,0); 
							floodFill(x, y);
					
							roi_label=call("ij.plugin.frame.RoiManager.getName", iRoi);
							print("Deleting Cell: "+roi_label);
							roiManager("delete");
					
							roiManager("deselect");
							makeOval(0,0,0,0); //Para des-seleccionar
							
							reColorRois(prev_slice);
							setSlice(prev_slice);	
							roiManager("Show All");
						}
					
				
					/******* <right>: SPLIT CELL ********/	
					}else{
					
						setBatchMode(true);
					
						//If (x,y) is not in a cell, look around
						for(ri=0; ri<r_cut; ri++){
							foundSeed=0;
							for(i=0; i<8; i++){
								theta=i*2*PI/8;
								xi=floor(x+ri*sin(theta));
								yi=floor(y+ri*cos(theta));
								iRoi=getRoiIndex(xi,yi,z);
								//print("... ("+xi+","+yi+"): "+iRoi);
								if(iRoi>-1){
									//print("("+x+","+y+") : ("+xi+","+yi+") : "+iRoi);
									x=xi;
									y=yi;
									foundSeed=1;
									break;
								}
							}
							if(foundSeed==1)
								break;
						}				
						
						if(iRoi>-1){
							roi_label=call("ij.plugin.frame.RoiManager.getName", iRoi);
							print("Splitting Cell: "+roi_label);
	
							//First determine best angle
							best_theta=0;
							min_l=2*r_cut;
							for(i=0; i<num_thetas; i++){
								theta=i*PI/num_thetas;
								makeLine(x-r_cut*sin(theta), y-r_cut*cos(theta),x+r_cut*sin(theta),y+r_cut*cos(theta));
								profile = getProfile();
								l=0;
								for (j=r_cut; j<profile.length; j++)
									 if(profile[j]>0)
										l=l+1;
									 else
										break;
								for (j=r_cut; j>=0; j--)
									 if(profile[j]>0)
										l=l+1;
									 else
										break;
					
								if(l<=min_l){
									min_l=l;
									best_theta=theta;
								}
							}
	
							//print("__ best_theta:"+best_theta);
							makeLine(x-r_cut*sin(best_theta), y-r_cut*cos(best_theta),x+r_cut*sin(best_theta),y+r_cut*cos(best_theta));
							profile = getProfile();
							setForegroundColor(0,0,0);
							for (j=r_cut; j<profile.length; j++)
								if(profile[j]==255){
									setLineWidth(2);
									drawLine(x+(j-r_cut)*sin(best_theta),y+(j-r_cut)*cos(best_theta),x+(j-r_cut+1)*sin(best_theta),y+(j-r_cut+1)*cos(best_theta));
								}else
									break;	
							for (j=r_cut-1; j>=0; j--)
								if(profile[j]==255){
									setLineWidth(2);
									drawLine(x+(j-r_cut)*sin(best_theta),y+(j-r_cut)*cos(best_theta),x+(j-r_cut-1)*sin(best_theta),y+(j-r_cut-1)*cos(best_theta));
								}else
									break;
					
							setBatchMode(false);
						
							run("Select None");
							roiManager("deselect");
							roiManager("delete");
							n=detect_all();
							print("  ... ("+n+" cells)");
							
							reColorRois(prev_slice);
							setSlice(prev_slice);	
							roiManager("Show All");
							roiManager("deselect");	
						
						
						}
						
					}
				}
				
			}
			wait(100);
			
			if(doAll==1 && currentFrame<=numFrames){
				print("Auto Next...");
				setPixel(0, 1, 255);//Next
			}
			
		}
		x2=x; y2=y; z2=z; flags2=flags;
	
	} //While correcting
	
	print("****");
	setOption("DisablePopupMenu", false);
	uJCloseFrame();
	setBatchMode(false);
} //for pos


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
/************************************************/	

// Finds Maxima in mask
// Load positions from Results table 
// Filter and save cells in ROI manager
// Returns number of cells detected
function detect_all(){
			imageID=getImageID();
			run("Find Maxima...", "noise=10 output=List");
			
			headings = split(String.getResultsHeadings); 
			n=0;
			for (row=0; row<nResults; row++) { 
				x=getResult(headings[0],row); 
				y=getResult(headings[1],row);
				
				detectCell(x, y);  
				roiManager("select",roiManager("count")-1);
				roi_label="_n"+n+"_x"+x+" y"+y;
				roi_label="roi_f"+currentFrame_reindexed+"_n"+n+""+"_x"+x+"_y"+y;
				//roiManager("rename",expeLabel+" "+list_pos[pos]+" x"+x+" y"+y+" z"+list_frames[z]);
				roiManager("rename",roi_label);
						
				//Filter too-small or too-large cells
				getStatistics(area, mean);
				if(area>max_area_cell || area<min_area_cell){
			//		print("---------------->deleting cell: "+roi_label+"\tarea="+area);
					roiManager("delete");
					setForegroundColor(0,0,0); 
					floodFill(x, y);
				}else{
					n=n+1;
				}
			}  //for rows in results

			if(n==0){
			setForegroundColor(255,255,255);
			y=512/2;
			x=640/2;
			
			fillOval(x,y , min_area_cell/4, min_area_cell/4);
			//n=detect_all();  
			print("  ... Added one fake cell");
			roiManager("Show All");
			roiManager("deselect");	
			run("Find Maxima...", "noise=10 output=List");
			
			headings = split(String.getResultsHeadings); 
			
			for (row=0; row<nResults; row++) { 
				x=getResult(headings[0],row); 
				y=getResult(headings[1],row);
				
				detectCell(x, y);  
				roiManager("select",roiManager("count")-1);
				roi_label="_n"+n+"_x"+x+" y"+y;
				roi_label="roi_f"+currentFrame_reindexed+"_n"+n+""+"_x"+x+"_y"+y;
				//roiManager("rename",expeLabel+" "+list_pos[pos]+" x"+x+" y"+y+" z"+list_frames[z]);
				roiManager("rename",roi_label);
						
				//Filter too-small or too-large cells
				getStatistics(area, mean);
			
				if(area>max_area_cell || area<min_area_cell){
					//print("---------------->deleting cell: "+roi_label+"\tarea="+area);
					roiManager("delete");
					setForegroundColor(0,0,0); 
					floodFill(x, y);
				}else{
					n=n+1;
				}
			}  //for rows in results




			
			}
		


			
			//setSlice(2);
			//reColorRois(2);
			roiManager("Show All");
			roiManager("deselect");	
			selectImage(imageID);
			return n;
}
/************************************************/			
function detectCell(x, y){
	if( getPixel(x,y) > 255/2 ) { //Only if not black
		//doWand(x, y, 5.0, "Legacy");
		doWand(x, y, 25.0, "Legacy");
		roiManager("add");
		roiManager("select",roiManager("count")-1);
		return roiManager("count")-1;
	}else{
		return -1;
	}
}
/************************************************/
function reColorRois(this_slice){
	n = roiManager("count");
	for (i=0; i<n; i++) {
		roiManager("select", i);
		c=floor(roiColors.length*random()); 
		if(c<1)
			c=1;
		strokeColor="#FF"+roiColors[c];
		if(this_slice==1){ //MASK
			fillColor="#FF"+roiColors[c];
			roiManager("Set Fill Color", fillColor);
		}else{
			fillColor="#66"+roiColors[c];
			roiManager("Set Color", strokeColor);
		}
	}
	roiManager("Show All");
}
/************************************************/
//Returns the index of an array
function array_index(a, value) { 
      for (i=0; i<a.length; i++) 
          if (a[i]==value) return i; 
      return -1; 
  } 

/************************************************/
// Funtion that loads (and prepares) time-lapse image
function loadImage(fileName){
	open(fileName);
	run("Remove Overlay"); // cleans up any previous overlays -- may be deleted or commented out
	//run("RGB Color");
	run("8-bit");
	
	setPixel(0,0,0);
	setPixel(0,1,0);
	setPixel(1,0,0);
	setPixel(1,1,0);
}

/************************************************/
// Function that returns the index of a ROI from
// the RoiManager given its position
function getRoiIndex(x,y,z){
	setBatchMode(true);
	n = roiManager("count");
	for (i=0; i<n; i++) {
		label=call("ij.plugin.frame.RoiManager.getName", i);
	  		roiManager("select", i);
			//run("Enlarge...", "enlarge=10");
			inside = Roi.contains(x, y); 
			if (inside > 0 ) {
				setBatchMode(false);
				return i;
			}
	}
	setBatchMode(false);
	return -1;
}

/************************************************/
function loadColorsfromLut(lut) {
    //path = getDirectory("luts")+lut;
    path=lut;
    list = File.openAsRawString(path);
    rgbColor = split(list,"\n");
    if(rgbColor.length<255 || rgbColor.length>257)
        verboseExit("Error reading "+path, "Reason: Found unexpected number of columns");
    hexColor = newArray(256);
    firstStrg = substring(rgbColor[0], 0, 1);
    if (isNaN(firstStrg))
        k = 1;    // the lut file has a header (Index Red Green Blue)
     else
        k = 0;    // there is no header. First row is index 0
    for(i=k; i<256; i++) {
        hex = rgbToHex(rgbColor[i]);    
        hexColor[i] = hex;
    }
    return  hexColor;  
}
/************************************************/
function rgbToHex(color) {
    color1 = split(color,"\t");
    if(color1.length==1)
        color1 = split(color," ");
    if(color1.length==1)
        verboseExit("Chosen LUT does not seem to be either a tab-delimited",
                    "or a space-delimited text file");
    if(color1.length==4)
        i = 1;    // first column of the lut file is the index number
    else
        i = 0;    // first column of the lut file is the red value
    r = color1[0+i]; g = color1[1+i]; b = color1[2+i];
    return ""+toHex(r)+""+toHex(g)+""+toHex(b);
}
/************************************************/
function toHex(n) {
    n = parseInt(n);
    if(n==0 || isNaN(n))
        return "00";
    n = maxOf(0,n); n = minOf(n,255); n = round(n);
    hex = ""+substring("0123456789ABCDEF",((n-n%16)/16),((n-n%16)/16)+1)+
          substring("0123456789ABCDEF",(n%16),(n%16)+1);
    return hex;
}

/************************************************/
//This function returns an array of sorted 3-digits-numbers corresponding to frames/times with respect to priority 
// 18 corresponds to 18 frames 
// 12 corresponds to 12 frames 
function makeNewIndexes(list_frames){
	first_frame=list_frames[0];
	last_frame=list_frames[list_frames.length-1];
	
	times_18=newArray();
	times_12=newArray();
	times_6=newArray();
	times_3=newArray();
	times_2=newArray();
	times_1=newArray();
	new_indexes=newArray();
	times_x=Array.getSequence(last_frame+1); 
	for(i=first_frame,i2=first_frame,i3=first_frame,i6=first_frame,i12=first_frame,i18=first_frame;i<=last_frame;i+=1){
		if(times_x[i18]!=0){times_18=Array.concat(times_18,i18);times_x[i18]=0;}
		if(times_x[i12]!=0){times_12=Array.concat(times_12,i12);times_x[i12]=0;}	
		if(times_x[i6]!=0){times_6=Array.concat(times_6,i6);times_x[i6]=0;}
		if(times_x[i3]!=0){times_3=Array.concat(times_3,i3);times_x[i3]=0;}
		if(times_x[i2]!=0){times_2=Array.concat(times_2,i2);times_x[i2]=0;}
		if(times_x[i]!=0){times_1=Array.concat(times_1,i);times_x[i]=0;}
	
		if(i18<=last_frame-18){i18+=18;}
		if(i12<=last_frame-12){i12+=12;}	
		if(i6<=last_frame-6 ){i6+=6;}
		if(i3<=last_frame-3 ){i3+=3;}
		if(i2<=last_frame-2 ){i2+=2;}
		}

	list_frames_reindexed=Array.concat(times_18,times_12,times_6,times_3,times_2,times_1);
	return list_frames_reindexed;
}

/***********************************/
function uJLoadFrame(currentFrame, currentPos, currentChannel, dirSEGMENTABLE, dirRAWMASKS, dirMASKS){

	setBatchMode(true);
	thisFrame=IJ.pad(currentFrame,3);   
	if(File.exists(dirMASKS+currentPos+"/"+currentChannel+"/"+expeLabel+"_"+currentPos+"_"+thisFrame+".tif")==1){
		print("Loading existing mask");
		loadImage(dirMASKS+currentPos+"/"+currentChannel+"/"+expeLabel+"_"+currentPos+"_"+thisFrame+".tif");  //LOADS MASK
	}else{
		print("Loading existing raw mask");
		loadImage(dirRAWMASKS+currentPos+"/"+currentChannel+"/"+expeLabel+"_"+currentPos+"_"+currentChannel+"_"+thisFrame+".tif");  //LOADS RAW MASK
	}
	
	//LOADS SEGMENTABLE IMAGE
	if(File.exists(dirSEGMENTABLE+currentPos+"/"+currentChannel+"/"+expeLabel+"_"+currentPos+"_"+currentChannel+"_"+thisFrame+".tif")){
		loadImage(dirSEGMENTABLE+currentPos+"/"+currentChannel+"/"+expeLabel+"_"+currentPos+"_"+currentChannel+"_"+thisFrame+".tif");  
	}else{  //Segmentable image not found. Load mask again
		loadImage(dirRAWMASKS+currentPos+"/"+currentChannel+"/"+expeLabel+"_"+currentPos+"_"+currentChannel+"_"+thisFrame+".tif");  //LOADS RAW MASK

	}
	
	run("Images to Stack", "name="+expeLabel+"_"+currentPos+"_"+thisFrame+".tif");
	
	setOption("DisablePopupMenu", true);
	

	
	if (isOpen("ROI Manager")) {
		selectWindow("ROI Manager");
		run("Close");
		run("Clear Results");
	}

	
	setBatchMode(false);
	//run("Set... ", "zoom=150");
	
	
	
	run("ROI Manager...");
	//print("zooming");
	
	while (!isOpen("ROI Manager")){
		print("ROI Manager not open at load\n\tTrying to re-open... ");
		run("ROI Manager...");
		wait(1000);
		}
	
	roiManager("deselect");	
//	setBatchMode(false);

	run("uJ keyListener");
	
}

/***********************************/
function uJSaveFrame(currentFrame, currentPos, currentChannel, dirMASKS, dirROIS){
	//setBatchMode(true);
	//wait(500);		
	imageID=getImageID();
	this_slice=getSliceNumber();
	
	if(currentFrame==0)
		currentFrame=1;
	
	thisFrame=IJ.pad(currentFrame,3);   //Here we should prioritize
	fileNAME=expeLabel+"_"+currentPos+"_"+thisFrame+".tif";
	roiNAME=expeLabel+"_"+currentPos+"_"+thisFrame+".zip";
	print("Saving frame "+currentFrame+": "+fileNAME);
	
	//Save mask
	//setSlice(1);  //Go to mask
	//roiManager("Show None");
	run("Make Substack...", "  slices=1");
	saveAs("Tiff",dirMASKS+currentPos+"/"+currentChannel+"/"+fileNAME);
	run("Close");
	//close(fileNAME);
	//Save roi
	if (!isOpen("ROI Manager")) {
	print("ROI Manager not open at save frame...");	
		}
	
	
	selectWindow("ROI Manager");
	
	roiManager("Save", dirROIS+currentPos+"/"+currentChannel+"/"+roiNAME);   
	
	//Return to original slice
	//setSlice(this_slice);   
	roiManager("deselect");	
	//roiManager("Show All");
	
	setBatchMode(false);
	selectImage(imageID);
	
}

/***********************************/
function uJCloseFrame(){
	
	run("Close All");
	
	selectWindow("ROI Manager");
	
	run("Close");
			
			
}

exit;