//BUG: not working in imageJ.  When data is imported to ResultsTable, strings are imported as NaN


print("\n> Data_2_Overlay");

setupFile=getArgument();
print("setupFile:"+setupFile);

args=split(getArgument(),",");
Array.print(args);

setupFile=args[0];
thisVariable=args[1];
maskChannel=args[2];
bgChannel=args[3];
lutFile=args[4];  

print("bgChannel: "+bgChannel);
print("Variable: "+thisVariable);


/********************************************************************/
// 				USER-DEFINED PARAMETERS
/********************************************************************/

min_val=-0.75;  //LUT[0]
max_val=0.75;  //LUT[255]

shift_x=-4;
shift_y=0;

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
dirROOT=getParam(setupFile, "dirROOT");
pathDATA=getParam(setupFile, "pathDATA");
dirROIS=pathDATA+getParam(setupFile, "dirROIS"); //Input
dirDATA=pathDATA+getParam(setupFile, "dirDATA"); //Input
dirMASKS=pathDATA+getParam(setupFile, "dirMASKS"); //Input
dirTIF=pathDATA+getParam(setupFile, "dirTIF"); //Input

dirBG="";
strBG="";
if(bgChannel!=""){
	dirBG=dirTIF;
	bgChannel=bgChannel;
	strBG=bgChannel+"+";
}

/********************************************************************/
//							MACRO
/********************************************************************/

lutPath=dirROOT+"uJ_src/luts/"+lutFile;
roiColors = loadColorsfromLut(lutPath);

setBatchMode(true); 
for(p=0; p<numPos; p++){
	thisPos=list_pos[p];  
	
	File.makeDirectory(dirTIF+"/"+thisPos+"/"+strBG+""+thisVariable);  //Create output directory (level2)

	for(f=0; f<list_frames.length; f++){
		currentFrame=list_frames[f];
		
		thisFrame=IJ.pad(currentFrame,3);
		
		fileNameData=expeLabel+"_"+thisPos+"_"+thisFrame+".txt";
		fileNameMask=expeLabel+"_"+thisPos+"_"+thisFrame+".tif";
		fileNameRoi=expeLabel+"_"+thisPos+"_"+thisFrame+".zip";
		fileNameBG=expeLabel+"_"+thisPos+"_"+bgChannel+"_"+thisFrame+".tif";
		
		print(dirDATA+""+thisVariable+"/"+thisPos+"/"+fileNameData);

		existData=File.exists(dirDATA+""+thisVariable+"/"+thisPos+"/"+fileNameData);
		existMask=File.exists(dirMASKS+""+thisPos+"/"+maskChannel+"/"+fileNameMask);
		existRoi=File.exists(dirROIS+""+thisPos+"/"+maskChannel+"/"+fileNameRoi);
		if((existData*existMask*existRoi)==1){
			//thisPos=substring(list_rois[i],indexOf(list_rois[i], "_xy")+1,indexOf(list_rois[i], "_xy")+5);
			//thisFrame=substring(list_rois[i],indexOf(list_rois[i], ".zip")-3,indexOf(list_rois[i], ".zip"));	
			print("\n"+i+": Overlaying "+" pos="+thisPos+", frame="+thisFrame);

			//Open Background
			if(dirBG!=""){
				print(dirBG+""+thisPos+"/"+bgChannel+"/"+fileNameBG);
				open(dirBG+""+thisPos+"/"+bgChannel+"/"+fileNameBG);
				run("Enhance Contrast...", "saturated=0.5");
				makeRectangle(0, 0, getWidth(), getHeight()); 
  				run("Copy"); 
  				makeRectangle(shift_x, shift_y, getWidth(), getHeight()); 
  				run("Paste"); 
				rename("BG");
			}else{ //Open mask
				open(dirMASKS+""+thisPos+"/"+maskChannel+"/"+fileNameMask);
				rename("BG");
			}
	
			//Loads ROI Manager
			roiManager("reset");
			roiManager("Open",dirROIS+""+thisPos+"/"+maskChannel+"/"+fileNameRoi);
	
			roiManager("Show All");
	
			//Imports Data into Results Table
			  lineseparator = "\n";
			  cellseparator = ",\t";

			  // copies the whole table to an array of lines
			  print("Importing data from "+dirDATA+""+thisVariable+"/"+thisPos+"/"+fileNameData);
			  lines=split(File.openAsString(dirDATA+""+thisVariable+"/"+thisPos+"/"+fileNameData), lineseparator);
			  if (lines.length==0) return;
			  path = File.directory + File.name;

			  // get the columns headers
			  labels=split(lines[0], cellseparator);
			  if (labels.length==1)
				 exit("This is not a tab or comma delimited text file.");
			  if (labels[0]==" ")
				 k=1; // it is an ImageJ Results table, skip first column
			  else
				 k=0; // it is not a Results table, load all columns

			  // is this a Results table?
			  if (k==1 || lines.length<2)
				 {importResults(); exit;}
			  items = split(lines[1]);
			  nonNumeric = false;
			  for (i=0; i<items.length; i++)
				 if (isNaN(parseFloat(items[i]))) nonNumeric=true;
				  importResults();


			//Now loop results table 
			selectWindow("Results");
			for (row = 0; row < nResults; row++) { 
		
				val=getResult(thisVariable, row);
				roi_label=getResultString("roi_label", row);
				//print(roi_label+" : "+val);
		
				//Now loop roi manager and color roi based on val
				n = roiManager("count");
				for (iRoi=0; iRoi<n; iRoi++) {
					  roiManager("select", iRoi);
			  
					  if(Roi.getName==roi_label){
							//icolor=floor(roiColors.length*random()); 
					
							icolor=floor(roiColors.length*((val-min_val)/(max_val-min_val))); 
							
							if(icolor<1)
								icolor=1;
							if(icolor>255)
								icolor=255;
		
							strokeColor="#FF"+roiColors[icolor];
							if(dirBG!="")
								fillColor="#99"+roiColors[icolor];
							else
								fillColor="#FF"+roiColors[icolor];
							roiManager("Set Fill Color", fillColor);
					  }
				}
		
		
	
		
			}
	
			run("Flatten");
			fileNameOverlay=expeLabel+"_"+strBG+""+thisVariable+"_"+thisPos+"_"+thisFrame+".tif";
			print("Saving: "+dirTIF+""+thisPos+"/"+strBG+""+thisVariable+"/"+fileNameOverlay);
			
			saveAs("Tiff",dirTIF+"/"+thisPos+"/"+strBG+""+thisVariable+"/"+fileNameOverlay);
			
			close();
			selectWindow("BG");
			close();
		} //If file exists
    }  // For frames
    
} //For pos
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
        	if (endsWith(list[i], ".zip"))
  				str_rois = str_rois+","+list[i];
        }
     }
     
     return str_rois;
  }
  
  
/********************************************************************/
  
  function importResults() {
      for (j=k; j<labels.length; j++)
        setResult(labels[j],0,0);
     run("Clear Results");
     if (k==1)
         call("ij.plugin.filter.Analyzer.setDefaultHeadings");
     for (i=1; i<lines.length; i++) {
     	items=split(lines[i], cellseparator);
        for (j=k; j<items.length; j++)
           setResult(labels[j],i-1,items[j]);
     }
     updateResults();
  }
 
  function importTable() {
      name = "["+File.name+"]";
      if (!isOpen(File.name))
          run("New... ", "name="+name+" type=Table");
      f = name;
      print(f, "\\Headings:"+lines[0]);
      for (i=1; i<lines.length; i++){
         print(f, lines[i]);
      }
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
  
print("Exit");
exit;