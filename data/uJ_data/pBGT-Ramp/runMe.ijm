//requires("1.51");
print("\\Clear"); 

/*******************/
// User-defined parameters

dirROOT="/Volumes/GoogleDrive/My Drive/SYNC_Projects/pBGT/ms/data/pBGT-Ramp/"; //getInfo("macro.filepath") ;
dirROOT = substring(dirROOT, 0, lastIndexOf(dirROOT, "/")+1);
print("dirROOT="+dirROOT);

pathDATA="/Volumes/GoogleDrive/My Drive/SYNC_Projects/pBGT/ms/data/pBGT-Ramp/"; //getInfo("macro.filepath") ; 
pathDATA = substring(pathDATA, 0, lastIndexOf(pathDATA, "/"));
pathDATA = substring(pathDATA, 0, lastIndexOf(pathDATA, "/")+1);
print("pathDATA="+pathDATA);

pathUJ="/Applications/Fiji.app/macros/"; //Path to uJ_src
setupFile=dirROOT+"setup.txt"; //Path to setup file

/****Auxiliar***/
// Run this macro if you have slpitted time-lapses with different set of traps
// Requires a position file defining absolute postitions for each time-lapse
//runMacro(pathUJ+"Pos-filler_2_N-Eclipses.ijm", setupFile);

	

/****Auxiliar***/
//  Run this macro if you have a splitted time-lapses
//  It merges them into one based on the signal file
//  *requires the signal file
//runMacro(pathUJ+"N_Eclipses_2_Raw.ijm", setupFile);


/*******************/
//Run Macro: Eclipse_2_Raw

//runMacro(pathUJ+"Eclipse_2_Raw.ijm", setupFile);

/*******************/
//Run Macro: Raw_2_TIF  
//runMacro(pathUJ+"Raw_2_TIF.ijm", setupFile); ///este no

//exec("ln -s "+pathDATA+"data_raw/ "+pathDATA+"data_tif");

/*******   Auxilliar ************/
//Run Macro: TIF_2_Montage  (DsRed+GFP)
//runMacro(pathUJ+"TIF_2_Montage.ijm", setupFile);
//runMacro(pathUJ+"TIF_2_Montage-SingleChan.ijm", setupFile);

/*******************/

/*******************/
//Run Macro: TIF_2_Segmentable

//runMacro(pathUJ+"TIF_2_Segmentable.ijm", setupFile);

/****************/
//Run Macro: Segmentable_2_DeepCell
 
//runMacro(pathUJ+"Segmentable_2_DeepCell.ijm", setupFile);

/*******************/
//
//Here we run DeepCell
//
/*******************/
//Run Macro: DeepCell_to_Raw_Maks
//runMacro(pathUJ+"DeepCell_2_RawMasks.ijm", setupFile);

/*******************/
//Run Macro: Correct Tracking

//bgChannel="GFP";
//args=setupFile+","+dirROOT+","+pathDATA+","+bgChannel;
//runMacro("uJ/RawMask_2_Masks.ijm", args);
//runMacro("uJ/TrackingCorrector.ijm", args);

/*******************/
//Run Macro: Overlay Data
/*
currentPos="xy46";
currentChannel="DIC";
currentLineages="all"; //analysis";
currentLayer=7;  //0:mask, 1,: highlight, 2:trackID, 3:Length, 4:GFP, 5:DsRed, 6:Division, 7:Dead

args=setupFile+","+dirROOT+","+pathDATA+","+currentPos+","+currentChannel+","+currentLineages+","+currentLayer;
runMacro("uJ/Data_2_Overlay.ijm", args);
*/

/*********************************************/
/*********************************************/
/*********************************************/
/*******************/
//Run Macro: TIF_2_Data

//runMacro(pathUJ+"uJ/TIF_2_Data.ijm", setupFile);


/*******************/
//Run Macro: TIF_2_Composite  (DsRed+GFP)
//runMacro(pathUJ+"TIF_2_Composite.ijm", setupFile);


/*******************/
//Run Macro: makeMovie
/*
list_pos=split("xy07"); //,xy08,xy09,xy12,xy13,xy17,xy26,xy27,xy32,xy34",",");  //Should be from setup.txt
thisChannel="relativeIntensity";
for(i=0; i<list_pos.length; i++){
	thisPos=list_pos[i];
	args=setupFile+","+thisPos+","+thisChannel;
	runMacro(pathUJ+"makeMovie.ijm", args);
}
*/

