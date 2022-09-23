print("> Pos-filler_2_N-Eclipses");

setupFile=getArgument();
print("   setupFile:"+setupFile);

/********************************************************************/
//					PARSE PARAMETERS
/********************************************************************/

//Parse user-defined Parameters
expeLabel=getParam(setupFile, "expeLabel");
list_pos=split(getParam(setupFile, "list_pos"),',');
list_frames=getParam(setupFile, "list_frames");
list_channel=split(getParam(setupFile, "list_channel"),',');
str_channel=split(getParam(setupFile, "str_channel"),',');
signal_file=getParam(setupFile, "signal_file");
pos_file=getParam(setupFile, "pos_file");

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
dirPrePosFix=pathDATA+"data_eclipse_pre_posfix/";
dirPosFix=pathDATA+"data_eclipse_posfix/";
File.makeDirectory(dirPosFix); 



/********************************************************************/
//							MACRO
/********************************************************************/

setBatchMode(true);
signal_file=pathDATA+""+signal_file;
pos_file=pathDATA+pos_file;


print(signal_file);

signal=File.openAsString(signal_file);
signal=split(signal,"\n");

pos_lines=File.openAsString(pos_file);

pos_lines=split(pos_lines,"\n");

folders=newArray();
traps=newArray();
for(pi=1;pi<pos_lines.length;pi++){
	this_line=split(pos_lines[pi],"=");
	folder=this_line[0];
	this_traps=this_line[1];
	this_traps=split(this_traps,",");
	traps=Array.concat(traps,this_traps);
	folders=Array.concat(folders,folder);
	
	}

//Array.print(traps);
traps=Array.sort(traps);
//Array.print(traps);

n_common=pos_lines.length-1;
commons=newArray();
for(i=0;i<traps.length;i++){
	it=traps[i];
	cont=1;
	for(j=i+1;j<traps.length;j++){
		jt=traps[j];
		if(it==jt){
			cont++;
			}
		else{
			j+=traps.length;
			}
	
	}

	if(cont==n_common){
		commons=Array.concat(commons,it);
		}
	
	}
print("Common traps are:");
Array.print(commons);

for(c=0;c<commons.length;c++){
//for(c=0;c<2;c++){
	common=commons[c];
	for(f=0;f<folders.length;f++){
	//for(f=0;f<1;f++){
		folder=folders[f];
		this_eclipse_dir=dirPosFix+folder+"/";
		File.makeDirectory(this_eclipse_dir); 
		traps=split(getParam(pos_file,folder),',');
		
		for(t=0;t<traps.length;t++){
			trap=traps[t];
			opentrap=IJ.pad(t+1,2);
			opentrap="xy"+opentrap;
			if(common==trap){
				for(chan=0;chan<list_channel.length;chan++){
//				for(chan=0;chan<1;chan++){
					channel=list_channel[chan];
					print(folder,trap,channel,opentrap);
					
					run("Image Sequence...", "open="+dirPrePosFix+folder+" file="+opentrap+channel+" sort");
					run("Image Sequence... ", "format=TIFF name="+expeLabel+"_"+trap+channel+"t start=1 digits=3 save="+this_eclipse_dir);

					close(); 
					}
				}
			}
		}
	}






/*
//Close all windows
 		while (nImages>0) { 
          selectImage(nImages); 
          close(); 
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