/*
run("Open...");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/data_raw/xy34/GFP sort");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/data_raw/xy34/DsRed sort");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/data_raw/xy34/DIC sort");

selectWindow("DIC");
run("Enhance Contrast...", "saturated=0.3 process_all");
selectWindow("GFP");
run("Enhance Contrast...", "saturated=0.3 normalize process_all");
selectWindow("DsRed");
run("Enhance Contrast...", "saturated=0.3 normalize process_all");

selectWindow("GFP");
run("Merge Channels...", "c2=GFP");
saveAs("Tiff", "/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/data_raw/xy34/HT13-CAZ_GFP.tif");

selectWindow("DsRed");
run("Merge Channels...", "c1=DsRed");
saveAs("Tiff", "/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/data_raw/xy34/HT13-CAZ_DsRed.tif");

selectWindow("DIC");
run("RGB Color");
saveAs("Tiff", "/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/data_raw/xy34/HT13-CAZ_DIC.tif");
*/

//run("Combine...", "stack1=HT13-CAZ_DIC stack2=HT13-CAZ_GFP.tif");
//selectWindow("HT13-CAZ_GFP.tif");

/*
selectWindow("DIC");
run("RGB Color");
selectWindow("RGB");
selectWindow("DIC");
run("Combine...");
selectWindow("RGB");
selectWindow("DIC");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/images_cells_tracked/xy34/DsRed sort");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/images_cells_tracked/xy34/GFP sort");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/images_cells_tracked/xy34/RelInt sort");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/images_cells_tracked/xy34/Tracking sort");
run("Image Sequence...", "open=/Users/ESB/RPM_SYNC/RPM_Work/uJ_master/data_sample/HT13-CAZ/images_cells_tracked/xy34/divisions sort");
close();
selectWindow("Tracking");
makeRectangle(165, 210, 747, 592);
makeRectangle(167, 259, 745, 543);
selectWindow("RelInt");
*/

makeRectangle(180, 260, 640, 512);
run("Crop");