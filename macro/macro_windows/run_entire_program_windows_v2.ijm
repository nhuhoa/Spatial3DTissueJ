//// Installation

//// Preparation
//// Create a folder for each islet, and put composite image of a given islet into this folder. 
//// Copy the color mapping file into this folder: cell_type_colormap.lut

//======================================================================================
//// STEP1: Splitting composite images into different channels
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
composite_image_fn="H1536_islet1.tiff";

open(dir+composite_image_fn);

selectWindow(composite_image_fn+" - C=0");
saveAs("Tiff", dir+"C1-delta.tif");

selectWindow(composite_image_fn+" - C=1");
saveAs("Tiff", dir+"C2-beta.tif");

selectWindow(composite_image_fn+" - C=2");
saveAs("Tiff", dir+"C4-dapi.tif");

selectWindow(composite_image_fn+" - C=3");
saveAs("Tiff", dir+"C3-alpha.tif");
run("Close All"); 





//======================================================================================
////  STEP2: Hysteresis Threshold for background cut off 
//// I know the way we define the threshold is a bit manual here 
//// but if you have same microscopy setting, you may need to define it once 
//// and reuse the set of parameters many times

dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";

open(dir+"C1-delta.tif");
open(dir+"C2-beta.tif");
open(dir+"C3-alpha.tif");


////Processing BETA marker channels, cut off background, background and intra tissue environment=0
low_thrs=70;
high_thrs=100;
observed_marker_fn="C2-beta.tif"
observed_marker = substring(observed_marker_fn,0,lastIndexOf(observed_marker_fn,".tif"));
print("Removing background from image: "+observed_marker_fn)

selectWindow(observed_marker_fn);
run("Duplicate...", "title=["+observed_marker+"_hysteresis_thresh]");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
////if you have same microscopy setting, you may need to define it once and reuse the set of parameters many times
run("3D Hysteresis Thresholding", "high="+high_thrs+" low="+low_thrs+" connectivity"); 

selectWindow(observed_marker_fn);
run("Median 3D...", "x=2 y=2 z=1");

imageCalculator("AND create stack", observed_marker_fn, observed_marker+"_hysteresis_thresh");
selectWindow("Result of "+observed_marker_fn);
saveAs("Tiff", dir+observed_marker+"_filtered.tif");
run("Close All"); 



////Processing DELTA marker channels, cut off background, background and intra tissue environment=0
low_thrs=30;
high_thrs=40;
observed_marker_fn="C1-delta.tif"
observed_marker = substring(observed_marker_fn,0,lastIndexOf(observed_marker_fn,".tif"));
print("Removing background from image: "+observed_marker_fn)

selectWindow(observed_marker_fn);
run("Duplicate...", "title=["+observed_marker+"_hysteresis_thresh]");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
////if you have same microscopy setting, you may need to define it once and reuse the set of parameters many times
//// delta cells marker staining usually lower intensity compared to beta cells, so the threshold here is also lower
run("3D Hysteresis Thresholding", "high="+high_thrs+" low="+low_thrs+" connectivity"); 

selectWindow(observed_marker_fn);
run("Median 3D...", "x=2 y=2 z=1");
imageCalculator("AND create stack", observed_marker_fn, observed_marker+"_hysteresis_thresh");
selectWindow("Result of "+observed_marker_fn);
saveAs("Tiff", dir+observed_marker+"_filtered.tif");
run("Close All"); 



////Processing ALPHA marker channels, cut off background, background and intra tissue environment=0
low_thrs=60;
high_thrs=90;
observed_marker_fn="C1-delta.tif"
observed_marker = substring(observed_marker_fn,0,lastIndexOf(observed_marker_fn,".tif"));
print("Removing background from image: "+observed_marker_fn)
selectWindow(observed_marker_fn);
run("Duplicate...", "title=["+observed_marker+"_hysteresis_thresh]");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
////if you have same microscopy setting, you may need to define it once and reuse the set of parameters many times
//// alpha cells marker staining usually sightly lower intensity compared to beta cells, so the threshold here is also lower, also from my image observation
run("3D Hysteresis Thresholding", "high="+high_thrs+" low="+low_thrs+" connectivity"); 

selectWindow(observed_marker_fn);
run("Median 3D...", "x=2 y=2 z=1");
imageCalculator("AND create stack", observed_marker_fn, observed_marker+"_hysteresis_thresh");
selectWindow("Result of "+observed_marker_fn);
saveAs("Tiff", dir+observed_marker+"_filtered.tif");
run("Close All"); 





//======================================================================================
//// STEP3: Nucleus segmentation
//// In 16 bits images, the intensity values are from 0 to 65536, so above 29000 can consider as signal intensity for maximal local intensity values
//// Similar to in 8 bits, ex: values above 120 from the image intensity range of [0, 255]
//// You can manually enter a threshold here or use automatic estimation. I suggest to use automatic seeds threshold estimation because some nucleus images have very low range of intensity, you can not detect object with high threshold above. 
//// For automatic threshold mode, program will calculate the mean value of image, and define the seeds threshold = mean + sd, I define a sd here is 500
//// You can run automatic mode and observe the results, if as not you expected --> increase or decrease threshold and use manual threshold mode. See the log file to have an idea here. 


small_nucleus_diameter_thrs=15;
large_nucleus_diameter_thrs=36;
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";

open(dir+"C4-dapi.tif");

//// Using automatic seeds selection mode
run("NUCLEI SEGMENTATION", "save_dir=["+dir+"] dapi=C4-dapi.tif maximal="+large_nucleus_diameter_thrs+" minimal="+small_nucleus_diameter_thrs+" seed=29000 automatic");
selectWindow("dapi-seg.tif"); //automatic save results into folder
run("Enhance Contrast", "saturated=0.35");
run("3-3-2 RGB"); //color map for better visualization

run("3D Draw Rois", "raw=C4-dapi seg=dapi-seg");
selectWindow("DUP_C4-dapi.tif"); //TO DO DOUBLE CHECK
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"DUP_C4-dapi_demo_only.tif"); // very useful for validation, and quick check accuracy of nucleus object detection
run("Close All"); 







//======================================================================================
//// STEP4: cell zone estimation
//// Estimating a cell zone for each cell, from each nucleus, extend into space a radius R, ex: R=5 here, then obtained region will be a cell zone. I observe image of dapi and markers to define a max radius for region growing. 
//// Then later, program will look into each cell zone and check the amount of marker that cover this cell zone. 


max_radius_extension_wat=5; //change this threshold if your cell zone area in the image is smaller, depend on image resolution, in general 3 to 5 is a good one. 
open(dir+"C4-dapi.tif");
open(dir+"dapi-seg.tif");
run("CELL ZONE ESTIMATION", "save_dir=["+dir+"] nuclei=dapi-seg.tif radius_max="+max_radius_extension_wat+" save");

selectWindow("dapi-seg-wat.tif");
run("3-3-2 RGB");
saveAs("Tiff", dir+"dapi-seg-wat.tif"); //color map for better visualization
run("Close All"); 






//======================================================================================
//// STEP5: Marker clustering

dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
open(dir+"C1-delta_filtered.tif");
open(dir+"C2-beta_filtered.tif");
open(dir+"C3-alpha_filtered.tif");

//// Note: the position of delta, alpha, beta is not important here, just need 3 channels. 
run("SLIC 3D CLUSTERING", "save_dir=["+dir+"] marker_1=C1-delta_filtered.tif marker_2=C2-beta_filtered.tif marker_3=C3-alpha_filtered.tif marker_4=*None* min_size=80 max_size=300 nb_iterations=10");
selectWindow("creator");
close();

selectWindow("combined_markers_SLIC.tif");
run("SLIC 3D 3 CHANNELS", "save_dir=["+dir+"] alpha_image=C3-alpha_filtered.tif beta_image=C2-beta_filtered.tif delta_image=C1-delta_filtered.tif combined_slic_image=combined_markers_SLIC.tif intensity_threshold=0");











//======================================================================================
//// STEP6: cell type detection

dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";

run("CELL TYPE DETECTION", "save_dir=["+dir+"] slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha_filtered.tif filtered_0=C2-beta_filtered.tif filtered_1=C1-delta_filtered.tif region_observed=[WATERSHED REGION] percent_marker_coverage=0.2 min_distance=0 max_distance_inside=2 max_distance_outside=3 save show");
selectWindow("TYPE_NUC");


/// Just for better visualization and publication, demo, not important computation part
selectWindow("TYPE_NUC");
open(dir+"cell_type_colormap.lut"); // I provide a color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"TYPE_NUC.tif"); 

selectWindow("composite_label.tif");
open(dir+"cell_type_colormap.lut"); // I provide a color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"composite_label.tif"); 


//// Different cell type detection methods, just testing here
//// percent_marker_coverage=0.25: at least 20% of a given cell zone area should be covered by a marker  
percent_marker_coverage=0.25;
run("CELL TYPE DETECTION", "save_dir=["+dir+"] slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha_filtered.tif filtered_0=C2-beta_filtered.tif filtered_1=C1-delta_filtered.tif region_observed=[INSIDE_OUTSIDE NUCLEUS] percent_marker_coverage="+percent_marker_coverage+" min_distance=0 max_distance_inside=1 max_distance_outside=3.000 save show");

selectWindow("TYPE_WAT");
open(dir+"cell_type_colormap.lut"); // I provide a predefined color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"TYPE_WAT.tif"); 

selectWindow("TYPE_NUC");
open(dir+"cell_type_colormap.lut"); // I provide a predefined color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"TYPE_NUC.tif"); 
run("Close All"); 


//======================================================================================
//// STEP7: cell to cell interaction analysis

run("CELLS INTERACTION ANALYSIS", "input="+dir+" save="+dir+" watershed_image=dapi-seg-wat.tif cell-cell_contact_analysis layer_contact histogram_contact cells_network_visualization save="+dir+"raw_data_stat/Log_Cellular_Interaction_Analysis.txt");
run("CELLS LAYER CONTACT", "input="+dir+" watershed_image=dapi-seg-wat.tif dapi_raw_image=C4-dapi.tif target_cell=BETA source_cell=DELTA");
selectWindow("LAYER_DISTANCE_BETA_DELTA.tif");
run("Fire");
saveAs("Tiff", save_dir_ct+"layer_distance/LAYER_DISTANCE_BETA_DELTA.tif");
run("CELLS LAYER CONTACT", "input="+dir+" watershed_image=dapi-seg-wat.tif dapi_raw_image=C4-dapi.tif target_cell=BETA source_cell=ALPHA");
selectWindow("LAYER_DISTANCE_BETA_DELTA.tif");
run("Fire");
saveAs("Tiff", save_dir_ct+"layer_distance/LAYER_DISTANCE_BETA_DELTA.tif");



//======================================================================================
//// STEP8: delta cell protrusion
open(dir+"dapi-seg-wat.tif");
open(dir+"C4-dapi.tif");
run("DELTA CELL PROTRUSION", "save_dir=["+dir+"] watershed_image=dapi-seg-wat.tif dapi_image=C4-dapi.tif min_distance=0 max_distance=35 region_value=0 target_cell=BETA protrusion_cell=DELTA range_observe=ALL_CELLS_PROTRUSION region_observe=CELL_REGION save=["+dir+"/protrusion_DELTA/Log_protrusion_DELTA_ALL_CELLS_PROTRUSION.txt]");




