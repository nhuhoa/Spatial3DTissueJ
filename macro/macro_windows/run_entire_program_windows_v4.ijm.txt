//======================================================================================
//// Installation
//// See at https://github.com/nhuhoa/Spatial3DTissueJ/tree/master/ImageJ_plugins_jar/list_plugin_readme.txt
//// First, you can copy these 3 folders into Fiji/plugins/ or ImageJ/plugins/ and quick check if there is any duplicate plugin, ex: 3DViewerXXXvXXX//// .jar, Fiji_PluginXXX.jar. You can keep only one version for each plugin. 
//// To update plugins, you can delele the folder 3D_suite/ and replace by the most updated folder from my github, most updated plugin: Spatial3DTissueJ_v22_windows.jar

//======================================================================================
//// Preparation
//// Create a folder for each islet, ex: H1536_islet1/ and put composite image of a given islet into this folder. 
//// Copy the color mapping file into this folder H1536_islet1/: cell_type_colormap.lut
When you download github files, you have this file in the github folder: https://github.com/nhuhoa/Spatial3DTissueJ/tree/master/macro/macro_windows/cell_type_colormap.lut
//// Noted: cell_type_colormap.lut is just a color map for visualization, if macro has an issue with this file, you can just remove the line 

//======================================================================================
//// How to run a macro here
//// I divide macro into many steps, you can run one by one step here
//// And then when you don't see any error, you can try to run entire pipeline
//// I tested macro in Windows systems and for tissue H1536_islet1





//======================================================================================
//// STEP1: Splitting composite images into different channels
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";

composite_image_fn="H1536_islet1.tiff";

//// Gajana, Marcus
////Need to choose an option split channels from biocimage tools automatic popup window
//// in order to split images into different channels.
//// I don't know how to write macro in order to split this composite image yet
open(dir+composite_image_fn); 



//// From this part, I can run all steps in one run
//// Noted: I name the images as: "C1-delta.tif", "C2-beta.tif", "C4-dapi.tif", "C3-alpha.tif"
//// Because in the results images, the values of intensity are: delta-1, beta-2, and alpha-3
//// just a convention, easy to remember it. 

composite_image_fn="H1536_islet1.tiff";
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
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


//======================================================================================
////  STEP21: Hysteresis Threshold for background cut off 

////Processing BETA marker channels, cut off background, background and intra tissue environment=0
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
low_thrs=70;
high_thrs=100;
observed_marker_fn1="C2-beta.tif"
open(dir + observed_marker_fn1);  //open(dir+"C2-beta.tif");
observed_marker1 = substring(observed_marker_fn1,0,lastIndexOf(observed_marker_fn1,".tif"));
print("Removing background from image: "+observed_marker_fn1)

selectWindow(observed_marker_fn1);
run("Duplicate...", "duplicate title=["+observed_marker1+"_hysteresis_thresh]");
selectWindow(observed_marker1+"_hysteresis_thresh");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
////if you have same microscopy setting, you may need to define it once and reuse the set of parameters many times
run("3D Hysteresis Thresholding", "high="+high_thrs+" low="+low_thrs+" connectivity"); 

selectWindow(observed_marker_fn1);
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
imageCalculator("AND create stack", observed_marker_fn1, observed_marker1+"_hysteresis_thresh");
selectWindow("Result of "+observed_marker_fn1);
saveAs("Tiff", dir+observed_marker1+"_filtered.tif");
selectWindow(observed_marker_fn1);
saveAs("Tiff", dir+observed_marker_fn1);
run("Close All"); 


//======================================================================================
////  STEP22: Hysteresis Threshold for background cut off 

////Processing DELTA marker channels, cut off background, background and intra tissue environment=0
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
low_thrs=30;
high_thrs=40;
observed_marker_fn2="C1-delta.tif";

open(dir + observed_marker_fn2);
observed_marker2 = substring(observed_marker_fn2,0,lastIndexOf(observed_marker_fn2,".tif"));
print("Removing background from image: "+observed_marker_fn2)

selectWindow(observed_marker_fn2);
run("Duplicate...", "duplicate title=["+observed_marker2+"_hysteresis_thresh]");
selectWindow(observed_marker2+"_hysteresis_thresh");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
////if you have same microscopy setting, you may need to define it once and reuse the set of parameters many times
//// delta cells marker staining usually lower intensity compared to beta cells, so the threshold here is also lower
run("3D Hysteresis Thresholding", "high="+high_thrs+" low="+low_thrs+" connectivity"); 

selectWindow(observed_marker_fn2);
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
imageCalculator("AND create stack", observed_marker_fn2, observed_marker2+"_hysteresis_thresh");
selectWindow("Result of "+observed_marker_fn2);
saveAs("Tiff", dir+observed_marker2+"_filtered.tif");
selectWindow(observed_marker_fn2);
saveAs("Tiff", dir+observed_marker_fn2);
run("Close All"); 


//======================================================================================
////  STEP23: Hysteresis Threshold for background cut off 

////Processing ALPHA marker channels, cut off background, background and intra tissue environment=0
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
low_thrs=60;
high_thrs=90;
observed_marker_fn3="C3-alpha.tif"
open(dir + observed_marker_fn3);  // open(dir+"C3-alpha.tif");
observed_marker3 = substring(observed_marker_fn3,0,lastIndexOf(observed_marker_fn3,".tif"));
print("Removing background from image: "+observed_marker_fn3)
selectWindow(observed_marker_fn3);
run("Duplicate...", "duplicate title=["+observed_marker3+"_hysteresis_thresh]");
selectWindow(observed_marker3+"_hysteresis_thresh");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
////if you have same microscopy setting, you may need to define it once and reuse the set of parameters many times
//// alpha cells marker staining usually sightly lower intensity compared to beta cells, so the threshold here is also lower, also from my image observation
run("3D Hysteresis Thresholding", "high="+high_thrs+" low="+low_thrs+" connectivity"); 

selectWindow(observed_marker_fn3);
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
imageCalculator("AND create stack", observed_marker_fn3, observed_marker3+"_hysteresis_thresh");
selectWindow("Result of "+observed_marker_fn3);
saveAs("Tiff", dir+observed_marker3+"_filtered.tif");
selectWindow(observed_marker_fn3);
saveAs("Tiff", dir+observed_marker_fn3);
run("Close All"); 





//======================================================================================
//// STEP3: Nucleus segmentation
//// In 16 bits images, the intensity values are from 0 to 65536, so above 29000 can consider as signal intensity for maximal local intensity values
//// Similar to in 8 bits, ex: values above 120 from the image intensity range of [0, 255]
//// You can manually enter a threshold here or use automatic estimation. I suggest to use automatic seeds threshold estimation because some nucleus images have very low range of intensity, you can not detect object with high threshold above. 
//// For automatic threshold mode, program will calculate the mean value of image, and define the seeds threshold = mean + sd, I define a sd here is 500
//// You can run automatic mode and observe the results, if as not you expected --> increase or decrease threshold and use manual threshold mode. See the log file to have an idea here. 
//// It take ~ 3 mins for this step

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
//// It take ~ 1 mins for this step
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
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
//// It take ~ 3 mins for this step for one small testing image 
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


//// Not important, just for visualization
selectWindow("composite_label.tif");
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
open(dir+"cell_type_colormap.lut"); // I provide a color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"composite_label.tif"); //color map for better visualization
run("Close All"); 










//======================================================================================
//// STEP6: cell type detection

dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
percent_marker_coverage=0.2;

open(dir+"C1-delta_filtered.tif");
open(dir+"C2-beta_filtered.tif");
open(dir+"C3-alpha_filtered.tif");
open(dir+"dapi-seg.tif");
open(dir+"dapi-seg-wat.tif");
open(dir+"composite_label.tif");
run("CELL TYPE DETECTION", "save_dir=["+dir+"] slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha_filtered.tif filtered_0=C2-beta_filtered.tif filtered_1=C1-delta_filtered.tif region_observed=[WATERSHED REGION] percent_marker_coverage="+percent_marker_coverage+" min_distance=0 max_distance_inside=2 max_distance_outside=3 save show");
/// Just for better visualization and publication, demo, not important computation part
selectWindow("TYPE_NUC");
open(dir+"cell_type_colormap.lut"); // I provide a color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"TYPE_NUC.tif"); 

selectWindow("composite_label.tif");
open(dir+"cell_type_colormap.lut"); // I provide a color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"composite_label.tif"); 

selectWindow("TYPE_WAT");
open(dir+"cell_type_colormap.lut"); // I provide a color map for good visualization of different channels here
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", dir+"TYPE_WAT.tif"); 
selectWindow("Unlabelled");
saveAs("Tiff", dir+"Unlabelled.tif"); 
run("Close All");

//// NOTED: save all results into Nuclei.zip and Regions.zip,
//// for all downstream analysis, data will be read from these files, 
//// so please point dir variable to this file as my example here, see more at Log window  
//// the easy way is keeping all results for each islet in one folder.  


//// Different cell type detection methods, just testing here
//// percent_marker_coverage=0.25: at least 20% of a given cell zone area should be covered by a marker  
//// Uncomment the part below if you want to use other cell type detection methods
//percent_marker_coverage=0.25;
//run("CELL TYPE DETECTION", "save_dir=["+dir+"] slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha_filtered.tif filtered_0=C2-beta_filtered.tif filtered_1=C1-delta_filtered.tif region_observed=[INSIDE_OUTSIDE NUCLEUS] percent_marker_coverage="+percent_marker_coverage+" min_distance=0 max_distance_inside=2 max_distance_outside=3 save show");
//selectWindow("TYPE_WAT");
//open(dir+"cell_type_colormap.lut"); // I provide a predefined color map for good visualization of different channels here
//run("Enhance Contrast", "saturated=0.35");
//saveAs("Tiff", dir+"TYPE_WAT.tif"); 

//selectWindow("TYPE_NUC");
//open(dir+"cell_type_colormap.lut"); // I provide a predefined color map for good visualization of different channels here
//run("Enhance Contrast", "saturated=0.35");
//saveAs("Tiff", dir+"TYPE_NUC.tif"); 
//run("Close All"); 


//======================================================================================
//// STEP7: cell to cell interaction analysis
//// NOTED: save all results into Nuclei.zip and Regions.zip,
//// for all downstream analysis, data will be read from these files, 
//// so please point dir variable to this file as my example here, see more at Log window  
//// the easy way is keeping all results for each islet in one folder.  
//// ~3 mins for this step

dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
//// for all downstream analysis, data will be read from Nuclei.zip and Regions.zip files, 
open(dir+"dapi-seg-wat.tif");
open(dir+"C4-dapi.tif");
run("CELLS INTERACTION ANALYSIS", "input="+dir+" save="+dir+" watershed_image=dapi-seg-wat.tif cell-cell_contact_analysis layer_contact histogram_contact cells_network_visualization save="+dir+"raw_data_stat/Log_Cellular_Interaction_Analysis.txt");
run("CELLS LAYER CONTACT", "input="+dir+" watershed_image=dapi-seg-wat.tif dapi_raw_image=C4-dapi.tif target_cell=BETA source_cell=DELTA");
selectWindow("LAYER_DISTANCE_BETA_DELTA.tif");
run("Fire");
saveAs("Tiff", dir+"layer_distance/LAYER_DISTANCE_BETA_DELTA.tif");
run("CELLS LAYER CONTACT", "input="+dir+" watershed_image=dapi-seg-wat.tif dapi_raw_image=C4-dapi.tif target_cell=BETA source_cell=ALPHA");
selectWindow("LAYER_DISTANCE_BETA_ALPHA.tif");
run("Fire");
saveAs("Tiff", dir+"layer_distance/LAYER_DISTANCE_BETA_ALPHA.tif");
run("Close All");


//======================================================================================
//// STEP8: delta cell protrusion
//// NOTED: save all results into Nuclei.zip and Regions.zip,
//// for all downstream analysis, data will be read from these files, 
//// so please point dir variable to this file as my example here, see more at Log window  
//// ~3 mins for this step
//// Noted: if you use graphical interface for this function, and tick on visualization options
//// too many windows of images will pop up. So I use default setting without visualization of delta cell protrusion images display. 
dir="C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/";
open(dir+"dapi-seg-wat.tif");
open(dir+"C4-dapi.tif");
run("DELTA CELL PROTRUSION", "save_dir=["+dir+"] watershed_image=dapi-seg-wat.tif dapi_image=C4-dapi.tif min_distance=0 max_distance=35 region_value=0 target_cell=BETA protrusion_cell=DELTA range_observe=ALL_CELLS_PROTRUSION region_observe=CELL_REGION save=["+dir+"/protrusion_DELTA/Log_protrusion_DELTA_ALL_CELLS_PROTRUSION.txt]");
run("Close All");



