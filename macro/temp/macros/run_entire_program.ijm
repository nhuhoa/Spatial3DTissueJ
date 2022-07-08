

dir="/Users/htran/Documents/storage_tmp/Spatial3DTissueJ/testing_images/";
//dir="yourdir/testing_images/";


print("Splitting composite images: ");
save_dir=dir+"filtered/";
open(dir+"raw_islet/islet_composite.tif");
run("PRE-FILTERING", "save="+save_dir+" composite=islet_composite.tif label_1=BETA label_2=ALPHA label_3=DELTA label_4=DAPI");
run("Close All"); 


dir="/Users/htran/Documents/storage_tmp/Spatial3DTissueJ/testing_images/";
save_dir_ct=dir+"cell_type/";
print("Nuc segmentation");
open(dir+"filtered/C4-dapi.tif");
run("NUCLEI SEGMENTATION", "save="+save_dir_ct+" dapi=C4-dapi.tif maximal=28 minimal=18 seed=29000 automatic");
selectWindow("dapi-seg.tif");
run("3-3-2 RGB");
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", save_dir_ct+"dapi-seg.tif");

print("Cell zone estimation");
run("CELL ZONE ESTIMATION", "save="+save_dir_ct+" nuclei=dapi-seg.tif radius_max=3 save_0");
selectWindow("dapi-seg-wat.tif");
run("3-3-2 RGB");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", save_dir_ct+"dapi-seg-wat.tif");
//run("Close All"); 


print("Markers clustering SLIC");
dir="/Users/htran/Documents/storage_tmp/Spatial3DTissueJ/testing_images/";
save_dir=dir+"filtered/";
save_dir_ct=dir+"cell_type/";
open(save_dir+"C1-delta.tif");
open(save_dir+"C2-beta.tif");
open(save_dir+"C3-alpha.tif");

run("SLIC 3D CLUSTERING", "save="+save_dir_ct+" marker_1=C1-delta.tif marker_2=C2-beta.tif marker_3=C3-alpha.tif marker_4=*None* min_size=80 max_size=300 nb_iterations=10");
run("SLIC 3D 3 CHANNELS", "save="+save_dir_ct+" alpha_image=C3-alpha.tif beta_image=C2-beta.tif delta_image=C1-delta.tif combined_slic_image=combined_markers_SLIC.tif intensity_threshold=1");

//open(save_dir_ct+"dapi-seg-wat.tif");
//open(save_dir_ct+"dapi-seg.tif");
print("Cell type detection");
run("CELL TYPE DETECTION", "save="+save_dir_ct+" slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha.tif filtered_0=C2-beta.tif filtered_1=C1-delta.tif region_observed=[WATERSHED REGION] percent_marker_coverage=0.1 min_distance=0 max_distance_inside=1 max_distance_outside=3.000 save_0 show");
selectWindow("TYPE_NUC");
open(dir+"raw_islet/composite_label.lut");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", save_dir_ct+"TYPE_NUC.tif");
close();
selectWindow("TYPE_WAT");
open(dir+"raw_islet/composite_label.lut");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", save_dir_ct+"TYPE_WAT.tif");
close();

selectWindow("Unlabelled");
saveAs("Tiff", save_dir_ct+"Unlabelled.tif");
close();



run("CELLS INTERACTION ANALYSIS", "input="+save_dir_ct+" save="+save_dir_ct+" watershed_image=dapi-seg-wat.tif cell-cell_contact_analysis layer_contact histogram_contact cells_network_visualization save="+save_dir_ct+"raw_data_stat/Log_Cellular_Interaction_Analysis.txt");
run("CELLS LAYER CONTACT", "input="+save_dir_ct+" watershed_image=dapi-seg-wat.tif dapi_raw_image=C4-dapi.tif target_cell=BETA source_cell=DELTA");
selectWindow("LAYER_DISTANCE_BETA_DELTA.tif");
run("Fire");
saveAs("Tiff", save_dir_ct+"layer_distance/LAYER_DISTANCE_BETA_DELTA.tif");

run("CELLS LAYER CONTACT", "input="+save_dir_ct+" watershed_image=dapi-seg-wat.tif dapi_raw_image=C4-dapi.tif target_cell=BETA source_cell=ALPHA");
selectWindow("LAYER_DISTANCE_BETA_ALPHA.tif");
run("Fire");
saveAs("Tiff", save_dir_ct+"layer_distance/LAYER_DISTANCE_BETA_ALPHA.tif");

selectWindow("CellsNetwork.tif");
open(dir+"raw_islet/composite_label.lut");
run("Enhance Contrast", "saturated=0.35");
setMinAndMax(0, 148);
saveAs("Tiff", save_dir_ct+"raw_data_stat/CellsNetwork.tif");
close();
close();

run("RANDOM DISTRIBUTION", "save="+save_dir_ct+" watershed_image=dapi-seg-wat.tif observed_cell_type=ALPHA source_cell_type_observed=BETA nb_randoms=100 save="+save_dir_ct+"random_distributions_targetCT_ALPHA_sourceCT_BETA/Log_random_distribution.txt");
run("DELTA CELL PROTRUSION", "input="+save_dir_ct+" watershed_image=dapi-seg-wat.tif dapi_image=C4-dapi.tif min_distance=0 max_distance=35 region_value=0 target_cell=BETA protrusion_cell=DELTA range_observe=ALL_CELLS_PROTRUSION region_observe=CELL_REGION save="+save_dir_ct+"protrusion_DELTA/Log_protrusion_DELTA_ALL_CELLS_PROTRUSION.txt");



selectWindow("Log");
saveAs("Text", save_dir_ct+ "run_entire_program_LOG.txt");
run("Close All"); 


