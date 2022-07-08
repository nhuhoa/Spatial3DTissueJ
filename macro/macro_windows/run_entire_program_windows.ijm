

open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/H1536_islet1.tiff");
selectWindow("H1536_islet1.tiff - C=2");
selectWindow("H1536_islet1.tiff - C=1");
selectWindow("H1536_islet1.tiff - C=0");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C1-delta.tif");
selectWindow("H1536_islet1.tiff - C=1");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C2-beta.tif");
selectWindow("C1-delta.tif");
selectWindow("H1536_islet1.tiff - C=2");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C4-dapi.tif");
selectWindow("H1536_islet1.tiff - C=3");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C3-alpha.tif");
run("Median 3D...", "x=2 y=2 z=1");






selectWindow("C2-beta.tif");
run("Duplicate...", "duplicate");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
run("3D Hysteresis Thresholding", "high=100 low=70 connectivity");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C2-beta_hysteresis_thresholding.tif");

imageCalculator("AND create stack", "C2-beta.tif","C2-beta_hysteresis_thresholding.tif");
selectWindow("Result of C2-beta.tif");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C2-beta_filtered.tif");
close();
close();




open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C1-delta.tif");
run("Duplicate...", "duplicate");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
run("3D Hysteresis Thresholding", "high=40 low=30 connectivity");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C1-delta_hysteresis_threshold.tif");
imageCalculator("AND create stack", "C1-delta.tif","C1-delta_hysteresis_threshold.tif");
selectWindow("Result of C1-delta.tif");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C1-delta_filtered.tif");
close();
close();





open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C3-alpha.tif");
run("8-bit");
run("Median 3D...", "x=2 y=2 z=1");
run("Enhance Contrast", "saturated=0.35");
run("3D Hysteresis Thresholding", "high=90 low=60 connectivity");
run("Image Calculator...");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C3-alpha_hysteresis_thresholding.tif");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C3-alpha.tif");
run("8-bit");
imageCalculator("AND create stack", "C3-alpha.tif","C3-alpha_hysteresis_thresholding.tif");
selectWindow("Result of C3-alpha.tif");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C3-alpha_filtered.tif");
close();
close();
selectWindow("C3-alpha_hysteresis_thresholding.tif");
close();




run("NUCLEI SEGMENTATION", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1] dapi=C4-dapi.tif maximal=37 minimal=17 seed=29000 automatic");
run("Enhance Contrast", "saturated=0.35");
run("3-3-2 RGB");
run("Save");
run("3D Draw Rois", "raw=C4-dapi seg=dapi-seg");
run("Enhance Contrast", "saturated=0.35");
close();
run("3D Draw Rois", "raw=C4-dapi-1 seg=dapi-seg");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/DUP_C4-dapi_ROI_demo_only.tif");
close();
selectWindow("C4-dapi-1.tif");
selectWindow("dapi-seg.tif");
close();
close();
close();
selectWindow("C4-dapi.tif");
close();





open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C1-delta_filtered.tif");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C2-beta_filtered.tif");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C3-alpha_filtered.tif");
selectWindow("C2-beta_filtered.tif");
selectWindow("C1-delta_filtered.tif");
selectWindow("C3-alpha_filtered.tif");
makeLine(337, 404, 337, 404);
selectWindow("C1-delta_filtered.tif");
selectWindow("C2-beta_filtered.tif");
selectWindow("C1-delta_filtered.tif");
makeLine(607, 94, 607, 94);
selectWindow("C3-alpha_filtered.tif");
run("SLIC 3D CLUSTERING", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] marker_1=C1-delta_filtered.tif marker_2=C2-beta_filtered.tif marker_3=C3-alpha_filtered.tif marker_4=*None* min_size=80 max_size=300 nb_iterations=10");

selectWindow("creator");
close();
selectWindow("C2-beta_filtered.tif");
selectWindow("C1-delta_filtered.tif");
selectWindow("C3-alpha_filtered.tif");
selectWindow("combined_markers_SLIC.tif");

run("SLIC 3D 3 CHANNELS", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] alpha_image=C3-alpha_filtered.tif beta_image=C2-beta_filtered.tif delta_image=C1-delta_filtered.tif combined_slic_image=combined_markers_SLIC.tif intensity_threshold=0");

open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C4-dapi.tif");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/dapi-seg.tif");
run("CELL ZONE ESTIMATION", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] nuclei=dapi-seg.tif radius_max=3 save");
selectWindow("C4-dapi.tif");
selectWindow("dapi-seg.tif");
run("3-3-2 RGB");
selectWindow("dapi-seg-wat");
close();
close();






selectWindow("composite_label.tif");
run("CELL ZONE ESTIMATION", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] nuclei=dapi-seg.tif radius_max=5 save");
run("3-3-2 RGB");
selectWindow("dapi-seg.tif");
makeLine(596, 107, 596, 107);
selectWindow("dapi-seg-wat.tif");
selectWindow("dapi-seg.tif");





run("CELL TYPE DETECTION", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha_filtered.tif filtered_0=C2-beta_filtered.tif filtered_1=C1-delta_filtered.tif region_observed=[WATERSHED REGION] percent_marker_coverage=0.1 min_distance=0 max_distance_inside=1 max_distance_outside=3.000 save show");
selectWindow("TYPE_NUC");
run("Edit LUT...");
run("Enhance Contrast", "saturated=0.35");
selectWindow("Unlabelled");
selectWindow("TYPE_WAT");
selectWindow("TYPE_NUC");
selectWindow("TYPE_WAT");
selectWindow("Unlabelled");
selectWindow("TYPE_WAT");
selectWindow("Unlabelled");
selectWindow("TYPE_WAT");


open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/cell_type_colormap.lut");
run("Enhance Contrast", "saturated=0.35");
run("Histogram", "stack");
close();
run("3D Manager");
selectWindow("TYPE_NUC");
Ext.Manager3D_AddImage();
selectWindow("C2-beta_filtered.tif");
selectWindow("composite_label.tif");
run("Merge Channels...", "c1=C2-beta_filtered.tif c2=C3-alpha_filtered.tif c3=C1-delta_filtered.tif c5=C4-dapi.tif create keep");
run("Enhance Contrast", "saturated=0.35");
run("Enhance Contrast", "saturated=0.35");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/DUP_C4-dapi_ROI_demo_only.tif");
run("Enhance Contrast", "saturated=0.35");
selectWindow("composite_label.tif");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/cell_type_colormap.lut");
run("Enhance Contrast", "saturated=0.35");



// Different cell type detection methods, just testing here
run("CELL TYPE DETECTION", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] slic_composite_img=composite_label.tif watershed_cellzone_img=dapi-seg-wat.tif nuclei_segmented_img=dapi-seg.tif filtered=C3-alpha_filtered.tif filtered_0=C2-beta_filtered.tif filtered_1=C1-delta_filtered.tif region_observed=[INSIDE_OUTSIDE NUCLEUS] percent_marker_coverage=0.15 min_distance=0 max_distance_inside=1 max_distance_outside=3.000 save show");
close();
close();
run("Save");
selectWindow("TYPE_WAT");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/cell_type_colormap.lut");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/TYPE_WAT.tif");
selectWindow("TYPE_NUC");
run("Save");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/cell_type_colormap.lut");
run("Enhance Contrast", "saturated=0.35");
saveAs("Tiff", "C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/TYPE_NUC.tif");
close();







open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/dapi-seg-wat.tif");
open("C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/C4-dapi.tif");

run("DELTA CELL PROTRUSION", "save_dir=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1/] watershed_image=dapi-seg-wat.tif dapi_image=C4-dapi.tif min_distance=0 max_distance=35 region_value=0 target_cell=BETA protrusion_cell=DELTA range_observe=ALL_CELLS_PROTRUSION region_observe=CELL_REGION save=[C:/Users/SALAB VR/Documents/Hoa/Spatial3DTissueJ-master/H1536_islet1//protrusion_DELTA/Log_protrusion_DELTA_ALL_CELLS_PROTRUSION.txt]");