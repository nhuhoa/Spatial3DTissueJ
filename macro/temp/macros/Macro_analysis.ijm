open("/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/seg/dapi-seg-wat.tif");
run("CELLS LAYER CONTACT", "save=/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/cell_type/ watershed_image=dapi-seg-wat.tif dapi_image=C4-dapi.tif target_cell=BETA source_cell=DELTA");
close();
run("CELLS INTERACTION ANALYSIS");
run("RANDOM DISTRIBUTION", "save=/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/cell_type/ watershed_image=dapi-seg-wat.tif observed_cell_type=ALPHA source_cell_type_observed=BETA nb_randoms=100");
run("CLUSTERS ANALYSIS", "save=/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/cell_type/ watershed_image=dapi-seg-wat.tif type_observed_cell=ALPHA clusters_spatial_statistic save=/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/cell_type/ALPHA_clusters_analysis/F_function_Euclidean_Distance.png");
selectWindow("dapi-seg-wat.tif");
selectWindow("Euclidean-Distance-G-function");
selectWindow("dapi-seg-wat.tif");
