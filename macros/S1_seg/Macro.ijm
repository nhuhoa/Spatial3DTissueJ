open("/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/filtered/C4-dapi.tif");
run("NUCLEI SEGMENTATION", "save=/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/seg/ dapi=C4-dapi.tif maximal=28 minimal=18 seed=29000 automatic");
close();
close();
selectWindow("C4-dapi.tif");
close();
