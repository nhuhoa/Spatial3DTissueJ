function AutoThresholdDark () {
	//setOption("BlackBackground", false);
	partofHisto=0;histo1=0;
	run("Select None");imageinit = getImageID();
	depth=bitDepth;
	run("Duplicate...", "title=AutoThreshold");
	workimage = getImageID();
	if (depth > 8) {run("8-bit");}
	run("Grays");
	resetThreshold();
	setOption("ScaleConversions", true);


	// get the auto threshold values
	//setAutoThreshold(); //Default method
	setAutoThreshold("Otsu dark"); //Using Otsu method here
	getThreshold(lower, upper );
	print("\\Clear");
	resetThreshold();
	// analyse of the histogram:
	getStatistics(area, mean, min, max, std, histogram);
	vol=newArray(3);
	vol[0]=volumeHisto (histogram,min,lower);
	vol[1]=volumeHisto (histogram,lower,upper);
	vol[2]=volumeHisto (histogram,upper,max);
	//print("\"Volume\" 1 between min (" + min + ") and lower ("+ lower +") = "+ vol[0]);
	//print("\"Volume\" 2 between lower (" + lower + ") and upper ("+ upper +") = "+ vol[1]);
	//print("\"Volume\" 3 between upper (" + upper + ") and max ("+ max +") = "+ vol[2]);
	maxvol=maxOf(vol[1], vol[2]);
	for (a=1; a<3; a++) {if (maxvol == vol[a]) partofHisto = a;}
	if (partofHisto == 1) histo1=lower;
	if (partofHisto == 2) histo1=upper;
   	setThreshold(histo1, max,"black & white");
	//print("Black (dark pixels) obtained by thresholding between " + histo1 +" and " + max);
	print("Objects obtained by thresholding between " + histo1 +" and " + max);
	run("Convert to Mask"); 
	//run("Invert");
	run("Grays");
	mask1 = getTitle;
}

function volumeHisto (histo,mini,maxi) {

	volhisto=0;
	if ((maxi-mini) > 0) {
		for (i=mini; i<= maxi; i++) {
			volhisto=volhisto +(i* histo[i]);
		}
	}
	return volhisto;
}


	setBatchMode(true);
	AutoThresholdDark ();
	setBatchMode("exit and display");


short_name = substring(name,0,lastIndexOf(name,suffixe));  // TO DO: cut suffixe instead of using indexOf function here, issue raised when there are 2 '.tif' words in a file name					
print("_______________________________________\n");
print(name);
