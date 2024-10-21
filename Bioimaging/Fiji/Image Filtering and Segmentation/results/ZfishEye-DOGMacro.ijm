open("C:/Files/Code/Bioinformatics/Bioimaging/Fiji/Image Filtering and Segmentation/data/zfish_eye.tif");
selectImage("zfish_eye.tif");
run("Gaussian Blur...", "sigma=8 stack");
open("C:/Files/Code/Bioinformatics/Bioimaging/Fiji/Image Filtering and Segmentation/data/zfish_eye.tif");
selectImage("zfish_eye-1.tif");
run("Gaussian Blur...", "sigma=2 stack");

imageCalculator("Subtract create stack", "zfish_eye-1.tif","zfish_eye.tif");
selectImage("Result of zfish_eye-1.tif");
setMinAndMax(3139, 57638);
setAutoThreshold("Yen dark no-reset");
//run("Threshold...");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Yen background=Dark calculate black");
run("Close");

//Split Channel to subtract the perimeter from blue fluorescent image
run("Split Channels");
selectImage("C1-Result of zfish_eye-1.tif");
run("Open");
run("Erode");
run("Open");
run("Open");
imageCalculator("Subtract create", "C1-Result of zfish_eye-1.tif","C3-Result of zfish_eye-1.tif");
selectImage("Result of C1-Result of zfish_eye-1.tif");

//Clearing the noise again
run("Open");
run("Fill Holes");
run("Dilate");
run("Fill Holes");

//Measure it and add is as overlay to original
open("C:/Files/Code/Bioinformatics/Bioimaging/Fiji/Image Filtering and Segmentation/data/zfish_eye.tif");
run("Analyze Particles...", "display clear summarize add");
selectImage("zfish_eye-2.tif");
roiManager("Measure");
roiManager("Show None");
roiManager("Show All");
