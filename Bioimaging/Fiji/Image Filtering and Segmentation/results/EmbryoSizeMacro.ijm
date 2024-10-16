open("C:/Users/Ricardo/Downloads/lecture_applied_bioimage_analysis_2020-master/02_Image_Filtering_and_Segmentation/examples/embryos.tif");
selectImage("embryos.tif");
run("Split Channels");
run("Gaussian Blur...", "sigma=3");
setAutoThreshold("Default no-reset");
//run("Threshold...");
//setThreshold(0, 58);
run("Convert to Mask");
run("Close");
run("Fill Holes");
run("Watershed");
run("Analyze Particles...", "add");
open("C:/Users/Ricardo/Downloads/lecture_applied_bioimage_analysis_2020-master/02_Image_Filtering_and_Segmentation/examples/embryos.tif");
selectImage("embryos.tif");
roiManager("Measure");
saveAs("Results", "C:/Users/Ricardo/Downloads/lecture_applied_bioimage_analysis_2020-master/02_Image_Filtering_and_Segmentation/results/EmbryoSize.csv");

