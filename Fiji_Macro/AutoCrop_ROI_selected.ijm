//This script is to automatically crop and save img based on saved ROI.
//Load just one image and corresponding ROI.

//First Merge all channels to be croped.
Stack.setChannel(1);
run("Magenta");
Stack.setChannel(2);
run("Green");
run("Subtract Background...", "rolling=50 stack");
run("Enhance Contrast", "saturated=0.35");
Property.set("CompositeProjection", "Sum");
Stack.setDisplayMode("composite");
Stack.setActiveChannels("11");

//Crop and Save Tiff
Savedir = "U:/VAMP/....../";
title = getTitle();
//idx=indexOf(title, ".vsi");
//sub = title.substring(0,idx);
print(title);
numROIs = roiManager("count");
print(numROIs);

for (i = 0; i < numROIs; i++){
nuc = i;
selectWindow(title);
roiManager("Select", nuc);
//currentZ = getSliceNumber();
//run("Duplicate...", "duplicate range=currentZ-currentZ use");
run("Duplicate...", " ");
}

for (i = 0; i < numROIs; i++){
img = i+1;
subtitle = title + "-" + img;
print(subtitle);
selectWindow(subtitle);
saveAs("Tiff", Savedir + subtitle);
}
close("*");
