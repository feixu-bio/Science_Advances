run("Close All");
//Change dir path. Please do not use spaces in the file name (use underscores or hyphens instead).
inputDir = "U:/VAMP/........../";
list = getFileList(inputDir);
for (image = 0; image < list.length; image++){
    run("Bio-Formats Windowless Importer", "open="+inputDir+list[image]);
    title = getTitle();
    print(title);
    selectWindow(title);
    idx=indexOf(title, ".tif");
    sub = title.substring(0,idx);
    addtitle = sub + "-1.tif";    
    run("Split Channels");
    selectImage("C1-" + title);
    run("Median...", "radius=4");
    run("Gaussian Blur...", "sigma=4");
    setMinAndMax(0, 300);
    setAutoThreshold("Mean dark 16-bit no-reset");
    setOption("BlackBackground", true);
    run("Convert to Mask");
    //dilate 1 pixel to include faint DAPI edge
    run("Dilate");
    selectImage("C1-" + title);
    run("Duplicate...", " ");
    selectImage("C1-" + addtitle);
    run("Distance Map");
    setAutoThreshold("Default dark 16-bit no-reset");
    //include a rouph 375 nm thickness of 11 pixels
    setThreshold(11, 1000, "raw");
    setOption("BlackBackground", true);
    run("Convert to Mask");
    imageCalculator("Subtract create", "C1-" + title,"C1-" + addtitle);
    selectImage("Result of C1-"+ title);
    run("Create Selection");
    selectImage("C4-"+ title);
    run("Translate...", "x=-3 y=-1 interpolation=None");
    run("Restore Selection");
    run("Measure");
    //count the whole nuclues intensity
    selectImage("C1-" + title);
    run("Create Selection");
    selectImage("C4-"+ title);
    run("Restore Selection");
    run("Measure");
    run("Close All");}
selectWindow("Log");
saveAs("Text", inputDir+"ImgList.txt")
selectWindow("Results");
saveAs("Results",inputDir+"Results.csv");
