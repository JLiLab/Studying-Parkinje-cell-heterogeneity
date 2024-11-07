// ImageJ Macro Script for Processing Cyc-Immunofluorescence Images
// This script processes .czi images named in the format "Ab1_Ab2_Ab3_Ab4_DAPI-ID-TR1.czi" for each scene.
// Images from individual slides should be placed in subfolders labeled as "ID_data".
// The script will recursively process all subfolders, automatically extracting antibody names from filenames 
// and using the DAPI channel to register multichannel images across different rounds of cyc IF.
// Note: This version includes optimized StarDist processing for faster execution.

// Author: James Li

// Get the current folder and a list of all subdirectories
inputFolder = getDirectory("images");

// list of tile positions
positions = newArray();
for (i = 1; i <= 50; i++) {
    positions[i - 1] = "TR" + i;
}

// make selection of mehtods for registration: "Rigid", "Similarity",or "Affine"
Dialog.create("Set Variables");
Dialog.addString("Set alignment method:", "Rigid"); // change this to either "Rigid", "Similarity", or "Affine"
Dialog.addString("Enter the maximal number of channels:", 6); // Add an input field to the dialog window
Dialog.addString("Enter the magnification:", "20X"); // Add an input field to the dialog window
Dialog.addString("Enter DAPI selection:", "all"); // select specific slice like "-1-3-"
Dialog.addString("Segmentation:", "YES"); 
Dialog.show();

mode = Dialog.getString();
maxChannels = Dialog.getString();
maxChannels = parseInt(maxChannels); 
magnification = Dialog.getString();
slices_to_keep_string = Dialog.getString();
segmentation = Dialog.getString();

print(segmentation);

roiManager("show none");
setBatchMode(true);
//setBatchMode(false);
run("Close All");

fileList = getFileList(inputFolder);
subdirectories = newArray();

// generate a list of subfolders
for (dir = 0; dir < fileList.length; dir++) {
    if (File.isDirectory(inputFolder + fileList[dir])) {
        subdirectories[subdirectories.length] = fileList[dir];
    }
}

// process images in each subfolder
for (k = 0; k < subdirectories.length; k++) {
    parts = split(subdirectories[k], "_");

    // Check if the subfolder name contains an underscore
    if (parts.length > 1) {
        ID = parts[0];
    } else {
        // Handle the case where there is no underscore
        // For example, use the entire folder name as ID
        // or skip the folder by continuing to the next iteration
        ID = subdirectories[k]; // or use 'continue;' to skip
    }

    part1 = split(subdirectories[k], "/");
    path = inputFolder + part1[0];
    print("Processing folder: " + path);

    outFolder = path + File.separator + "aligned_" + ID;
    File.makeDirectory(outFolder);
    
    list = getFileList(path);
	
	for (tile = 0; tile < positions.length; tile++) {
		suffix = positions[tile]+".czi";
		Ab_name = "";
		
		// open all images end with a certain suffix
		for (i = 0; i < list.length; i++) {
			if(endsWith(list[i], suffix)){
				print(list[i]);
				imagePath = path + File.separator + list[i];
				
				run("Bio-Formats Importer", "open=[" + imagePath + "] color_mode=Default view=Hyperstack stack_order=XYCZT autoscale ");
				
				// perform backgroun subtraction
				run("Subtract Background...", "rolling=50 stack");
				
				// enhance contrast to facilitate registration
				setSlice(nSlices);
				run("Enhance Contrast...", "saturated=0.35");
				
				// extract marker names
				Image_name = list[i];
				parts = split(Image_name, "-");
				string = parts[0];
				
				// get the current image and its number of channels
	  			getDimensions(w, h, channels, slices, frames);
	  
	 			 // If the current image has fewer channels than the maximal number of channels, insert blank channels
				if (channels < maxChannels) {
	  				add_name = "";
	    			for (n = channels+1; n <= maxChannels; n++) {
	      			run("Add Slice", "add=channel prepend");
	      			add_name = "blank_"+add_name;
	    			}
	    			string = add_name + "_" + string;	
	    			Ab_name = Ab_name + "_" + string;
				} else {
				  	Ab_name = Ab_name + "_" + string;
				}
			}	
		}
		
		// register all open imaegs 
		print(Ab_name);
		if (nImages > 1) {
	 	labels = split(Ab_name, "_");
		getPixelSize(unit, pw, ph, pd);
		
	  	run("Concatenate...", "all_open title=Hyperstack open");
		run("Linear Stack Alignment with SIFT MultiChannel", "registration_channel="+maxChannels+" initial_gaussian_blur=1.60 steps_per_scale_octave=3 minimum_image_size=64 maximum_image_size=1024 feature_descriptor_size=4 feature_descriptor_orientation_bins=8 closest/next_closest_ratio=0.92 maximal_alignment_error=25 inlier_ratio=0.05 expected_transformation=" + mode+ " interpolate");
		close("Hyperstack");			
		
		// add antibody label to each channel of the aligned image
		selectWindow("Aligned_Hyperstack");
		run("Stack to Images");
		
		ids=newArray(nImages); 
		for (m=0; m<nImages; m++) { 
		selectImage(m+1); 
		ids[m]=getImageID;
		rename(labels[m]);
		}			
		
		// delete blank channels that were added as fillers
		close("blank"); 
		
		// extract DAPI slices
		run("Images to Stack", "name=DAPI title=DAPI use");	
		run("Images to Stack", "name=Multiplexed use");
		
		// scale 10x images for Stardist seqmentation
		if (magnification == "10X") {
			selectWindow("DAPI");
			run("Scale...", "x=2 y=2 z=1.0 width=7932 height=8600 depth=7 interpolation=Bilinear average process create");
			close("DAPI");
			
			selectWindow("DAPI-1");
			rename("DAPI");
			
			selectWindow("Multiplexed");
			run("Scale...", "x=2 y=2 z=1.0 width=7932 height=8600 depth=7 interpolation=Bilinear average process create");
			close("Multiplexed");
			selectWindow("Multiplexed-1");
			rename("Multiplexed");
		}
		
		run("Grays");
		
		selectWindow("DAPI");
		if (slices_to_keep_string == "all") {
			// make orthogonal projection of all DAPI images
			run("Z Project...", "projection=[Sum Slices]");
			run("Enhance Contrast...", "saturated=0.35");
		} else {
			run("Duplicate...", "title=tempStack duplicate");
			total_slices = nSlices;
			// Start from the end of the stack and move to the beginning
			for (i = total_slices; i >= 1; i--) {
			    slice_str = "-" + i + "-";  // Convert current slice number to the formatted string
			    if (indexOf(slices_to_keep_string, slice_str) == -1) {  // Check if the slice exists in the slices_to_keep_string
			        setSlice(i);
			        run("Delete Slice");
			    }
			}
			// Perform Maximal Intensity Projection
			total_slices = nSlices;
			if (total_slices > 1){
				run("Z Project...", "projection=[Sum Slices]");
			}
			run("Enhance Contrast...", "saturated=0.35");
			rename("SUM_DAPI");
				    
			// Close the temporary stack
			close("tempStack");
			}
		
		if (segmentation == "YES") {
			// perform nuclei segmentation using the Stardist plugin
			//run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'SUM_DAPI', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.5', 'nmsThresh':'0.4', 'outputType':'Both', 'nTiles':'4', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
		
			selectWindow("SUM_DAPI");
			run("StarDist2D... (Advanced)", "model_path=C:\\Users\\zeiss\\Documents\\Fiji.app\\models\\2D_versatile_fluo x_tiles=4 y_tiles=4 min_norm=1.0 max_norm=99.8 prob_thresh=0.5 nms_thresh=0.4");
			selectWindow("SUM_DAPI-stardist2D");
			rename("Label Image");
			run("glasbey inverted");
			run("Label image to ROIs");
		
			selectWindow("Label Image");
			run("8-bit");
			saveAs("Tiff", outFolder + File.separator + "Mask_" + ID + "_" + positions[tile]+ ".tif");
			close();
			
			selectWindow("DAPI");
			run("16-bit");
			saveAs("Tiff", outFolder + File.separator + "DAPI_" + ID + "_" + positions[tile]+ ".tif");
			close();
		
			// Add othorgonal projected DAPI to the multiplex images		
			// save segmentation results
			selectWindow("SUM_DAPI");
			run("Duplicate...", "title=DAPI");
			
			selectWindow("SUM_DAPI");
			roiManager("Show All");
			run("Flatten");
			saveAs("Tiff", outFolder + File.separator + "Label_" + ID + "_" + positions[tile]+ ".tif");
			close();
			close("SUM_DAPI");
			
					
			selectWindow("Multiplexed");
			run("Stack to Images");  
		    run("Images to Stack", "name=Multiplexed use");
		    
		    //measure each frames (expect the first frame)
			selectWindow("Multiplexed");	
			//run("Set Measurements...", "area mean min perimeter shape integrated skewness stack display redirect=None decimal=2");
			run("Set Measurements...", "area mean min centroid perimeter shape integrated skewness stack display redirect=None decimal=2");
			for (ch=0; ch < nSlices; ch++){	
				roiManager("Deselect");
				roiManager("Measure");
				run("Next Slice [>]");
			}
			
			// save results
			saveAs("Results", outFolder + File.separator + "Measurement_" + ID + "_" + positions[tile]+ ".csv");
				
			getDimensions(width,height,channels,slices,frames);
			run("Properties...", "channels="+slices+" slices="+channels+" frames="+frames+" unit="+unit+" pixel_width="+pw+" pixel_height="+ph+" voxel_depth="+pd);
			saveAs("Tiff", outFolder + File.separator + "Multiplexed_" + ID + "_" + positions[tile]+ ".tif");
			
			// Save the ROIs
			//roiManager("Save", outFolder + File.separator + "roiManager_" + ID + "_" + positions[tile]+ ".zip");
			
		} else {
			selectWindow("SUM_DAPI");
			rename("DAPI");
			
			//measure each frames (expect the first frame)
			selectWindow("Multiplexed");
			run("Stack to Images");  
		    run("Images to Stack", "use");
			
			getDimensions(width,height,channels,slices,frames);
			run("Properties...", "channels="+slices+" slices="+channels+" frames="+frames+" unit="+unit+" pixel_width="+pw+" pixel_height="+ph+" voxel_depth="+pd);
			saveAs("Tiff", outFolder + File.separator + "Multiplexed_" + ID + "_" + positions[tile]+ ".tif");	
			}
		run("Clear Results");
		run("Close All");	
		} 
		
	}
}	

setBatchMode(false);