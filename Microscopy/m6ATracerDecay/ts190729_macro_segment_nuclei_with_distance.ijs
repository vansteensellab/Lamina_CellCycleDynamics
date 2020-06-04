macro 'segment-nuclei' {

    /*
    * - Segment nuclei
    * This macro will do a few things:
    * 
    *   given a .lif file, for every set of images:
    *     segment the nuclei
    *     for every nuclei:
    *       duplicate the window for a nuclei
    *       create a DAPI distance mask
    *       create a lamina distance mask
    *       save the raw intensities as pnm
    *       save the distance masks as pnm
    * 
    * Tom, 2019
    * 
    */
    
    run("Bio-Formats Macro Extensions");
    
    // Global parameters
    ID1="DAPI";
    ID2="tracer";
    ID3="lamina";

    target="DAPI"; // To use for segmentation

    // Input paramaters
    dir1 = getDirectory("Choose folder with lif files ");
    dir2 = getDirectory("Choose Destination Directory ");
    
    // List the files
    list = getFileList(dir1);
    add_pixels = 5;
        
    setBatchMode(true);
    
    for (k=0; k<list.length; k++) {
        showProgress(k+1, list.length);
        print("processing ... "+k+1+"/"+list.length+"\n         "+list[k]);
        path=dir1+list[k];
        
        if (! matches(list[k], ".*pADam.*")) {
            print("No pA-DamID experiment");
            continue;
        }
        if (! endsWith(path, ".lif")) {
            print("No lif file experiment");
            continue; 
        }

        //how many series in this lif file?
        Ext.setId(path);//-- Initializes the given path (filename).
        Ext.getSeriesCount(seriesCount); //-- Gets the number of image series in the active dataset.
    
        for (j=1; j<=seriesCount; j++) {
        
        
            run("Bio-Formats", "open=path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+j);
            
            // Get the title, prepare the output directory and change the names
            title = getTitle();

            // if (! matches(title, ".*LMNB.*")) {
            //     continue;
            // }

            experiment = replace(title, "\\.lif.*", "");
            experiment_dir = dir2 + "/" + experiment + "_analysis";
            if (!File.exists(experiment_dir)) {
                File.makeDirectory(experiment_dir);
            }

            title_new = replace(title, ".*lif - ", "");
            print("Processing: " + title_new);

            rename(title_new);

            // Split channels + get names
            run("Split Channels");
            C1 = "C1-" + title_new;
            C2 = "C2-" + title_new;
            C3 = "C3-" + title_new;


            // Segment target and create ROIs
            if (target == "DAPI") {
               selectWindow(C1);
            } else if (target == "lamina") {
               if (ID2 == "lamina") {
                  selectWindow(C2);
               } else {
                  selectWindow(C3);
               }
            }

            run("Duplicate...", "duplicate");
            rename("segmentation");

            run("Gaussian Blur...", "sigma=2");
            run("Enhance Contrast...", "saturated=0 normalize");
            run("Convert to Mask", "method=Li");
            
            if (target == "DAPI") {
               run("Outline");
            } else {
               run("Skeletonize");
            }
            
            run("Fill Holes");
            run("Erode");
            run("Dilate");
            run("Dilate");
            // run("Dilate");
            
            // run("Convert to Mask", "method=Li");
            // run("Dilate");
            // run("Erode");
            // run("Fill Holes");
            
            run("Duplicate...", "duplicate");
            rename("distance_mask");

            // Define ROIs
            run("3D Objects Counter", "threshold=1 min.=2000 max.=100000000 exclude_objects_on_edges objects");

            // Add this to the ROI manager
            run("3D Manager");
            Ext.Manager3D_AddImage();


            // Create distance maps for DAPI and lamina
            selectWindow("distance_mask");
            run("Distance Map");

            // For every ROI, create a duplicate and save all as pcg


            // Get the total counts - note that the last one is rubbish
            Ext.Manager3D_Count(nb_obj);
            nb_obj = nb_obj - 1;

            print("number of objects", nb_obj);

            // Get all individual cells and save these
            for (i=0; i<nb_obj; ++i) {

                print(i);
                Ext.Manager3D_Bounding3D(i, x0, x1, y0, y1, z0, z1);


                // Select the window of interest
                // 1) DAPI
                selectWindow(C1);
                run("Duplicate...", "title=crop");
                makeRectangle(x0-add_pixels, y0-add_pixels, x1-x0+2*add_pixels, y1-y0+2*add_pixels);
                run("Crop");
                saveAs("pgm", experiment_dir + "/" + title_new + "_" + (i+1) + "_DAPI.pgm");
                close();


                // 2) Lamina
                if (ID2 == "lamina") {
                   selectWindow(C2);
                } else {
                   selectWindow(C3);
                }
                
                run("Duplicate...", "title=crop");
                makeRectangle(x0-add_pixels, y0-add_pixels, x1-x0+2*add_pixels, y1-y0+2*add_pixels);
                run("Crop");
                saveAs("pgm", experiment_dir + "/" + title_new + "_" + (i+1) + "_lamina.pgm");
                close();


                // 3) Tracer
                if (ID2 == "lamina") {
                   selectWindow(C3);
                } else {
                   selectWindow(C2);
                }
                
                run("Duplicate...", "title=crop");
                makeRectangle(x0-add_pixels, y0-add_pixels, x1-x0+2*add_pixels, y1-y0+2*add_pixels);
                run("Crop");
                saveAs("pgm", experiment_dir + "/" + title_new + "_" + (i+1) + "_tracer.pgm");
                close();


                // 4) Distance mask
                selectWindow("distance_mask");
                run("Duplicate...", "title=crop");
                makeRectangle(x0-add_pixels, y0-add_pixels, x1-x0+2*add_pixels, y1-y0+2*add_pixels);
                run("Crop");
                saveAs("pgm", experiment_dir + "/" + title_new + "_" + (i+1) + "_distancemap.pgm");
                close();
                
                // 5) Segmentation
                selectWindow("segmentation");
                run("Duplicate...", "title=crop");
                makeRectangle(x0-add_pixels, y0-add_pixels, x1-x0+2*add_pixels, y1-y0+2*add_pixels);
                run("Crop");
                saveAs("pgm", experiment_dir + "/" + title_new + "_" + (i+1) + "_segmentation.pgm");
                close();
    

            }
            

            // Cleanup
            run("Close All");
            run("Collect Garbage");

            // list = getList("window.titles");
            // print(list.length);
            // for (i=0; i<list.length; i++){
            //     winame = list[i];
            //     selectWindow(winame);
            //     run("Close");
            // }
        
        }

    }
    showMessage(" -- finished --");    
    run("Close All");
    setBatchMode(false);

} // macro