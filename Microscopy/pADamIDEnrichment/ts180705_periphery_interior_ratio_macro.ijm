macro 'pA-DamID enrichments' {

    /*
    * - pA-DamID enrichments
    * 
    * Tom, 2018
    * 
    */
    
    run("Bio-Formats Macro Extensions");
        
    dir1 = getDirectory("Choose folder with lif files ");
    list = getFileList(dir1);
        
    setBatchMode(true);
    
    for (i=0; i<list.length; i++) {
        showProgress(i+1, list.length);
        print("processing ... "+i+1+"/"+list.length+"\n         "+list[i]);
        path=dir1+list[i];
        
        if (! endsWith(path, ".lif")) {
            continue; 
        }

        //how many series in this lif file?
        Ext.setId(path);//-- Initializes the given path (filename).
        Ext.getSeriesCount(seriesCount); //-- Gets the number of image series in the active dataset.
    
        for (j=1; j<=seriesCount; j++) {
        
            run("Bio-Formats", "open=path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+j);
        
            getDimensions(width, height, channels, slices, frames);
        
            // Global parameters
            working_dir = dir1 + "/analysis";
            if (!File.exists(working_dir)) {
                File.makeDirectory(working_dir);
            }


            // Get the title, prepare the output directory and change the names
            title = getTitle();

            experiment = replace(title, "\\.lif.*", "");
            experiment_dir = working_dir + "/" + experiment + "_analysis";
            if (!File.exists(experiment_dir)) {
                File.makeDirectory(experiment_dir);
            }

            title_new = replace(title, ".*lif - ", "");

            rename(title_new);

            // Split channels + get names
            run("Split Channels");
            C1 = "C1-" + title_new;
            C2 = "C2-" + title_new;
            C3 = "C3-" + title_new;


            // 1) DAPI segmentation
            selectWindow(C1);
            run("Duplicate...", "title=nuclei");

            run("Gaussian Blur...", "sigma=1");
            run("Enhance Contrast...", "saturated=0 normalize");
            setOption("BlackBackground", false);

            run("Make Binary");
            run("Fill Holes");

            // Create file of DAPI segmentation?


            // 2) Peripheral segmentation
            run("Duplicate...", "title=periphery");
            run("Erode");
            //run("Erode");
            run("Outline");
            //run("Dilate");
            run("Dilate");
            run("Dilate");
            run("Dilate");

            // imageCalculator("AND create", "periphery", C2);
            // selectWindow("Result of periphery");
            // rename("periphery - GFP");
            //
            // imageCalculator("AND create", "periphery", C3);
            // selectWindow("Result of periphery");
            // rename("periphery - antibody");

            // Create file of periphery


            // 3) Interior segmentation
            imageCalculator("Subtract create", "nuclei","periphery");
            selectWindow("Result of nuclei");
            rename("interior");

            // imageCalculator("AND create", "interior", C2);
            // selectWindow("Result of interior");
            // rename("interior - GFP");
            //
            // imageCalculator("AND create", "interior", C3);
            // selectWindow("Result of interior");
            // rename("interior - antibody");


            // 4) Generate statistics
            // 4a) GFP channel
            if (channels == 2) {
                
                run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=10 show_numbers white_numbers store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=[" + C2 + "]");

                selectWindow("periphery");
                run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=1048576 exclude_objects_on_edges statistics");
                IJ.renameResults("Statistics for periphery redirect to " + C2, "Result");
                selectWindow("Result");
                saveAs("Results",  experiment_dir + "/" + title_new + ",periphery,antibody.csv");
                run("Close");

                selectWindow("interior");
                run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=1048576 exclude_objects_on_edges statistics");
                IJ.renameResults("Statistics for interior redirect to " + C2, "Result");
                selectWindow("Result");
                saveAs("Results",  experiment_dir + "/" + title_new + ",interior,antibody.csv");
                run("Close");
                
            } else {
                
                run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=10 show_numbers white_numbers store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=[" + C2 + "]");

                selectWindow("periphery");
                run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=1048576 exclude_objects_on_edges statistics");
                IJ.renameResults("Statistics for periphery redirect to " + C2, "Result");
                selectWindow("Result");
                saveAs("Results",  experiment_dir + "/" + title_new + ",periphery,tracer.csv");
                run("Close");

                selectWindow("interior");
                run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=1048576 exclude_objects_on_edges statistics");
                IJ.renameResults("Statistics for interior redirect to " + C2, "Result");
                selectWindow("Result");
                saveAs("Results",  experiment_dir + "/" + title_new + ",interior,tracer.csv");
                run("Close");

                // 4b) antibody channel
                run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=10 show_numbers white_numbers store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=[" + C3 + "]");

                selectWindow("periphery");
                run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=1048576 exclude_objects_on_edges statistics");
                IJ.renameResults("Statistics for periphery redirect to " + C3, "Result");
                selectWindow("Result");
                saveAs("Results",  experiment_dir + "/" + title_new + ",periphery,antibody.csv");
                run("Close");

                selectWindow("interior");
                run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=1048576 exclude_objects_on_edges statistics");
                IJ.renameResults("Statistics for interior redirect to " + C3, "Result");
                selectWindow("Result");
                saveAs("Results",  experiment_dir + "/" + title_new + ",interior,antibody.csv");
                run("Close");
                
            }
            


            // Cleanup
            // selectWindow("interior");
            // run("Close");
            // selectWindow("periphery");
            // run("Close");
            // selectWindow("nuclei");
            // run("Close");
            // selectWindow(C1);
            // run("Close");
            // selectWindow(C2);
            // run("Close");
            // selectWindow(C3);
            // run("Close");
        
        
        
            run("Close All");
            run("Collect Garbage");
        
        }

    }
    showMessage(" -- finished --");    
    run("Close All");
    setBatchMode(false);

} // macro