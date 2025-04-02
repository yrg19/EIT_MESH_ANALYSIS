This repository contains the final versions of the data and code used for all results found in the upgrade report. 

In the code you will find a few different functions: 

1) Functions relating to mesh generation:

   - segmentation_to_elec_coords uses the mask from seg3d CT scan to identify the electrode coordinates
   - this can be used as an input to create_elec_mask which uses mk_cyl, rotate_cyls, as well as the pixel size and electrode coordinates to make a 3D mesh mask to add to seg3d
   - the next step would be to create an inr file from the segmentation and put this through mesher. the grid search function is also here (the only shell script).(Data is also in data, mesh_segmentations) 
   - the mesh convergence analysis function is also here along with find perturbation indices to simulate activity. mk_stim_kai is the function used to create protocols

2) Functions for scalp electrodes 
   - make_32_elecs generates the scalp electrode coodinates from the fudicials.
   - scalp_electrodes does all the analysis comparing depth only with depth plus scalp

3) Functions for human EIT data analysis
   - clean_EIT_EEG_final does that, it is the same as the Clean_EIT_EEG_Nov24 function locally which is the last version used.
   - to reconstruct the timeseries you can use EIT_reconstruction_patient which makes all the VTKs for that timeseries.

Data is also here:

1) mesh optimisation are the segmentations for the mesh as well as the electrode coordiantes used to generate the meshes (always in pixels) but also in mm.
   - for the first patient, E05, i have also included the older version of grid search so you have the final names of the files used to make the meshes for this patient as the combi_mesh_convergence
     has the names for the second patient

2) the data in EIT001 and EIT002 is just the final EIT timeseries and EEG timeseries per spike.

Results: 
1) The raw mesh files and raw data are on your hard drive. 
