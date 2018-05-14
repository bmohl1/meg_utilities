# meg_utilities
This set of scripts is used to process MEG data by calling methods from SPM, Fieldtrip, and other packages.

The main function that will run the other modules is the "grand_function.m". In short, the function loops participants through each of the steps that are defined in the other batch files.

Final stats can be run from the grand function or separately from the stats scripts called in grand function. All stats require text files with MNI coordinate locations and common, neurological names of the peaks of interest.  

e.g. label_file:  
lstg  
rstg  
e.g. corresponding coords_file (MNI space):  
48 -14 -2  
-44 -64 -12  
