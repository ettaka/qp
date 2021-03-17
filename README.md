Script files for characterizing MQXF mechanical behaviour


** Coil_size program:

*) The size (and gap) template files are available in the Data/Test folder. For new coils they need to be modified using a text editor. Be careful no to modify the file structure and format.


** Cross-sect to ANSYS program:

*) The code needs to be launched with the data of four coils. Both the .txt file and the .size file are required. Examples are found in the Data/Test folder.

*) For the reading of the files, the path can be changed in the code if needed (Line 153, where the program replaces the .txt by .size and looks for the file one level up.). 
