// imagej macro for batched vesselness computation using ... plugin
// search for image files (*.zip) within the given directory (for other formats replace "zip" with the ext)
// apply the set of parameters when processing each image
// java terminal command for running the macro:
// java -Xmx4g -jar ~/ImageJ/ij.jar -ijpath ~/ImageJ/plugins/ -batch vess.ijm given_directory_path
// java -Xmx4g -jar ~/ImageJ/ij.jar -ijpath ~/ImageJ/plugins/ -batch ~/vess/macro/vess.ijm ~/vess.test/a/

setBatchMode(true);

sigmas = "2,4,6";

// 3d
alpha = 0.5;
beta = 0.5;
c = 500;
zdist = 2.0;
threaded3d = false;
if (threaded3d) 		  threaded3d_arg 			= " threaded3d";	
else 					  threaded3d_arg            = "";

// 2d
betaone = 0.5;
betatwo = 15;

dark_foreground = false;
if (dark_foreground)      dark_foreground_arg      = " darkforeground";
else    			      dark_foreground_arg      = "";

calculate_directions = true;
if (calculate_directions) calculate_directions_arg = " directions";
else                      calculate_directions_arg = "";

save_output = true;
if (save_output)    	  save_output_arg 		   = " saveoutput";
else                      save_output_arg 		   = "";


print(calculate_directions_arg);
 
dir = getArgument;
if (dir=="") exit("Command needs an argument at the end: path to the directory with the images to process.");
if (!File.isDirectory(dir)) exit(dir + " is not a directory.");
if (!endsWith(dir, "/")) dir += "/";

t_start = getTime();

files = getFileList(dir);

cnt = 0;
for (i=0; i<files.length; i++) {

	if (endsWith(files[i], ".zip")) { // alternative paterns "[0-9][0-9][0-9].ext" or matches(files[i], "*.ext")
		
		open(dir + files[i]);
		getDimensions(null, null, null, nr_slices, null);
		file_name = File.nameWithoutExtension;
		close();
		run("Close All");
		
		// run macro command depending on the number of the slices (2D or 3D image, different parameter list)
		if (nr_slices==1) {
			run("Vess", "select="+dir+files[i]+" sigmas="+sigmas+dark_foreground_arg+calculate_directions_arg+save_output_arg+" betaone="+betaone+" betatwo="+betatwo);
		}	
		else {
			run("Vess", "select="+dir+files[i]+" sigmas="+sigmas+dark_foreground_arg+calculate_directions_arg+save_output_arg+" alpha="+alpha+" beta="+beta+" c="+c+" zdist="+zdist+threaded3d_arg);
		}

	}
}

t_end = getTime();

print("done. elapsed " + ((t_end-t_start)/60000) + " min.");