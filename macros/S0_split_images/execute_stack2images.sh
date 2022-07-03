#!/bin/sh

# batch macro input_filename
# Ref: page 20: https://imagej.nih.gov/ij/docs/macro_reference_guide.pdf


## Configuration
imagej_dir="/Applications/ImageJ.app" ## MacOS ImageJ java run folder
macro_dir="/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/macros/"
project_dir="/Users/hoatran/Documents/python_workspace/Spatial3DTissueJ/testing_images/"


# imagej_exe_file=${imagej_dir}jars/ij-1.53q.jar  ## Fiji file location, in general file name is ij.jar, sometimes include version here
imagej_exe_file=${imagej_dir}/Contents/Java/ij.jar ## ImageJ exe file location, in general file name is ij.jar, sometimes include version here

macro_fn="${macro_dir}S0_split_images/stack2images.txt"
sub_folder="raw_islet" # a folder where you keep all large tiff files that are extracted from imc_converter
task="splitChannels" # Stack to images
log_file="${macro_dir}S0_split_images/${task}.log"
echo $log_file
exec >> $log_file 2>&1 && tail $log_file

input_dir="${project_dir}${sub_folder}/"
echo "__________________________________\n"
echo "Input directory is: \n"
echo $input_dir
echo "Splitting raw composite images into channels: \n"
## You can change memory amount, ex: 20000m to 30000m or more, so program will run faster.
## MacOS background mode here
java -Xmx20000m -jar ${imagej_exe_file} -ijpath $imagej_dir/ -batch $macro_fn $input_dir 

## Linux background mode here, using xvfb-run in case you run in server due to lack of graphical environment
## If you don't have right permission in server, you can install conda first, and install xvfb package on top of that
# xvfb-run -a java -Xmx15000m -jar $imagej_dir/ij.jar -ijpath $imagej_dir/ -batch $macro_fn $input_dir 

## In Linux local computer, java command is sufficient
#java -Xmx15000m -jar $imagej_dir/ij.jar -ijpath $imagej_dir/ -batch $macro_fn $input_dir 

echo "Completed!"
echo "__________________________________\n"





