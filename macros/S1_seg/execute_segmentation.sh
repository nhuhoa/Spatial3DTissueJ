#!/bin/sh

# batch macro input_filename
# Ref: page 20: https://imagej.nih.gov/ij/docs/macro_reference_guide.pdf

## Configuration MacOS
imagej_dir="/Applications/ImageJ.app" ## MacOS ImageJ java run folder
macro_dir="/Users/hoatran/Documents/BCCRC_projects/IMC/imc_pipeline/quality_control_pipeline/XP1487_IMCpipeline/script/"
project_dir="/Users/hoatran/Documents/BCCRC_projects/IMC/XP1487/testing/XP1487_IMCpipeline/"

series="raw/220202_Cambridge_142309_PreStained_s0_a1_ac.ome_rawdata_NUC/" 
task="NUC_SEG" # Nuc segmentation

## ImageJ/ Fiji can be installed to Applications, or put into any folder in your drive. 
imagej_dir="/Applications/ImageJ.app" ## MacOS ImageJ java run folder
# imagej_dir="/Applications/Fiji.app" ## MacOS Fiji java run folder
# imagej_dir="/Users/htran/Downloads/Fiji.app/" ## MacOS Fiji java run folder - from my computer


# imagej_exe_file=${imagej_dir}jars/ij-1.53q.jar  ## Fiji file location, in general file name is ij.jar, sometimes include version here
imagej_exe_file=${imagej_dir}/Contents/Java/ij.jar ## ImageJ exe file location, in general file name is ij.jar, sometimes include version here


macro_fn="${macro_dir}S1_seg/nuc_seg.txt"



log_file="${macro_dir}S1_seg/${task}.log"
echo $log_file
exec >> $log_file 2>&1 && tail $log_file

input_dir="${project_dir}${series}"
echo "__________________________________\n"
echo "Input directory is: \n"
echo $input_dir
echo "Segment images \n"
## You can change memory amount, ex: 20000m to 30000m so program will run faster.
## MacOS background mode here

## If you have java environment jre installed in your computer
java -Xmx20000m -jar ${imagej_exe_file} -ijpath $imagej_dir/ -batch $macro_fn $input_dir

## Otherwise using existing jre env from Fiji here
# /Users/htran/Downloads/Fiji.app/java/macosx/adoptopenjdk-8.jdk/jre/Contents/Home/bin/java -Xmx20000m -jar $imagej_dir/Contents/Java/ij.jar -ijpath $imagej_dir/ -batch $macro_fn $input_dir
# ${imagej_dir}java/macosx/adoptopenjdk-8.jdk/jre/Contents/Home/bin/java -Xmx20000m -jar ${imagej_exe_file} -ijpath $imagej_dir/ -batch $macro_fn $input_dir

## Linux background mode here, using xvfb-run in case you run in server (graphical env), in local computer, java command is enough
# xvfb-run -a java -Xmx15000m -jar $imagej_dir/ij.jar -ijpath $imagej_dir/ -batch $macro_fn $input_dir 

## In Linux local computer, java command is sufficient
# java -Xmx15000m -jar $imagej_dir/ij.jar -ijpath $imagej_dir/ -batch $macro_fn $input_dir 

echo "Nucleus Segmentation Completed!"
echo "__________________________________\n"





