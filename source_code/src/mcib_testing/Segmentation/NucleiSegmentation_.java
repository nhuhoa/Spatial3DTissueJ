/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Segmentation;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import fiji.util.gui.GenericDialogPlus;
import ij.plugin.PlugIn;
import ij.process.ImageStatistics;
import java.io.File;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageInt;
/**
 *
 * @author tranhoa
 */
public class NucleiSegmentation_ implements PlugIn
{
    String dir = null;
    int max_nuc_radius = 28, min_nuc_radius = 18, seed_threshold = 29000;
    public String save_dir = "";
    public void run(String arg) 
    {	
        if (IJ.versionLessThan("1.37f")) return;
        int[] wList = WindowManager.getIDList();

        if (wList==null) {
                IJ.showMessage("Nuclei Segmentation", "There must be at least one image open.");
                return;
        }
        String[] titles = new String[wList.length];
        for (int i=0, k=0; i<wList.length; i++) {
                ImagePlus imp = WindowManager.getImage(wList[i]);
                if (null !=imp)
                        titles[k++] = imp.getTitle();
        }
        ImagePlus temp = WindowManager.getImage(wList[wList.length-1]);
        if (null !=temp){
            WindowManager.setTempCurrentImage(temp);
            save_dir = IJ.getDirectory("image");
        }
        
        boolean prefilter = false, automatic_seed=false; 
        GenericDialogPlus gd = new GenericDialogPlus("3D Nuclei Segmentation");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("See and quote reference:\nA novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : DAPI nuclei image");
        gd.addMessage("Output: segmented nuclei image");
        gd.addMessage(" \n");
        gd.addMessage("Option 1: automatic seed threshold - recommended\n Option 2: predefined a threshold");
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("DAPI Nuclei Image :", titles, titles[0]);
//        gd.addCheckbox("Pre-filtering (3D Median 4-4-2)", false);
        gd.addNumericField("Maximal Nucleus Size", max_nuc_radius, 0, 10, "");
        gd.addNumericField("Minimal Nucleus Size", min_nuc_radius, 0, 10, "");
        gd.addCheckbox("Automatic Seed Threshold", true);
        gd.addNumericField("Seed Threshold: ", seed_threshold, 0, 10, "");
        
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
        
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return;
        }
        
        int i1Index = gd.getNextChoiceIndex();
//        prefilter  = gd.getNextBoolean();
        max_nuc_radius = (int) gd.getNextNumber();
        min_nuc_radius = (int) gd.getNextNumber();
        automatic_seed = gd.getNextBoolean();
        seed_threshold = (int) gd.getNextNumber();
        
        ImagePlus imp1 = WindowManager.getImage(wList[i1Index]);
//        WindowManager.setTempCurrentImage(imp1);
//        dir = IJ.getDirectory("image");
        
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        IJ.log("====================================================================");
        IJ.log("Initializing...");
        IJ.log("Max nucleus radius: "+max_nuc_radius+"  min nucleus radius: "+min_nuc_radius+"  seed:"+seed_threshold);
        segment(imp1, prefilter, automatic_seed);
        IJ.log("Finished.");
        IJ.log("====================================================================");
    }
    
    public void segment(ImagePlus imp1, boolean prefilter, boolean automatic_seed)
    {
        if(save_dir!=null)
        {
            //            String subdir = null;
//            File fl = new File(dir+"SEG");
            File fl = new File(save_dir);
            if (!fl.exists()) 
            {
                fl.mkdir();
//                if (fl.mkdir()) {
//                    subdir = dir + "SEG/";
//                } else {
//                    subdir = dir;
//                }
            }
//            else{
//                subdir = dir + "SEG/";
//            }

        }
        ImagePlus imp2 = null;
//        if(prefilter)
//        {
//            IJ.log("Pre-filtering...");
//            IJ.run(imp1, "Median 3D...", "x=4 y=4 z=2");
////            IJ.log("Duplicate...");
//            imp2 = imp1.duplicate();
//            imp2.setTitle("C4-dapi");
////            if(save_dir!=null)
////            {
////                IJ.selectWindow(imp2.getTitle());
////                IJ.saveAs("Tiff", save_dir+imp2.getTitle()+".tif");
////                IJ.log("Save the nuclei image as C4-dapi.tif in "+save_dir);
////            }    
//            
//        }
//        else{
//            IJ.log("Duplicate...");
            imp2 = imp1.duplicate();
            imp2.setTitle("C4-dapi");
//            if(dir!=null)
//            {
//                IJ.saveAs("Tiff", dir+"C4-dapi.tif");
//                IJ.log("Save the nuclei image as C4-dapi.tif in "+dir);
//            }  
//        }
        
        //band pass filtering
        IJ.run(imp2, "32-bit", "");
        IJ.log("Computing Band Pass...");
        IJ.run(imp2, "DF_Bandpass", "maximal_feature_size="+max_nuc_radius+" minimal_feature_size="+min_nuc_radius+" z/x_aspect_ratio=3.000");
        IJ.selectWindow("BP_C4-dapi");
        if(save_dir!=null)
        {
            IJ.selectWindow("BP_C4-dapi");
            IJ.saveAs("Tiff", save_dir+"BP_C4-dapi.tif");
            IJ.log("Save the band pass filtered image as BP_C4-dapi.tif in "+save_dir);
            ImagePlus imp3 = IJ.openImage(save_dir+"BP_C4-dapi.tif");
//            imp3.setSlice(imp3.getNSlices()/2);
            IJ.run(imp3, "Enhance Contrast...", "saturated=0"); 
            IJ.run(imp3, "16-bit", "");
            
            if(automatic_seed){
                
                IJ.log("Local maximal seeds are defined automatically here: ");
//                ImageStatistics s = imp3.getStatistics(area, mean, min, max);
//                imp3.getStatistics();
//                ImageInt draw = new ImageInt(imp3);
                ImageStatistics s = imp3.getRawStatistics();
                int sd = 500;
                seed_threshold = (int)s.mean + sd;
                IJ.log("\n Automatic Seed Threshold estimated is: "+seed_threshold);
                IJ.log("\n Please increase seed threshold and re-run this function using pre-defined threshold if number of achieved objects is less than expected!");
                IJ.log("\n Please decrease seed threshold and re-run this function using pre-defined threshold if number of achieved objects is more than expected! \n");
            }

            //nuclei segmentation
            IJ.log("Segmenting...");
            String seg_option = "seeds_threshold="+seed_threshold+" local_background=0 radius_0=2 radius_1=4 radius_2=6 weigth=0.50 radius_max=20 sd_value=1.8 local_threshold=[Gaussian fit] seg_spot=Block watershed volume_min=100 volume_max=1000000 seeds=Automatic spots=BP_C4-dapi radius_for_seeds=2 output=[Label Image]";
            IJ.run("3D Spot Segmentation",seg_option);
            IJ.selectWindow("seg");
            IJ.saveAs("Tiff", save_dir+"dapi-seg.tif");
            IJ.log("Save the segmented nuclei image as dapi-seg.tif in "+save_dir);
        }
        else
        {
            ImagePlus imp3 = WindowManager.getImage("BP_C4-dapi");
            IJ.log(imp3.getTitle()+" data here");
    //        ImagePlus imp3 = IJ.openImage(dir+"BP_C4-dapi.tif");
//            imp3.setSlice(imp3.getNSlices()/2);
            IJ.run(imp3, "Enhance Contrast...", "saturated=0"); 
            IJ.run(imp3, "16-bit", "");

            //nuclei segmentation
            IJ.log("Segmenting...");
            String seg_option = "seeds_threshold="+seed_threshold+" local_background=0 radius_0=2 radius_1=4 radius_2=6 weigth=0.50 radius_max=20 sd_value=1.8 local_threshold=[Gaussian fit] seg_spot=Block watershed volume_min=100 volume_max=1000000 seeds=Automatic spots=BP_C4-dapi radius_for_seeds=2 output=[Label Image]";
            IJ.run("3D Spot Segmentation",seg_option);
            IJ.selectWindow("seg");
            IJ.log("Program can not read saving directory folder, please manual save image");
        }
          
        
    }        
}
