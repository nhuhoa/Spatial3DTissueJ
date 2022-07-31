/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.PreProcessing;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import java.io.File;

/**
 *
 * @author tranhoa
 */
public class Pre_Filtering_V2 implements PlugIn
{
//    private final int NONE = 0, DELTA = 1, BETA = 2, ALPHA = 3, DAPI = 4;
    String dir = null, subdir = null;
    public void run(String arg) 
    {
        String[] label = {"*NONE*", "M1", "M2","M3", "DAPI"};
        int def = 0;
        if (IJ.versionLessThan("1.37f")) return;
        int[] wList = WindowManager.getIDList();

        if (wList==null) {
                IJ.showMessage("Pre-filtering", "There must be at least one image open.");
                return;
        }
        String[] titles = new String[wList.length];
        for (int i=0, k=0; i<wList.length; i++) {
                ImagePlus imp = WindowManager.getImage(wList[i]);
                if (null !=imp)
                        titles[k++] = imp.getTitle();
        }
        
        boolean prefilter = false; 
        GenericDialog gd = new GenericDialog("3D Filtering");
        gd.addMessage("3D Tissue Spatial Analysis");
        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
        "organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : composite image of many channels");
        gd.addMessage("Output: filtered images");
        gd.addChoice("Composite Image :", titles, titles[0]);
        gd.addChoice("Label_1 : ", label, label[def]);
        gd.addChoice("Label_2 : ", label, label[def]);
        gd.addChoice("Label_3 : ", label, label[def]);
        gd.addChoice("Label_4: ", label, label[def]);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
        int[] idx = new int[4];
        int i1Index = gd.getNextChoiceIndex();
        idx[0] = gd.getNextChoiceIndex();
        idx[1] = gd.getNextChoiceIndex();
        idx[2] = gd.getNextChoiceIndex();
        idx[3] = gd.getNextChoiceIndex();
        
        ImagePlus imp1 = WindowManager.getImage(wList[i1Index]);
        WindowManager.setTempCurrentImage(imp1);
        dir = IJ.getDirectory("image");
        
        IJ.log("Initializing...");
        File fl = new File(dir+"filtered");
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = dir + "filtered/";
            } else {
                subdir = dir;
            }
        }
        else{
            subdir = dir + "filtered/";
        }
        preFiltering(imp1, idx);
        IJ.log("Finished.");
        
    }
    private String getTitle(int label)
    {
        String titleImg = "UNKNOWN";
        switch (label) {
            case 0:  titleImg = null;
                     break;
            case 1:  titleImg = "C1";
                     break;
            case 2:  titleImg = "C2";
                     break;
            case 3:  titleImg = "C3";
                     break;
            case 4:  titleImg = "C4-dapi";
                     break;
            default: titleImg = "UNKNOWN";
                     break;
        }
        return titleImg;
    }        
    private void preFiltering(ImagePlus imp1, int[] idx)
    {
        ImagePlus[] channels = ChannelSplitter.split(imp1);
        imp1.changes = false;
        imp1.setIgnoreFlush(true);
        imp1.close();
        if(channels.length > 0)
        {
            IJ.log("Number of channels: " + channels.length);
//            ImagePlus[] imp = new ImagePlus[channels.length];
            for(int i=0; i<channels.length; i++)
            {
                channels[i].show();
//                imp[i] = channels[i];
//                imp[i].show();
//                IJ.log("Title: " + imp[i].getTitle());
                IJ.run(channels[i], "Median 3D...", "x=4 y=4 z=2");
                IJ.log("Save output into folder: "+subdir);
                if(i<idx.length)
                {
                    String title = getTitle(idx[i]);
                    
                    if(title!=null)
                    {
                        IJ.selectWindow(title);
                        IJ.saveAs("Tiff", subdir+title+".tif");
                        IJ.log("Filtered and save channel_"+(i+1)+" as "+title+".tif");
                    }    
                }
            }    
        }    
        
    }        
    
}
