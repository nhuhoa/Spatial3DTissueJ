/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.PreProcessing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import java.io.File;

/**
 *
 * @author tranhoa
 */
public class Pre_Filtering_windows implements PlugIn
{
    private final int NONE = 0, DELTA = 1, BETA = 2, ALPHA = 3, DAPI = 4;
//    String dir = null, subdir = null;
    public String save_dir = "";
    public void run(String arg) 
    {
        String[] label = {"*NONE*", "DELTA", "BETA","ALPHA", "DAPI"};
        int def = 0;
        if (IJ.versionLessThan("1.37f")) return;
        int[] wList = WindowManager.getIDList();

        if (wList==null) {
                IJ.showMessage("Pre-filtering raw composite images", "There must be at least one image open.");
                return;
        }
        String[] titles = new String[wList.length];
        
        for (int i=0, k=0; i<wList.length; i++) {
                ImagePlus imp = WindowManager.getImage(wList[i]);
                if (null !=imp){
                        titles[k++] = imp.getTitle();
//                        temp = imp;
                }        
        }
        ImagePlus temp = WindowManager.getImage(wList[wList.length-1]);
        if (null !=temp){
            WindowManager.setTempCurrentImage(temp);
            save_dir = IJ.getDirectory("image");
        }else{
            save_dir = "";
        }
        

        boolean prefilter = false; 
        GenericDialogPlus gd = new GenericDialogPlus("3D Filtering");
//        GenericDialog gd = new GenericDialog("3D Filtering");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        "organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : composite raw image");
        gd.addMessage("Output: filtered images");
        gd.addMessage(" \n");
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addChoice("Composite_Image :", titles, titles[0]);
        gd.addChoice("Label_1 : ", label, label[def]);
        gd.addChoice("Label_2 : ", label, label[def]);
        gd.addChoice("Label_3 : ", label, label[def]);
        gd.addChoice("Label_4: ", label, label[def]);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return;
        }
        
        int[] idx = new int[4];
        int i1Index = gd.getNextChoiceIndex();
        idx[0] = gd.getNextChoiceIndex();
        idx[1] = gd.getNextChoiceIndex();
        idx[2] = gd.getNextChoiceIndex();
        idx[3] = gd.getNextChoiceIndex();
        
//        IJ.log("Idx 0: "+idx[0]);
//        IJ.log("Idx 1: "+idx[1]);
//        IJ.log("Idx 2: "+idx[2]);
//        IJ.log("Idx 3: "+idx[3]);
//        IJ.log("Predefined intensity values: NONE = 0, DELTA = 1, BETA = 2, ALPHA = 3, DAPI = 4");
        ImagePlus imp1 = WindowManager.getImage(wList[i1Index]);
//        WindowManager.setTempCurrentImage(imp1);
//        dir = IJ.getDirectory("image");
//        IJ.log("Automatic detecting image folder: ");
//        IJ.log(dir);
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        
        
        IJ.log("Initializing...");
        preFiltering(imp1, idx);
        IJ.log("Completed!");
        IJ.selectWindow("Log");
        IJ.saveAs("Text", save_dir+ "Log_Filtering.txt");
        
    }
    private String getTitle(int label)
    {
//        String titleImg = null;
        String[] label_desc = {null, "C1-delta", "C2-beta","C3-alpha", "C4-dapi"};
//        switch (label) {
//            case 0:  titleImg = null;
//                     break;
//            case 1:  titleImg = "C1-delta";
//                     break;
//            case 2:  titleImg = "C2-beta";
//                     break;
//            case 3:  titleImg = "C3-alpha";
//                     break;
//            case 4:  titleImg = "C4-dapi";
//                     break;
//            default: titleImg = null;
//                     break;
//        }
        return label_desc[label];
    }        
    private void preFiltering(ImagePlus imp1, int[] idx)
    {
        ImagePlus[] channels = ChannelSplitter.split(imp1);
        imp1.changes = false;
        imp1.setIgnoreFlush(true);
        imp1.close();
        IJ.log("DEBUG: ");
        
        if(channels.length > 0)
        {
            IJ.log("Number of channels as input in composite images: " + channels.length);
//            ImagePlus[] imp = new ImagePlus[channels.length];
            for(int i=0; i<channels.length; i++)
            {
                channels[i].show();
//                imp[i] = channels[i];
//                imp[i].show();
//                IJ.log("Title: " + imp[i].getTitle());
                IJ.run(channels[i], "Median 3D...", "x=4 y=4 z=2");
                IJ.log("Save output into folder: " + save_dir);
                if(i<idx.length)
                {
                    IJ.log("Index: "+idx[i]);
                    String title = getTitle(idx[i]);
                    IJ.log("Title: "+title);
                    if(title!=null)
                    {
                        channels[i].setTitle(title);
//                        IJ.selectWindow(title);
                        IJ.saveAs("Tiff", save_dir+title+".tif");
                        IJ.log("Filtered and save image channel_"+(i+1)+" as "+title+".tif");
                    }    
                }
            }    
        }    
        
    }        
    
}
