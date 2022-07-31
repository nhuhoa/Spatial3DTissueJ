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
//import ij.gui.GenericDialog;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import mcib3d.image3d.ImageHandler;
import mcib_testing.CellType.VoxS;
//import mcib3d.geom.Objects3DPopulation;
//import ij.measure.Calibration;
//import mcib3d.image3d.ImageInt;

/**
 *
 * @author tranhoa
 */
public class Normalise_Imgs_ implements ij.plugin.PlugIn 
{
//    String dir = null, subdir = null;
    public String save_dir = "";
    double thresh = 5;
    ImageHandler img = null;
//    ImageHandler resImg = null;
    public void run(String arg) 
    {
        if (Dialog())
        {
            ImageHandler res = normalizeValuePixels(img, thresh);
            res.show();
            IJ.saveAs("Tiff", save_dir+res.getTitle()+".tif");
            IJ.log("Save the normalized image as "+res.getTitle()+" into folder : "+save_dir);
        }    
            
        
//        int nbima = WindowManager.getImageCount();
//        double minDist = 5;
//        ImageHandler[] imgs = new ImageHandler[nbima];
//        GenericDialog gd = new GenericDialog("NORMALIZE IMAGE");
//        gd.addNumericField("Threhold", minDist, 0);
//        gd.showDialog();
//        if (gd.wasCanceled()) {
//            return;
//        }
//        minDist = (int) gd.getNextNumber();
//        for (int i = 0; i < nbima; i++) {
//            imgs[i] = ImageHandler.wrap(WindowManager.getImage(i + 1));
//            IJ.log("Beginning, normalise the value of pixels of image: " + imgs[i].getTitle());
//            
//            ImageHandler resImg1 = normalizeValuePixels(imgs[i], minDist);
//            resImg1.show();
//        }

        IJ.log("Finished");
    }
       
    private boolean Dialog() 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return false;
        }
        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) 
        {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        ImagePlus temp = WindowManager.getImage(wList[wList.length-1]);;
        if (null !=temp){
            WindowManager.setTempCurrentImage(temp);
            save_dir = IJ.getDirectory("image");
        }else{
            save_dir = "";
        }
        GenericDialogPlus gd = new GenericDialogPlus("Normalizing image");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("3D Tissue Spatial Analysis");
//        gd.addMessage("See and quote reference:\nA novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans\n");
        gd.addMessage("Input : marker image");
        gd.addMessage("Output: normalized image");
        gd.addMessage(" ");
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("Input_Image : ", titles, titles[0]);
        gd.addNumericField("Threshold : ", thresh, 0);
        gd.showDialog();
        
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return false;
        }
        
        int idxImg = gd.getNextChoiceIndex();
        thresh = (int) gd.getNextNumber();
        img = ImageHandler.wrap(WindowManager.getImage(wList[idxImg]));
//        WindowManager.setTempCurrentImage(WindowManager.getImage(wList[idxImg]));
//        dir = IJ.getDirectory("image");
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        
        
        IJ.log("Input image: "+img.getTitle()+" threshold: "+thresh);
        return gd.wasOKed();
    }
    public ImageHandler normalizeValuePixels(ImageHandler imgOriginal, double minDist) 
    {
        ImageHandler imgOrig = imgOriginal.duplicate();
        imgOrig.setTitle("normalized_"+imgOriginal.getTitle());
        // int count = 0;
        ArrayList<VoxS> idxList = new ArrayList<VoxS>();
        //VoxEVF[] idx = new VoxEVF[mask.countMaskVolume()];
        //double volume = idx.length;
        for (int x = 0; x < imgOrig.sizeX; x++) {
            for (int y = 0; y < imgOrig.sizeY; y++) {
                for (int z = 0; z < imgOrig.sizeZ; z++) {
                    if (imgOrig.getPixel(x, y, z) > minDist) 
                    {
                        idxList.add(new VoxS(imgOrig.getPixel(x, y, z), x, y, z));
                    } else 
                    {
                        imgOrig.setPixel(x, y, z, 0);
                    }
                }    
            }
        }
//        IJ.log("nb of pixels: "+idxList.size());
        if (idxList.isEmpty()) {
            return null;
        }
        VoxS[] idx = new VoxS[idxList.size()];
        idx = (VoxS[]) idxList.toArray(idx);
//        double volume = idx.length;
        Arrays.sort(idx);
//        IJ.log("nb of pixels after: "+idx.length);
        double maxV = idx[idx.length - 1].distance;
        IJ.log("Max value is : " + maxV);
        for (int i = 0; i < idx.length; i++) 
        {
            float valRe = (float) ((float)255.0 * (float)(idx[i].distance / maxV));
            //IJ.log(idx[i].distance + " ----- "+ valRe);
            imgOrig.setPixel(idx[i].x, idx[i].y, idx[i].z, valRe);
            
        }
        return imgOrig;
    }
   
    
}
