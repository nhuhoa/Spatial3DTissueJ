/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.SLIC;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
//import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;

/**
 *
 * @author thomasb
 */
public class SLIC3D_3channels_V2 implements PlugIn {

//    String dir = null;
    int thresh = 3;
    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
//    public String dir = null;
    public String save_dir = "";
//    @Override
    public void run(String string) 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return;
        }
        if(wList.length<3)
        {
            IJ.error("At least 3 imgs : marker1.tif, marker2.tif and SLIC clustering label slic.tif should be open.");
            return;
        }    
        
        String[] titles = new String[wList.length+1];
        titles[0] = "*None*";
        ImagePlus temp = null;
        int count=0;
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
//            IJ.log("Idx i: "+i+"  named: "+wList[i]);
            titles[i+1] = imp!=null?imp.getTitle():"";
            if (null !=imp){
                count++;
                if(count==1){
                    temp = imp;
                    WindowManager.setTempCurrentImage(temp);
                    save_dir = IJ.getDirectory("image");
                }
                
            }    
        }
        
        
        GenericDialogPlus gd = new GenericDialogPlus("3D SLIC 3 Channels");
//        gd.addMessage("3D Tissue Spatial Analysis");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        "organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("    Spatial3DTissueJ  ");
        gd.addMessage("Input : markers images");
        gd.addMessage("Input : combined markers SLIC");
        gd.addMessage("Output: SLIC composite label ");
        gd.addMessage(" \n");
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("ALPHA_image: ", titles, titles[0]);
        gd.addChoice("BETA_image: ", titles, titles[0]);
        gd.addChoice("DELTA_image: ", titles, titles[0]);
        gd.addChoice("combined_SLIC_image: ", titles, titles[0]);
        gd.addNumericField("Intensity_threshold: ", thresh, 0, 10, "");
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                IJ.log("Need to input save directory!!!");
                return;
        }
        int[] index = new int[4];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        index[3] = gd.getNextChoiceIndex();
        thresh = (int) gd.getNextNumber();

        IJ.log("SLIC Multi-SuperVoxels Clustering ...");
        IJ.log("Noted: require alpha,beta,SLIC images, or : alpha,beta,delta,SLIC images as input");
        IJ.log("Input Images: ");
        if(index[0]!=0 && index[1]!=0 && index[2]!=0 && index[3]!=0)
        {
            ImageHandler al = ImageHandler.wrap(WindowManager.getImage(wList[index[0]-1]));
            ImageHandler be = ImageHandler.wrap(WindowManager.getImage(wList[index[1]-1]));
            ImageHandler del = ImageHandler.wrap(WindowManager.getImage(wList[index[2]-1]));
            ImageInt label = ImageInt.wrap(WindowManager.getImage(wList[index[3]-1]));
//            WindowManager.setTempCurrentImage(WindowManager.getImage(wList[index[3]-1]));
//            dir = IJ.getDirectory("image");
            IJ.log(" Marker 1 :"+al.getTitle()+"   Marker 2 :"+be.getTitle()+"   Marker 3 :"+del.getTitle()+"   SLIC label :"+label.getTitle());
            ImageHandler labelF = getTypeRegion(al, be, del, label);
            labelF.show();
            IJ.selectWindow(labelF.getTitle());
            IJ.saveAs("Tiff", save_dir + labelF.getTitle()+".tif");
            IJ.log("Save SLIC clustering result as " + labelF.getTitle()+".tif into the folder " + save_dir);
            IJ.log("Output Image: "+ labelF.getTitle());
        } 
        else if(index[0]!=0 && index[1]!=0 && index[2]==0 && index[3]!=0)
        {
            ImageHandler al = ImageHandler.wrap(WindowManager.getImage(wList[index[0]-1]));
            ImageHandler be = ImageHandler.wrap(WindowManager.getImage(wList[index[1]-1]));
            ImageInt label = ImageInt.wrap(WindowManager.getImage(wList[index[3]-1]));
//            WindowManager.setTempCurrentImage(WindowManager.getImage(wList[index[3]-1]));
//            dir = IJ.getDirectory("image");
            IJ.log(" Marker 1 :"+al.getTitle()+"   Marker 2 :"+be.getTitle()+"   SLIC label :"+label.getTitle());
            IJ.log("Noted: do not exist delta channel from this islet, take into account only alpha, beta channels");
            ImageHandler labelF = getTypeRegion(al, be, label);
            labelF.show();
            IJ.saveAs("Tiff", save_dir + labelF.getTitle()+".tif");
            IJ.log("Save SLIC clustering result as " + labelF.getTitle()+".tif into the folder "+save_dir);
            IJ.log("Output Image: "+ labelF.getTitle());
        }
        else{
            IJ.error("Need to choose at least 3 input images: (alpha,beta,delta,SLIC) or (alpha,beta,SLIC) as input");
        }
        IJ.log("Index Value in output Image : UNLABELED = 0, DELTA = 1, BETA = 2, ALPHA = 3");
        IJ.log("Finished");
//        dir = IJ.getDirectory("image");
//        ImageInt label = ImageInt.wrap(WindowManager.getImage("label.tif"));
        
    }
    private ImageHandler getTypeRegion(ImageHandler al, ImageHandler be, ImageHandler del, ImageInt label)
    {
        Objects3DPopulation pop = new Objects3DPopulation(label);
        ImageHandler labelF = new ImageShort("composite_label", label.sizeX, label.sizeY, label.sizeZ);

        for (int i = 0; i < pop.getNbObjects(); i++) 
        {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            double valC1 = O.getPixMedianValue(al);
            double valC2 = O.getPixMedianValue(be);
            double valC3 = O.getPixMedianValue(del);
 
            if ((valC2 > thresh) && (valC2 > valC1) && (valC2 >= valC3)) 
            {
                O.draw(labelF, BETA);
            }
            if ((valC1 > thresh) && (valC1 >= valC2) && (valC1 >= valC3)) 
            {
                O.draw(labelF, ALPHA);             
            }
            
            if ((valC3 > thresh) && (valC3 >= 0.8*valC1) && (valC3 >= 0.8*valC2)) 
            {
                O.draw(labelF, DELTA);
            }
        }
        return labelF;
    }    
    private ImageHandler getTypeRegion(ImageHandler al, ImageHandler be, ImageInt label)
    {
        Objects3DPopulation pop = new Objects3DPopulation(label);
        ImageHandler labelF = new ImageShort("composite_label", label.sizeX, label.sizeY, label.sizeZ);

        for (int i = 0; i < pop.getNbObjects(); i++) 
        {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            double valC1 = O.getPixMedianValue(al);
            double valC2 = O.getPixMedianValue(be);
            if ((valC2 > thresh) && (valC2 > valC1)) 
            {
                O.draw(labelF, BETA);
            }
            if ((valC1 > thresh) && (valC1 >= 0.8*valC2)) 
            {
                O.draw(labelF, ALPHA);             
            }
        }
        return labelF;
    }        
}
