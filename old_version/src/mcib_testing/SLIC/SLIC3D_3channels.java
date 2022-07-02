/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.SLIC;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
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
public class SLIC3D_3channels implements PlugIn {

//    String dir = null;
    int thresh = 3;
    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    public String dir = null;
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
            IJ.error("At least 3 imgs : insulin, glucagon, and SLIC clustering label should be open.");
            return;
        }    
        
        String[] titles = new String[wList.length+1];
        titles[0] = "*None*";
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
//            IJ.log("Idx i: "+i+"  named: "+wList[i]);
            titles[i+1] = imp!=null?imp.getTitle():"";
        }
        GenericDialog gd = new GenericDialog("3D SLIC 3 Channels");
        gd.addMessage("3D Tissue Spatial Analysis");
        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
        "organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : different channels and SLIC label");
        gd.addMessage("Output: SLIC composite label.");
        gd.addMessage(" ");
        gd.addChoice("alpha_img ", titles, titles[0]);
        gd.addChoice("beta_img: ", titles, titles[0]);
        gd.addChoice("delta_image: ", titles, titles[0]);
        gd.addChoice("slic_label: ", titles, titles[0]);
        gd.addNumericField("Threshold: ", thresh, 0, 10, "");
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        int[] index = new int[4];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        index[3] = gd.getNextChoiceIndex();
        thresh = (int) gd.getNextNumber();

        IJ.log("SLIC Multi-SuperVoxels Clustering ...");
        IJ.log("Input Images: ");
        if(index[0]!=0 && index[1]!=0 && index[2]!=0 && index[3]!=0)
        {
            ImageHandler al = ImageHandler.wrap(WindowManager.getImage(wList[index[0]-1]));
            ImageHandler be = ImageHandler.wrap(WindowManager.getImage(wList[index[1]-1]));
            ImageHandler del = ImageHandler.wrap(WindowManager.getImage(wList[index[2]-1]));
            ImageInt label = ImageInt.wrap(WindowManager.getImage(wList[index[3]-1]));
            WindowManager.setTempCurrentImage(WindowManager.getImage(wList[index[3]-1]));
            dir = IJ.getDirectory("image");
            IJ.log(" al :"+al.getTitle()+"   beta :"+be.getTitle()+"   del :"+del.getTitle()+"   label :"+label.getTitle());
            ImageHandler labelF = getTypeRegion(al, be, del, label);
            labelF.show();
            IJ.selectWindow(labelF.getTitle());
            IJ.saveAs("Tiff", dir + labelF.getTitle()+".tif");
            IJ.log("Save SLIC clustering result as " + labelF.getTitle()+".tif into the folder "+dir);
            IJ.log("Output Image: "+ labelF.getTitle());
        } 
        else if(index[0]!=0 && index[1]!=0 && index[2]==0 && index[3]!=0)
        {
            ImageHandler al = ImageHandler.wrap(WindowManager.getImage(wList[index[0]-1]));
            ImageHandler be = ImageHandler.wrap(WindowManager.getImage(wList[index[1]-1]));
            ImageInt label = ImageInt.wrap(WindowManager.getImage(wList[index[3]-1]));
            WindowManager.setTempCurrentImage(WindowManager.getImage(wList[index[3]-1]));
            dir = IJ.getDirectory("image");
            IJ.log(" al :"+al.getTitle()+"   beta :"+be.getTitle()+"   label :"+label.getTitle());
            ImageHandler labelF = getTypeRegion(al, be, label);
            labelF.show();
            IJ.saveAs("Tiff", dir + labelF.getTitle()+".tif");
            IJ.log("Save SLIC clustering result as " + labelF.getTitle()+".tif into the folder "+dir);
            IJ.log("Output Image: "+ labelF.getTitle());
        }
        else{
            IJ.error("Need to choose at least 3 input images: (al, be, del, label) or (al, be, label) as input");
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
 

            if ((valC1 > thresh) && (valC1 >= valC2) && (valC1 >= valC3)) 
            {
                O.draw(labelF, ALPHA);             
            }
            if ((valC2 > thresh) && (valC2 > valC1) && (valC2 >= valC3)) 
            {
                O.draw(labelF, BETA);
            }
            if ((valC3 > thresh) && (valC3 > valC1) && (valC3 > valC2)) 
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
 

            if ((valC1 > thresh) && (valC1 >= valC2)) 
            {
                O.draw(labelF, ALPHA);             
            }
            if ((valC2 > thresh) && (valC2 > valC1)) 
            {
                O.draw(labelF, BETA);
            }
        }
        return labelF;
    }        
}
