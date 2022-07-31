/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Evaluation;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import mcib3d.image3d.ImageHandler;
import mcib_testing.CellType.VoxS;

/**
 *
 * @author tranhoa
 */
public class Validation_NucSeg_ implements PlugIn
{
    
    public double segVal = 1, gtVal = 10;
    public void run(String string) 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return;
        }

        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
//        String none = "*None*";
//        titles[wList.length] = none;

        GenericDialog gd = new GenericDialog("Evaluation Segmentation");
        gd.addChoice("Segmentation Result ", titles, titles[0]);
        gd.addChoice("Ground truth ", titles, titles[1]);
        gd.addNumericField("Seg Obj Val: ", segVal, 1);
        gd.addNumericField("GT Obj Val: ", gtVal, 10);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        int[] index = new int[2];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        segVal = (int) gd.getNextNumber();
        gtVal = (int) gd.getNextNumber();
        ImageHandler imgSeg = ImageHandler.wrap(WindowManager.getImage(wList[index[0]]));
        ImageHandler imgGT = ImageHandler.wrap(WindowManager.getImage(wList[index[1]]));
        IJ.log("Initialing ...");
        IJ.log("----------------------------------------------------------");
        validate(imgSeg, imgGT);
        IJ.log("----------------------------------------------------------");
        IJ.log("Finished ...");
    }
    
    private void validate(ImageHandler imgSeg, ImageHandler imgGT)
    {
        int TP=0, TN=0, FP=0, FN=0;
        for (int x = 0; x < imgGT.sizeX; x++) {
            for (int y = 0; y < imgGT.sizeY; y++) {
                for (int z = 0; z < imgGT.sizeZ; z++) {
                    if (imgGT.getPixel(x, y, z) == gtVal) 
                    {
                        if(imgSeg.getPixel(x, y, z)==segVal)
                        {
                            TP++;
                        }
                        else{
                            FN++;
                        }
                    } 
                    else 
                    {
                        if(imgSeg.getPixel(x, y, z)==segVal)
                        {
                            FP++;
                        }
                        else{
                            TN++;
                        }
                        
                    }
                }    
            }
        }
        float precision=0, recall=0, accuracy=0;
        precision = TP/(float)(TP+FP);
        recall = TP/(float)(TP+FN);
        accuracy = (TP+TN)/(float)(TP+TN+FP+FN);
        IJ.log("Precision: "+precision+" recall: "+recall+" accuracy: "+accuracy);
    }        
    
}
