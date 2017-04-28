/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;

/**
 *
 * @author tranhoa
 */
public class TestDeltaRegion_ implements ij.plugin.PlugIn {
    
    public void run(String arg) 
    {
        int nbima = WindowManager.getImageCount();
        if (nbima < 2) {
            IJ.error("Needs two images");
        }
        String[] namesA = new String[nbima];
        String[] namesB = new String[nbima];
        for (int i = 0; i < nbima; i++) {
            namesA[i] = WindowManager.getImage(i + 1).getShortTitle();
            namesB[i] = WindowManager.getImage(i + 1).getShortTitle();
        }

        GenericDialog dia = new GenericDialog("Extract Region");
        dia.addChoice("Image_A", namesA, namesA[0]);
        dia.addChoice("Image_B", namesB, namesB[1]);


        dia.showDialog();
        if (dia.wasOKed()) {
            int idxA = dia.getNextChoiceIndex();
            int idxB = dia.getNextChoiceIndex();

            ImagePlus plusA = WindowManager.getImage(idxA + 1);
            ImagePlus plusB = WindowManager.getImage(idxB + 1);
            ImageInt imgA = ImageInt.wrap(plusA);
            ImageInt imgB = ImageInt.wrap(plusB);
            IJ.log("Compute the difference between two images A and B");
            ImageHandler labelSimilar = new ImageShort("labelNuclei", imgA.sizeX, imgA.sizeY, imgA.sizeZ);
            int n = 0, count = 0, count1=0;
            for (int z = 0; z < imgA.sizeZ; z++) {
                for (int x = 0; x < imgA.sizeX; x++) {
                    for (int y = 0; y < imgA.sizeY; y++) {
                        //double val = imgA.getPixel(x, y, z);
                        if (imgA.getPixel(x, y, z) > 0 && imgB.getPixel(x, y, z) > 0) 
                        {   
                            labelSimilar.setPixel(x, y, z, imgB.getPixel(x, y, z));
                        }
                    }
                }
            }
            
            labelSimilar.show();
        }    
    }    
    
}
