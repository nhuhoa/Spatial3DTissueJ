/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.SLIC;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.filter.PlugInFilter;
import static ij.plugin.filter.PlugInFilter.DOES_16;
import static ij.plugin.filter.PlugInFilter.DOES_8G;
import ij.process.ImageProcessor;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;

/**
 *
 * @author thomasb
 */
public class SLIC3D_4channels implements PlugInFilter {

    ImagePlus plus;
    String dir = null;
    int thresh = 0;
    @Override
    public int setup(String string, ImagePlus ip) {
        plus = ip;

        return DOES_16 + DOES_8G;
    }

    @Override
    public void run(ImageProcessor ip) 
    {
        ImageHandler img1 = ImageHandler.wrap(WindowManager.getImage("C1-delta.tif"));
        ImageHandler img2 = ImageHandler.wrap(WindowManager.getImage("C2-beta.tif"));
        ImageHandler img3 = ImageHandler.wrap(WindowManager.getImage("C1-alpha.tif"));
        ImageHandler img4 = ImageHandler.wrap(WindowManager.getImage("C4-dapi.tif"));
        dir = IJ.getDirectory("image");
        ImageInt label = ImageInt.wrap(WindowManager.getImage("label.tif"));
        Objects3DPopulation pop = new Objects3DPopulation(label);
        ImageHandler labelF = new ImageShort("labelF", img1.sizeX, img1.sizeY, img1.sizeZ);

        for (int i = 0; i < pop.getNbObjects(); i++) {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            double valC1 = O.getPixMedianValue(img1);
            double valC2 = O.getPixMedianValue(img2);
            double valC3 = O.getPixMedianValue(img3);
            double valC4 = O.getPixMedianValue(img4);

            if ((valC1 > thresh) && (valC1 >= valC2) && (valC1 >= valC3)&& (valC1 >= valC4)) 
            {
                O.draw(labelF, 1);             
            }
            if ((valC2 > thresh) && (valC2 > valC1) && (valC2 > valC3)&& (valC2 > valC4)) 
            {
                O.draw(labelF, 2);
            }
            if ((valC3 > thresh) && (valC3 > valC1) && (valC3 > valC2)&& (valC3 > valC4)) 
            {
                O.draw(labelF, 3);
            }
            if ((valC4 > thresh) && (valC4 > valC1) && (valC4 > valC2)&& (valC4 > valC3)) 
            {
                O.draw(labelF, 4);
            }
        }
        labelF.show();
    }
}
