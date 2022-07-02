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
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import mcib3d.geom.Object3D;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Vector3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;

/**
 *
 * @author tranhoa
 */
public class EllipsoidFitting_ implements PlugIn
{
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
        GenericDialog gd = new GenericDialog("Ellipsoid Fitting");
        gd.addChoice("Input Image: ", titles, titles[0]);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        int[] index = new int[1];
        index[0] = gd.getNextChoiceIndex();
        ImageHandler img = ImageHandler.wrap(WindowManager.getImage(wList[index[0]]));
        IJ.log("Initialing ...");
        IJ.log("----------------------------------------------------------");
        computeEllipsoidFitting(img);
        IJ.log("----------------------------------------------------------");
        IJ.log("Finished ...");
    }
    private void computeEllipsoidFitting(ImageHandler img)
    {
        Calibration calibration = img.getCalibration();
        if (calibration == null) {
            IJ.log("Image not calibrated");
            calibration = new Calibration();
            calibration.setUnit("pix");
            calibration.pixelWidth = 1;
            calibration.pixelHeight = 1;
            calibration.pixelDepth = 1;
        }

        ImageInt inImage = (ImageInt) img;
        ImageInt segImage;
        if (inImage.isBinary(0)) {
            IJ.log("Segmenting image...");
            inImage = inImage.threshold(0, false, true);
            ImageLabeller labels = new ImageLabeller(false);
            segImage = labels.getLabels(inImage);
            segImage.show("Labelled Image");
        } else {
            segImage = (ImageInt) inImage.duplicate();
        }
        segImage.setCalibration(calibration);
        Objects3DPopulation pop = new Objects3DPopulation();
        pop.addImage(segImage, calibration);
        pop.setCalibration(calibration);
        ObjectCreator3D objs;
        objs = new ObjectCreator3D(img.sizeX, img.sizeY, img.sizeZ);
        for (int i = 0; i < pop.getNbObjects(); i++) 
        {
            Object3D obj = pop.getObject(i);
            Vector3D center = obj.getCenterAsVector();
            double dis = obj.getDistCenterMax();
            objs.createEllipsoid(center, dis, dis, dis, 255);
        }    
        ImagePlus imgBoundingBox = new ImagePlus("BoundingBox_", objs.getStack());
        imgBoundingBox.show();
    }        
}
