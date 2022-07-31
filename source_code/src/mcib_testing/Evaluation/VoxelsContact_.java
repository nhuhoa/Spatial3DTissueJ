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
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import java.util.ArrayList;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DPoint;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.regionGrowing.Watershed3DVoronoi;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;

/**
 *
 * @author tranhoa
 */
public class VoxelsContact_ implements PlugIn
{
    // seeds
    ImagePlus seedPlus = null;
    public int min = 5, max = 15;
    boolean save = true;
    String dir = null;
    private Objects3DPopulation popRegions = null;
    @Override
    public void run(String arg) 
    {
        
        if (Dialog())
        {
            getNbContacts();
        }
            
            
    }
    private boolean Dialog() 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return false;
        }
        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        ImagePlus img = WindowManager.getImage(wList[0]);
        Calibration calibration = img.getCalibration();
        String unit = "pixel";
        if (calibration != null) unit = calibration.getUnits();
        GenericDialog gd = new GenericDialog("3D Cell Zone");
        gd.addMessage("3D Tissue Spatial Analysis");
        gd.addMessage(" ");
        gd.addChoice("Nuclei Segmented Image : ", titles, titles[0]);
        gd.addNumericField("Distance min : ", min,  0, 10, "");
        gd.addNumericField("Distance max : ", max,  0, 10, "");
        gd.showDialog();
        int idxImg = gd.getNextChoiceIndex();
        min = (int) gd.getNextNumber();
        max = (int) gd.getNextNumber();
        seedPlus = WindowManager.getImage(wList[idxImg]);
        WindowManager.setTempCurrentImage(seedPlus);
        dir = IJ.getDirectory("image");
        return gd.wasOKed();
    }
    private ImageHandler WatershedVoronoi(float radMax) 
    {
        Watershed3DVoronoi watershed3DVoronoi = new Watershed3DVoronoi(ImageInt.wrap(seedPlus), radMax);
        watershed3DVoronoi.setLabelSeeds(false);
        ImageInt wat = watershed3DVoronoi.getVoronoiZones(false);
        Objects3DPopulation popWat = new Objects3DPopulation(wat, 1);  //exclude value 1 used by border
//        IJ.log("Nb wat regions: " + popWat.getNbObjects());
        ImageHandler drawWat = wat.createSameDimensions();
        popWat.draw(drawWat);
        return drawWat;
    }
        
    
    private void getNbContacts() 
    {
        ResultsTable rtContact = new ResultsTable();
        int c = 0; 
        IJ.log("min: "+ min+" max: "+max);
        for(int k = min; k<max; k+=1)
        {
             
            ImageHandler imgWat = WatershedVoronoi(k);
//            imgWat.show("wat_"+k);
            Objects3DPopulation popRegions = new Objects3DPopulation((ImageInt) imgWat); 
            Objects3DPopulation popNuclei = new Objects3DPopulation(seedPlus);
            ArrayList<Object3D> popObjs = new ArrayList<Object3D>(popRegions.getNbObjects());

            HashMap<Integer, Object3D> region2Cell = new HashMap<Integer, Object3D>(popRegions.getNbObjects());

            // get nucleus label for each region
            for (int i = 0; i < popRegions.getNbObjects(); i++) 
            {
                Object3D reg = popRegions.getObject(i);
                Object3D nuc = popNuclei.getObjectByValue(reg.getValue());
                if (nuc == null) {
                    continue;
                }
                region2Cell.put(reg.getValue(), reg);
            }

            ArrayList<Integer>[] nei = computeAsso(imgWat, 0);

            rtContact.incrementCounter();
                 
            // index refers to region label
            int nbcont = 0;
            int count = 0;
            for (int i = 0; i < nei.length; i++) {
                if (nei[i].isEmpty()) {
                    continue;
                }
                Object3D C = region2Cell.get(i);
                if (C != null) {
                    count++;
                    ArrayList<Object3D> nei1 = new ArrayList<Object3D>();
                    for (int j : nei[i]) {
                        Object3D C1 = region2Cell.get(j);
                        if ((C1 != null) && (C1 != C)) {
                            nei1.add(C1);
                        }
                    }
//                    IJ.log("Cell id "+C.getValue()+" nb conts: "+nei1.size());
                    rtContact.setValue("obj_"+C.getValue(), c, nei1.size());
                    nbcont += nei1.size();
                }
            }
            rtContact.setValue("distance", c, k);
            rtContact.setValue("nbcontacts", c, nbcont);
            rtContact.setValue("nbcontacts_average", c, (float)nbcont/(float)count);
            c++;
        }    
        rtContact.show("Number of contacts");
    }
    private ArrayList<Integer>[] computeAsso(ImageHandler img, int BorderValue) {
        //ImagePlus imp = WindowManager.getCurrentImage();
        //ImageHandler img = ImageHandler.wrap(imp);
        int max = (int) img.getMax();

        ArrayList<Integer>[] nei = new ArrayList[max + 1];
        for (int i = 0; i < nei.length; i++) {
            nei[i] = new ArrayList();
        }

        for (int x = 0; x < img.sizeX; x++) {
            for (int y = 0; y < img.sizeY; y++) {
                for (int z = 0; z < img.sizeZ; z++) {
                    if (((BorderValue >= 0) && (img.getPixel(x, y, z) == BorderValue)) || (BorderValue < 0)) {
                        ArrayUtil tab = img.getNeighborhood3x3x3(x, y, z);
                        tab = tab.distinctValues();
                        for (int j = 0; j < tab.getSize(); j++) {
                            int val = (int) tab.getValue(j);
                            if (!tab.hasOnlyValuesInt(nei[val])) {
                                for (int i = 0; i < tab.getSize(); i++) {
                                    if (!nei[val].contains((int) tab.getValue(i))) {
                                        nei[val].add((int) tab.getValue(i));
//                                        if (val == 9) {
//                                            IJ.log(" " + tab);
//                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return nei;
    }
    
}
