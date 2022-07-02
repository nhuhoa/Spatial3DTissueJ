/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.organisation;

import mcib_testing.Utils.Cell;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.ThreadUtil;
import mcib3d.utils.exceptionPrinter;
/**
 *
 * @author tranhoa
 */
public class VoxelSurfaceContact_ implements ij.plugin.PlugIn
{
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg=null, imgWat=null;
    ImageHandler img;
    String dir = null;
    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    
    public void run(String arg) 
    {
        int nbima = WindowManager.getImageCount();

        if (nbima < 1) {
            IJ.showMessage("No image opened !");
            return;
        }
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        img = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        double dilateDistance = 1;
        double minVoxContact = 20;
        GenericDialog gd = new GenericDialog("Voxel Surface Contact");
        gd.addNumericField("Dilate Distance: ", dilateDistance, 1);
        gd.addNumericField("minimum vox contact: ", minVoxContact, 20);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        dilateDistance = (double)gd.getNextNumber();
        minVoxContact = (double)gd.getNextNumber();
        IJ.log("Analysing ...");
        IJ.log("Reading data from " + dir);
        popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        popRegions.loadObjects(dir + "Regions.zip");
        popNuclei.loadObjects(dir + "Nuclei.zip");
        
        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        
        // init cells 
        initCells2();
        IJ.log("size is" + popCells.size());
        IJ.log("Association1 ...");
        IJ.log("Dilate distance: "+dilateDistance + " min voxels contact: "+minVoxContact);
        computeAssoCells(dilateDistance, minVoxContact);
        IJ.log("Finished");
        
    } 
    private void initCells2() {

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer,Cell>(popNuclei.getNbObjects());

        // get nucleus label for each region
        int c = 1;
        for (int i = 0; i < popRegions.getNbObjects(); i++) {
            Object3D reg = popRegions.getObject(i);
            System.out.println("reg=" + reg + " " + reg.getValue());
            Object3D nuc = popNuclei.getObjectByValue(reg.getValue());
            System.out.println("nuc=" + nuc);
            if (nuc == null) {
                continue;
            }
            Cell cell = new Cell();
            cell.region = reg;
            cell.nucleus = nuc;
            popCells.add(cell);
            cell.id = c++;
            region2Cell.put(reg.getValue(), cell);
            nucleus2Cell.put(nuc.getValue(), cell);
            cell.type = (byte) nuc.getType();
        }
    }
    
    private void computeAssoCells(double dilateDistance, double minNbVoxContact)
    {
//        for(Cell C : popCells)
        for(int i=0; i<popCells.size()/10; i++)
        {
            Cell C = popCells.get(i);
            ImageHandler label = new ImageShort("reg_", imgWat.sizeX, imgWat.sizeY, imgWat.sizeZ);
            Object3D reg = C.region;
            int regVal = reg.getValue();
            reg.draw(label, regVal);
//            label.show();
            ImageByte imgDilate = dilateRegion(dilateDistance, reg);
            imgDilate.substractMask((ImageInt)label);
//            imgDilate.show("substract_img");
            ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
            for(int x=0; x < imgDilate.sizeX;x++)
            {
                for(int y=0; y<imgDilate.sizeY;y++)
                {
                    for(int z=0; z<imgDilate.sizeZ;z++)
                    {
                        if(imgDilate.getPixel(x, y, z) > 0)
                        {
                            Voxel3D v = new Voxel3D(x, y, z, regVal);
                            arr.add(v);
                        }
                    }    

                }    
            } 
            Object3DVoxels oo = new Object3DVoxels(arr);
            ArrayUtil ac = oo.listValues(imgWat); 
//            IJ.log("Cell_"+(i+1));
            HashMap<Integer, Integer> map = countNbVoxContact(ac);
            Iterator<Integer> keySetIterator = map.keySet().iterator();
            C.nei1 = new ArrayList<Cell>();
            while(keySetIterator.hasNext())
            {
                Integer key = keySetIterator.next();
                if(map.get(key)>minNbVoxContact)
                {
//                    IJ.log("val: "+key+" "+map.get(key));
                    Cell C1 = region2Cell.get(key);
                    if ((C1 != null) && (C1 != C)) 
                    {
                        C.nei1.add(C1);
                    }
                }    
            }
            IJ.log("Cell_"+(i+1)+"  nb of neighbors: "+ C.nei1.size());

        }    
        
    } 
    public ImageByte dilateRegion(double microDistance, Object3D region)
    {
        ImageHandler label = new ImageShort("dilate_", img.sizeX, img.sizeY, img.sizeZ);
        region.draw(label, region.getValue());
        try {
            int nbCPUs = ThreadUtil.getNbCpus();
            float rDilate = (float) microDistance;
            float radiusZ = (float) rDilate/2;  //minimum should be 1
            ImageFloat edm = EDT.run(label, 0, 1, rDilate / radiusZ, true, nbCPUs);
            ImageByte temp = edm.threshold(rDilate, true, false);
            edm.flush();
            edm = null;
            temp.setOffset(label);
            temp.setScale(label);
            for(int xy=0; xy < temp.sizeXY;xy++)
            {
                for(int z=0; z<temp.sizeZ;z++)
                {
                    if(temp.getPixel(xy, z) > 0)
                    {
                        temp.setPixel(xy, z, region.getValue());
                    }    
                }    
            }    
            return temp;
        } catch (Exception e) {
            exceptionPrinter.print(e, null, true);
        }
        return null;
    }
    public HashMap countNbVoxContact(ArrayUtil ac) {
        ac.sort();
        HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
        int s = 0;
        double tmp = ac.getValue(0);
        s++;
        int p = 1;
        int c = 1;
        int si = ac.getSize();
        while (p < si) {
            while ((p < si) && (ac.getValue(p) == tmp)) {
                p++;
                c++;
            }
            map.put((int)tmp, c);
            c = 0;
            if (p < si) {
                tmp = ac.getValue(p);
                s++;
                p++;
                c++;
            }
        }
        return map;
//        IJ.log(MapUtils.toString(map));
        
    }
}
