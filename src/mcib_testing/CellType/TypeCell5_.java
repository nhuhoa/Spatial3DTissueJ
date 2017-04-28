/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;

/**
 *
 * @author ttnhoa
 */
public class TypeCell5_ implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;
    private final int COLOCAB = 23;
    private final int COLOCAD = 13;
    private final int COLOCBD = 12;
   
    
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    
    ImageInt imgLabel;
    String dir = null;
    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    ImageInt imgWat=null, imgSeg=null;
    ImageFloat imgExpand =null;
    @Override
    public void run(String arg) 
    {
        imgLabel = ImageInt.wrap(WindowManager.getImage("labelF.tif"));
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        //IJ.log("dir is: "+dir);
        double ratio = 0.1;
        GenericDialog gd = new GenericDialog("SLIC");
        gd.addNumericField(" Ratio ", ratio, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        ratio = (double) gd.getNextNumber();
        imgWat = ImageInt.wrap(plus);
        imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
  
        IJ.log("Initialization Tissue Analysis ...");
        initCells(imgSeg, imgWat);
        //testExpand();
        //imgExpand = (ImageFloat) ImageHandler.wrap(WindowManager.getImage("final.tif"));
        IJ.log("Type ...");
        this.computeTypeCellsNucleus(imgLabel, ratio);
        
//        IJ.log("Save the result");
//        popNuclei.saveObjects(dir + "Nuclei.zip");
//        popRegions.saveObjects(dir + "Regions.zip");
//        
        
//        // draw 
        //drawCellTypes(false).show("TYPE");
        IJ.log("Finished");

    }

    private void initCells(ImageInt nucLabel, ImageInt regionLabel) 
    {
        popNuclei = new Objects3DPopulation(nucLabel);
        popRegions = new Objects3DPopulation(regionLabel, 1); // exclude value 1 used by borders

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer, Cell>(popNuclei.getNbObjects());

        // get nucleus label for each region
        int c = 1;
        //int count = 0;
        for (Object3D region : popRegions.getObjectsList()) 
        {
            int nuc = (int) region.getPixModeNonZero(nucLabel);
            //IJ.log("nuc " + nuc);
            if(nuc==-1){ continue;}
            else
            {
                Cell cell = new Cell();
                cell.region = region;
                cell.nucleus = popNuclei.getObjectByValue(nuc);
                popCells.add(cell);
                cell.id = c++;
                region2Cell.put(region.getValue(), cell);
                nucleus2Cell.put(nuc, cell);
            }
            
                  
        }
        //IJ.log("number of error: " + count);
        IJ.log("Number of cells is: "+popCells.size());
    }

    public void normalizeDistanceMap(ImageFloat distanceMap, ImageInt mask, float ratioUse) 
    {
        int count = 0;
        VoxS[] idx = new VoxS[mask.countMaskVolume()];
        double volume = idx.length;
        for (int z = 0; z < distanceMap.sizeZ; z++) {
            for (int y = 0; y < distanceMap.sizeY; y++) {
                for (int x = 0; x < distanceMap.sizeX; x++) {
                    if (mask.getPixelInt(x, y, z) != 0) {
                        idx[count] = new VoxS(distanceMap.getPixel(x, y, z), x, y, z);
                        count++;
                    }
                }    
            }
        }
        Arrays.sort(idx);

        for (int i = 0; i < idx.length - 1; i++) {
            // gestion des repetitions
            if (idx[i + 1].distance == idx[i].distance) {
                int j = i + 1;
                while (j < (idx.length - 1) && idx[i].distance == idx[j].distance) {
                    j++;
                }
                double median = (i + j) / 2d;
                for (int k = i; k <= j; k++) {
                    idx[k].index = median;
                }
                i = j;
            } else {
                idx[i].index = i;
            }
        }
        if (idx[idx.length - 1].index == 0) {
            idx[idx.length - 1].index = idx.length - 1;
        }
        for (VoxS idx1 : idx) 
        {
            if((float) (idx1.index / volume) <= ratioUse)
            {
                distanceMap.setPixel(idx1.x, idx1.y, idx1.z, (float) (idx1.index / volume) * (float)1000.0);
            }
            else{
                distanceMap.setPixel(idx1.x, idx1.y, idx1.z, 0);
            }
            
        }
    }
    public ImageFloat getDomain(ImageHandler img)
    {
        boolean inverse = true;
        int threshold = 1;
        ImageFloat r = EDT.run(img, threshold, inverse, Runtime.getRuntime().availableProcessors());
        if (r != null) 
        {   
            ImageFloat r2 = r.duplicate(); 
            normalizeDistanceMap(r2, img.threshold(threshold, inverse, true), (float) 0.05);
            return r2;
        }
        return null;

    } 
//    public void testExpand()
//    {
//        ImageFloat imgBorder = getDomain(imgSeg);
//        imgBorder.show("after EDT 10");
//        for(Cell C : popCells)
//        {
//            Object3D nuc = C.nucleus;
//            nuc.draw(imgBorder, 200);
//        }
//        ImageFloat imgBorder1 = imgBorder.duplicate();
//        imgBorder1.show("result contour");
//        imgBorder.intersectMask(imgWat);
//        imgBorder.show("result final");
//    }        
    public int getCellType2(int vA, int vB, int vD, double vol, double ratio)
    {
        //double ratioDel = 0.25;
        //double ratioDel = 0.6;
        if(vA >= vB && vA > vD && vA/vol >= ratio)
        {
            return ALPHA;
//            if(vD >= vA * 0.5 && vD/vol >= ratio)
//            {
//                return DELTA;
//            }
//            else{
//                return ALPHA;
//            }
        } 
        else if(vB > vA && vB > vD && vB/vol >= ratio)
        {
            if(vD >= vB * 0.9 && vD/vol >= ratio)  //0.5
            {
                return DELTA;
            }
            else if(vA>=vB*0.8 && vA>vD && vA/vol >= ratio)
            { 
                return ALPHA;
            }
            else{
                return BETA;
            }
        }
        else if(vD >= vA && vD >= vB && vD/vol >= ratio)
        {
            return DELTA;
        }
//        else if(vA==vB && vA>vD && vA/vol >= ratio)
//        { 
//            return ALPHA;
//            //return BETA;
//        }
        else
        {
            return UNLABELED;
        }  
        
    }
//    public int getCellType(int vA, int vB, int vD)
//    {
//        if(vA > vB && vA > vD)
//        {
//            if(vD >= vA * 0.5)
//            {
//                return DELTA;
//            }
//            else if(vB >= vA*0.95)
//            {
//                return COLOCAB;
//            }
//            else if(vD >= vA*0.7)
//            {
//                return COLOCAD;
//            }
//            else{
//                return ALPHA;
//            }
//        }
//        else if(vB > vA && vB > vD)
//        {
//            if(vD >= vB * 0.45)
//            {
//                return DELTA;
//            }
//            else if(vA >= vB*0.8)
//            {
//                return ALPHA;
//            }
//            else if(vA >= vB*0.6)
//            {
//                return COLOCAB;
//            }
//            else if(vD >= vB*0.35)
//            {
//                return COLOCBD;
//            }
//            else{
//                return BETA;
//            }
//        }
//        else if(vD > vA && vD > vB)
//        {
//            return DELTA;
//        }
//        else if(vA==vB && vA>vD)
//        { 
//            return COLOCAB;
//        }
//        else if((vA==vD && vA>vB))
//        {
//            return COLOCAD;
//        }
//        else if(vD==vB && vD>vA)
//        {
//            return COLOCBD;
//        } 
//        else
//        {
//            return UNLABELED;
//        }  
//        
//    }      
    private void computeTypeCellsNucleus(ImageInt label, double ratio) 
    {
        int cA=0, cB=0, cD=0, cUnlabelled=0, cColocAB=0, cColocAD=0, cColocBD=0;
        ResultsTable cellM = new ResultsTable();
        imgExpand = getDomain(imgSeg);
//        ImageFloat imgExpand = getDomain(imgSeg);
//        for(Cell C : popCells)
//        {
//            Object3D nuc = C.nucleus;
//            nuc.draw(imgExpand, 200);
//        }
//        
//        imgExpand.intersectMask(imgWat);
//        imgExpand.show("result final");
            for (Cell C : popCells) 
            {
                Object3D nucleus = C.nucleus;
                if(nucleus==null) continue;
                double vol = nucleus.getVolumeUnit();
                Object3D region = C.region;
                int reVal = region.getValue();
                ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
                int xmin0 = region.getXmin();
                int ymin0 = region.getYmin();
                int zmin0 = region.getZmin();
                int xmax0 = region.getXmax();
                int ymax0 = region.getYmax();
                int zmax0 = region.getZmax();
                for (int k1 = zmin0; k1 <= zmax0; k1++) {
                    for (int j1 = ymin0; j1 <= ymax0; j1++) {
                        for (int i1 = xmin0; i1 <= xmax0; i1++) {
                            if (imgExpand.contains(i1, j1, k1) && imgExpand.getPixel(i1, j1, k1) > 0) 
                            {
                                Voxel3D v1 = new Voxel3D(i1, j1, k1, reVal);
                                arr.add(v1);
                            }
                        }
                    }
                }
                Object3DVoxels oo = new Object3DVoxels(arr);
                //ArrayUtil ac1 = region.listValues(label); 
                ArrayUtil ac1 = oo.listValues(label); 
                
                cellM.incrementCounter();
                cellM.setValue("cell", C.id, reVal);
                cellM.setValue("id", C.id, C.id);
                
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<ac1.getSize();i++)
                {
                    if(ac1.getValue(i)==ALPHA)nbA++;
                    if(ac1.getValue(i)==BETA)nbB++;
                    if(ac1.getValue(i)==DELTA)nbD++;                
                }
                C.type = UNLABELED;
                //double ratio = 0.6;
                int typeCell = getCellType2(nbA, nbB, nbD, vol, ratio);
                C.type = (byte)typeCell;
                cellM.setValue("method5", C.id, typeCell);
                if(typeCell==UNLABELED){cUnlabelled++;}
                if(typeCell==ALPHA){
                    cA++;
                    //IJ.log("al: "+100 * (nbA/vol));
                }
                if(typeCell==BETA){
                    cB++;
                    //IJ.log("be: "+100 * (nbB/vol));
                }
                if(typeCell==DELTA){
                    cD++;
                    //IJ.log("del: "+100 * (nbD/vol));
                }
                if(typeCell==COLOCAB){cColocAB++;}
                if(typeCell==COLOCAD){cColocAD++;}
                if(typeCell==COLOCBD){cColocBD++;}
                String content = nbA +","+nbB +","+nbD +","+typeCell +"\n";
                //out.write(content);
                nucleus.setType(C.type);
                nucleus.setName("Nuc" + C.id);
                nucleus.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
            }
//            out.close();
//        }
//        catch (IOException ex) 
//        {
//            Logger.getLogger(TypeCell1_.class.getName()).log(Level.SEVERE, null, ex);
//        }
        //cellM.show("Method5");
        
//        IJ.log("Nb of each cells type:  A: "+cA+"   B: "+cB+"   D:"+cD+"   Unlabelled: "+cUnlabelled);
//        ImageHandler unLabelledCell;
//        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);
//
//        for (Cell C : popCells) {
//            if(C.type==UNLABELED){
//               Object3D nucleus = C.nucleus;
//               nucleus.draw(unLabelledCell, 255);
//            }     
//        }
//        unLabelledCell.show("Unlabelled5");
        IJ.log("Nb of each cells type: A: "+cA+"   B: "+cB+"   D:"+cD+"   Unlabelled: "+cUnlabelled);
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);
        ImageHandler typeImg;
        typeImg = new ImageShort("TYPE", label.sizeX, label.sizeY, label.sizeZ);
        for (Cell C : popCells) {
            Object3D nucleus = C.nucleus;
            if(C.type==UNLABELED)
            {  
               nucleus.draw(unLabelledCell, 255);
            }
            else
            {
                nucleus.draw(typeImg, C.type);
            }
        }
        unLabelledCell.show();
        typeImg.show();
    }

    
    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);
//        ImageHandler drawDelta = new ImageShort("DELTA", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);
//        ImageHandler drawDeltaCenter = new ImageShort("DELTA CENTER", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);
        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.type);
            } else {
//                if(C.type==DELTA)
//                {
//                    C.region.draw(drawDelta, C.type);
//                    Point3D center = C.nucleus.getCenterAsPoint();
//                    Voxel3D v = new Voxel3D(center, 255);
//                    arr.add(v);
//                }
                C.region.draw(draw, C.type);
            }
        }
//        Object3DVoxels obj = new Object3DVoxels(arr);
//        obj.draw(drawDeltaCenter);
//        drawDelta.show();
//        drawDeltaCenter.show();
        return draw;
    }

}
