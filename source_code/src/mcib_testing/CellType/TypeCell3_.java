/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;

//import mcib_testing.RafaelV2.*;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
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
import mcib_testing.Utils.Cell;

/**
 *
 * @author tranhoa
 */
/**
 * Using the percentage to compute the percentage of cell type. 
 * @author tranhoa
 */
public class TypeCell3_ implements ij.plugin.PlugIn
{
    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int COLOCAB = 23;
    private final int COLOCAD = 13;
    private final int COLOCBD = 12;
    private final int COLOCABD = 123;
    ImageInt imgSeg, label, imgWat;
    int thresh = 0;
    String dir;
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;
    ImageHandler dapi, alpha, beta, delta;
    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    ArrayUtil idxCR2 = null;
    ArrayUtil idxRe = null;
    public void run(String arg) 
    {
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        imgWat = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        label = ImageInt.wrap(WindowManager.getImage("label.tif"));
        imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
        
        delta = ImageHandler.wrap(WindowManager.getImage("C1-delta.tif"));
        beta = ImageHandler.wrap(WindowManager.getImage("C2-beta.tif"));
        alpha = ImageHandler.wrap(WindowManager.getImage("C1-alpha.tif"));
        dapi = ImageHandler.wrap(WindowManager.getImage("C4-dapi.tif"));
      

        IJ.log("Initialization Type Cell Analysis ...");
        initCells(imgSeg, imgWat);
        mappingRC();
        IJ.log("Saving the result");
        popNuclei.saveObjects(dir + "method3/Nuclei.zip");
        popRegions.saveObjects(dir + "method3/Regions.zip");
        IJ.log("Drawing type cell");
        drawCellTypes(false).show("TYPE");
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
    }
    
    
    
    /**
     * For each of region, find the type of this region
        One of 5 type : alpha, beta, delta, dapi, unlabelled
        in case of nb of alpha=beta || alpha=delta || beta=delta not verify
        * mean value extracted to find the posibility of each signal.
     * @param feature
     * @return
     */
    public ArrayUtil getTypeRegion()
    {
        //ImageInt labelR = new ImageByte("label", label.sizeX, label.sizeY, label.sizeZ);
        //Objects3DPopulation pop = new Objects3DPopulation(label);
        Objects3DPopulation pop = new Objects3DPopulation();
        pop.addImage(label, -1, new Calibration());
        int nbRegion = pop.getNbObjects();
        ArrayUtil idx = new ArrayUtil(nbRegion);
        //int index = 0;
        //ResultsTable rt = new ResultsTable();
        for (int i = 0; i < nbRegion; i++) 
        {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
//            ArrayList<Voxel3D> cont = O.getContours();
//            if (cont.isEmpty()) {
//                IJ.log("Empty cont " + O);
//            }
//            Object3DVoxels vox = new Object3DVoxels();
//            vox.addVoxels(cont);
//            vox.draw(labelR, 255);
            int value = O.getValue();
            //int valType = 0;
            
            double valC1 = O.getPixMedianValue(delta);
            double valC2 = O.getPixMedianValue(beta);
            double valC3 = O.getPixMedianValue(alpha);
            double valC4 = O.getPixMedianValue(dapi);
            
            int typeRegion = getTypeR(valC1, valC2, valC3, valC4);
            idx.putValue(value, typeRegion);
            //O.draw(labelR, typeRegion);
//            rt.incrementCounter();
//            rt.setValue("value", i+1, value);
//            rt.setValue("m delta", i+1, valC1);
//            rt.setValue("m beta", i+1, valC2);
//            rt.setValue("m alpha", i+1, valC3);
//            rt.setValue("m dapi", i+1, valC4);
//            rt.setValue("type", i+1, typeRegion);
//            //rt.setValue("vol pixel", i+1, O.getVolumePixels());
//            rt.setValue("vol unit", i+1, O.getVolumeUnit());
//            rt.setValue("diameter", i+1, O.getDistCenterMean());
//            rt.setValue("elongation", i+1, O.getMainElongation());
        }
//        for (Cell C : popCells)
//        {
//            Object3D region = C.region;
//            region.setLabelImage(imgWat);
//            ArrayList<Voxel3D> cont = region.getContours();
//            Object3DVoxels vos = new Object3DVoxels();
//            vos.addVoxels(cont);
//            vos.draw(labelR, 255);
//            region.setLabelImage(null);
//        } 
        //labelR.show("label region");
        //rt.show("Region Attribute");
        return idx;
    }
    
    
    public int getTypeR(double valC1, double valC2, double valC3, double valC4)
    {
        int valType = 0;
        if ((valC1 > thresh) && (valC1 > valC2) && (valC1 > valC3) && (valC1 > valC4)) 
        {
            valType = 1;
        }
        else if ((valC2 > thresh) && (valC2 > valC1) && (valC2 > valC3) && (valC2 > valC4)) 
        {
            valType = 2;
        }
        else if ((valC3 > thresh) && (valC3 > valC1) && (valC3 > valC2) && (valC3 > valC4)) 
        {
            valType = 3;
        }
        else if ((valC4 > thresh) && (valC4 > valC1) && (valC4 > valC2) && (valC4 > valC3)) 
        {
            valType = 4;
        }
        else if ((valC4 > thresh) && (valC4 == valC1) && (valC4 > valC3) && (valC4 > valC2)) 
        {
            valType = 41;
        }
        else if ((valC4 > thresh) && (valC4 == valC2) && (valC4 > valC3) && (valC4 > valC1)) 
        {
            valType = 42;
            }
        else if ((valC4 > thresh) && (valC4 == valC3) && (valC4 > valC2) && (valC4 > valC1)) 
        {
            valType = 43;
        }
        else if ((valC1 > thresh) && (valC1 == valC2) && (valC1 > valC3) && (valC1 > valC4)) 
        {
            valType = 12;
        }
        else if ((valC1 > thresh) && (valC1 == valC3) && (valC1 > valC2) && (valC1 > valC4)) 
        {
            valType = 13;
        }
        else if ((valC2 > thresh) && (valC2 == valC3) && (valC2 > valC1) && (valC2 > valC4)) 
        {
            valType = 23;
        }
        else
        {
            valType = 0;
        }
        
        return valType;
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
            normalizeDistanceMap(r2, img.threshold(threshold, inverse, true), (float) 0.1);
            return r2;
        }
        return null;

    }  
        
    public void mappingRC()
    {
        ImageFloat imgBorder = getDomain(imgSeg);
//        for (Cell C : popCells)
//        {
//            Object3D region = C.region;
//            region.setLabelImage(imgWat);
//            ArrayList<Voxel3D> cont = region.getContours();
//            Object3DVoxels vos = new Object3DVoxels();
//            vos.addVoxels(cont);
//            vos.draw(imgBorder, 255);
//            region.setLabelImage(null);
//        } 
        for(Cell C : popCells)
        {
            Object3D nuc = C.nucleus;
            nuc.draw(imgBorder, 200);
        }
        ImageFloat imgBorder1 = imgBorder.duplicate();
        imgBorder1.show("result contour");
        imgBorder.intersectMask(imgWat);
        imgBorder.show("result final");
        
        
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        //ArrayUtil idxRe = getTypeRegion();
        idxRe = getTypeRegion();
        IJ.log("Get type region, size of array is: " + idxRe.getSize());
        
        ResultsTable mapR = new ResultsTable();
        ResultsTable cellM = new ResultsTable();
        try {
            FileWriter ryt = new FileWriter(dir+"method3/method3.arff");
            BufferedWriter out=new BufferedWriter(ryt);
            out.write("@relation cell\n");
            out.write("@attribute nbA numeric\n");
            out.write("@attribute nbB numeric\n");
            out.write("@attribute nbD numeric\n");
            out.write("@attribute vA numeric\n");
            out.write("@attribute vB numeric\n");
            out.write("@attribute vD numeric\n");
            out.write("@attribute class {0,1,2,3,23,13,12}\n");
            out.write("@data\n");
            
            for (int k = 0; k < popCells.size(); k++) 
            {
                Cell C = popCells.get(k);
                Object3D nuc = C.nucleus;
                if(nuc==null) continue;
                Object3D region = C.region;
                //Object3DVoxels sRegion = new Object3DVoxels();
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
                            if (imgBorder.contains(i1, j1, k1) && imgBorder.getPixel(i1, j1, k1)!=0.0) 
                            {
                                Voxel3D v1 = new Voxel3D(i1, j1, k1, region.getValue());
                                arr.add(v1);
                            }
                        }
                    }
                }
                Object3DVoxels oo = new Object3DVoxels(arr);
                //ArrayUtil ac1 = region.listValues(label); 
                ArrayUtil ac1 = oo.listValues(label); 
                ac1 = ac1.distinctValues();
                int reVal = region.getValue();
                mapR.incrementCounter();
                mapR.setValue("cell", C.id, reVal);
                cellM.incrementCounter();
                cellM.setValue("cell", C.id, reVal);
                
                int nbA=0, nbB=0, nbD=0, nbNuc=0;
                int vA=0, vB=0, vD=0, vNuc=0;
                for(int j=1; j<=20; j++){
                    mapR.setValue("idx"+j, C.id, 0);
                }
                int count = 0;
                for (int i = 0; i < ac1.getSize(); i++) 
                {
                    int indexRegion = ac1.getValueInt(i);
                       
                    int type = idxRe.getValueInt(indexRegion);
                    //IJ.log("Index region is: " + indexRegion+" type is: "+ type);
                    Object3DVoxels intersectR = new Object3DVoxels();
                    //double vol = 0;
                    if(type!=0)
                    {    
                        if(popR.getObjectByValue(indexRegion)!=null)
                        {
                            Object3D re = popR.getObjectByValue(indexRegion);
                            //double volObj = re.getVolumeUnit();
                            intersectR.addVoxelsIntersection(re, region);
                            double vol = intersectR.getVolumeUnit();
                            if(vol >= 300)
                            {
                                count++;
                                mapR.setValue("idx"+count, C.id, vol);
                                if(type==1){nbD++; vD += vol;}
                                if(type==2){nbB++; vB += vol;}
                                if(type==3){nbA++; vA += vol;}
                                if(type==4){nbNuc++; vNuc +=vol;}
                                if(type==41){nbD++; nbNuc++; vD += vol; vNuc +=vol;}
                                if(type==42){nbB++; nbNuc++; vB += vol; vNuc +=vol;}
                                if(type==43){nbA++; nbNuc++; vA += vol; vNuc +=vol;}
                                if(type==12){nbD++; nbB++; vD += vol; vB += vol;}
                                if(type==13){nbD++; nbA++; vD += vol; vA += vol;}
                                if(type==23){nbB++; nbA++; vA += vol; vB += vol;}
                            }    
                            
                                //if(type==0){nbUnlabelled++;}
                        }
                    }
                     
                }
                mapR.setValue("nbNuc", C.id, nbNuc);
                mapR.setValue("nbA", C.id, nbA);
                mapR.setValue("nbB", C.id, nbB);
                mapR.setValue("nbD", C.id, nbD);
                mapR.setValue("vA", C.id, vA);
                mapR.setValue("vB", C.id, vB);
                mapR.setValue("vD", C.id, vD);
                mapR.setValue("vNuc", C.id, vNuc);
                
                int typeCell = getCellType(vA, vB, vD);
                C.type = (byte) typeCell;
                cellM.setValue("method3", C.id, typeCell);
                
                nuc.setType(C.type);
                nuc.setName("Nuc" + C.id);
                nuc.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
                String content = nbA +","+nbB +","+nbD +","+vA +","+vB +","+vD +","+typeCell +"\n";
                out.write(content);
            }   
            
            out.close();
        }   
        catch (IOException ex) {
//            Logger.getLogger(RegionSLIC_V2_.class.getName()).log(Level.SEVERE, null, ex);
            Logger.getLogger("Issue here").log(Level.SEVERE, null, ex);
        }
        mapR.show("Mapping RC");
        cellM.show("Method3");
        ImageHandler unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);

        for (Cell C : popCells) {
            if(C.type==UNLABELED){
               Object3D nucleus = C.nucleus;
               nucleus.draw(unLabelledCell, 255);
            }     
        }
        unLabelledCell.show("Unlabelled");
        
    }
    
    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageByte("TYPE", label.sizeX, label.sizeY, label.sizeZ);

        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.type);
            } else {
                C.region.draw(draw, C.type);
            }
        }
        return draw;
    }
    public int getCellType(int vA, int vB, int vD)
    {
        if((vA < 100) && (vB < 100) && (vD < 100))
        {
            return UNLABELED;
        }
        else if(vA > vB && vA > vD)
        {
            if(vD >= vA * 0.8)
            {
                return DELTA;
            }
            else if(vB >= vA*0.95)
            {
                return COLOCAB;
            }
            else if(vD >= vA*0.7){
                return COLOCAD;
            }
            else{
                return ALPHA;
            }
        }
        else if(vB > vA && vB > vD)
        {
            if(vD >= vB * 0.45)
            {
                return DELTA;
            }
            else if(vA >= vB*0.8)
            {
                return ALPHA;
            }
            else if(vA >= vB*0.6)
            {
                return COLOCAB;
            }
            else if(vD >= vB*0.35)
            {
                return COLOCBD;
            }
            else{
                return BETA;
            }
        }
        else if(vD > vA && vD > vB)
        {
            return DELTA;
        }
        else if(vA==vB && vA>vD)
        { 
            return COLOCAB;
        }
        else if((vA==vD && vA>vB))
        {
            return COLOCAD;
        }
        else if(vD==vB && vD>vA)
        {
            return COLOCBD;
        }
        else if((vD==vB && vA==vB && vA>0))
        {
            return COLOCABD;
        }
        else
        {
            return UNLABELED;
        }  
        
    }      
    
}
