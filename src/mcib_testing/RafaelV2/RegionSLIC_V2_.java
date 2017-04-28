/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.ResultsTable;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */
public class RegionSLIC_V2_ implements ij.plugin.PlugIn
{
    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;
    private final int COLOC = 5;
    ImageInt imgSeg, label, imgWat;
    int thresh = 0;
    String dir;
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;
    ImageHandler dapi, alpha, beta, delta;
    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
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
      

        IJ.log("Initialization Region Analysis ...");
        /*IJ.log("Reading data from " + dir);
        popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        popRegions.loadObjects(dir + "Regions.zip");
        popNuclei.loadObjects(dir + "Nuclei.zip");
        initCells2();
        */
        initCells(imgSeg, imgWat);
        //ArrayUtil idxRe = getTypeRegion();
        //IJ.log("Size of array is: " + idxRe.getSize());
        
//        ArrayUtil idxCR = mappingCR();
//        IJ.log("mapping region to cell region, size: "+idxCR.getSize());
        //mappingCRV2();
        mappingRCV2();
        //mappingRC();
        //test();
        //extractRegion();
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
    private void initCells2() {

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer, Cell>(popNuclei.getNbObjects());

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
    
    /**
     * For each of region, find the type of this region
     * One of 5 type : alpha, beta, delta, dapi, unlabelled
     * in case of nb of alpha=beta || alpha=delta || beta=delta not verify
     * @return 
     */
    public ArrayUtil getTypeRegion()
    {
        Objects3DPopulation pop = new Objects3DPopulation(label);
        int nbRegion = pop.getNbObjects();
        ArrayUtil idx = new ArrayUtil(nbRegion);
        //int index = 0;
        ResultsTable rt = new ResultsTable();
        for (int i = 0; i < nbRegion; i++) 
        {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            int value = O.getValue();
            int valType = 0;
            double valC1 = O.getPixMedianValue(delta);
            double valC2 = O.getPixMedianValue(beta);
            double valC3 = O.getPixMedianValue(alpha);
            double valC4 = O.getPixMedianValue(dapi);
            
//            double valC1 = O.getPixMeanValue(delta);
//            double valC2 = O.getPixMeanValue(beta);
//            double valC3 = O.getPixMeanValue(alpha);
//            double valC4 = O.getPixMeanValue(dapi);
            int typeRegion = getTypeR(valC1, valC2, valC3, valC4);
            idx.putValue(value, typeRegion);
            rt.incrementCounter();
            rt.setValue("value", i+1, value);
            rt.setValue("m delta", i+1, valC1);
            rt.setValue("m beta", i+1, valC2);
            rt.setValue("m alpha", i+1, valC3);
            rt.setValue("m dapi", i+1, valC4);
            rt.setValue("type", i+1, typeRegion);
            //rt.setValue("vol pixel", i+1, O.getVolumePixels());
            rt.setValue("vol unit", i+1, O.getVolumeUnit());
            rt.setValue("diameter", i+1, O.getDistCenterMean());
            rt.setValue("elongation", i+1, O.getMainElongation());
        }
        rt.show("Region Attribute");
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
            //IJ.log("test dapi: c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4);
        }
        else if ((valC4 > thresh) && (valC4 == valC2) && (valC4 > valC3) && (valC4 > valC1)) 
        {
            valType = 42;
            //IJ.log("test dapi: c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4);
        }
        else if ((valC4 > thresh) && (valC4 == valC3) && (valC4 > valC2) && (valC4 > valC1)) 
        {
            valType = 43;
            //IJ.log("test: dapi c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4);
        }
        else if ((valC1 > thresh) && (valC1 == valC2) && (valC1 > valC3) && (valC1 > valC4)) 
        {
            valType = 12;
            //IJ.log("-----------------------------------");
            //IJ.log("test: coloc c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4);
        }
        else if ((valC1 > thresh) && (valC1 == valC3) && (valC1 > valC2) && (valC1 > valC4)) 
        {
            valType = 13;
            //IJ.log("-----------------------------------");
            //IJ.log("test: coloc c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4);
        }
        else if ((valC2 > thresh) && (valC2 == valC3) && (valC2 > valC1) && (valC2 > valC4)) 
        {
            valType = 23;
            //IJ.log("-----------------------------------");
            //IJ.log("test: coloc c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4);
        }

        else
        {
            valType = 0;
            //IJ.log("-----------------------------------");
            //IJ.log("test background: c1: "+valC1+" c2: "+valC2+" c3: "+valC3+" c4: "+valC4 + " c41: "+valC41);

        }
        
        return valType;
    }   
    
    public ArrayUtil mappingCRV2()
    {
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        int nbRegion = popR.getNbObjects();
        ArrayUtil idxCR = new ArrayUtil(nbRegion);
        ArrayUtil idxCR2 = new ArrayUtil(nbRegion);
        ImageInt labelCR = new ImageByte("labelRC", label.sizeX, label.sizeY, label.sizeZ);
        try {
            FileWriter ryt = new FileWriter("/home/tranhoa/Pictures/map.txt");
            BufferedWriter out=new BufferedWriter(ryt);
            out.write(nbRegion);
            
            
            for (int i = 0; i < nbRegion; i++) 
            {
                if(popR.getObject(i)!=null)
                {
                    Object3DVoxels O = (Object3DVoxels) popR.getObject(i);

                    int value = O.getValue();
                    ArrayUtil ac1 = O.listValues(imgWat); 
                    ac1 = ac1.distinctValues();
                    int count = 0, idx = 0;
                    double maxVol = 0, maxVol2 = 0;
                    for (int k = 0; k < ac1.getSize(); k++) 
                    {
                        count++;
                        int idxCell = ac1.getValueInt(k);
                        //mapR.setValue("idx"+count, i+1, ac1.getValueInt(k));
                        Object3D reCell = popRegions.getObjectByValue(idxCell);
                        if(reCell!=null)
                        {
                            Object3DVoxels intersectR = new Object3DVoxels();
                            intersectR.addVoxelsIntersection(popR.getObject(i), reCell);
                            double vol = 0;
                            vol = intersectR.getVolumeUnit();
                            if(maxVol < vol){ maxVol = vol; idx = idxCell;}
                            //mapR.setValue("vol"+count, i+1, vol);
                        }

                    }   
                    idxCR.putValue(value, idx);
                    String content = value +" "+idx+"\n";
                    out.write(content);
                    O.draw(labelCR, idx);
                }   
            } 
            out.close();
            labelCR.show("map-res");
            
        } catch (IOException ex) {
            Logger.getLogger(RegionSLIC_V2_.class.getName()).log(Level.SEVERE, null, ex);
        }
        return idxCR;
        
    } 
    public ArrayUtil mappingCR()
    {
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        //ResultsTable mapR = new ResultsTable();
        int nbRegion = popR.getNbObjects();
        ArrayUtil idxCR = new ArrayUtil(nbRegion);
        //ImageInt labelCR = label.duplicate();
        ImageInt labelCR = new ImageByte("labelRC", label.sizeX, label.sizeY, label.sizeZ);
        for (int i = 0; i < nbRegion; i++) 
        {
            if(popR.getObject(i)!=null)
            {
                Object3DVoxels O = (Object3DVoxels) popR.getObject(i);
         
                int value = O.getValue();
                ArrayUtil ac1 = O.listValues(imgWat); 
                ac1 = ac1.distinctValues();
                //mapR.incrementCounter();
                //mapR.setValue("region", i+1, value);
//                for(int j=1; j<=5; j++){
//                    mapR.setValue("idx"+j, i+1, 0);
//                }
//                for(int j=1; j<=5; j++){
//                    mapR.setValue("vol"+j, i+1, 0);
//                }
                int count = 0, idx = 0;
                double maxVol = 0;
                for (int k = 0; k < ac1.getSize(); k++) 
                {
                    count++;
                    int idxCell = ac1.getValueInt(k);
                    //mapR.setValue("idx"+count, i+1, ac1.getValueInt(k));
                    Object3D reCell = popRegions.getObjectByValue(idxCell);
                    if(reCell!=null)
                    {
                        Object3DVoxels intersectR = new Object3DVoxels();
                        intersectR.addVoxelsIntersection(popR.getObject(i), reCell);
                        double vol = 0;
                        vol = intersectR.getVolumeUnit();
                        if(maxVol < vol){ maxVol = vol; idx = idxCell;}
                        //mapR.setValue("vol"+count, i+1, vol);
                    }
                    
                }    
                //mapR.setValue("max vol", i+1, maxVol);
                //mapR.setValue("idx max", i+1, idx);
                idxCR.putValue(value, idx);
                O.draw(labelCR, idx);
            }    
            
            
        } 
        for (Cell C : popCells)
        {
            Object3D region = C.region;
            region.setLabelImage(imgWat);
            ArrayList<Voxel3D> cont = region.getContours();
            Object3DVoxels vos = new Object3DVoxels();
            vos.addVoxels(cont);
            vos.draw(labelCR, 255);
            region.setLabelImage(null);
        } 
        //mapR.show("Mapping CR");
        labelCR.show("result-mapping");
        return idxCR;
    }
    
    public void mappingRCV2()
    {
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        ArrayUtil idxRe = getTypeRegion();
        IJ.log("Get type region, size of array is: " + idxRe.getSize());
        ArrayUtil idxCR = mappingCRV2();
        IJ.log("mapping region to cell region, size: "+idxCR.getSize());
        int cA=0, cB=0, cD=0, cColoc=0, cUnlabelled=0;
        
        ResultsTable mapR = new ResultsTable();
        //for (Cell C : popCells) 
        
        try {
            FileWriter ryt = new FileWriter("/home/tranhoa/Pictures/mapping2.arff");
            BufferedWriter out=new BufferedWriter(ryt);
            out.write("@relation cell\n");
            out.write("@attribute nbA numeric\n");
            out.write("@attribute nbB numeric\n");
            out.write("@attribute nbD numeric\n");
            out.write("@attribute vA numeric\n");
            out.write("@attribute vB numeric\n");
            out.write("@attribute vD numeric\n");
            out.write("@attribute class {0,1,2,3}\n");
            out.write("@data\n");
            
                for (int k = 0; k < popCells.size(); k++) 
                {
                    Cell C = popCells.get(k);
                    Object3D nuc = C.nucleus;
                    if(nuc==null) continue;
                    Object3D region = C.region;
                    //for(int i = 0; i < idxCR.getSize(); i++)
                    ArrayUtil ac1 = region.listValues(label); 
                    ac1 = ac1.distinctValues();
                    int reVal = region.getValue();
                    mapR.incrementCounter();
                    mapR.setValue("cell ", C.id, reVal);

                    int nbA=0, nbB=0, nbD=0, nbNuc=0, nbUnlabelled=0, nbColoc=0;
                    int vA=0, vB=0, vD=0, vNuc=0;
                    for(int j=1; j<=15; j++){
                        mapR.setValue("idx"+j, C.id, 0);
                    }
                    int count = 0;
                    for (int i = 0; i < ac1.getSize(); i++) 
                    {
                        int indexRegion = ac1.getValueInt(i);
                        if(idxCR.getValue(indexRegion)==reVal)
                        {    
                            int type = idxRe.getValueInt(indexRegion);
                            //IJ.log("Index region is: " + indexRegion+" type is: "+ type);
                            Object3DVoxels intersectR = new Object3DVoxels();
                            double vol = 0;
                            if(type!=0)
                            {    
                                if(popR.getObjectByValue(indexRegion)!=null)
                                {
                                    Object3D re = popR.getObjectByValue(indexRegion);
                                    //double volObj = re.getVolumeUnit();
                                    intersectR.addVoxelsIntersection(re, region);
                                    vol = intersectR.getVolumeUnit();

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
                                        //if(type==0){nbUnlabelled++;}
                                }
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
                    int type = 0;
                    if(vA>vB && vA>vD)
                    {
                        type = 3;
                    }
                    else if(vB>vA && vB>vD)
                    {
                        type = 2;
                    }
                    else if(vD>vA && vD>vB)
                    {
                        type = 1;
                    }
                    else if((vA==vB && vA>vD) || (vA==vD && vA>vB) || (vA==vB && vA>vD))
                    { 
                        type=4;
                    }
                    else{
                        type=0;
                    }
                    String content = nbA +","+nbB +","+nbD +","+vA +","+vB +","+vD +","+type +"\n";
                    out.write(content);
                }   
            
            out.close();
        }   
        catch (IOException ex) {
            Logger.getLogger(RegionSLIC_V2_.class.getName()).log(Level.SEVERE, null, ex);
        }
        mapR.show("Mapping RC");
    
    }
    public void mappingRC()
    {
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        ArrayUtil idxRe = getTypeRegion();
        IJ.log("Size of array is: " + idxRe.getSize());
        int cA=0, cB=0, cD=0, cColoc=0, cUnlabelled=0;
        
        ResultsTable mapR = new ResultsTable();
        //for (Cell C : popCells) 
        for (int k = 0; k < popCells.size()/20; k++) 
        {
            Cell C = popCells.get(k);
            Object3D nuc = C.nucleus;
            if(nuc==null) continue;
            Object3D region = C.region;
            ArrayUtil ac1 = region.listValues(label); 
            ac1 = ac1.distinctValues();
            
            mapR.incrementCounter();
            mapR.setValue("cell ", C.id, region.getValue());
            
            int nbA=0, nbB=0, nbD=0, nbNuc=0, nbUnlabelled=0, nbColoc=0;
            int vA=0, vB=0, vD=0, vNuc=0;
            for(int j=1; j<=10; j++){
                mapR.setValue("idx"+j, C.id, 0);
            }
            int count = 0;
            for (int i = 0; i < ac1.getSize(); i++) 
            {
                int indexRegion = ac1.getValueInt(i);
                int type = idxRe.getValueInt(indexRegion);
                //IJ.log("Index region is: " + indexRegion+" type is: "+ type);
                Object3DVoxels intersectR = new Object3DVoxels();
                double vol = 0;
                if(type!=0)
                {    
                    if(popR.getObjectByValue(indexRegion)!=null)
                    {
                        Object3D re = popR.getObjectByValue(indexRegion);
                        double volObj = re.getVolumeUnit();
                        intersectR.addVoxelsIntersection(re, region);
                        vol = intersectR.getVolumeUnit();
                        if(vol >= volObj*0.33)
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
                            //if(type==0){nbUnlabelled++;}
                        }
                        if((vol > volObj*0.1 && vol < volObj*0.33) && (type==4 || type==41 || type==42 || type==43)){
                            nbNuc++; vNuc +=vol;
                        }

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
        }   
        mapR.show("Mapping RC");
    
    }
    public void extractRegion()
    {
        
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        ImageInt newLabel = label.duplicate();
        //ImageHandler imgReg = new ImageShort("cell-region", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        ImageHandler imgUnlabelled = new ImageShort("cell-unlabelled", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        ArrayUtil idxRe = getTypeRegion();
        IJ.log("Size of array is: " + idxRe.getSize());
        int cA=0, cB=0, cD=0, cColoc=0, cUnlabelled=0;
        for (Cell C : popCells) 
        {
            Object3D nuc = C.nucleus;
            if(nuc==null) continue;
            Object3D region = C.region;
            ArrayUtil ac1 = region.listValues(label); 
            ac1 = ac1.distinctValues();
            int nbA=0, nbB=0, nbD=0, nbNuc=0, nbUnlabelled=0, nbColoc=0;
            int vA=0, vB=0, vD=0;
            for (int i = 0; i < ac1.getSize(); i++) 
            {
                int indexRegion = ac1.getValueInt(i);
                int type = idxRe.getValueInt(indexRegion);
                //IJ.log("Index region is: " + indexRegion+" type is: "+ type);
                Object3DVoxels obj = null;
                if(popR.getObjectByValue(indexRegion)!=null)
                {
                    obj = (Object3DVoxels) popR.getObjectByValue(indexRegion);
                    
                    if(type < 5)
                    {
                        obj.setValue(type*400);
                    }
                    else{
                        obj.setValue(2000);
                    }
                    
//                    if(type==4) {obj.setValue(1426);}
//                    if(type==3) {obj.setValue(166);}
//                    if(type==2) {obj.setValue(1619);}
//                    if(type==3) {obj.setValue(74);}
//                    if(type==0) {obj.setValue(0);}
                    obj.draw(newLabel);
//                    obj.setLabelImage(label);
//                    ArrayList<Voxel3D> contourRegion = obj.getContours();
//                    Object3DVoxels vox = new Object3DVoxels();
//                    vox.addVoxels(contourRegion);
//                    vox.draw(newLabel, 2020);
//                    obj.setLabelImage(null);
                }
                if(type==1){nbD++; vD += obj.getVolumePixels();}
                if(type==2){nbB++; vB += obj.getVolumePixels();}
                if(type==3){nbA++; vA += obj.getVolumePixels();}
                if(type==4){nbNuc++;}
                //if(type==5){nbColoc++;}
                if(type==41){nbD++; nbNuc++; vD += obj.getVolumePixels();}
                if(type==42){nbB++; nbNuc++; vB += obj.getVolumePixels();}
                if(type==43){nbA++; nbNuc++; vA += obj.getVolumePixels();}
                if(type==12){nbD++; nbB++; vD += obj.getVolumePixels();}
                if(type==13){nbD++; nbA++; vD += obj.getVolumePixels();}
                if(type==23){nbB++; nbA++; vA += obj.getVolumePixels(); vB += obj.getVolumePixels();}
                if(type==0){nbUnlabelled++;}
                //IJ.log(" Value of region is: " + indexRegion + " & type of region is: " + type);
            }
//            ArrayList<Voxel3D> cont = region.getContours();
//            Object3DVoxels vos = new Object3DVoxels();
//            vos.addVoxels(cont);
//            vos.draw(newLabel, 2040);
            
            if(nbD>0){
                IJ.log("*********************test delta cell ********** ");
                IJ.log("Nb of al: " + nbA + "  beta: " + nbB + " del: " + nbD + " nuc: " + nbNuc 
                        + " unlabelled: " + nbUnlabelled + " nbColoc: " + nbColoc
                        + " vA: "+vA+" vB: "+vB+" vD: "+vD);
            }
            int typeCell = getCellType(nbA, nbB, nbD, nbNuc, nbColoc, vA, vB, vD);
            
            if(typeCell==UNLABELED){
                IJ.log("****************cell is unlabelled***********");
                IJ.log("Nb of al: " + nbA + "  beta: " + nbB + " del: " + nbD + " nuc: " + nbNuc 
                        + " unlabelled: " + nbUnlabelled + " nbColoc: " + nbColoc);
                C.type = UNLABELED;
                cUnlabelled++; 
            }
            if(typeCell==ALPHA)
            {
                //IJ.log("cell is alpha");
                C.type = ALPHA;
                cA++;
            }
            if(typeCell==BETA)
            {
                //IJ.log("cell is beta"); 
                C.type = BETA;
                cB++;
            }
            if(typeCell==DELTA)
            {
                //IJ.log("cell is delta"); 
                C.type = DELTA;
                cD++;
            }
            if(typeCell==COLOC)
            {
                //IJ.log("cell is colocalise"); 
                C.type = COLOC;
                IJ.log("*********************coloc********** ");
                IJ.log("Nb of al: " + nbA + "  beta: " + nbB + " del: " + nbD + " nuc: " + nbNuc 
                        + " unlabelled: " + nbUnlabelled + " nbColoc: " + nbColoc);
                cColoc++;
            }
            //IJ.log("Nb of al: " + nbA + "  beta: " + nbB + " del: " + nbD + " nuc: " + nbNuc + " unlabelled: " + nbUnlabelled);
        } 
        int countCell = 0;
        //IJ.log("Number of region: " + popCells.size());
        for (Cell C : popCells)
        {
            if(C.nucleus!=null)
            {
                countCell++;
            }
            if(C.type==UNLABELED){
                Object3D nuc1 = C.nucleus;
                nuc1.draw(imgUnlabelled, 2000);
                double valMean = nuc1.getPixMeanValue(dapi);
                IJ.log("*********************unlabelled********** ");
                IJ.log("Value mean is: " + valMean);
            }
            Object3D nuc = C.nucleus;
            if (nuc!=null) {
                C.nucleus.draw(newLabel,1600);
            }
            Object3D region = C.region;
            region.setLabelImage(imgWat);
            ArrayList<Voxel3D> cont = region.getContours();
            Object3DVoxels vos = new Object3DVoxels();
            vos.addVoxels(cont);
            vos.draw(newLabel, 2047);
            region.setLabelImage(null);
        } 
        IJ.log("Nb of cell exact is: "+countCell+" nb of region not cell: " + (popCells.size()- countCell));
        newLabel.show("New Label");
        imgUnlabelled.show();
        IJ.log("Result of analyse cell type is: al: "+cA+" beta: "+cB+" del: "+cD+" coloc:"+cColoc+" unlabelled: "+cUnlabelled);
    }        
    
    public int getCellType(int nbA, int nbB, int nbD, int nbNuc, int nbColoc, int vA, int vB, int vD)
    {
          
        if(nbA > nbB && nbA > nbD)
        {
            if(nbD>0 && vD > vA * 0.4)
            {
                return DELTA;
            }
            else{
                return ALPHA;
            }
            
        }   
        else if(nbB > nbA && nbB > nbD)
        {
            if(nbD > 0 && vD > vB * 0.4)
            {
                return DELTA;
            }
            else{
                return BETA;
            }
        }
        else if(nbD > nbA && nbD > nbB)
        {
            return DELTA;
        }
        else if((nbA==nbB && nbA>nbD))
        {
            return COLOC;
        }
        else if((nbA==nbD && nbA>nbB) || (nbB==nbD && nbB>nbA))
        {
            return DELTA;
        } 
        else{
            return UNLABELED;
        }
        
    }        
    private class Cell {

        int id;
        Object3D nucleus;
        Object3D region;
        byte type;
        int layer;

        ArrayList<Cell> nei1 = null;
        ArrayList<Cell> nei2 = null;

        @Override
        public String toString() {
            return "(" + nucleus.getValue() + ", " + region.getValue() + ", " + type + ")";
        }

        private int[] computeNeiType(ArrayList<Cell> nei) {
            int[] res = new int[3];
            for (Cell C : nei) {
                if (C.type > 0) {
                    res[C.type - 1]++;
                }
            }
            return res;
        }

        public int[] computeNei1Type() {
            return computeNeiType(nei1);
        }

        public int[] computeNei2Type() {
            return computeNeiType(nei2);
        }

        public int[] computeNeiRangeType(double dist) {
            return computeNeiType(getCellsRange(dist));
        }

        private boolean hasContact(int type, ArrayList<Cell> nei) {
            for (Cell N : nei) {
                if (N.type == type) {
                    return true;
                }
            }
            return false;
        }

        public boolean hasContact1(int type) {
            return hasContact(type, nei1);
        }

        public boolean hasContact2(int type) {
            return hasContact(type, nei2);
        }

        private ArrayList<Cell> getCellsRange(double dist) {
            ArrayList<Cell> res = new ArrayList<Cell>();
            for (Cell C : popCells) {
                double d = this.nucleus.distCenterUnit(C.nucleus);
                if ((d > 0) && (d < dist)) {
                    res.add(C);
                }
            }
            return res;
        }
    }

}
