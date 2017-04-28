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
import java.util.ArrayList;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */
public class Test_ implements ij.plugin.PlugIn
{
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
        IJ.log("Testing......");
        initCells(imgSeg, imgWat);
        fusionTable();
        //getTypeRegionV2();
        IJ.log("Finished");
        
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
    public ArrayUtil getTypeRegion()
    {
        ImageInt labelR = new ImageByte("label", label.sizeX, label.sizeY, label.sizeZ);
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
            
            int typeRegion = getTypeR(valC1, valC2, valC3, valC4);
            idx.putValue(value, typeRegion);
            O.draw(labelR, typeRegion * 20);
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
        for (Cell C : popCells)
        {
            Object3D region = C.region;
            region.setLabelImage(imgWat);
            ArrayList<Voxel3D> cont = region.getContours();
            Object3DVoxels vos = new Object3DVoxels();
            vos.addVoxels(cont);
            vos.draw(labelR, 255);
            region.setLabelImage(null);
        } 
        //labelR.show("label region");
        //rt.show("Region Attribute");
        return idx;
    }
    
    public void getTypeRegionV2()
    {
        ImageInt labelR = new ImageByte("label", label.sizeX, label.sizeY, label.sizeZ);
        Objects3DPopulation pop = new Objects3DPopulation(label);
        int nbRegion = pop.getNbObjects();
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
            if(typeRegion==4 || typeRegion==41 || typeRegion==42 || typeRegion==43)
            {
                O.draw(labelR, 160);
            }    

        }
        for (Cell C : popCells)
        {
            Object3D nuc = C.nucleus;
            //Object3D region = C.region;
            //region.setLabelImage(imgWat);
            nuc.setLabelImage(imgSeg);
            ArrayList<Voxel3D> cont = nuc.getContours();
            Object3DVoxels vos = new Object3DVoxels();
            vos.addVoxels(cont);
            vos.draw(labelR, 255);
//            ArrayList<Voxel3D> contR = region.getContours();
//            Object3DVoxels vosR = new Object3DVoxels();
//            vosR.addVoxels(contR);
//            vosR.draw(labelR, 190);
//            region.setLabelImage(null);
            nuc.setLabelImage(null);
        } 
        labelR.show("label region");
        //rt.show("Region Attribute");
        //return idx;
    }
    
    public void fusionTable()
    {
        ResultsTable total = new ResultsTable();
        ResultsTable r1 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/Islet1-3B-small/SEG/v7/Method1.csv");
        ResultsTable r2 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/Islet1-3B-small/SEG/v6/Method2.csv");
        ResultsTable r3 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/Islet1-3B-small/SEG/v6/Method2Copy.csv");
        for(int i=0; i<r1.size(); i++)
        {
            double val1 = r1.getValue("cell", i);
            double val2 = r2.getValue("cell", i);
            double val3 = r3.getValue("cell", i);
            IJ.log("val1: "+val1 +" val2: "+val2+" val3:"+val3);
            if(val1==val2 && val1==val3)
            {
                double type1 = r1.getValue("method1", i);
                double type2 = r2.getValue("method2", i);
                double type3 = r3.getValue("method3", i);
                IJ.log("t1: "+type1 +" t2: "+type2+" t3:"+type3);
                total.incrementCounter();
                total.setValue("cell", i, val1);
                total.setValue("method1", i, type1);
                total.setValue("method2", i, type2);
                total.setValue("method3", i, type3);
            }
        }
        total.show("Result");
    }
    /**
     * Found the unlabelled cell.
     */
    public void readTable()
    {
        ImageInt labelCR = new ImageByte("labelRC", label.sizeX, label.sizeY, label.sizeZ);
        Objects3DPopulation popR = new Objects3DPopulation((ImageInt)label);
        ArrayUtil idxRe = getTypeRegion();
        IJ.log("Get type region, size of array is: " + idxRe.getSize());
        ResultsTable rt = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/Islet1-3B-small/SEG/v5/MappingRC.csv");
        for(int i=0; i<rt.size(); i++)
        {
            double al = rt.getValue("vA", i);
            double bt = rt.getValue("vB", i);
            double dl = rt.getValue("vD", i);
            double nuc = rt.getValue("vNuc", i);
            if(al==0 && bt==0 && dl==0)
            {
                double idxCell = rt.getValue("cell", i);
                Object3D region = popRegions.getObjectByValue((int)idxCell);
                if(region!=null)
                {
                    //IJ.log("holala " + idxCell);
                    //region.draw(labelCR, 255);
                    ArrayUtil ac1 = region.listValues(label); 
                    ac1 = ac1.distinctValues();
                    int nbA=0, nbB=0, nbD=0, nbNuc=0, nbUnlabelled=0, nbColoc=0;
                    int vA=0, vB=0, vD=0, vNuc=0;
                    for (int k = 0; k < ac1.getSize(); k++) 
                    {
                        int indexRegion = ac1.getValueInt(k);
                        Object3DVoxels intersectR = new Object3DVoxels();
                        //double vol = 0;
                        if(popR.getObjectByValue(indexRegion)!=null)
                        {
                            Object3D re = popR.getObjectByValue(indexRegion);
                            int type = idxRe.getValueInt(indexRegion);
                            intersectR.addVoxelsIntersection(re, region);
                            double volR = intersectR.getVolumeUnit();
                            if(type==1){nbD++; vD += volR;}
                            if(type==2){nbB++; vB += volR;}
                            if(type==3){nbA++; vA += volR;}
                            if(type==4){nbNuc++; vNuc +=volR;}
                            if(type==41){nbD++; nbNuc++; vD += volR; vNuc +=volR;}
                            if(type==42){nbB++; nbNuc++; vB += volR; vNuc +=volR;}
                            if(type==43){nbA++; nbNuc++; vA += volR; vNuc +=volR;}
                            if(type==12){nbD++; nbB++; vD += volR; vB += volR;}
                            if(type==13){nbD++; nbA++; vD += volR; vA += volR;}
                            if(type==23){nbB++; nbA++; vA += volR; vB += volR;}

                            //if(type==4 || type==41 || type==42 || type==43){
                                intersectR.draw(labelCR, type);
                            //}
                            
                        }
                        
                    }  
                    
                    Cell c1 = region2Cell.get(region.getValue());
                    if(c1.nucleus!=null)
                    {
                        IJ.log("holala " + idxCell);
                        IJ.log("volume: vA: "+vA+" vB: "+vB+" vD: "+vD+" vNuc:"+vNuc);
                        IJ.log("data vNuc: "+nuc);
                        //region.draw(labelCR, 255);
                        region.setLabelImage(imgWat);
                        ArrayList<Voxel3D> cont = region.getContours();
                        Object3DVoxels vos = new Object3DVoxels();
                        vos.addVoxels(cont);
                        vos.draw(labelCR, 255);
                        region.setLabelImage(null);
                    }
//                    if(vA>0 || vB>0 || vD>0)
//                    {
//                        
//                        Cell c1 = region2Cell.get(region.getValue());
//                        if(c1.nucleus!=null)
//                        {
//                            IJ.log("holala " + idxCell);
//                            IJ.log("volume: vA: "+vA+" vB: "+vB+" vD: "+vD+" vNuc:"+vNuc);
//                            //region.draw(labelCR, 255);
//                            region.setLabelImage(imgWat);
//                            ArrayList<Voxel3D> cont = region.getContours();
//                            Object3DVoxels vos = new Object3DVoxels();
//                            vos.addVoxels(cont);
//                            vos.draw(labelCR, 255);
//                            region.setLabelImage(null);
//                        }
//                    }
//                    region.setLabelImage(imgWat);
//                    ArrayList<Voxel3D> cont = region.getContours();
//                    Object3DVoxels vos = new Object3DVoxels();
//                    vos.addVoxels(cont);
//                    vos.draw(labelCR, 255);
//                    region.setLabelImage(null);
                }
                
                 
                
                
            }
        }
        labelCR.show();
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
