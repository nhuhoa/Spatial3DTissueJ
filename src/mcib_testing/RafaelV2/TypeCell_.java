/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

import mcib_testing.Rafael.*;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.ResultsTable;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author thomasb
 */
public class TypeCell_ implements ij.plugin.PlugIn {

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
    
    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;

    @Override
    public void run(String arg) {
        imgLabel = ImageInt.wrap(WindowManager.getImage("labelF.tif"));
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        
        WindowManager.setTempCurrentImage(plus);
        String dir = IJ.getDirectory("image");
        ImageInt imgWat = ImageInt.wrap(plus);
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
  
        IJ.log("Initialization Tissue Analysis ...");
        initCells(imgSeg, imgWat);
        IJ.log("Type ...");
        
        this.computeTypeCellsNucleus(imgLabel);
        popNuclei.saveObjects(dir + "/v7/Nuclei.zip");
        popRegions.saveObjects(dir + "/v7/Regions.zip");
        IJ.log("save the result");
        
        // draw 
        drawCellTypes(false).show("TYPE");
        IJ.log("Finished");

    }

    private void computeResults() {
        ResultsTable rtALPHA = new ResultsTable();
        ResultsTable rtBETA = new ResultsTable();
        ResultsTable rtDELTA = new ResultsTable();
        int rA = -1, rB = -1, rD = -1;

        // DELTA population
        Objects3DPopulation popD = createPopulationType(DELTA);
        popD.createKDTreeCenters();

        for (Cell C : popCells) {
            Object3D obj = C.nucleus;
            ResultsTable rt = null;
            int c = 0;
            if (C.type == ALPHA) {
                rt = rtALPHA;
                rA++;
                c = rA;
            } else if (C.type == BETA) {
                rt = rtBETA;
                rB++;
                c = rB;
            } else if (C.type == DELTA) {
                rt = rtDELTA;
                rD++;
                c = rD;
            }
            if (rt != null) {
                rt.incrementCounter();
                rt.setValue("type", c, C.type);
                rt.setValue("val", c, C.id);
                rt.setValue("vol", c, obj.getVolumeUnit());
                rt.setValue("diam", c, obj.getDistCenterMean());
                rt.setValue("elon", c, obj.getMainElongation());
                int[] nei = C.computeNei1Type();
                rt.setValue("Nei1_A", c, nei[0]);
                rt.setValue("Nei1_B", c, nei[1]);
                rt.setValue("Nei1_D", c, nei[2]);
                nei = C.computeNei2Type();
                rt.setValue("Nei2_A", c, nei[0]);
                rt.setValue("Nei2_B", c, nei[1]);
                rt.setValue("Nei2_D", c, nei[2]);
                // nb in range distance
                double dist = 20;
                nei = C.computeNeiRangeType(dist);
                rt.setValue("NeiR_A", c, nei[0]);
                rt.setValue("NeiR_B", c, nei[1]);
                rt.setValue("NeiR_D", c, nei[2]);

                // closest D 
                Object3D cloD = popD.closestCenter(obj, 0);
                if (cloD == null) {
                    rt.setValue("cloD", c, -1);
                    rt.setValue("cloDi", c, -1);
                } else {
                    rt.setValue("cloD", c, obj.distCenterUnit(cloD));
                    rt.setValue("cloDi", c, nucleus2Cell.get(cloD.getValue()).id);
                }
            }
        }
        rtALPHA.show("Alpha");
        rtBETA.show("Beta");
        rtDELTA.show("Delta");
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
                    if (img.getPixel(x, y, z) == BorderValue) {
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

    private Objects3DPopulation createPopulationType(int type) {
        Objects3DPopulation pop = new Objects3DPopulation();
        pop.setCalibration(popNuclei.getCalibration());

        for (Cell C : popCells) {
            if (C.type == type) {
                pop.addObject(C.nucleus);
            }
        }
        return pop;
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

    public int getCellType(int vA, int vB, int vD)
    {
        if(vA > vB && vA > vD)
        {
            if(vD >= vA * 0.8)
            {
                return DELTA;
            }
            else if(vB >= vA*0.95)
            {
                return COLOCAB;
            }
            else if(vD >= vA*0.7)
            {
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
        else
        {
            return UNLABELED;
        }  
        
    }      
    private void computeTypeCellsNucleus(ImageInt label) 
    {
        int cA=0, cB=0, cD=0, cUnlabelled=0, cColocAB=0, cColocAD=0, cColocBD=0;
        ResultsTable cellM = new ResultsTable();
        try {
            FileWriter ryt = new FileWriter("/home/tranhoa/Raphael/data_test/testV2/Islet1-3B-small/SEG/v7/method1.arff");
            BufferedWriter out=new BufferedWriter(ryt);
            out.write("@relation cell\n");
            out.write("@attribute nbA numeric\n");
            out.write("@attribute nbB numeric\n");
            out.write("@attribute nbD numeric\n");
            out.write("@attribute class {0,1,2,3,23,13,12}\n");
            out.write("@data\n");
            for (Cell C : popCells) 
            {
                Object3D nucleus = C.nucleus;
                Object3D region = C.region;
                int reVal = region.getValue();
                cellM.incrementCounter();
                cellM.setValue("cell", C.id, reVal);
                if(nucleus==null) continue;
                ArrayUtil list = nucleus.listValues(label);
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<list.getSize();i++)
                {
                    if(list.getValue(i)==ALPHA)nbA++;
                    if(list.getValue(i)==BETA)nbB++;
                    if(list.getValue(i)==DELTA)nbD++;                
                }
                C.type = UNLABELED;
                int typeCell = getCellType(nbA, nbB, nbD);
                C.type = (byte)typeCell;
                cellM.setValue("method1", C.id, typeCell);
                if(typeCell==UNLABELED){cUnlabelled++;}
                if(typeCell==ALPHA){cA++;}
                if(typeCell==BETA){cB++;}
                if(typeCell==DELTA){cD++;}
                if(typeCell==COLOCAB){cColocAB++;}
                if(typeCell==COLOCAD){cColocAD++;}
                if(typeCell==COLOCBD){cColocBD++;}
                String content = nbA +","+nbB +","+nbD +","+typeCell +"\n";
                out.write(content);
                nucleus.setType(C.type);
                nucleus.setName("Nuc" + C.id);
                nucleus.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
            }
            out.close();
        }
        catch (IOException ex) 
        {
            Logger.getLogger(TypeCell_.class.getName()).log(Level.SEVERE, null, ex);
        }
        cellM.show("Method1");
        IJ.log("Nb of each cells type: A: "+cA+" B: "+cB+" D:"+cD+" Unlabelled: "+cUnlabelled
                +" AB: "+cColocAB+" AD: "+cColocAD+" BD: "+cColocBD);
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);

        for (Cell C : popCells) {
            if(C.type==UNLABELED){
               Object3D nucleus = C.nucleus;
               nucleus.draw(unLabelledCell, 255);
            }     
        }
        unLabelledCell.show("Unlabelled");
    }

    private void computeLayer(int type) {
        // init
        for (Cell C : popCells) {
            C.layer = -1;
        }
        // layer --> ARRAY LIST
        ArrayList<Cell> layer0 = new ArrayList<Cell>();
        ArrayList<Cell> layer1;

        for (Cell C : popCells) {
            if (C.type == type) {
                C.layer = 0;
                layer0.add(C);
            }
        }

        boolean stop = false;
        while (!stop) {
            stop = true;
            layer1 = new ArrayList<Cell>();
            for (Cell C : layer0) {
                for (Cell C1 : C.nei1) {
                    if ((C1.type > 0) && (C1.layer == -1)) {
                        C1.layer = C.layer + 1;
                        layer1.add(C1);
                        stop = false;
                    }
                }
            }
            layer0 = layer1;
        }
    }

    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);

        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.type);
            } else {
                C.region.draw(draw, C.type);
            }
        }

        return draw;
    }

    private ImageHandler drawCellLayer(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);

        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.layer + 1);
            } else {
                C.region.draw(draw, C.layer + 1);
            }
        }

        return draw;
    }

    private ImageHandler drawCellNuclei() {
        ImageHandler draw = new ImageShort("NUC", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);

        for (Cell C : popCells) {
            C.nucleus.draw(draw, C.id);
        }

        return draw;
    }

    private ImageHandler drawCellRegions() {
        ImageHandler draw = new ImageShort("REG", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);

        for (Cell C : popCells) {
            C.region.draw(draw, C.id);
        }

        return draw;
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
