/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Rafael;


import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.ResultsTable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DPoint;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */

public class Rafael_CorrectNuclei_ implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;

    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    
    ImageInt imgLabel;
    String dir = null;
    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;

    @Override
    public void run(String arg) 
    {
        //ImagePlus plus = WindowManager.getImage("dapi-seg-wat.tif");
        
        //ImageInt imgWat = ImageInt.wrap(plus);
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg-SLIC-25000.tif"));
        //ImageInt imgNucleiOrigin = ImageInt.wrap(WindowManager.getImage("C4-dapi.tif"));
        
        //ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg-SLIC-25000.tif"));
        dir = IJ.getDirectory("image");
        imgLabel = ImageInt.wrap(WindowManager.getImage("labelF.tif"));
        
        //IJ.log("Initialization Tissue Analysis ...");
        //initCells(imgSeg, imgWat);
        //IJ.log("Type ...");
        
        IJ.log("Correct the result of segmentation");
        this.correctNucleusSLIC(imgSeg, imgLabel);
        //this.computeTypeCellsNucleusSLIC(imgLabel);

        
        //ImageHandler nucleiLabelDraw = null;
        //nucleiLabelDraw.setTitle("result-extract");
        //nucleiLabelDraw = extractRegion_DAPI_SLIC(imgNucleiOrigin, imgLabel);
        //nucleiLabelDraw.show("SLIC");
        //nucleiLabelDraw.save(dir+"extraction-SLIC.tif");
        
        
//        // draw 
        //drawCellTypes(false).show("TYPE");
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

    private ImageHandler extractRegion_DAPI_SLIC(ImageInt imgNucleiOrigin, ImageInt label)
    {
        IJ.log("Start, extract the region that appear in img nuclei and label at the same time");
        ImageHandler nucleiLabelDraw = null;
        nucleiLabelDraw = new ImageShort("Nuclei DAPI", imgNucleiOrigin.sizeX, imgNucleiOrigin.sizeY, imgNucleiOrigin.sizeZ);
        /*if((imgNucleiOrigin.sizeX != label.sizeX) && (imgNucleiOrigin.sizeY != label.sizeY) 
           && (imgNucleiOrigin.sizeZ != label.sizeZ))
        { return null; }
        */
        //ImageInt imgNuclei = imgNucleiOrigin.duplicate();
        int count = 0;
        for (int z = 0; z < imgNucleiOrigin.sizeZ; z++) {
            for (int y = 0; y < imgNucleiOrigin.sizeY; y++) {
                for (int x = 0; x < imgNucleiOrigin.sizeX; x++) {
                    int value = imgNucleiOrigin.getPixelInt(x, y, z);
                    if((label.getPixelInt(x,y,z)==4) && value > 0)
                    {
                        Object3D o = new Object3DPoint(value, x, y, z);
                        o.draw(nucleiLabelDraw, value);
                    }
                }
            }
        }
        
        return nucleiLabelDraw;
    }
    private void correctNucleusSLIC(ImageInt nucLabel, ImageInt label)
    {
        IJ.log("Correct the result of segmentation using image labelF, result of SLIC");
        //int[] histype = {0};
        int rm = 0; 
        int ac = 0;
        popNuclei = new Objects3DPopulation(nucLabel);
        IJ.log("Number of cells detected at the first time : " + popNuclei.getNbObjects());
        
        //for(Object3D nuc : popNuclei.getObjectsList())
        for (int i = 0; i < popNuclei.getNbObjects(); i++)
        {
            Object3D nuc = popNuclei.getObject(i);
            int nb_dapi = 0;
            int nbA=0; 
            int nbB=0; 
            int nbD=0;
            int nbNoise = 0;
            int nbColoc = 0; 
            ArrayUtil list = nuc.listValues(label);
            int nbTotal = list.getSize();
            for(int k=0; k<nbTotal; k++)
            {
                if(list.getValue(k)!=0)
                {
                    nbColoc++;
                }
                if(list.getValue(k)==DAPI)
                {
                    nb_dapi++;
                }   
                /*else if(list.getValue(k)==ALPHA)
                {
                    nbA++;
                }
                else if(list.getValue(k)==BETA)
                {
                    nbB++;
                }
                else if(list.getValue(k)==DELTA)
                {
                    nbD++;
                }
                else
                {
                    nbNoise++;
                }*/
            }
            if((nb_dapi > nbTotal * 0.25) && (nbColoc >= nbTotal * 0.75))
            {
                //IJ.log("Percentage of exact segment: " + nb_dapi +" " + nbTotal);
                ac++;
                
            }
            else
            {
                //IJ.log("Error of segmentation, not exact a nucleus");
                rm++;
                popNuclei.removeObject(i);
           
            }
        }  
        IJ.log("Number of error of segmentation is: " + rm);
        ImageHandler nucLabelDraw;
        nucLabelDraw = new ImageShort("DAPI Correct", nucLabel.sizeX, nucLabel.sizeY, nucLabel.sizeZ);

        for (Object3D nuc : popNuclei.getObjectsList()) {
            //C.nucleus.draw(draw, C.id);
            if(nuc==null) continue; 
            nuc.draw(nucLabelDraw, nuc.getValue());
        }
        nucLabelDraw.show("NUC Correct");
        nucLabelDraw.setTitle(nucLabel.getTitle()+"-correct");
        //nucLabelDraw.save(dir);
        
    }  
    
    private void initCells(ImageInt nucLabel, ImageInt regionLabel) {
        popNuclei = new Objects3DPopulation(nucLabel);
        popRegions = new Objects3DPopulation(regionLabel, 1); // exclude value 1 used by borders

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer, Cell>(popNuclei.getNbObjects());

        // get nucleus label for each region
        int c = 1;
        for (Object3D region : popRegions.getObjectsList()) {
            int nuc = (int) region.getPixModeNonZero(nucLabel);
            Cell cell = new Cell();
            cell.region = region;
            cell.nucleus = popNuclei.getObjectByValue(nuc);
            popCells.add(cell);
            cell.id = c++;
            region2Cell.put(region.getValue(), cell);
            nucleus2Cell.put(nuc, cell);
        }
    }



    private void computeTypeCellsNucleusSLIC(ImageInt label) {
        int[] histype = {0, 0, 0, 0};
        for (Cell C : popCells) {
            Object3D nucleus = C.nucleus;
            // get list of regions inside nucleus
            //IJ.log("Nuc="+nucleus+" label="+label);
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
            if((nbA>nbB)&&(nbA>nbD)) C.type=ALPHA;
            if((nbB>nbA)&&(nbB>nbD)) C.type=BETA;
            if((nbD>nbA)&&(nbD>nbB)) C.type=DELTA;
            
            histype[C.type]++;
        }
        IJ.log("NB : " + histype[UNLABELED] + " " + histype[ALPHA] + " " + histype[BETA] + " " + histype[DELTA]);
    }
    
          



    private void computeAssoCells(ImageHandler img, int borderValue) {

        ArrayList<Integer>[] nei = computeAsso(img, borderValue);

        // index refers to region label
        for (int i = 0; i < nei.length; i++) {
            if (nei[i].isEmpty()) {
                continue;
            }
            Cell C = region2Cell.get(i);
            if (C != null) {
                C.nei1 = new ArrayList<Cell>();
                for (int j : nei[i]) {
                    Cell C1 = region2Cell.get(j);
                    if ((C1 != null) && (C1 != C)) {
                        C.nei1.add(C1);
                    }
                }
            }
        }
    }

    private void computeAssoCells2() {
        for (Cell C : popCells) {
            C.nei2 = new ArrayList<Cell>();
            for (Cell C1 : C.nei1) {
                for (Cell C2 : C1.nei1) {
                    if ((C2 != C) && (!C.nei2.contains(C2))) {
                        C.nei2.add(C2);
                    }
                }
            }
        }
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
