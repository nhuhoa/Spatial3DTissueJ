/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Rafael;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.ResultsTable;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DPoint;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author thomasb
 */
public class Rafael_TypeSLIC_V2 implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;

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
        plus = WindowManager.getImage("dapi-seg-SLIC-25000-wat.tif");
        
        WindowManager.setTempCurrentImage(plus);
        String dir = IJ.getDirectory("image");
        ImageInt imgWat = ImageInt.wrap(plus);
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg-SLIC-25000-correct.tif"));
        
        
        //ImageInt imgAlphaFillLabel = ImageInt.wrap(WindowManager.getImage("alpha-fill-holes.tif"));
//        ImageInt imgAlpha = ImageInt.wrap(WindowManager.getImage("alpha-bin.tif"));
//        if (imgAlpha.getMax() == 255) {
//            imgAlpha.divideByValue(255);// 0-1
//        }
//        ImageInt imgBeta = ImageInt.wrap(WindowManager.getImage("beta-bin.tif"));
//        if (imgBeta.getMax() == 255) {
//            imgBeta.divideByValue(255);// 0-1
//        }
//        ImageInt imgDelta = ImageInt.wrap(WindowManager.getImage("delta-bin.tif"));
//        if (imgDelta.getMax() == 255) {
//            imgDelta.divideByValue(255);// 0-1
//        }
        // signals
        //signals = new ImageHandler[]{imgAlpha, imgBeta, imgDelta};
        

        IJ.log("Initialization Tissue Analysis ...");
        initCells(imgSeg, imgWat);
        IJ.log("Type ...");
        
        this.computeTypeCellsNucleusSLIC(imgLabel);
                
        //computeTypeCellsLayer(10);
        popNuclei.saveObjects(dir + "Nuclei.zip");
        popRegions.saveObjects(dir + "Regions.zip");
        IJ.log("save the result");
        Manager3DType manager = new Manager3DType();
        manager.directory=dir;
        //ArrayList<Object3D> sortNuc = popNuclei.getObjectsList();
        //Collections.sort(sortNuc, new CompareNucleus());
        //Objects3DPopulation popNucleiSort=new Objects3DPopulation();
        //popNucleiSort.addObjects(sortNuc);
        manager.addObjects3DPopulation(popNuclei);
        manager.save();
        // THE  REST SHOULD BE DONE AFTER VALIDATING TYPES
//        computeTypeCellsNucleus(0.1);
        // randomizing
        //randomizeTypeCell();
//        IJ.log("Association1 ...");
//        computeAssoCells(imgWat, 1);
//        IJ.log("Association2 ...");
//        computeAssoCells2();
//
//        // draw 
          drawCellTypes(false).show("TYPE");
//        drawCellNuclei().show("NUC");
//        drawCellRegions().show("REG");
//        drawBetaDeltaContacts().show("BD");
//
//        computeResults();
//
//        // TEST LAYER
//        computeLayer(DELTA);
//        drawCellLayer(false).show("LAYER");
//
//        // TEST Manager type
//        Manager3DType type = new Manager3DType();
//        type.addObjects3DPopulation(popNuclei);
          IJ.log("Finished");

    }

//    private ImageHandler drawBetaDeltaContacts() {
//        ImageHandler draw = new ImageShort("NUC", signals[0].sizeX, signals[0].sizeY, signals[0].sizeZ);
//
//        for (Cell C : popCells) {
//            int col;
//            if (C.type == DELTA) {
//                col = 4;
//                C.region.draw(draw, col);
//            } else if (C.type == BETA) {
//                boolean co = C.hasContact1(DELTA);
//                if (co) {
//                    col = 1;
//                } else {
//                    co = C.hasContact2(DELTA);
//                    if (co) {
//                        col = 2;
//                    } else {
//                        col = 3;
//                    }
//                }
//                C.region.draw(draw, col);
//            }
//        }
//
//        return draw;
//    }

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

    private void randomizeTypeCell() {
        ArrayList<Integer> okcells = new ArrayList();
        for (int i = 0; i < popCells.size(); i++) {
            int ty = popCells.get(i).type;
            if (ty != UNLABELED) {
                okcells.add(i);
            }
        }

        ArrayList<Integer> okcellsorig = new ArrayList();
        okcellsorig.addAll(okcells);
        Collections.shuffle(okcells);

        for (int i = 0; i < okcells.size(); i++) {
            popCells.get(okcells.get(i)).type = popCells.get(okcellsorig.get(i)).type;
        }
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

//    private void computeTypeCellsNucleus(double threshold) {
//        int[] histype = {0, 0, 0, 0};
//        for (Cell C : popCells) {
//            C.type = UNLABELED;
//            Object3D nucleus = C.nucleus;
//            ArrayUtil ids = new ArrayUtil(signals.length);
//            for (byte i = 0; i < signals.length; i++) {
//                double IDN = nucleus.getIntegratedDensity(signals[i]); // = volume since all voxels = 1 
//                ids.addValue(i, IDN);
//
//            }
//            double volN = nucleus.getVolumePixels();
//            int[] idx = ids.sortIndexShellMeitzner();
//            C.type = (byte) (idx[idx.length - 1] + 1);
//            if (ids.getValue(idx[idx.length - 1]) < 2 * ids.getValue(idx[idx.length - 2]) && (ids.getValue(idx[idx.length - 2]) / volN > threshold)) {
//                //C.type = UNLABELED;
//                IJ.log("coloc " + C.nucleus + " " + C.type);
//            }
//            if (ids.getValue(idx[idx.length - 1]) / volN < threshold) {
//                C.type = UNLABELED;
//            }
//            histype[C.type]++;
//        }
//
//        Object3D nuc0 = popCells.get(popCells.size() / 2).nucleus;
//        Object3D reg0 = popCells.get(popCells.size() / 2).region;
//        ImageHandler test = signals[0].createSameDimensions();
//        reg0.draw(test, 255);
//        nuc0.draw(test, 0);
//        Object3DVoxels layer0 = nuc0.getLayerEVFObject((ImageInt) test, 0.1f);
//        test.draw(layer0, 128);
//        test.show("test");
//
//        IJ.log("NB : " + histype[0] + " " + histype[1] + " " + histype[2] + " " + histype[3]);
//    }

    
    private void computeTypeCellsNucleusSLIC(ImageInt label) 
    {
        double[] histype = {0, 0, 0, 0};
        for (Cell C : popCells) 
        {
            Object3D nucleus = C.nucleus;
            Object3D region = C.region;
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
            nucleus.setType(C.type);
            nucleus.setName("Nuc" + C.id);
            nucleus.setValue(C.id);
            region.setType(C.type);
            region.setName("Reg" + C.id);
            region.setValue(C.id);
        }
        IJ.log("Number of cell type is : " + histype[UNLABELED] + " " + histype[ALPHA] + " " + histype[BETA] + " " + histype[DELTA]);
        
        
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);

        for (Cell C : popCells) {
            if(C.type==UNLABELED){
               Object3D nucleus = C.nucleus;
               nucleus.draw(unLabelledCell, 255);
            }     
        }
        unLabelledCell.show("Unlabelled");
        unLabelledCell.setTitle("unlabelled");
        //unLabelledCell.save(dir);


        double[] y = {0, 1, 2, 3};
        PlotWindow.noGridLines = false; // draw grid lines
        Plot plot = new Plot("Plot","X Axis","Y Axis", y, histype);
        plot.setLimits(0, 5, 5, 900);
        plot.setLineWidth(1);
        plot.changeFont(new Font("Helvetica", Font.PLAIN, 16));
        plot.setColor(Color.blue);
        plot.show();
        IJ.log("NB : unlabelled : " + histype[UNLABELED] + " alpha: " + histype[ALPHA] + " beta: " + histype[BETA] 
                + " delta: " + histype[DELTA]);
    }

//    private void computeTypeCellsLayer(double threshold) {
//        int[] histype = {0, 0, 0, 0};
//
//        // create EDT for all nuclei
//        ImageHandler dup = signals[0].createSameDimensions();
//        dup.setCalibration(signals[0].getCalibration());
//        for (Cell C : popCells) {
//            C.nucleus.draw(dup, 255);
//        }
//        ImageFloat edt = EDT.run(dup, 128, (float) dup.getCalibration().pixelWidth, (float) dup.getCalibration().pixelDepth, true, 0);
//        //edt.show("edt layer");
//
//        ImageHandler test = signals[0].createSameDimensions();
//        for (Cell C : popCells) {
//            IJ.showStatus("Cell type " + C.id + "              ");
//            C.type = UNLABELED;
//            Object3D nucleus = C.nucleus;
//            Object3D region = C.region;
//            region.draw(test, 255);
//            nucleus.draw(test, 0);
//            ImageFloat edtdup = edt.duplicate();
//            EDT.normalizeDistanceMap(edtdup, (ImageInt) test, true);
//            ImageByte bin = edtdup.thresholdRangeExclusive(0, 0.1f);
//            if (C.id == popCells.size() / 2) {
//                //edtdup.show("EVF");
//                //bin.show("bin evf");
//            }
//
//            ArrayUtil ids = new ArrayUtil(signals.length);
//            for (byte i = 0; i < signals.length; i++) {
//                double IDN = region.getIntegratedDensity(signals[i], bin); // = volume since all voxels = 1 
//                ids.addValue(i, IDN);
//            }
//            int[] idx = ids.sortIndexShellMeitzner();
//            C.type = (byte) (idx[idx.length - 1] + 1);
//            if (ids.getValue(idx[idx.length - 1]) < 2 * ids.getValue(idx[idx.length - 2]) && (ids.getValue(idx[idx.length - 2]) > threshold)) {
//                //C.type = UNLABELED;
//                nucleus.setComment("*");
//                IJ.log("coloc " + C.id + " : " + (idx[idx.length - 1] + 1) + " " + ids.getValue(idx[idx.length - 1]) + " " + (idx[idx.length - 2] + 1) + " " + ids.getValue(idx[idx.length - 2]));
//            }
//            if (ids.getValue(idx[idx.length - 1]) < threshold) {
//                C.type = UNLABELED;
//            }
//            histype[C.type]++;
//            // update type in Object3D
//            nucleus.setType(C.type);
//            nucleus.setName("Nuc" + C.id);
//            nucleus.setValue(C.id);
//            region.setType(C.type);
//            region.setName("Reg" + C.id);
//            region.setValue(C.id);
//            // erase in image
//            region.draw(test, 0);
//        }
//
//        IJ.log("NB : " + histype[0] + " " + histype[1] + " " + histype[2] + " " + histype[3]);
//    }

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
