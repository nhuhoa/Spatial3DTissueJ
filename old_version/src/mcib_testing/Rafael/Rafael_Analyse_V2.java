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
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author thomasb
 */
public class Rafael_Analyse_V2 implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 1;
    private final int BETA = 2;
    private final int DELTA = 3;

    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler img;

    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;

    @Override
    public void run(String arg) {
//        ImagePlus plus = WindowManager.getImage("dapi-seg-wat.tif");
//        String dir = IJ.getDirectory("image");
//        ImageInt imgWat = ImageInt.wrap(plus);
//        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
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
//        // signals
//        signals = new ImageHandler[]{imgAlpha, imgBeta, imgDelta};
//        
//        IJ.log("Initialization Tissue Analysis ...");
//        initCells(imgSeg, imgWat);
//        IJ.log("Type ...");
//        computeTypeCellsLayer(10);
//        popNuclei.saveObjects(dir + "Nuclei.zip");
//        popRegions.saveObjects(dir + "Regions.zip");

//        Manager3DType manager = new Manager3DType();
//        manager.addObjects3DPopulation(popNuclei);
        ImagePlus plus = WindowManager.getImage("dapi-seg-SLIC-25000-wat.tif");
        img = ImageHandler.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        String directory = IJ.getDirectory("image");
        IJ.log("test dir is : " + directory);
        IJ.log("Reading data from " + directory);
        popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        popRegions.loadObjects(directory + "Regions.zip");
        popNuclei.loadObjects(directory + "Nuclei.zip");
        
        // read file and update from user
        Manager3DType manager = new Manager3DType();
        manager.directory = directory;
        manager.addObjects3DPopulation(popNuclei);
        manager.load();

        ImageHandler imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
        imgWat.show("Reg");
        ImageHandler imgWat2 = img.createSameDimensions();
        popNuclei.draw(imgWat2);
        imgWat2.show("Nuc");
                    
        // init cells 
        initCells2();
        // THE  REST SHOULD BE DONE AFTER VALIDATING TYPES
        //computeTypeCellsNucleus(0.1);        
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);// 0 for thomas´s cell seg; -1 for farsight
        IJ.log("Association2 ...");
        computeAssoCells2();
        computeAllAssociations();
                
        // draw 
        drawCellTypes(false).show("TYPE");
        drawCellNuclei().show("NUC");
        drawCellRegions().show("REG");
        drawBetaDeltaContacts().show("BD");
        
        // TEST LAYER
        computeLayer(ALPHA);
        drawCellLayer(false).show("LAYER ALPHA");
        computeLayerBorder();
        drawCellLayer(false).show("LAYER BORDER");
        computeResults("");
        
        // randomizing
        randomizeTypeCell();
        computeAllAssociations();
        computeLayer(ALPHA);
        computeLayerBorder();
        computeResults("_rand");
        drawCellLayer(false).show("LAYER ALPHA RAND");
        
        IJ.log("Finished");
               
    }

    private ImageHandler drawBetaDeltaContacts() {
        ImageHandler draw = new ImageShort("NUC", img.sizeX, img.sizeY, img.sizeZ);

        for (Cell C : popCells) {
            int col;
            if (C.type == DELTA) {
                col = 4;
                C.region.draw(draw, col);
            } else if (C.type == BETA) {
                boolean co = C.hasContact1(DELTA);
                if (co) {
                    col = 1;
                } else {
                    co = C.hasContact2(DELTA);
                    if (co) {
                        col = 2;
                    } else {
                        col = 3;
                    }
                }
                C.region.draw(draw, col);
            }
        }

        return draw;
    }

    private void computeResults(String name) 
    {
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
                int[] nei = C.computeNei1TypeHisto();
                rt.setValue("Nei1_A", c, nei[1]);
                rt.setValue("Nei1_B", c, nei[2]);
                rt.setValue("Nei1_D", c, nei[3]);
                nei = C.computeNei2TypeHisto();
                rt.setValue("Nei2_A", c, nei[1]);
                rt.setValue("Nei2_B", c, nei[2]);
                rt.setValue("Nei2_D", c, nei[3]);
                // nb in range distance
                double dist = 20;
                nei = C.computeNeiRangeTypeHisto(dist);
                rt.setValue("NeiR_A", c, nei[1]);
                rt.setValue("NeiR_B", c, nei[2]);
                rt.setValue("NeiR_D", c, nei[3]);
                rt.setValue("Layer_A", c, C.info1);
                rt.setValue("Layer_P", c, C.info2);

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
        rtALPHA.show("Alpha" + name);
        rtBETA.show("Beta" + name);
        rtDELTA.show("Delta" + name);
    }

    private void computeAllAssociations() {
        float aa = 0, bb = 0, dd = 0, ab = 0, ad = 0, bd = 0;
        float na = 0, nb = 0, nd = 0;

        for (Cell C : popCells) {
            if (C.type == ALPHA) {
                na++;
                int[] ne = C.computeNei1TypeHisto();
                aa += ne[ALPHA];
                ab += ne[BETA];
                ad += ne[DELTA];
            } else if (C.type == BETA) {
                nb++;
                int[] ne = C.computeNei1TypeHisto();
                ab += ne[ALPHA];
                bb += ne[BETA];
                bd += ne[DELTA];
            } else if (C.type == DELTA) {
                nd++;
                int[] ne = C.computeNei1TypeHisto();
                ad += ne[ALPHA];
                bd += ne[BETA];
                dd += ne[DELTA];
            }
        }
        aa /= 2;
        bb /= 2;
        dd /= 2;
        ab /= 2;
        ad /= 2;
        bd /= 2;
        float sum = aa + bb + dd + ab + ad + bd;
        float sumnb = na + nb + nd;

        IJ.log("Nb : na=" + na + " nb=" + nb + " nd=" + nd + " sum=" + sumnb);
        IJ.log("Nb : pa=" + na / sumnb + " pb=" + nb / sumnb + " pd=" + nd / sumnb);
        IJ.log("Assos : aa=" + aa + " bb=" + bb + " dd=" + dd + " ab=" + ab + " ad=" + ad + " bd=" + bd + " all=" + sum);
        IJ.log("Assos : aa=" + aa / sum + " bb=" + bb / sum + " dd=" + dd / sum + " ab=" + ab / sum + " ad=" + ad / sum + " bd=" + bd / sum);

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
                    if (((BorderValue>=0)&&(img.getPixel(x, y, z) == BorderValue))||(BorderValue<0)) {
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

        byte[] typeShuffle = new byte[okcellsorig.size()];
        for (int i = 0; i < okcellsorig.size(); i++) {
            typeShuffle[i] = popCells.get(okcellsorig.get(i)).type;
        }
        for (int i = 0; i < okcells.size(); i++) {
            popCells.get(okcells.get(i)).type = typeShuffle[i];
        }

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

        for (Cell C : popCells) {
            C.info1 = C.layer;
        }

        layer0 = null;
        layer1 = null;

    }

    private void computeLayerBorder() {
        // init
        for (Cell C : popCells) {
            C.layer = -1;
        }
        // layer --> ARRAY LIST
        ArrayList<Cell> layer0 = new ArrayList<Cell>();
        ArrayList<Cell> layer1;

        for (Cell C : popCells) {
            if (C.type == UNLABELED) {
                if (C.computeTouchBorderXYImageRegion(img)) {
                    C.layer = 0;
                    layer0.add(C);
                }
            }
        }

        boolean stop = false;
        while (!stop) {
            stop = true;
            layer1 = new ArrayList<Cell>();
            for (Cell C : layer0) {
                for (Cell C1 : C.nei1) {
                    if ((C1.type == UNLABELED) && (C1.layer == -1) && (!C1.hasContact1(ALPHA)) && (!C1.hasContact1(BETA)) && (!C1.hasContact1(DELTA))) {
                        C1.layer = C.layer + 1;
                        layer1.add(C1);
                        stop = false;
                    }
                }
            }
            layer0 = layer1;
        }

        for (Cell C : popCells) {
            if ((C.layer >= 0) && (C.layer < 10)) {
                C.border = 1;
                for (Cell C1 : C.nei1) {
                    if ((C1.type == UNLABELED) && (C1.layer == -1)) {
                        C1.layer = 10;
                        C1.border = 1;
                    }
                }
            }
        }

        // init
        for (Cell C : popCells) {
            C.layer = -1;
        }
        // layer --> ARRAY LIST
        layer0 = new ArrayList<Cell>();

        for (Cell C : popCells) {
            if (C.border == 1) {
                C.layer = 0;
                layer0.add(C);
            }
        }

        stop = false;
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

        for (Cell C : popCells) {
            C.info2 = C.layer;
        }

        layer0 = null;
        layer1 = null;

    }

    private ImageHandler drawBorderCells(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", img.sizeX, img.sizeY, img.sizeZ);
//        int val;
//        for (Cell C : popCells) {
//            if (C.type == UNLABELED) {
//                val = 1;
//                if (C.hasContact1(ALPHA) || C.hasContact1(BETA) || C.hasContact1(DELTA)) {
//                    val = 2;
//                }
//                if (nuc) {
//                    C.nucleus.draw(draw, val);
//                } else {
//                    C.region.draw(draw, val);
//                }
//            }
//        }

        int val;
        for (Cell C : popCells) {
            if (C.type == UNLABELED) {
                val = 1;
                if (C.computeTouchBorderXYImageRegion(img)) {
                    val = 2;
                }

                if (nuc) {
                    C.nucleus.draw(draw, val);
                } else {
                    C.region.draw(draw, val);
                }
            }
        }

        return draw;

    }

    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", img.sizeX, img.sizeY, img.sizeZ);

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
        ImageHandler draw = new ImageShort("TYPE", img.sizeX, img.sizeY, img.sizeZ);

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
        ImageHandler draw = new ImageShort("NUC", img.sizeX, img.sizeY, img.sizeZ);

        for (Cell C : popCells) {
            C.nucleus.draw(draw, C.id);
        }

        return draw;
    }

    private ImageHandler drawCellRegions() {
        ImageHandler draw = new ImageShort("REG", img.sizeX, img.sizeY, img.sizeZ);

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
        int border = -1;
        // info
        int info1 = -1;
        int info2 = -1;

        ArrayList<Cell> nei1 = null;
        ArrayList<Cell> nei2 = null;

        @Override
        public String toString() {
            return "(" + nucleus.getValue() + ", " + region.getValue() + ", " + type + ")";
        }

        private int[] computeNeiTypeHisto(ArrayList<Cell> nei) {
            int[] res = new int[4];
            // al types , includin unlabelled
            for (Cell C : nei) {
                //if (C.type > 0) {
                res[C.type]++;
                //}
            }
            return res;
        }

        private void computeBorder() {
            int[] ne = computeNei1TypeHisto();
            border = ne[0];
        }

        public int nbBorder() {
            computeBorder();
            return border;
        }

        public int[] computeNei1TypeHisto() {
            return computeNeiTypeHisto(nei1);
        }

        public int[] computeNei2TypeHisto() {
            return computeNeiTypeHisto(nei2);
        }

        public int[] computeNeiRangeTypeHisto(double dist) {
            return computeNeiTypeHisto(getCellsRange(dist));
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

        public boolean computeTouchBorderXYImageRegion(ImageHandler img) {
            if (region.getXmin() <= 0) {
                return true;
            }
            if (region.getXmax() >= img.sizeX - 1) {
                return true;
            }
            if (region.getYmin() <= 0) {
                return true;
            }
            if (region.getYmax() >= img.sizeY - 1) {
                return true;
            }

            return false;
        }
    }
}
