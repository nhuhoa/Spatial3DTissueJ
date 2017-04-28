/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.organisation;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.ThreadUtil;
import mcib3d.utils.exceptionPrinter;
import mcib_testing.Utils.Cell;
import mcib_testing.Utils.Objects3DPopulationGrid;

/**
 *
 * @author tranhoa
 */
public class Protrusion_V2_ implements PlugIn
{
    
    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
//    private final int COLOCAB = 23;
//    private final int COLOCAD = 13;
//    private final int COLOCBD = 12;

    private Objects3DPopulationGrid popRegions = null;
    private Objects3DPopulationGrid popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg = null, imgWat = null;
    ImageHandler wat, dapi;

    ArrayList<Cell> popCells = null;
    ArrayList<Cell> popInputCells = null;
    ArrayList<Cell> popS = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    HashMap<Integer, Integer> cacheId = null;
    String dir = null, subdir = null;
    ResultsTable rtLayer = new ResultsTable();
    ResultsTable rtInf = new ResultsTable();
    ResultsTable rtContact = new ResultsTable();
    public boolean verbose = false, disp = false;
    public double ratioMicro2Px=0, ratioZ2XY=0;
    int minDistance = 0, maxDistance = 35;
    public byte typeCellDef = 0;
    public byte typeInfDef = 0;
    public int op = 0, reg = 0;
    @Override
    public void run(String arg) 
    {
        if (Dialog())
        {
            IJ.log("Analysing ...");
            IJ.log("Reading data Nuclei.zip & Regions.zip from " + dir);
            popRegions = new Objects3DPopulationGrid();
            popNuclei = new Objects3DPopulationGrid();
            popRegions.loadObjects(dir + "Regions.zip");
            popNuclei.loadObjects(dir + "Nuclei.zip");

            imgWat = wat.createSameDimensions();
            popRegions.draw(imgWat);
            //imgWat.show("Reg");
            imgSeg = wat.createSameDimensions();
            popNuclei.draw(imgSeg);
            //imgSeg.show("Nuc");

            // init cells 
            initCells2();
            IJ.log("Nb cells in this tissue : " + popCells.size());
            IJ.log("Association 1 ...");
            computeAssoCells(imgWat, 0);
            IJ.log("**************PROTRUSION-EXPANDED**************************");
            byte typeInf = DELTA;
    //        
            Calibration cal = dapi.getCalibration();
            ratioMicro2Px = getRatioMicro2Pixel(cal);
            ratioZ2XY = getRatio_Z_XY(cal);
            IJ.log("Ratio micro to pixel: "+ratioMicro2Px);
            IJ.log("Ratio XY to Z: "+ratioZ2XY);
            if(op==0)
            {
                getTotalContacted(typeCellDef, typeInfDef, minDistance, maxDistance, reg);
            }
            else{
                getEachCellProContacted(typeCellDef, typeInf, minDistance, maxDistance, reg);
            }
            IJ.log("Finished");
        }    
    }
    private boolean Dialog() 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open");
            return false;
        }
        String[] titles = new String[wList.length];
        if (wList.length<2) {
            IJ.error("Open at least 2 images: watershed image and nuclei dapi image");
            return false;
        }
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        String[] typeObserved = {"BETA","ALPHA"};
        String[] typeInf = {"DELTA", "BETA","ALPHA"};
        String[] option = {"ALL DELTA CELLS","EACH DELTA CELL"};
        String[] regionObs = {"NUCLEUS REGION","CELL REGION"};
        
        
        
        String unit = "microns";
        GenericDialog gd = new GenericDialog("3D Tissue Spatial Analysis - Protrusion");
        gd.addChoice("Watershed_Image :", titles, titles[0]);
        gd.addChoice("DAPI_Image :", titles, titles[1]);
        gd.addNumericField("Min_Distance: ", minDistance, 0 , 10, unit);
        gd.addNumericField("Max_Distance: ", maxDistance, 0 , 10, unit);
        gd.addChoice("Type_Cell_Observed : ", typeObserved, typeObserved[typeCellDef]);
        gd.addChoice("Type_Cell_Inf : ", typeInf, typeInf[typeInfDef]);
        gd.addChoice("Range_Observed : ", option, option[op]);
        gd.addChoice("Region_Observed : ", regionObs, regionObs[reg]);
        gd.addCheckbox("Show_Protrusion_Image", verbose);
        gd.addCheckbox("Show_Observed_Cells", disp);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        int[] index = new int[2];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        minDistance = (int) gd.getNextNumber();
        maxDistance = (int) gd.getNextNumber();
        typeCellDef = (byte)(gd.getNextChoiceIndex()+2);
        typeInfDef = (byte)(gd.getNextChoiceIndex()+1);
        op =  (int)(gd.getNextChoiceIndex());
        reg =  (int)(gd.getNextChoiceIndex());
        verbose = gd.getNextBoolean();
        disp = gd.getNextBoolean();
        if(gd.wasOKed())
        {
            IJ.log("---------------------------------------------------");
            IJ.log("watershed img: "+wList[index[0]]+" dapi: "+wList[index[1]]+"\n Min distance: "+minDistance
                   +" microns   &   max distance: "+maxDistance + " microns \n & type cell observed: "+typeObserved[typeCellDef-2]
                   +" & type cell inf: "+typeInf[typeInfDef-1]+"\n&  range(0-all, 1-each): "+op+"  & region (0-nuc, 1-reg): " + reg);
            IJ.log("---------------------------------------------------");
        } 
        ImagePlus plus;
        plus = WindowManager.getImage(wList[index[0]]);
        wat = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        ImagePlus plus1;
        plus1 = WindowManager.getImage(wList[index[1]]);
        dapi = ImageInt.wrap(plus1);
        dir = IJ.getDirectory("image");
        File fl = new File(dir+"protrusion");
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = dir + "protrusion/";
                IJ.log("Create folder protrusion in : "+dir);
            } else {
                subdir = dir;
            }
        }
        else{
            subdir = dir + "protrusion/";
        }
        return gd.wasOKed();
    }
    private void initCells2() {

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer, Cell>(popNuclei.getNbObjects());

        // get nucleus label for each region
        int c = 1;
        for (int i = 0; i < popRegions.getNbObjects(); i++) 
        {
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

    private double getRatioMicro2Pixel(Calibration cal)
    {
        return (1 * 1.0f / cal.pixelWidth);
    }        
    private double getRatio_Z_XY(Calibration cal)
    {
        return (cal.pixelWidth * 1.0f / cal.pixelDepth);
    }
    private void getTotalContacted(byte typeObserved, byte typeInf, int minDistance, int maxDistance, int reg)
    {
        ArrayList<Cell> popTmp = new ArrayList<Cell>();
        String idStr = " Id Delta Cell Before: ";
        for (Cell C : popCells) 
        {
            if (C.type == typeInf) 
            {
                idStr += "  "+C.id;
                popTmp.add(C);
            }
        }
        IJ.log(idStr);
        int st = 0;
        //testing
//        ImageByte res = protrusion(minDistance, typeInf, popTmp, reg, st);
//        res.show();
        computeContacted(typeObserved, typeInf, minDistance, maxDistance, popTmp, reg);
    }
    private void getEachCellProContacted(byte typeObserved, byte typeInf, int minDistance, int maxDistance, int reg)
    {
        int st = 0;
        for (Cell C : popCells) 
        {
            if (C.type == typeInf && st < 2) 
            {
                st++;
                ArrayList<Cell> popTmp = new ArrayList<Cell>();
                popTmp.add(C);
                IJ.log("------------------ "+st+" ----------------------");
                computeContacted(typeObserved, typeInf, minDistance, maxDistance, popTmp, reg);
            }
        }
    }        
    public void computeContacted(byte typeObserved, byte typeInf, int minDistance, int maxDistance, ArrayList<Cell> popTmp, int reg)
    {
        IJ.log("Nb delta cells: "+popTmp.size());
        ArrayList<Integer> idDelLs = new ArrayList<Integer>();
        ArrayList<Cell> lsCells = new ArrayList<Cell>();
        for (Cell C : popTmp) 
        {
            idDelLs.add(C.id);
            for (Cell Ci : C.nei1) 
            {
                if (Ci.type == typeObserved && !lsCells.contains(Ci)) 
                {
                    lsCells.add(Ci);
                }
            }
        }
        String title = "", title1 = "";
          
        if(typeInf== DELTA)
        {
            title += "_Del_";
        }    
        if(typeObserved==BETA)
        {
            title1 += "Beta"; 
            title += title1;
        }    
        if(typeObserved==ALPHA)
        {
            title1 += "Al";
            title += title1;
        }
        title += "_Contacted";
        
        //compute xaxis 
        int nbCells = popCells.size();
        int sizeDis = maxDistance - minDistance + 1;
//        ArrayUtil xAxis = new ArrayUtil(sizeDis);
//        ArrayUtil nbCellInf = new ArrayUtil(sizeDis);
        
        ResultsTable rtCellsPro = new ResultsTable();
        IJ.log("----------------BEFORE---------------------");
        int count = 0;
        for (Cell C : popCells) 
        {
            if(C.type == typeObserved)
            {
                count++;
            }    
//            if (C.type == typeInf) 
//            {
//                for (Cell Ci : C.nei1) 
//                {
//                    if (Ci.type == typeObserved && !lsCells.contains(Ci)) 
//                    {
//                        lsCells.add(Ci);
//                    }
//                }
////                nbInf += nbLabelled;
//            }
        }
        int nbInf = lsCells.size();
//        nbCellInf.putValue(0, lsCells.size());
        rtCellsPro.incrementCounter();
        rtCellsPro.setValue("MicroM", 0, 0);
        rtCellsPro.setValue("nbCells"+title, 0, nbInf);
        rtCellsPro.setValue("nbCells"+title1, 0, count);
        rtCellsPro.setValue("percent_Beta_Total"+title1, 0,  nbInf * 100.0f / count);
        rtCellsPro.setValue("nbCells_ABD", 0, nbCells);
        rtCellsPro.setValue("percent_"+title1+"_ABD", 0, nbInf * 100.0f / nbCells);
        
        IJ.log("----------------AFTER PROTRUSION---------------------");
        //compute nb of influenced
        int nb = 1;
        for(int d=minDistance; d<=maxDistance; d++)
        {
            IJ.log("-------------------------------------------");
            IJ.log("Protrusion with distance: "+d+" micro m");
            ArrayList<Cell> lsCells2 = new ArrayList<Cell>();
            ImageByte res = protrusion(d, typeInf, popTmp, reg);
            if(verbose)
            {
                res.show();
            }    
            Objects3DPopulationGrid popObjs = new Objects3DPopulationGrid((ImageInt)res);
            for (int i = 0; i < popObjs.getNbObjects(); i++) 
            {
                Object3D o = popObjs.getObject(i);
                ArrayUtil ac = o.listValues(imgWat);
                ac = ac.distinctValues();
                for (int j = 0; j < ac.getSize(); j++) 
                {
                    if (!idDelLs.contains(ac.getValueInt(j)) && ac.getValueInt(j) != 0) 
                    {
                        int idx = ac.getValueInt(j);
                        Cell C1 = region2Cell.get(idx);
                        if ((C1 != null) && C1.type==typeObserved && !lsCells2.contains(C1)) 
                        {
                            lsCells2.add(C1);
                        }
                    }
                }
            }
            
            int nbInfAfter = lsCells2.size();
            IJ.log("Distance:  "+d+ "   inscreased: "+ lsCells2.size());
            rtCellsPro.incrementCounter();
//            nbCellInf.putValue(nb, nbInf+lsCells2.size());
            rtCellsPro.setValue("MicroM", nb, d);
            rtCellsPro.setValue("nbCells"+title, nb, nbInfAfter);
            rtCellsPro.setValue("nbCells"+title1, nb, count);
            rtCellsPro.setValue("percent_Beta_Total"+title1, nb, (nbInfAfter)* 100.0f/count);
            rtCellsPro.setValue("nbCells_ABD", nb, nbCells);
            rtCellsPro.setValue("percent_"+title1+"_ABD", nb, (nbInfAfter)* 100.0f/nbCells);
            nb++;
        }
        if(popTmp.size()>1)
        {
            title += "_All_";
        }  
        else{
            title += "_objId_"+popTmp.get(0).id;
        }
//        rtCellsPro.show(title);
        try {
            rtCellsPro.saveAs(subdir+title+".csv");
            IJ.log("Saving the output table in the "+subdir);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
//        Plot layerPlot = PlotDraw.createPlot("Level 1 Influenced","distance (micro m)","nb cells influenced by delta cell", xAxis, nbCellInf);
//        PlotWindow.noGridLines = true;
//        layerPlot.draw();
//        PlotWindow plotLayer = layerPlot.show();
    }
    public ImageByte protrusion(double microDistance, byte typeInf, ArrayList<Cell> popTmp, int reg)
    {
        IJ.log("Nb objs of seg image: "+ popTmp.size());
        String title = "dilated_"+microDistance+"_"+reg;
        if(popTmp.size()>1)
        {
            title += "_all_objs";
        }    
        else{
            title += "_idObj_"+popTmp.get(0).id;
        }
        ImageHandler label = new ImageShort("seg_", wat.sizeX, wat.sizeY, wat.sizeZ);
//        ImageHandler dilateImg = new ImageShort(title, img.sizeX, img.sizeY, img.sizeZ);
        for (Cell C : popTmp) 
        {
            Object3D region = C.region;
            Object3D nuc = C.nucleus;
            if(reg==0)  //observed nucleus region
            {
                nuc.draw(label, region.getValue());
            }
            else{        //observed cell region
                region.draw(label, region.getValue());
            }
        }
        if(disp)
        {
            label.show();
        }
        try {
            int nbCPUs = ThreadUtil.getNbCpus();
            float rDilate = (float) (microDistance * ratioMicro2Px);
            float radiusZ = (float) (rDilate*ratioZ2XY);
            ImageFloat edm = EDT.run(label, 0, 1, rDilate / radiusZ, true, nbCPUs);
            ImageByte temp = edm.threshold(rDilate, true, false);
            edm.flush();
            edm = null;
            temp.setOffset(label);
            temp.setScale(label);
            temp.setTitle(title);
            return temp;
        } catch (Exception e) {
            exceptionPrinter.print(e, null, true);
        }
        return null;
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
