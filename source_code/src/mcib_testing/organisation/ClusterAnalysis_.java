/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.organisation;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;
import mcib_testing.Utils.Objects3DPopulationGrid;
import mcib_testing.Utils.SpatialTissueAnalysis;
import stdlib.MathStat;
import stdlib.Queue;

/**
 *
 * @author tranhoa
 */
public class ClusterAnalysis_ implements ij.plugin.PlugIn
{
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int UNLABELLED = 0;
    private Objects3DPopulationGrid popRegions = null;
    private Objects3DPopulationGrid popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg=null, imgWat=null;
    ImageHandler img;

    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    String subdir = null;
    public String save_dir = " ";
    public void run(String arg) 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return;
        }
        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        boolean clusterStat = false;
        
        String[] typeObserved = {" ", "DELTA", "BETA","ALPHA"};
        byte typeCellObs = ALPHA;
        GenericDialogPlus gd = new GenericDialogPlus("Cluster Analysis ");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("3D Tissue Spatial Analysis");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : watershed image");
        gd.addMessage("Input : type observed cell");
        gd.addMessage("Output: clusters analysis ");
        gd.addMessage(" ");
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("Watershed_Image ", titles, titles[0]);
        gd.addChoice("Type_Observed_Cell : ", typeObserved, typeObserved[typeCellObs]);
        gd.addCheckbox("Clusters_Spatial_Statistic", clusterStat);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        save_dir = gd.getNextString();
        if (" ".equals(save_dir)) 
        {
                return;
        }
        
        int idxImg = gd.getNextChoiceIndex();
        typeCellObs = (byte)gd.getNextChoiceIndex();
        clusterStat = gd.getNextBoolean();
        
        if(typeCellObs==0)
        {
            IJ.log(" ");
            IJ.log("Error...");
            IJ.log("Please choose observed cell type...");
            return;
        }   
        ImagePlus plus;
        plus = WindowManager.getImage(wList[idxImg]);
        img = ImageInt.wrap(plus);
//        WindowManager.setTempCurrentImage(plus);
//        dir = IJ.getDirectory("image");
        String subStr = typeObserved[typeCellObs]+"_clusters_analysis";
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        File fl = new File(save_dir + subStr);
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = save_dir + subStr +"/";
            } else {
                subdir = save_dir;
            }
        }
        else{
            subdir = save_dir + subStr +"/";
        }
        IJ.log("====================================================================");
        IJ.log("Analysing ...");
        IJ.log("Reading data from " + save_dir);
        IJ.log("Load list objects "+" Regions.zip"+ " and Nuclei.zip from "+save_dir);
        popRegions = new Objects3DPopulationGrid();
        popNuclei = new Objects3DPopulationGrid();
        popRegions.loadObjects(save_dir + "Regions.zip");
        popNuclei.loadObjects(save_dir + "Nuclei.zip");

        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        
        // init cells 
        initCells2();
        IJ.log("Nb cells: " + popCells.size());
        //drawCellTypes(false).show("TYPE BEFORE");
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);  
        
        getMaxCluster(typeCellObs);
        if(clusterStat)
        {
            calculClusterEuclidStatistic(typeCellObs);
        }
        IJ.log("Saved the output into directory: "+subdir);
        IJ.log("Finished");
        IJ.log("====================================================================");
        IJ.selectWindow("Log");
        IJ.saveAs("Text", subdir+ "Log_"+subStr+"_Euclidean_distance.txt");
    }
    public ArrayList<Point3D> getNbCluster(byte typeCellObserved)
    {
        ArrayList<Point3D> plist = new ArrayList<Point3D>();
        //Objects3DPopulation pop = new Objects3DPopulation();
        for (Cell Ci : popCells) 
        {
            Ci.marked = false;
            Ci.connected = false;
        }
        int val = 0, maxSize = 0;
        ArrayList<Voxel3D> arrObs = new ArrayList<Voxel3D>();
        for (Cell Cii : popCells) 
        {
            if((Cii.type==typeCellObserved) && (!Cii.marked))
            {
                ArrayList<Cell> cls = BFS(Cii, typeCellObserved);
                
                if(cls.size() > 0)
                {
                    val++;
                    ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
                    for(Cell c : cls)
                    {
                        Object3DVoxels nuc = (Object3DVoxels)c.nucleus;
                        arr.addAll(nuc.getVoxels());
                    }
                    Object3DVoxels obj = new Object3DVoxels(arr);
                    Point3D cent = obj.getCenterAsPoint();
                    plist.add(cent);
//                    Voxel3D v = new Voxel3D(cent, val);
//                    ArrayList<Voxel3D> centroid = new ArrayList<Voxel3D>();
//                    centroid.add(v);
//                    Object3DVoxels objCent = new Object3DVoxels(centroid);
//                    //obj.setCalibration(calibration);
//                    obj.setName("Obj" + val);
//                    pop.addObject(obj);
                }
            }
             
        }
        return plist;
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
    private void initCells2() {

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer,Cell>(popNuclei.getNbObjects());

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
    /**
     * Breadth First Search
     * @param s start from one cell, initialisation
     * @param typeObserved
     * @return 
     */
    private ArrayList<Cell> BFS(Cell s, byte typeObserved) 
    {
        ArrayList<Cell> cls = new ArrayList<Cell>();
        Queue<Cell> q = new Queue<Cell>();
        s.marked = true;
        s.currLevel = 1;
        q.enqueue(s);
        cls.add(s);
        
        while (!q.isEmpty()) {
            Cell v = q.dequeue();
            ArrayList<Cell> nei1 = v.nei1;
            boolean flag = true;
            for(Cell ne : nei1)
            {
                if((!ne.marked) && (ne.type == typeObserved))
                {
                    ne.marked = true;
                    ne.currLevel = v.currLevel + 1;
                    cls.add(ne);
                    q.enqueue(ne);
                    flag = false;
                }    
            }
            if(flag)
            {
                v.lastnode = true;
            }    
        }
        return cls;
    }
    public void getMaxCluster(byte typeObserved)
    {
        String[] typeArr = {"", "DELTA", "BETA","ALPHA"};
        HashMap<Integer,ArrayList<Cell>> map = new HashMap<Integer,ArrayList<Cell>>();
        for (Cell Ci : popCells) 
        {
            Ci.marked = false;
            Ci.connected = false;
        }
        int val = 0, maxSize = 0; 
        for (Cell Cii : popCells) 
        {
            if((!Cii.marked) && (Cii.type==typeObserved))
            {
                ArrayList<Cell> cls = BFS(Cii, typeObserved);
                
                if(cls.size() > 0)
                {
                    val++;
                    IJ.log("Size of cluster "+val+" is: "+cls.size());
                    map.put(val, cls);
                    if(maxSize < cls.size())
                    {
                        maxSize = cls.size();
                    }    
                }
            }
             
        }
        IJ.log("Size of max cluster is: "+ maxSize);
        IJ.log("Drawing...");
        ImageHandler drawReg = new ImageShort("Clusters_Wat_"+typeArr[typeObserved], img.sizeX, img.sizeY, img.sizeZ);
        ImageHandler drawNuc = new ImageShort("Clusters_Nuc_"+typeArr[typeObserved], img.sizeX, img.sizeY, img.sizeZ);
        ResultsTable rs = new ResultsTable();
        int count = 0;
        double sC[] = new double[map.size()];
        for (Map.Entry<Integer, ArrayList<Cell>> ee : map.entrySet()) 
        {
            int key = ee.getKey();
            ArrayList<Cell> lsCells = ee.getValue();
            rs.setValue("cluster_"+key, 0, lsCells.size());
            sC[count] = lsCells.size();
            count++;
            for (Cell C : lsCells) 
            {
                C.connected = true;
                if(C.currLevel == -1)
                {
                    IJ.log("Have the problem");
                }
                else{
                    C.region.draw(drawReg, key);
                    C.nucleus.draw(drawNuc, key);
                }

            }
//            if(maxSize != 0 && lsCells.size()==maxSize)
//            {
//                for (Cell C : lsCells) 
//                {
//                    C.connected = true;
//                    if(C.currLevel == -1)
//                    {
//                        IJ.log("have the problem");
//                    }
//                    else{
//                        C.region.draw(drawReg, C.currLevel);
//                    }
//                    
//                }
//            }    
            
        }
        rs.setValue("nb_cluster", 0, map.size());
        rs.setValue("max_cluster", 0, maxSize);
        rs.setValue("aver_nbCells_in_cluster", 0, MathStat.getMean(sC));
//        rs.show("All_Clusters_"+typeArr[typeObserved]);
        ImagePlus re = drawReg.getImagePlus();
//        re.show("Clusters_Wat_"+typeArr[typeObserved]);
        ImagePlus nu = drawNuc.getImagePlus();
//        nu.show("Clusters_Nuc_"+typeArr[typeObserved]);
        try {
            rs.saveAs(subdir+"All_Clusters_"+typeArr[typeObserved]+".csv");
            IJ.log("Saving the output table All_Clusters_"+typeArr[typeObserved]+".csv"+" in the folder: "+subdir);
            IJ.saveAs(re, "Tiff", subdir+"Clusters_Wat_"+typeArr[typeObserved]+".tif");
            IJ.saveAs(nu, "Tiff", subdir+"Clusters_Nuc_"+typeArr[typeObserved]+".tif");
            IJ.log("Saving the output image "+"Clusters_Wat_"+typeArr[typeObserved]+".tif"+" in the folder: "+subdir);
            IJ.log("Saving the output image "+"Clusters_Nuc_"+typeArr[typeObserved]+".tif"+" in the folder: "+subdir);
//            re.close();
//            nu.close();
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
    }
    private void calculClusterEuclidStatistic(byte typeObserved) 
    {
        ImageHandler espaceObserved = new ImageShort("Espace Observed", img.sizeX, img.sizeY, img.sizeZ);
        for (int z = 0; z < espaceObserved.sizeZ; z++) 
        {
                for (int xy = 0; xy < espaceObserved.sizeXY; xy++) 
                {
                      espaceObserved.setPixel(xy, z, 255);
                }
        }
        ImageHandler imgM = espaceObserved.duplicate();
        ImagePlus imgMask = imgM.getImagePlus();
//        ImageHandler drawImg = new ImageShort("STATISTIC", img.sizeX, img.sizeY, img.sizeZ);
        Calibration calibration = imgMask.getCalibration();
        if (calibration == null) 
        {
            IJ.log("Image not calibrated");
            calibration = new Calibration();
            calibration.setUnit("pix");
            calibration.pixelWidth = 1;
            calibration.pixelHeight = 1;
            calibration.pixelDepth = 1;
        }
        Objects3DPopulationGrid pop = new Objects3DPopulationGrid();
        ArrayList<Point3D> pls = getNbCluster(typeObserved);
        int val = 0;
        for(Point3D p : pls)
        {
            val++;
            Voxel3D v = new Voxel3D(p, val);
            ArrayList<Voxel3D> centroid = new ArrayList<Voxel3D>();
            centroid.add(v);
            Object3DVoxels objCent = new Object3DVoxels(centroid);
            objCent.setCalibration(calibration);
            objCent.setName("Obj" + val);
            pop.addObject(objCent);
        }    
        
        IJ.log("Nb observed cells in the list: "+pop.getNbObjects());
        
        
        ArrayList<Point3D> lsCenterCells = new ArrayList<Point3D>();
        
//        ArrayList<Point3D> lsReferencePoints = new ArrayList<Point3D>();
        for (Cell cc : popCells) 
        {
            if(cc.type!=UNLABELLED && cc.type != typeObserved)
            {
                Object3D nuc = cc.nucleus;
                Point3D center = nuc.getCenterAsPoint();
                lsCenterCells.add(center);
            }   
        }
        lsCenterCells.addAll(pls);
        
//        int numPoints = 1000;
        String desc = "_clusters_stat";
        int numPoints = lsCenterCells.size();
        int numRandomSamples = 100;
        double distHardCore = 2;
        double env = 0.05;
        SpatialTissueAnalysis spa = new SpatialTissueAnalysis(numPoints, numRandomSamples, distHardCore, env);
        spa.process(pop, imgMask, lsCenterCells, "FG", true, true, true, subdir, desc);
//        rs.setValue("F_SDI_Eu"+typeObserved, 0, out.getValue(0));
//        rs.setValue("G_SDI_Eu"+typeObserved, 0, out.getValue(1));
//        rs.show("Statistic Function");
        
    }
    
}
