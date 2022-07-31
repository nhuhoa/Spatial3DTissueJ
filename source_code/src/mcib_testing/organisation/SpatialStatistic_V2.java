/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.organisation;

import mcib_testing.Utils.CellSpatialAnalysis;
import mcib_testing.Utils.Cell;
import stdlib.Vertex;
import stdlib.DirectedEdgeC;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cells3DPopulation;
import mcib_testing.Utils.Objects3DPopulationGrid;
import mcib_testing.Utils.SpatialTissueAnalysis;
import stdlib.MapUtils;
import stdlib.Queue;

/**
 *
 * @author tranhoa
 */
public class SpatialStatistic_V2 implements ij.plugin.PlugIn
{
    private final int M3 = 3;
    private final int M2 = 2;
    private final int M1 = 1;
    private final int UNLABELLED = 0;
//    private final int CLUSTER = 1;
//    private final int UNIFORM = 2;
//    private final int RANDOM = 3;
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg=null, imgWat=null;
    ImageHandler img;

    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    HashMap<Integer, Integer> cacheId = null;
    String dir = null, subdir = null;
    
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
        String[] typeObserved = {"*NONE*", "M1", "M2","M3"};
        String[] typeDistance = {"Euclidean_Distance", "Cell_Distance"};
        String[] typeStat = {"Raw_Data_Spatial_Statistic", "Random_Organization_Spatial_Statistic",
                             "Cluster_Spatial_Statistic"};
        byte typeCellObs = M3;
        byte typeCellSource = M2;
        int typeDis = 0, typeStatis = 0;
        boolean f = true, g = true;
        String func = "";
        GenericDialog gd = new GenericDialog("Spatial Statistic Analysis");
        gd.addMessage("3D Tissue Spatial Analysis");
        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : watershed image");
        gd.addMessage("Output: result of F-function, G-function.");
        gd.addMessage(" ");
        gd.addChoice("Watershed_Image ", titles, titles[0]);
        gd.addChoice("Type_Cell_Observed : ", typeObserved, typeObserved[typeCellObs]);
        gd.addChoice("Type_Cell_Source : ", typeObserved, typeObserved[typeCellSource]);
        gd.addChoice("Type_Distance : ", typeDistance, typeDistance[typeDis]);
        gd.addChoice("Spatial Statistic : ", typeStat, typeStat[typeStatis]);
        gd.addCheckbox("F-function statistic", true);
        gd.addCheckbox("G-function statistic", true);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        
        int idx = gd.getNextChoiceIndex();
        typeCellObs = (byte)gd.getNextChoiceIndex();
        typeCellSource = (byte)gd.getNextChoiceIndex();
        typeDis = (int)gd.getNextChoiceIndex();
        typeStatis = (int)gd.getNextChoiceIndex();
        f  = gd.getNextBoolean();
        g  = gd.getNextBoolean();
        if(typeCellObs==0 || typeCellSource==0)
        {
            IJ.log(" ");
            IJ.log(" ");
            IJ.log("Choose observed and source cell type...");
            return;
        }    
            
        if(f){ func += "F";}
        if(g){ func += "G";}
        ImagePlus plus;
        plus = WindowManager.getImage(wList[idx]);
        img = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        String subStr = typeObserved[typeCellObs]+"_"+typeStat[typeStatis]+"_"+typeDistance[typeDis];
        File fl = new File(dir + subStr);
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = dir + subStr +"/";
            } else {
                subdir = dir;
            }
        }
        else{
            subdir = dir + subStr +"/";
        }
        IJ.log("====================================================================");
        
        IJ.log("Compute spatial statistic of "+typeStat[typeStatis]+"    Observed distance: "+typeDistance[typeDis]);
        IJ.log("  ");
        IJ.log("Input watershed image : "+img.getTitle());
        IJ.log("Reading data from " + dir);
        IJ.log("Load list objects "+" Regions.zip"+ "  and Nuclei.zip from "+dir);
        
        popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        popRegions.loadObjects(dir + "Regions.zip");
        popNuclei.loadObjects(dir + "Nuclei.zip");
        
        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        
        // init cells 
        IJ.log("Initialing...");
        initCells2();
        IJ.log("Nb cells: " + popCells.size());
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);// 0 for thomas´s cell seg; -1 for farsight
        getMaxCluster();
        String str = typeObserved[typeCellObs]+"_";
        if(typeStatis==0)
        {
            str +="Raw_Data";
            if(typeDis==0)
            {
                calculEuclidStatistic(typeCellObs, func);
            }
            if(typeDis==1)
            {
                calculCellDistanceStatistic(typeCellSource, typeCellObs, func);
            }    
        } 
        else if(typeStatis==1)
        {
            str +="Random_Organization";
            createRandomDistribution(typeCellObs);
            if(typeDis==0)
            {
                calculEuclidStatistic(typeCellObs, func);
            }
            if(typeDis==1)
            {
                calculCellDistanceStatistic(typeCellSource, typeCellObs, func);
            }  
            ImageHandler cellTypeWat = drawCellTypes(false, imgSeg);
            cellTypeWat.show("Cell_Type_"+str);
            IJ.selectWindow(cellTypeWat.getTitle());
            IJ.saveAs("Tiff", subdir+cellTypeWat.getTitle()+".tif");
            cellTypeWat.closeImagePlus();
        }
        else{
            str +="Raw_Data";
            IJ.log("Please noted that: ");
            IJ.log("For cluster statistic analysis, we do not provide spatial statistic for cell distance");
            calculClusterEuclidStatistic(typeCellObs, func);
        }
        
        IJ.log("Saved the output into directory: "+subdir);
        IJ.log("Finished");
        IJ.selectWindow("Log");
        IJ.saveAs("Text", subdir+ "log_"+typeObserved[typeCellObs]+"_"
                +typeStat[typeStatis]+"_"+typeDistance[typeDis]+".txt");
        IJ.log("====================================================================");
    }
    private ImageHandler drawCellTypes(boolean nuc, ImageHandler imgSeg) {
        ImageHandler draw = new ImageShort("TYPE_", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.type);
            } else {
                C.region.draw(draw, C.type);
            }
        }
        return draw;
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
    
    private void calculCellDistanceStatistic(byte sourceType, byte desType, String func)
    {
        String[] typeObserved = {"*None*","DELTA", "BETA","ALPHA"};
        getMaxCluster();
        int sizeV = popCells.size();
        ArrayList<Vertex> vertexes = new ArrayList<>(sizeV);
        ArrayList<DirectedEdgeC> edges = new ArrayList<>();
        HashMap<Integer, Integer> cacheId = null;
        cacheId = new HashMap<Integer, Integer>();
        int idx = 0;
        for (Cell C : popCells) 
        {
            if (C.type != UNLABELLED && C.connected) 
            {
                //Vertex vC = new Vertex(C.id - 1, "Cell" + C.id, C.type);
                cacheId.put(idx, C.id);
                Vertex vC = new Vertex(idx, "Cell" + C.id, C.type);
                vertexes.add(vC);
                idx++;
            }
        }
        IJ.log("Nb of vertex is: " + vertexes.size());
        for (Cell C : popCells) 
        {
            if (C.type != UNLABELLED && C.connected) 
            {
                for (Cell C1 : C.nei1) 
                {
                    if (C1.type != UNLABELLED && C1.connected) 
                    {
                        int keyC = (int) MapUtils.getKeyFromValue(cacheId, C.id);
                        int keyC1 = (int) MapUtils.getKeyFromValue(cacheId, C1.id);
                        double weight = 1.0;
                        Vertex source = null, des = null;
                        for(Vertex ve : vertexes)
                        {
                            if(ve.getId()==keyC)
                            {
                                source = ve;
                            }
                            if(ve.getId()==keyC1)
                            {
                                des = ve; 
                            }
                        }
                        if(source!=null && des!=null)
                        {
                            DirectedEdgeC e = new DirectedEdgeC(source, des, weight);  //problem here
                            edges.add(e);
                        }    
                        
                    }
                }
            }
        }
//        for (Cell C : popCells) 
//        {
//            for (Cell Cnei : C.nei1) 
//            {
//                double weight = 1.0;
//                Vertex source = null, des = null;
//                for(Vertex ve : vertexes)
//                {
//                    if(ve.getId()==C.id - 1)
//                    {
//                        source = ve;
//                    }
//                    if(ve.getId()==Cnei.id - 1)
//                    {
//                        des = ve; 
//                    }
//                }
//                if(source!=null && des!=null)
//                {
//                    DirectedEdgeC e = new DirectedEdgeC(source, des, weight);  //problem here
//                    edges.add(e);
//                }  
//            }    
//            
//        }    
        IJ.log("Nb of edges is: " + edges.size());

        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : vertexes) 
        {
            if (v1.getType() == sourceType)
            {
                sourceLst.add(v1);
            }
            if (v1.getType() == desType) 
            {
                desLst.add(v1);
            }
        }
        String str = typeObserved[sourceType]+"_"+typeObserved[desType];
        IJ.log("Size source: beta: " + sourceLst.size() + " desLst: delta " + desLst.size());
        int numPoints = 1000;
        int numRandomSamples = 100;
        double env = 0.05;
        CellSpatialAnalysis spa = new CellSpatialAnalysis(numPoints, numRandomSamples, env, subdir);
//        spa.process(sourceLst, desLst, vertexes, edges, str);
        spa.process(desLst, vertexes, edges, func);
    }  
    
    private void calculEuclidStatistic(byte typeObserved, String func) 
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
        int idx = 0;
        for (Cell C : popCells) 
        {
            if (C.type == typeObserved) 
            {
                idx++;
                Object3D nuc = C.nucleus;
                Point3D center = nuc.getCenterAsPoint();
                Voxel3D v = new Voxel3D(center, idx);
                ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
                arr.add(v);
                Object3DVoxels obj = new Object3DVoxels(arr);
                obj.setCalibration(calibration);
                obj.setName("Obj" + idx);
                pop.addObject(obj);
            }
            
        }
        IJ.log("Nb observed cells in the list: "+pop.getNbObjects());
        
        
        ArrayList<Point3D> lsCenterCells = new ArrayList<Point3D>();  //list total cells
        for (Cell cc : popCells) 
        {
            if(cc.type != UNLABELLED)
            {
                Object3D nuc = cc.nucleus;
                Point3D center = nuc.getCenterAsPoint();
                lsCenterCells.add(center);
                
            }   
        }
        
//        int numPoints = 1000;
        String des = "_";
        int numPoints = lsCenterCells.size();
        int numRandomSamples = 100;
        double distHardCore = 2;
        double env = 0.05;
        SpatialTissueAnalysis spa = new SpatialTissueAnalysis(numPoints, numRandomSamples, distHardCore, env);
        spa.process(pop, imgMask, lsCenterCells, func, true, true, true, subdir, des);
//        rs.setValue("F_SDI_Eu"+typeObserved, 0, out.getValue(0));
//        rs.setValue("G_SDI_Eu"+typeObserved, 0, out.getValue(1));
//        rs.show("Statistic Function");
        
    }
    private void calculClusterEuclidStatistic(byte typeObserved, String func) 
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
        spa.process(pop, imgMask, lsCenterCells, func, true, true, true, subdir, desc);
//        rs.setValue("F_SDI_Eu"+typeObserved, 0, out.getValue(0));
//        rs.setValue("G_SDI_Eu"+typeObserved, 0, out.getValue(1));
//        rs.show("Statistic Function");
        
    }
    private void createRandomDistribution(byte typeObserved)
    {
        Cells3DPopulation popTotal = initCellPopulation();
        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();  //all cells
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();  //observed cells
        ArrayList<Vertex> vertexes = popTotal.getVertexes();
        for (Vertex v1 : vertexes) 
        {
            if (v1.getType() == typeObserved) 
            {
                desLst.add(v1);
            }
            else
            {
                sourceLst.add(v1);
            }
            
        }
        int count = 0;
        while(desLst.size() > 0)
        {
            int rDes = (int) (Math.random() * desLst.size());
            Vertex d = desLst.get(rDes);
            int rSource = (int) (Math.random() * sourceLst.size());
            Vertex s = sourceLst.get(rSource);
            byte typeC = s.getType();
            s.setType(d.getType());
            d.setType(typeC);
            desLst.remove(d);
            sourceLst.remove(s);
            sourceLst.add(d);
            count++;
        }
        for (Vertex vi : vertexes) 
        {
            int cellId = cacheId.get(vi.getId());
            for(Cell c : popCells)
            {
               if(c.id == cellId)
               {
                   c.type = vi.getType();
               }    
            }    
        }
    }
    private Cells3DPopulation initCellPopulation()
    {
        int sizeV = popCells.size();
        ArrayList<Vertex> vertexes = new ArrayList<>(sizeV);
        ArrayList<DirectedEdgeC> edges = new ArrayList<>();
        
        cacheId = new HashMap<Integer, Integer>();
        int idx = 0;
        for (Cell C : popCells) 
        {
            if (C.type != UNLABELLED && C.connected) 
            {
                //Vertex vC = new Vertex(C.id - 1, "Cell" + C.id, C.type);
                cacheId.put(idx, C.id);
                Vertex vC = new Vertex(idx, "Cell" + C.id, C.type);
                vertexes.add(vC);
                idx++;
            }
        }
//        IJ.log("Nb of vertex is: " + vertexes.size());
        for (Cell C : popCells) 
        {
            if (C.type != UNLABELLED && C.connected) 
            {
                for (Cell C1 : C.nei1) 
                {
                    if (C1.type != UNLABELLED && C1.connected) 
                    {
                        int keyC = (int) MapUtils.getKeyFromValue(cacheId, C.id);
                        int keyC1 = (int) MapUtils.getKeyFromValue(cacheId, C1.id);
                        double weight = 1.0;
                        Vertex source = null, des = null;
                        for(Vertex ve : vertexes)
                        {
                            if(ve.getId()==keyC)
                            {
                                source = ve;
                            }
                            if(ve.getId()==keyC1)
                            {
                                des = ve; 
                            }
                        }
                        if(source!=null && des!=null)
                        {
                            DirectedEdgeC e = new DirectedEdgeC(source, des, weight);  //problem here
                            edges.add(e);
                        }    
                        
                    }
                }
            }
        }
//        IJ.log("Nb of edges is: " + edges.size());
        Cells3DPopulation pop = new Cells3DPopulation(vertexes, edges);
        pop.setMask(vertexes, edges);
        return pop;
        
    }
    public void getMaxCluster()
    {
        HashMap<Integer,ArrayList<Cell>> map = new HashMap<Integer,ArrayList<Cell>>();
        for (Cell Ci : popCells) 
        {
            Ci.marked = false;
            Ci.connected = false;
        }
        int val = 0, maxSize = 0; 
        for (Cell Cii : popCells) 
        {
            if((!Cii.marked) && (Cii.type!=UNLABELLED))
            {
                ArrayList<Cell> cls = BFS(Cii);
                
                if(cls.size() > 0)
                {
                    val++;
//                    IJ.log("Size of cluster "+val+" is: "+cls.size());
                    map.put(val, cls);
                    if(maxSize < cls.size())
                    {
                        maxSize = cls.size();
                    }    
                }
            }
             
        }
//        IJ.log("Size of max cluster is: "+ maxSize);
//        IJ.log("Drawing...");
        for (Map.Entry<Integer, ArrayList<Cell>> ee : map.entrySet()) 
        {
            int key = ee.getKey();
            ArrayList<Cell> lsCells = ee.getValue();
            if(maxSize != 0 && lsCells.size()==maxSize)
            {
                for (Cell C : lsCells) 
                {
                    C.connected = true;
                }
            }    
            
        }
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
    private ArrayList<Cell> BFS(Cell s, byte typeCellObserved) 
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
                if((!ne.marked) && (ne.type ==typeCellObserved))
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
    /**
     * Breadth First Search
     * @param s
     * @return 
     */
    private ArrayList<Cell> BFS(Cell s) 
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
                if((!ne.marked) && (ne.type != UNLABELLED))
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
}
