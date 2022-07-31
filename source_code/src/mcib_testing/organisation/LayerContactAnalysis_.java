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
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;
import mcib_testing.Utils.Cells3DPopulation;
import stdlib.DirectedEdgeC;
import stdlib.MapUtils;
import stdlib.Queue;
import stdlib.Vertex;

/**
 *
 * @author tranhoa
 */
public class LayerContactAnalysis_ implements ij.plugin.PlugIn
{
//    private final int ALPHA = 3;
//    private final int BETA = 2;
//    private final int DELTA = 1;
    private final int UNLABELLED = 0;
//    private final int CLUSTER = 1;
//    private final int UNIFORM = 2;
//    private final int RANDOM = 3;
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg=null, imgWat=null;
    ImageHandler img, dapi;

    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    public String subdir = null;
    public boolean contactCal = true, neighborCal=true, hisContact=true, cellsNetwork=true;
    public boolean show = false;
    String[] typeObserved = {"DELTA", "BETA","ALPHA"};
    public byte targetCell = 1, sourceCell = 0;
    public String input_dir = "";
    public void run(String arg) 
    {
        
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open");
            return;
        }
        String[] titles = new String[wList.length];
        if (wList.length<2) {
            IJ.error("Open at least 2 images: watershed image and nuclei dapi image");
            return;
        }
        ImagePlus temp = null;
        int count=0;
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
            if (null !=imp){
                count++;
                if(count==1){
                    temp = imp;
                    WindowManager.setTempCurrentImage(temp);
                    input_dir = IJ.getDirectory("image");
                }
                
            }    
        }
        int cellId = 0;
        GenericDialogPlus gd = new GenericDialogPlus("Layer Contact Analysis");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("3D Tissue Analysis Framework");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans");
//        gd.addMessage("Input : watershed image");
//        gd.addMessage("Output: cell-cell contact.");
        gd.addMessage(" ");
//        gd.addStringField("Input Dir: ", input_dir);
        gd.addDirectoryField("Input_Dir: ", input_dir, 20);
        gd.addChoice("Watershed_Image ", titles, titles[0]);
        gd.addChoice("DAPI_Raw_Image :", titles, titles[1]);
        gd.addChoice("Target_Cell : ", typeObserved, typeObserved[targetCell]);
        gd.addChoice("Source_Cell : ", typeObserved, typeObserved[sourceCell]);
//        gd.addNumericField("Cell Id", cellId, 0);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        input_dir = gd.getNextString();
        if ("".equals(input_dir)) 
        {
                return;
        }
        
        int[] index = new int[2];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        targetCell = (byte)(gd.getNextChoiceIndex()+1);
        sourceCell = (byte)(gd.getNextChoiceIndex()+1);
//        cellId = (int)gd.getNextNumber();
        if (targetCell==sourceCell) {
            IJ.error("Error, target cell is the same source cell, please choose again.");
            return;
        }
        ImagePlus plus;
        plus = WindowManager.getImage(wList[index[0]]);
        img = ImageInt.wrap(plus);
//        WindowManager.setTempCurrentImage(plus);
        ImagePlus plus1;
        plus1 = WindowManager.getImage(wList[index[1]]);
        dapi = ImageInt.wrap(plus1);
//        dir = IJ.getDirectory("image");
        if (!input_dir.endsWith("/"))
            input_dir = input_dir + "/";
        IJ.log("Noted: input dir: where you keep list objects Regions.zip and Nuclei.zip");
        IJ.log("Noted: need raw dapi image as input to capture the real micrometer object scale");
        
        File wdir = new File(input_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        File fl = new File(input_dir+"layer_distance");
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = input_dir + "layer_distance/";
                IJ.log("Create a layer_distance folder in the directory: "+input_dir);
            } else {
                subdir = input_dir;
            }
        }
        else{
            subdir = input_dir + "layer_distance/";
        }
   
        IJ.log("====================================================================");
        IJ.log("Input watershed image : "+img.getTitle());
        IJ.log("Target cell: "+typeObserved[targetCell-1]+"  source cell observed: "+typeObserved[sourceCell-1]);
        IJ.log("Reading data from " + input_dir);
        IJ.log("Load list objects "+" Regions.zip"+ " and Nuclei.zip from "+input_dir);
//        IJ.log("cell id: "+cellId);
        popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        popRegions.loadObjects(input_dir + "Regions.zip");
        popNuclei.loadObjects(input_dir + "Nuclei.zip");
        
        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
//        imgWat = img;
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        
        // init cells 
        IJ.log("Initializing...");
        initCells2();
        IJ.log("size is" + popCells.size());
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);// 0 for thomas´s cell seg; -1 for farsight
//        showCellDistance(cellId);
//        getMaxCluster();
        calculHistogramContact(targetCell, sourceCell);
        
        IJ.log("Finished");
        IJ.log("====================================================================");
        
    }
    
    
    /**
     * breadth-first search from a single source
     * @param s
     * @return 
     */
    private ArrayList<Cell> BFS(Cell s) 
    {
        ArrayList<Cell> cls = new ArrayList<Cell>();
        Queue<Cell> q = new Queue<Cell>();
        s.marked = true;
        q.enqueue(s);
        cls.add(s);
        
        while (!q.isEmpty()) {
            Cell v = q.dequeue();
            for(Cell ne : v.nei1)
            {
                if((!ne.marked) && (ne.type != UNLABELLED))
                {
                    ne.marked = true;
                    cls.add(ne);
                    q.enqueue(ne);
                }    
            }  
        }
        return cls;
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
        ImageHandler drawConnected = new ImageShort("connected_component", img.sizeX, img.sizeY, img.sizeZ);
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
                    C.region.draw(drawConnected, C.id);
                }
                IJ.log("Size of max cluster is: "+lsCells.size());
            }    
            
        }
//        drawConnected.show();
    }
    private void showCellDistance(int cellId)
    {
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

        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : vertexes) 
        {
            int vexId = (int)MapUtils.getKeyFromValue(cacheId, cellId);
            if(v1.getId()==vexId)
            {
                desLst.add(v1);
            }
            else{
                sourceLst.add(v1);
            }
        }
        
//        for (Vertex v1 : vertexes) 
//        {
//            if(v1.getId()==(int)MapUtils.getKeyFromValue(cacheId, 475))
//            {
//                desLst.add(v1);
//            }    
//            else{
//                sourceLst.add(v1);
//            }
//        }
        
//        IJ.log("Size source list: " + sourceLst.size() + " size destionation list: " + desLst.size());
        //String str = typeObserved[sourceType-1]+"_"+typeObserved[desType-1];
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(vertexes);
        pop.setMask(vertexes, edges);
        ArrayUtil arr = pop.computeMinDistance(sourceLst, desLst);
//        spa.process(desLst, vertexes, edges, "G");
        
        ImageHandler drawLayer = new ImageShort("LAYER_DISTANCE_"+cellId, img.sizeX, img.sizeY, img.sizeZ);
        for (Vertex vd : desLst) 
        {
            int cId = cacheId.get(vd.getId());
            for(Cell c : popCells)
            {
                if(c.id==cId)
                {
                    c.region.draw(drawLayer, (int)arr.getMaximum()+2);
                }    
            }
        }
        for (Vertex vs : sourceLst) 
        {
            int cId = cacheId.get(vs.getId());
            if(vs.getShortestDistance()!=-1)
            {
                for(Cell c : popCells)
                {
                    if(c.id==cId)
                    {
                        c.region.draw(drawLayer, (int)vs.getShortestDistance()+1);
                    }    
                }
            }    
        }
        drawLayer.setScale(dapi);
        drawLayer.show();
//        IJ.selectWindow(drawLayer.getTitle());
//        IJ.saveAs(drawLayer.getImagePlus(), "Tiff", subdir+drawLayer.getTitle()+".tif");
//        IJ.log("Saving the output image "+drawLayer.getTitle() +".tif into the folder: "+subdir);
    }        
    private void calculHistogramContact(byte sourceType, byte desType)
    {
        
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
        
//        IJ.log("Size source list: " + sourceLst.size() + " size destionation list: " + desLst.size());
        String str = typeObserved[sourceType-1]+"_"+typeObserved[desType-1];
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(vertexes);
        pop.setMask(vertexes, edges);
        ArrayUtil arr = pop.computeMinDistance(sourceLst, desLst);
//        spa.process(desLst, vertexes, edges, "G");
        
        ImageHandler drawLayer = new ImageShort("LAYER_DISTANCE_"+str, img.sizeX, img.sizeY, img.sizeZ);
        
        for (Vertex vd : desLst) 
        {
            int cId = cacheId.get(vd.getId());
            for(Cell c : popCells)
            {
                if(c.id==cId)
                {
                    c.region.draw(drawLayer, 1);
                }    
            }
        }
        for (Vertex vs : sourceLst) 
        {
            int cId = cacheId.get(vs.getId());
            if(vs.getShortestDistance()!=-1)
            {
                for(Cell c : popCells)
                {
                    if(c.id==cId)
                    {
                        c.region.draw(drawLayer, (int)vs.getShortestDistance()+1);
                    }    
                }
            }    
        }
        drawLayer.setScale(dapi);
        drawLayer.show();
        IJ.selectWindow(drawLayer.getTitle());
        IJ.saveAs(drawLayer.getImagePlus(), "Tiff", subdir+drawLayer.getTitle()+".tif");
        IJ.log("Saving the output image "+drawLayer.getTitle() +".tif into the folder: "+subdir);
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
    
    private void initCells2() 
    {
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
    
    
}
