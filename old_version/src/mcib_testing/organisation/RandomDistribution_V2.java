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
import ij.measure.ResultsTable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;
import mcib_testing.Utils.Cells3DPopulation;
import mcib_testing.Utils.Objects3DPopulationGrid;
import stdlib.DirectedEdgeC;
import stdlib.MapUtils;
import stdlib.Queue;
import stdlib.Vertex;

/**
 *
 * @author tranhoa
 */
public class RandomDistribution_V2 implements ij.plugin.PlugIn
{
    private final int M3 = 3;
    private final int M2 = 2;
    private final int M1 = 1;
    private final int UNLABELLED = 0;
    private Objects3DPopulationGrid popRegions = null;
    private Objects3DPopulationGrid popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg=null, imgWat=null;
    ImageHandler img;

    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    public String dir = null, subdir = null;
    HashMap<Integer, Integer> cacheId = null;
    ResultsTable rtContact = null, rtLayer = null;
    int max_distance = 17;
    int[] sRes = null;
    public float taa = 0, tbb = 0, tdd = 0, tab = 0, tad = 0, tbd = 0, tsum = 0;
    public boolean show = false;
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
        String[] typeObserved = {"*None*","M1", "M2","M3"};
        byte typeCellObserved = (byte) M3;
        byte typeSource = (byte) M2;
        int nbRandoms = 100;
        GenericDialog gd = new GenericDialog("Random Distribution");
        gd.addMessage("3D Tissue Analysis Framework");
        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : watershed image, nb random models");
        gd.addMessage("Output: random distribution computation.");
        gd.addMessage(" ");
        gd.addChoice("Watershed_Image ", titles, titles[0]);
        gd.addChoice("Observed_Cell_Type: ", typeObserved, typeObserved[typeCellObserved]);
        gd.addChoice("Source_Cell_Type_Observed : ", typeObserved, typeObserved[typeSource]);
        gd.addNumericField("Nb_Randoms : ", nbRandoms,  0, 10, "");
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        int idxImg = gd.getNextChoiceIndex();
        typeCellObserved = (byte)gd.getNextChoiceIndex();
        typeSource = (byte)gd.getNextChoiceIndex();
        nbRandoms = (int) gd.getNextNumber();
        ImagePlus plus;
        plus = WindowManager.getImage(wList[idxImg]);
        img = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        File fl = new File(dir+"100_random_distributions");
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = dir + "100_random_distributions/";
            } else {
                subdir = dir;
            }
        }
        else{
            subdir = dir + "100_random_distributions/";
        }
        
        IJ.log("Analysing ...");
        IJ.log("Reading data from " + dir);
        popRegions = new Objects3DPopulationGrid();
        popNuclei = new Objects3DPopulationGrid();
        popRegions.loadObjects(dir + "Regions.zip");
        popNuclei.loadObjects(dir + "Nuclei.zip");

        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
//        imgWat = img;
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        
        // init cells 
        initCells2();
        IJ.log("size is" + popCells.size());
        //drawCellTypes(false).show("TYPE BEFORE");
        IJ.log("Association1 ...");
//        computeAssoCells(imgWat, 0); 
        computeAssoCells(imgWat, 0); 
        
        rtContact = new ResultsTable();
        rtLayer = new ResultsTable();
        sRes = new int[max_distance];
        getMaxCluster();
        IJ.log("Type observed : "+typeObserved[typeCellObserved]);
        IJ.log("Type source : "+typeObserved[typeSource]);
        IJ.log("Create random distribution with cells :"+typeObserved[typeCellObserved]+" nb random organizations: "+nbRandoms);
        createRandomDistribution(typeCellObserved, typeSource, nbRandoms);
        IJ.log("Finished");
    }
    private void createRandomDistribution(byte typeObserved, byte typeSource, int nbRandoms)
    {
        Cells3DPopulation popTotal = initCellPopulation();
        computeLayerRawData(typeSource, typeObserved, popTotal);  //verify source??
        for(int i=0; i< sRes.length; i++)
        {
            sRes[i] = 0;
        }
//        int nbRandoms = 100;
//        IJ.log("No problem here");
//        String str = MapUtils.toString(cacheId);
//        IJ.log("cache  :"+str);
        for(int k = 0 ; k < nbRandoms; k++)
        {
            IJ.log("Simulated Random Organization " + (k+1));
            ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();  //all cells
            ArrayList<Vertex> desLst = new ArrayList<Vertex>();  //observed cells
            ArrayList<Vertex> vertexes = popTotal.getVertexes();
            for (Vertex v1 : vertexes) 
            {
                if (v1.getType() != typeObserved)
                {
                    sourceLst.add(v1);
                }
                if (v1.getType() == typeObserved) 
                {
                    desLst.add(v1);
                }
            }
            int sizeObserved = desLst.size();
//            IJ.log("source: "+sourceLst.size()+"  des: "+desLst.size());
//            ArrayList<Vertex> popObserved = new ArrayList<Vertex>();
            
//            popObserved.add(s);
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
            computeLayer(typeSource, typeObserved, popTotal, k);
            computeAllContacts(k);
        }
        rtLayer.incrementCounter();
        for(int i=1; i< sRes.length; i++)
        {
            rtLayer.setValue("d_"+i, nbRandoms, sRes[i]);
        }
        rtLayer.incrementCounter();
        for(int i=1; i< sRes.length; i++)
        {
            rtLayer.setValue("d_"+i, nbRandoms+1, (double)sRes[i]/nbRandoms);
        }
        rtContact.incrementCounter();
        rtContact.setValue("aa", nbRandoms, taa);
        rtContact.setValue("bb", nbRandoms, tbb);
        rtContact.setValue("dd", nbRandoms, tdd);
        rtContact.setValue("ab", nbRandoms, tab);
        rtContact.setValue("ad", nbRandoms, tad);
        rtContact.setValue("bd", nbRandoms, tbd);
        rtContact.setValue("all", nbRandoms, tsum);
        rtContact.incrementCounter();
        rtContact.setValue("aa", nbRandoms+1, taa/tsum);
        rtContact.setValue("bb", nbRandoms+1, tbb/tsum);
        rtContact.setValue("dd", nbRandoms+1, tdd/tsum);
        rtContact.setValue("ab", nbRandoms+1, tab/tsum);
        rtContact.setValue("ad", nbRandoms+1, tad/tsum);
        rtContact.setValue("bd", nbRandoms+1, tbd/tsum);
        rtContact.setValue("all", nbRandoms+1, tsum/tsum);
        if(show)
        {
            rtLayer.show("LayerDistance_"+nbRandoms+"_randoms.csv");
            rtContact.show("AllContacts_"+nbRandoms+"_randoms.csv");
        }    
        
        try {
            rtLayer.saveAs(subdir+"LayerDistance_"+nbRandoms+"_randoms.csv");
            rtContact.saveAs(subdir+"AllContacts_"+nbRandoms+"_randoms.csv");
            IJ.log("Saving the output table "+ "LayerDistance_"+nbRandoms+"_randoms.csv into the folder: "+subdir);
            IJ.log("Saving the output table "+ "AllContacts_"+nbRandoms+"_randoms.csv into the folder: "+subdir);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
    }
    private void computeLayer(byte sourceType, byte destinationType, Cells3DPopulation popTotal, int st)
    {
        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : popTotal.getVertexes()) 
        {
            if (v1.getType() == sourceType)
            {
                sourceLst.add(v1);
            }
            if (v1.getType() == destinationType) 
            {
                desLst.add(v1);
            }
        }
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(popTotal.getVertexes());
        pop.setMask(popTotal.getVertexes(), popTotal.getEdges());
        int maxDistanceLayer = pop.computeLayerDistance(sourceLst, desLst);
        
        int count = 0;
        int distanceLength = maxDistanceLayer+1;
        int[] res = new int[distanceLength];
        
        rtLayer.incrementCounter();
        for (Vertex v : sourceLst) 
        {
            if(v.getDistanceLayer()!=-1)
            {
                res[(int)v.getDistanceLayer()]++;
            }    
            
        }
        for(int j=1; j< res.length; j++)
        {
            rtLayer.setValue("d_"+j, st, res[j]);
            if(j < max_distance)
            {
                sRes[j] += res[j];
            }    
            
        }
        
    }
    private void computeLayerRawData(byte sourceType, byte destinationType, Cells3DPopulation popTotal)
    {
        ResultsTable rtLayerOrigin = new ResultsTable();
        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : popTotal.getVertexes()) 
        {
            if (v1.getType() == sourceType)
            {
                sourceLst.add(v1);
            }
            if (v1.getType() == destinationType) 
            {
                desLst.add(v1);
            }
        }
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(popTotal.getVertexes());
        pop.setMask(popTotal.getVertexes(), popTotal.getEdges());
        int maxDistanceLayer = pop.computeLayerDistance(sourceLst, desLst);
        
        int count = 0;
        int distanceLength = maxDistanceLayer+1;
        int[] res = new int[distanceLength];
        
        rtLayerOrigin.incrementCounter();
        for (Vertex v : sourceLst) 
        {
            if(v.getDistanceLayer()!=-1)
            {
                
                res[(int)v.getDistanceLayer()]++;
            }    
            
        }
        for(int j=1; j< res.length; j++)
        {
            rtLayerOrigin.setValue("d_"+j, 0, res[j]);
        }
//        rtLayerOrigin.show("Observed Layer");
        try {
            rtLayerOrigin.saveAs(subdir+"Layer_Observed_Raw_Data.csv");
            IJ.log("Saving the output table Layer_Observed_Raw_Data.csv in the folder: "+subdir);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
        
    }
    private void computeAllContacts(int st)
    {
//        rtContact = new ResultsTable();
//        float taa = 0, tbb = 0, tdd = 0, tab = 0, tad = 0, tbd = 0, tsum = 0;
        float aa = 0, bb = 0, dd = 0, ab = 0, ad = 0, bd = 0;
        float na = 0, nb = 0, nd = 0;
        int nbUnlabelled = 0;
        for (Cell C : popCells) 
        {
            if (C.type == M3) 
            {
                na++;
                int[] ne = C.computeNei1TypeHisto();
                aa += ne[M3-1];
                ab += ne[M2-1];
                ad += ne[M1-1];
            } else if (C.type == M2) {
                nb++;
                int[] ne = C.computeNei1TypeHisto();
                ab += ne[M3-1];
                bb += ne[M2-1];
                bd += ne[M1-1];
            } else if (C.type == M1) {
                nd++;
                int[] ne = C.computeNei1TypeHisto();
                ad += ne[M3-1];
                bd += ne[M2-1];
                dd += ne[M1-1];
            }
            else{
                nbUnlabelled++;
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
        float pa = 0, pb = 0, pd = 0;
        pa = na / sumnb;
        pb = nb / sumnb;
        pd = nd / sumnb;
        rtContact.incrementCounter();
        
        rtContact.setValue("total_cells", st, popCells.size());
        rtContact.setValue("nb_unlabelled_cells", st, nbUnlabelled);
        rtContact.setValue("na", st, na);
        rtContact.setValue("nb", st, nb);
        rtContact.setValue("nd", st, nd);
        rtContact.setValue("totalabd", st, sumnb);
        
        rtContact.setValue("pa", st, pa);
        rtContact.setValue("pb", st, pb);
        rtContact.setValue("pd", st, pd);
        
        rtContact.setValue("aa", st, aa);
        rtContact.setValue("bb", st, bb);
        rtContact.setValue("dd", st, dd);
        rtContact.setValue("ab", st, ab);
        rtContact.setValue("ad", st, ad);
        rtContact.setValue("bd", st, bd);
        rtContact.setValue("all_contacts", st, sum);
        
        taa += aa;
        tbb += bb;
        tdd += dd;
        tab += ab;
        tad += ad;
        tbd += bd;
        tsum += sum;
        rtContact.setValue("paa", st, aa / sum);
        rtContact.setValue("pbb", st, bb / sum);
        rtContact.setValue("pdd", st, dd / sum);
        rtContact.setValue("pab", st, ab / sum);
        rtContact.setValue("pad", st, ad / sum);
        rtContact.setValue("pbd", st, bd / sum);
//        rtContact.show("Observed_raw_data");
//        try {
//            rtContact.saveAs(dir+"AllContact.csv");
//            IJ.log("Saving the output table AllContact.csv in the folder: "+dir);
//        } catch (IOException ex) {
//            IJ.log("Have the problem of saving data");
//        }
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
        IJ.log("Nb of vertexs is: " + vertexes.size());
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
        IJ.log("Nb of edges is: " + edges.size());
        Cells3DPopulation pop = new Cells3DPopulation(vertexes, edges);
        pop.setMask(vertexes, edges);
        return pop;
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
        ImageHandler drawBFS = new ImageShort("BFS", img.sizeX, img.sizeY, img.sizeZ);
        for (Map.Entry<Integer, ArrayList<Cell>> ee : map.entrySet()) 
        {
            int key = ee.getKey();
            ArrayList<Cell> lsCells = ee.getValue();
            
            if(maxSize != 0 && lsCells.size()==maxSize)
            {
                for (Cell C : lsCells) 
                {
                    C.connected = true;
                    if(C.currLevel == -1)
                    {
                        IJ.log("have the problem");
                    }
                    else{
                        C.region.draw(drawBFS, C.currLevel);
                    }
                    
                }
            }    
            
        }
//        drawBFS.show();
    }
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
}
