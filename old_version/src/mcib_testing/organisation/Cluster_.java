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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
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
import mcib_testing.Utils.CellSpatialAnalysis;
import mcib_testing.Utils.Cells3DPopulation;
import mcib_testing.Utils.SpatialTissueAnalysis;
import stdlib.CentroidC;
import stdlib.DirectedEdgeC;
import stdlib.MapUtils;
import stdlib.Queue;
import stdlib.Vertex;

/**
 *
 * @author tranhoa
 */
public class Cluster_ implements ij.plugin.PlugIn
{
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int UNLABELLED = 0;
    private final int CLUSTER = 1;
    private final int UNIFORM = 2;
    private final int RANDOM = 3;
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler imgSeg=null, imgWat=null;
    ImageHandler img;

    ArrayList<Cell> popCells = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    HashMap<Integer, Integer> cacheId = null;
    String dir = null;
    ResultsTable rs = null;
    ResultsTable rtLayer = null;
    ResultsTable rtContact = null;
    public float taa = 0, tbb = 0, tdd = 0, tab = 0, tad = 0, tbd = 0, tsum = 0;
    int dis = 17;
    int[] sRes = null;
    public void run(String arg) 
    {
        int nbima = WindowManager.getImageCount();

        if (nbima < 1) {
            IJ.showMessage("No image opened !");
            return;
        }
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        img = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        
        
        byte typeCellObserved = ALPHA;
        byte typeSource = (byte) BETA;
        int distributionType = CLUSTER;
        double ratioDis = 0.25;
        String[] typeObserved = {"", "DELTA", "BETA","ALPHA"};
        String[] type_distribution = {"","CLUSTER", "UNIFORM", "RANDOM"};
        GenericDialog gd = new GenericDialog("Spatial Distribution");
        gd.addNumericField("Ratio Distribution: ", ratioDis, 1);
        gd.addChoice("Type of Distribution: ", type_distribution, type_distribution[distributionType]);
        gd.addChoice("Type Cell Observed : ", typeObserved, typeObserved[typeCellObserved]);
        gd.addChoice("Type Source : ", typeObserved, typeObserved[typeSource]);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        ratioDis = (double)gd.getNextNumber();
        distributionType = gd.getNextChoiceIndex();
        typeCellObserved = (byte)gd.getNextChoiceIndex();
        typeSource = (byte)gd.getNextChoiceIndex();
        IJ.log("Rad dis: "+ratioDis+" kind of distribution: "
                + type_distribution[distributionType] +"  type cell observed: "+typeObserved[typeCellObserved]
                +"  type source: "+typeObserved[typeSource]);
        
        
        IJ.log("Analysing ...");
        IJ.log("Reading data from " + dir);
        popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        popRegions.loadObjects(dir + "Regions.zip");
        popNuclei.loadObjects(dir + "Nuclei.zip");
        
        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        
        // init cells 
        initCells2();
        IJ.log("size is" + popCells.size());
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);  //replace by 1
//        test(4);
        

//        getNbCluster(typeCellObserved);
        
        
        
        rs = new ResultsTable();
        rs.incrementCounter();
//        createDistribution(typeCellObserved, distributionType, (float)ratioDis, typeSource);
//        computeAllContacts("AFTER");
        
        calculStatistic_(typeCellObserved);
//        calculEuclidStatistic(typeCellObserved);
//        calculStatistic_();
        
        
//        drawCellTypes(true).show("TYPE_NUCLEI");
//        drawCellTypes(false).show("TYPE_WAT");
        
        
//        showDistanceMap(typeCellObserved, distributionType, (float)ratioDis, typeSource);

        IJ.log("Finished");
    }
    
     
    private void computeAssoCells(ImageHandler img, int borderValue) {

        ArrayList<Integer>[] nei = computeAsso(img, borderValue);
//        int nbVoxSurf = 4;
        // index refers to region label
        for (int i = 0; i < nei.length; i++) 
        {
            if (nei[i].isEmpty()) {
                continue;
            }
            Cell C = region2Cell.get(i);
            if (C != null) {
                Object3D reg = C.region;
                C.nei1 = new ArrayList<Cell>();
                for (int j : nei[i]) {
                    Cell C1 = region2Cell.get(j);
                    if ((C1 != null) && (C1 != C)) 
                    {
//                        int[] surf = reg.surfaceContact(C1.region, Math.sqrt(3));
                        C.nei1.add(C1);
                    }
                }
            }
        }
    }
    public void test(double microDistance)
    {
//        ObjectCreator3D draw = new ObjectCreator3D(imgWat.sizeX, imgWat.sizeY, imgWat.sizeZ);
        Cell c = popCells.get(10);
        ImageHandler label = new ImageShort("reg_", imgWat.sizeX, imgWat.sizeY, imgWat.sizeZ);
        Object3D reg = c.region;
        int regVal = reg.getValue();
        reg.draw(label, regVal);
        label.show();
//        ArrayList<Voxel3D> arr = reg.getContours();
//        reg.drawContours(draw, reg.getValue());
//        draw.getImageHandler().show("cell10_contours");
        
        ImageByte imgDilate = dilateRegion(microDistance, reg);
        imgDilate.substractMask((ImageInt)label);
        imgDilate.show("substract_img");
        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
        for(int x=0; x < imgDilate.sizeX;x++)
        {
            for(int y=0; y<imgDilate.sizeY;y++)
            {
                for(int z=0; z<imgDilate.sizeZ;z++)
                {
                    if(imgDilate.getPixel(x, y, z) > 0)
                    {
                        Voxel3D v = new Voxel3D(x, y, z, regVal);
                        arr.add(v);
                    }
                }    
                    
            }    
        } 
        Object3DVoxels oo = new Object3DVoxels(arr);
        ArrayUtil ac = oo.listValues(imgWat); 
        IJ.log("nei : "+ac.toString());
        countNbVoxContact(ac);
//        ArrayUtil ac1 = ac.distinctValues();
//        IJ.log("nei distinct :  "+ac1.toString());
//        for(int i=0;i<ac1.getSize();i++)
//        {
//            
//        }
          
    } 
    public void countNbVoxContact(ArrayUtil ac) {
        ac.sort();
        HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
        ArrayUtil V = new ArrayUtil(ac.getSize());
        ArrayUtil V1 = new ArrayUtil(ac.getSize());
        int s = 0;
        double tmp = ac.getValue(0);
        V.addValue(0, tmp);
        s++;
        int p = 1;
        int c = 1;
        int si = ac.getSize();
        while (p < si) {
            while ((p < si) && (ac.getValue(p) == tmp)) {
                p++;
                c++;
            }
            map.put((int)tmp, c);
            V1.addValue(s-1, c);
            c = 0;
            if (p < si) {
                tmp = ac.getValue(p);
                V.addValue(s, tmp);
                s++;
                p++;
                c++;
            }
        }
        IJ.log("vals : "+V.getSubTabUtil(0, s).toString());
        IJ.log("count : "+V1.getSubTabUtil(0, s).toString());
        IJ.log(MapUtils.toString(map));
        
    }
    public ImageByte dilateRegion(double microDistance, Object3D region)
    {
        ImageHandler label = new ImageShort("dilate_", img.sizeX, img.sizeY, img.sizeZ);
        region.draw(label, region.getValue());
        try {
            int nbCPUs = ThreadUtil.getNbCpus();
            float rDilate = (float) microDistance;
            float radiusZ = (float) rDilate/2;
            ImageFloat edm = EDT.run(label, 0, 1, rDilate / radiusZ, true, nbCPUs);
            ImageByte temp = edm.threshold(rDilate, true, false);
            edm.flush();
            edm = null;
            temp.setOffset(label);
            temp.setScale(label);
            for(int xy=0; xy < temp.sizeXY;xy++)
            {
                for(int z=0; z<temp.sizeZ;z++)
                {
                    if(temp.getPixel(xy, z) > 0)
                    {
                        temp.setPixel(xy, z, region.getValue());
                    }    
                }    
            }    
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
    
         
    
    private Cell getCenterShellCore()
    {
        ArrayList<Cell> popS = (ArrayList<Cell>) popCells.clone();
        IJ.log("size of pop centroid before : "+ popS.size());
        Iterator<Cell> it = popS.iterator();
        int maxLabelledDegree = 0;
        Cell tmp = null;
        while(it.hasNext())
        {
            Cell Ci = it.next();
            int cLabelled = 0, cUnLabelled = 0;
            //double ratio = 0;
            if (Ci.type != UNLABELLED) 
            {
                for (Cell C1 : Ci.nei1) 
                {
                    if ((C1.type > 0)) 
                    {
                        cLabelled++;
                    } else {
                        cUnLabelled++;
                    }
                }
            }
//            Ci.degreeLabelled = cLabelled;
//            Ci.degreeUnlabelled = cUnLabelled;
            if(Ci.type != UNLABELLED && maxLabelledDegree < cLabelled)
            {
                maxLabelledDegree = cLabelled;
                tmp = Ci;
            }    
        }
//        IJ.log("size of pop cell at the shell core : "+ popS.size());
//        IJ.log("max degree labelled: "+maxLabelledDegree);
        return tmp;
    }
    private ArrayList<Cell> getListShellCore(byte typeObserved)
    {
        ArrayList<Cell> popS = (ArrayList<Cell>) popCells.clone();
//        IJ.log("size of pop centroid before : "+ popS.size());
        Iterator<Cell> it = popS.iterator();
        int maxLabelledDegree = 0;
        while(it.hasNext())
        {
            Cell Ci = it.next();
            int cLabelled = 0, cUnLabelled = 0;
            //double ratio = 0;
            if (Ci.type > 0) 
            {
                for (Cell C1 : Ci.nei1) 
                {
                    if ((C1.type > 0)) 
                    {
                        cLabelled++;
                    } else {
                        cUnLabelled++;
                    }
                }
            }
            Ci.degreeLabelled = cLabelled;
            Ci.degreeUnlabelled = cUnLabelled;
            if(maxLabelledDegree < cLabelled)
            {
                maxLabelledDegree = cLabelled;
            }    
            if (Ci.degreeLabelled < 7) 
            {
                it.remove();
            }
        }
        IJ.log("size of pop cell at the shell core : "+ popS.size());
        IJ.log("max degree labelled: "+maxLabelledDegree);
        return popS;
    }
    private ArrayList<CentroidC> getCentroid(ArrayList<Cell> popInputCells, byte typeObserved)
    {
        //int lsSize = popInputCells.size();
        ArrayList<CentroidC> lsCells = new ArrayList<CentroidC>();
        for (Cell cc : popInputCells) 
        {
            if(cc.type != UNLABELLED && cc.type != typeObserved && cc.connected)
            {
                Object3D nuc = cc.nucleus;
                Point3D center = nuc.getCenterAsPoint();
                CentroidC cent = new CentroidC(cc.id, cc.type, center);
                lsCells.add(cent);
            }    
            
        }
        return lsCells;
    }
    
    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", img.sizeX, img.sizeY, img.sizeZ);

        for (Cell C : popCells) 
        {
            if (nuc) 
            {
                C.nucleus.draw(draw, C.type);
            } else 
            {
                C.region.draw(draw, C.type);
            }
        }

        return draw;
    }
//    private ArrayList<Cell> getListShellCoreV1(byte typeObserved)
//    {
//        ArrayList<Cell> popS = (ArrayList<Cell>) popCells.clone();
//        IJ.log("size of pop centroid before : "+ popS.size());
//        Iterator<Cell> it = popS.iterator();
//        int maxLabelledDegree = 0;
//        while(it.hasNext())
//        {
//            Cell Ci = it.next();
//            int cLabelled = 0, cUnLabelled = 0;
//            //double ratio = 0;
//            if (Ci.type > 0) 
//            {
//                for (Cell C1 : Ci.nei1) 
//                {
//                    if ((C1.type > 0)) 
//                    {
//                        cLabelled++;
//                    } else {
//                        cUnLabelled++;
//                    }
//                }
//            }
//            Ci.degreeLabelled = cLabelled;
//            Ci.degreeUnlabelled = cUnLabelled;
//            if(maxLabelledDegree < cLabelled)
//            {
//                maxLabelledDegree = cLabelled;
//            }    
//            if (Ci.degreeLabelled < 8 || Ci.type==typeObserved) 
//            {
//                it.remove();
//            }
//        }
//        IJ.log("size of pop cell at the shell core : "+ popS.size());
//        IJ.log("max degree labelled: "+maxLabelledDegree);
//        return popS;
//    }
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
        IJ.log("Nb of edges is: " + edges.size());
        Cells3DPopulation pop = new Cells3DPopulation(vertexes, edges);
        pop.setMask(vertexes, edges);
        return pop;
        
    }    
    private void calculStatistic_(byte typeObserved)
    {
        Cells3DPopulation popTotal = initCellPopulation();
//        ArrayList<Vertex> vertexes = popTotal.getVertexes();
//        computeLayer(typeSource, typeObserved, popTotal);
        calculCellDistanceStatistic(popTotal, typeObserved);
//        calculEuclidStatistic(typeObserved);
//        drawGraph();
//        createAdjMatrix();
    }    
    private void calculStatistic_()
    {
        Cells3DPopulation popTotal = initCellPopulation();
        ArrayList<Vertex> vertexes = popTotal.getVertexes();
        
        byte typeSource = BETA, typeObserved = ALPHA;
        IJ.log("Observed cell: "+typeObserved+" ALPHA"+ "compute distance to: "+typeSource+" BETA");
//        computeLayer(typeSource, typeObserved, popTotal);
//        calculCellDistanceStatistic(popTotal, typeObserved);
        calculEuclidStatistic(typeObserved);
        
        typeObserved = DELTA;
        IJ.log("Observed cell: "+typeObserved+" DELTA"+ "compute distance to: "+typeSource+" BETA");
//        computeLayer(typeSource, typeObserved, popTotal);
//        calculCellDistanceStatistic(popTotal, typeObserved);
        calculEuclidStatistic(typeObserved);
    } 
    
    private void createRandomDistribution(byte typeObserved, byte typeSource, int nbRandoms)
    {
        Cells3DPopulation popTotal = initCellPopulation();
        computeLayerV2(typeSource, typeObserved, popTotal);
        for(int i=1; i< sRes.length; i++)
        {
            sRes[i] = 0;
        }
//        int nbRandoms = 100;
        for(int k = 0 ; k < nbRandoms; k++)
        {
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
//            computeAllContacts(k);
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
//        rtLayer.show("Layer Distance_"+nbRandoms+"_times");
//        rtContact.show("AllContact_"+nbRandoms+"_times");
        try {
            rtLayer.saveAs(dir+"LayerDistance.csv");
            rtContact.saveAs(dir+"AllContact.csv");
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
    }        
    private void createDistribution(byte typeObserved, int typeDistribution, float ratioDis, byte typeSource)
    {
        Cells3DPopulation popTotal = initCellPopulation();
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
//        IJ.log("Size source: " + sourceLst.size() + " desLst: " + desLst.size());
        ArrayList<Vertex> popObserved = new ArrayList<Vertex>();
        
        //init, choose a first cell to swap type cell
        
//        int ran = (int) (Math.random() * sourceLst.size());
//        Vertex s = sourceLst.get(ran);
        ArrayList<Cell> lsCore = getListShellCore(typeObserved);
        int ran = (int) (Math.random() * lsCore.size());
        Cell cTmp = lsCore.get(ran);
        int keyC = (int) MapUtils.getKeyFromValue(cacheId, cTmp.id);
        Vertex s = null;
        for(Vertex ve : sourceLst)
        {
            if(ve.getId()==keyC)
            {
                s = ve;
            }
        }
        if(s==null)
        {
            int r = (int) (Math.random() * sourceLst.size());
            s = sourceLst.get(r);
        }    
        int ranInit = (int) (Math.random() * desLst.size());
        Vertex d = desLst.get(ranInit);
        byte typeC = s.getType();
        s.setType(d.getType());
        d.setType(typeC);
        desLst.remove(d);
        popObserved.add(s);
        
        //compute the distance map and swap objs
        int count = 1;
        while(desLst.size() > 0)
        {
            ArrayList<Vertex> sourceTmp = new ArrayList<Vertex>();
            for (Vertex v1 : vertexes) 
            {
                if (v1.getType() != typeObserved)
                {
                    sourceTmp.add(v1);
                }
            }
            popTotal.computeShortestDistance(sourceTmp, popObserved);
            ArrayList<Vertex> lsTmp = normalizeDistanceMap(sourceTmp, ratioDis, typeDistribution);
            IJ.log("nb possible :" + lsTmp.size());
            if(lsTmp.size()>0)
            {
                int randomC;
                ArrayList<Vertex> arrTmp = new ArrayList<Vertex>();
                Vertex randomCell = null;
                if(typeDistribution==CLUSTER)
                {
                    Vertex first = lsTmp.get(0);
                    arrTmp.add(first);
                    for(Vertex vC : lsTmp)
                    {
                        if(vC.getShortestDistance()==first.getShortestDistance())
                        {
                            arrTmp.add(vC);
                        }    
                    }
                    randomC = (int) (Math.random() * arrTmp.size());
                    randomCell = arrTmp.get(randomC);
                }
                else if(typeDistribution==UNIFORM)
                {
                    Vertex last = lsTmp.get(lsTmp.size() - 1);
                    for(Vertex vC : lsTmp)
                    {
                        if(vC.getShortestDistance()==last.getShortestDistance())
                        {
                            arrTmp.add(vC);
                        }    
                    }
                    randomC = (int) (Math.random() * arrTmp.size()); 
                    randomCell = arrTmp.get(randomC);
                }
                else{   //random
                    randomC = (int) (Math.random() * lsTmp.size());
                    randomCell = lsTmp.get(randomC);
                }
                int ranInput = (int) (Math.random() * desLst.size());
                Vertex d2 = desLst.get(ranInput);
                
//                IJ.log("Change : "+ d2.getId() + " type: "+d2.getType()
//                        + " cell: " + randomCell.getId() +" type: "+randomCell.getType());
                byte typeTmp = randomCell.getType();
                randomCell.setType(d2.getType());
                d2.setType(typeTmp);
                popObserved.add(randomCell);
                desLst.remove(d2);
                count++;
//                IJ.log("_____"+st+"_____");
            }
            else{
                IJ.log("not have cell to change the type, increase the ratio of distance map");
                break;
            }
            
        } 
        IJ.log("Nb of count : "+count + " nb cell observed: "+sizeObserved);
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
        
//        computeLayer(typeSource, typeObserved, popTotal);
        calculCellDistanceStatistic(popTotal, typeObserved);
        calculEuclidStatistic(typeObserved);
        
    }
    
    private void calculCellDistanceStatistic(Cells3DPopulation popTotal, byte typeObserved)
    {
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : popTotal.getVertexes()) 
        {
            if (v1.getType() == typeObserved) 
            {
                desLst.add(v1);
            }
        }
        
        IJ.log(" desLst:  "+ typeObserved+"  " + desLst.size());
        int numPoints = popTotal.getVertexes().size();
        int numRandomSamples = 100;
        double env = 0.05;
        CellSpatialAnalysis spa = new CellSpatialAnalysis(numPoints, numRandomSamples, env, dir);
        spa.process(desLst, popTotal, "FG");
        
    }
    private void calculEuclidStatistic(byte typeObserved) 
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
        Objects3DPopulation pop = new Objects3DPopulation();
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
        IJ.log("Nb alpha in pop: "+pop.getNbObjects());
        
        
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
        int numPoints = lsCenterCells.size();
        int numRandomSamples = 100;
        double distHardCore = 2;
        double env = 0.05;

//        ImagePlus imagePlus = drawImg.duplicate().getImagePlus();
        

        SpatialTissueAnalysis spa = new SpatialTissueAnalysis(numPoints, numRandomSamples, distHardCore, env);
//        spa.process(pop, imgMask, lsCenterCells, "FG", false, true, false, dir);
//        rs.setValue("F_SDI_Eu"+typeObserved, 0, out.getValue(0));
//        rs.setValue("G_SDI_Eu"+typeObserved, 0, out.getValue(1));
//        rs.show("Statistic Function");
        
    }
    private void showDistanceMap(byte typeObserved, int typeDistribution, float ratioDis, byte typeSource)
    {
        Cells3DPopulation popTotal = initCellPopulation();
        ImageFloat drawDMR = new ImageFloat("DISTANCE MAP WAT", img.sizeX, img.sizeY, img.sizeZ);
        ImageFloat drawDMN = new ImageFloat("DISTANCE MAP NUC", img.sizeX, img.sizeY, img.sizeZ);
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
//        IJ.log("Size source: " + sourceLst.size() + " desLst: " + desLst.size());
        ArrayList<Vertex> popObserved = new ArrayList<Vertex>();
        
        //init, choose a first cell to swap type cell
        
        int ran = (int) (Math.random() * sourceLst.size());
        Vertex s = sourceLst.get(ran);
//        ArrayList<Cell> lsCore = getListShellCore(typeObserved);
//        int ran = (int) (Math.random() * lsCore.size());
//        Cell cTmp = lsCore.get(ran);
//        int keyC = (int) MapUtils.getKeyFromValue(cacheId, cTmp.id);
//        Vertex s = null;
//        for(Vertex ve : sourceLst)
//        {
//            if(ve.getId()==keyC)
//            {
//                s = ve;
//            }
//        }
//        if(s==null)
//        {
//            int r = (int) (Math.random() * sourceLst.size());
//            s = sourceLst.get(r);
//        }    
        int ranInit = (int) (Math.random() * desLst.size());
        Vertex d = desLst.get(ranInit);
        byte typeC = s.getType();
        s.setType(d.getType());
        d.setType(typeC);
        desLst.remove(d);
        popObserved.add(s);
        int cellObs = cacheId.get(s.getId());
        for(Cell c : popCells)
        {
           if(c.id == cellObs)
           {
               drawDMR.draw(c.region, 1100);
               drawDMN.draw(c.nucleus,1100);
           }    
        }    
        ArrayList<Vertex> sourceTmp = new ArrayList<Vertex>();
        for (Vertex v1 : vertexes) 
        {
            if (v1.getType() != typeObserved)
            {
                sourceTmp.add(v1);
            }
        }
        popTotal.computeShortestDistance(sourceTmp, popObserved);
        ArrayList<Vertex> lsTmp = normalizeDistanceMap(sourceTmp, ratioDis, typeDistribution);
        IJ.log("nb possible :" + lsTmp.size());
        
        for (Vertex vi : lsTmp) 
        {
            int cellId = cacheId.get(vi.getId());
            for(Cell c : popCells)
            {
               if(c.id == cellId)
               {
                   drawDMR.draw(c.region, (float)vi.getShortestDistance());
                   drawDMN.draw(c.nucleus, (float)vi.getShortestDistance());
               }    
            }    
        }
        drawDMN.show();
        drawDMR.show();
    }
    private void computeLayerV2(byte sourceType, byte destinationType, Cells3DPopulation popTotal)
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
            rtLayerOrigin.saveAs(dir+"Observed_Layer.csv");
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
            if(j < dis)
            {
                sRes[j] += res[j];
            }    
            
        }
        
    }
    public ArrayList<Vertex> normalizeDistanceMap(ArrayList<Vertex> lsCent, float ratio, int typeDistribution) 
    {
        ArrayList<Vertex> arr = new ArrayList<Vertex>();
        ArrayList<Vertex> lst = new ArrayList<Vertex>();
        int c = 0;
        for(Vertex ce : lsCent)
        {
            if(ce.getShortestDistance()!=-1)
            {
                c++;
            }
        }
        int count = 0;
        Vertex[] idx = new Vertex[c];
        double volume = idx.length;
//        IJ.log("volume is: "+volume);
        for(Vertex ce : lsCent)
        {
            if(ce.getShortestDistance()!=-1)
            {
                idx[count] = ce;
                count++;
            }    
            
        }
        Arrays.sort(idx);

        for (int i = 0; i < idx.length - 1; i++) 
        {
            // gestion des repetitions
            if (idx[i + 1].getShortestDistance() == idx[i].getShortestDistance()) {
                int j = i + 1;
                while (j < (idx.length - 1) && idx[i].getShortestDistance() == idx[j].getShortestDistance()) {
                    j++;
                }
                double median = (i + j) / 2d;
                for (int k = i; k <= j; k++) {
                    idx[k].index = median;
                }
                i = j;
            } else {
                idx[i].index = i;
            }
        }
        if (idx[idx.length - 1].index == 0) {
            idx[idx.length - 1].index = idx.length - 1;
        }
        
        if(typeDistribution==CLUSTER)
        {
            for (Vertex idx1 : idx) 
            {
                if((float) (idx1.index / volume) <= ratio)
                {
                    //IJ.log("id: "+idx1.getId()+"  "+"distributionType : initial: "+idx1.getShortestDistance());
                    idx1.setShortestDistance((float) (idx1.index / volume) * (float)1000.0);
                    arr.add(idx1);
                    //IJ.log("_________: after: "+idx1.getShortestDistance());
                }
                else{
                    //IJ.log("----------------------------------------------------------");
                    //IJ.log("id: "+idx1.getId()+"  "+"distributionType : initial: "+idx1.getShortestDistance());
                    //IJ.log("----------------------------------------------------------");
                    idx1.setShortestDistance(-1);
                }
            }
        }
        if(typeDistribution==UNIFORM)
        {
            for (Vertex idx1 : idx) 
            {
                if((float) (idx1.index / volume) >= ratio)
                {
                    //IJ.log("id: "+idx1.getId()+"  "+"distributionType : initial: "+idx1.getShortestDistance());
                    idx1.setShortestDistance((float) (idx1.index / volume) * (float)1000.0);
                    arr.add(idx1);
                    //IJ.log("_________: after: "+idx1.getShortestDistance());
                }
                else{
                    //IJ.log("----------------------------------------------------------");
                    //IJ.log("id: "+idx1.getId()+"  "+"distributionType : initial: "+idx1.getShortestDistance());
                    //IJ.log("----------------------------------------------------------");
                    idx1.setShortestDistance(-1);
                }
            }
            
        }
        if(typeDistribution==RANDOM)
        {
//            for (CentroidC idx1 : idx) 
//            {
//                idx1.setShortestDistance((float) (idx1.index / volume) * (float)1000.0);
//                arr.add(idx1);
//            }
            for (Vertex idx1 : idx) 
            {
                idx1.setShortestDistance((float) (idx1.index / volume) * (float)1000.0);
                arr.add(idx1);
//                if((float) (idx1.index / volume) <= ratio)
//                {
//                    //IJ.log("id: "+idx1.getId()+"  "+"distributionType : initial: "+idx1.getShortestDistance());
//                    idx1.setShortestDistance((float) (idx1.index / volume) * (float)1000.0);
//                    arr.add(idx1);
//                    //IJ.log("_________: after: "+idx1.getShortestDistance());
//                }
//                else{
//                    //IJ.log("----------------------------------------------------------");
//                    //IJ.log("id: "+idx1.getId()+"  "+"distributionType : initial: "+idx1.getShortestDistance());
//                    //IJ.log("----------------------------------------------------------");
//                    idx1.setShortestDistance(-1);
//                }
            }
        }
        Collections.sort(arr);
//        String str = "normalize: ";
//        for(Vertex v : arr)
//        {
//            str += v.getShortestDistance()+"   ";
//        } 
//        IJ.log(str);
        return arr;
       
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
    public void getBoundary()
    {
        for (Cell Ci : popCells) 
        {
            Ci.marked = false;
            Ci.connected = false;
        }
        Cell center = getCenterShellCore();
//        ArrayList<Cell> popShell = getListShellCore((byte)ALPHA);
//        int ranIdx = (int) (Math.random() * popShell.size());
//        Cell center = popShell.get(ranIdx);
        ArrayList<Cell> cls = BFS(center);
        IJ.log("Size of cluster is: "+ cls.size());
//        ImageHandler drawBFS = new ImageShort("BFS", img.sizeX, img.sizeY, img.sizeZ);
//        ImageHandler drawBFSNUC = new ImageShort("BFS", img.sizeX, img.sizeY, img.sizeZ);
//        ImageHandler drawLastNode = new ImageShort("LAST NODE", img.sizeX, img.sizeY, img.sizeZ);
        int maxLevel = -1;
        for (Cell C : cls) 
        {
            C.connected = true;
//            if(C.currLevel != 1)
//            {
//                C.region.draw(drawBFS, C.currLevel-1);
//                C.nucleus.draw(drawBFSNUC, C.currLevel-1);
//            }    
            
            if(maxLevel < C.currLevel)
            {
                maxLevel = C.currLevel;
            }    
//            if(C.lastnode)
//            {
//                C.region.draw(drawLastNode, C.currLevel);
//            }    
        }
//        for (Cell C : cls) 
//        {
//            if(C.currLevel==1)
//            {
//                C.region.draw(drawBFS, maxLevel+1);
//                C.nucleus.draw(drawBFSNUC, maxLevel+1);
//            }    
//        }
//        drawBFS.show();
//        drawBFSNUC.show();
//        drawLastNode.show();
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
        drawBFS.show();
    }
    
    
}
