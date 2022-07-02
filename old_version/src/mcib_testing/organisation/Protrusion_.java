/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
Plugins>3D Tissue Spatial Analysis,"-"
Plugins>3D Tissue Spatial Analysis,"DEEP SEG WAT",mcib_testing.RafaelV2.DeepSegWat_
 */
package mcib_testing.organisation;


import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.ResultsTable;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Point3D;
import mcib3d.geom.Vector3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import stdlib.PlotDraw;

import ij.measure.Calibration;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ThreadUtil;
import mcib3d.utils.exceptionPrinter;
import mcib_testing.Utils.Objects3DPopulationGrid;
import stdlib.DijkstraSPC;
import stdlib.DirectedEdgeC;
import stdlib.EdgeWeightedDigraphC;
import stdlib.Vertex;
/**
 * Description: Method for create the cluster organization 
 * steps: 
 * 1. Open the image nuclei and watershed, obtain the list of cell and type correspond. 
 * 2. Find list of delta cell and the first cell delta. 
 * 3. Compute the distance map, obtain the float image. 
 * 4. Compute the list of cell with ratio minimum
 * or equal 0.3 : list1 
 * 5. Change the position of delta cell with list 1 
 * 6. Redo step 3 to 5 until all delta cell have new position.
 *
 * @author tranhoa
 */
public class Protrusion_ implements ij.plugin.PlugIn {

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
    ImageHandler img, dapiImg;

    ArrayList<Cell> popCells = null;
    ArrayList<Cell> popInputCells = null;
    ArrayList<Cell> popS = null;
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    HashMap<Integer, Integer> cacheId = null;
    String dir = null;
    ResultsTable rtLayer = new ResultsTable();
    ResultsTable rtInf = new ResultsTable();
    ResultsTable rtContact = new ResultsTable();
    public boolean verbose = false;
    public double ratioMicro2Px=0, ratioZ2XY=0;
    public void run(String arg) 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return;
        }
        String[] titles = new String[wList.length];
        if (wList.length<2) {
            IJ.error("Open DAPI image.");
            return;
        }
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        String[] typeObserved = {"BETA","ALPHA"};
        String[] option = {"ALL DELTA CELLS","EACH DELTA CELL"};
        String[] regionObs = {"NUCLEUS REGION","CELL REGION"};
        byte typeCell = 0;
        int op = 0, reg = 0;
        int minDistance = 0, maxDistance = 35;
        
        GenericDialog gd = new GenericDialog("Protrusion");
        gd.addChoice("Watershed Image :", titles, titles[0]);
        gd.addChoice("DAPI Image :", titles, titles[1]);
        gd.addNumericField("Min Distance (MicroM) : ", minDistance, 0);
        gd.addNumericField("Max Distance (MicroM) : ", maxDistance, 35);
        gd.addChoice("Type Cell Observed : ", typeObserved, typeObserved[typeCell]);
        gd.addChoice("Range Observed : ", option, option[op]);
        gd.addChoice("Region Observed : ", regionObs, regionObs[reg]);
        gd.addCheckbox(" Show Protrusion Image", verbose);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        int[] index = new int[2];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        minDistance = (int) gd.getNextNumber();
        maxDistance = (int) gd.getNextNumber();
        typeCell = (byte)(gd.getNextChoiceIndex()+2);
        op =  (int)(gd.getNextChoiceIndex());
        reg =  (int)(gd.getNextChoiceIndex());
        verbose = gd.getNextBoolean();
        if(gd.wasOKed())
        {
            IJ.log("img: "+wList[index[0]]+" "+wList[index[1]]+" min distance: "+minDistance
                   +" max distance: "+maxDistance + " type cell observed: "+typeCell
                   +"    range(0-all, 1-each): "+op+"   region (0-nuc, 1-reg): " + reg);
        } 
        ImagePlus plus;
        plus = WindowManager.getImage(wList[index[0]]);
        img = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        ImagePlus plus1;
        plus1 = WindowManager.getImage(wList[index[1]]);
        dapiImg = ImageInt.wrap(plus1);
        dir = IJ.getDirectory("image");
        IJ.log("Analysing ...");
        IJ.log("Reading data from " + dir);
        popRegions = new Objects3DPopulationGrid();
        popNuclei = new Objects3DPopulationGrid();
        popRegions.loadObjects(dir + "Regions.zip");
        popNuclei.loadObjects(dir + "Nuclei.zip");

        imgWat = img.createSameDimensions();
        popRegions.draw(imgWat);
        //imgWat.show("Reg");
        imgSeg = img.createSameDimensions();
        popNuclei.draw(imgSeg);
        //imgSeg.show("Nuc");

        // init cells 
        initCells2();
        IJ.log("size is: " + popCells.size());
        //drawCellTypes(false).show("TYPE BEFORE");
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);// 0 for thomas´s cell seg; -1 for farsight
//        computeAllContacts(0);
//        IJ.log("***************LAYER LEVEL BEFORE*************************");
//        byte stB = 1, optionB = 1;
//        computeShortestDistanceCell2Cell(stB, optionB);
        
        IJ.log("**************PROTRUSION-EXPANDED**************************");
        byte typeInf = DELTA;
//        
        Calibration cal = dapiImg.getCalibration();
        ratioMicro2Px = getRatioMicro2Pixel(cal);
        ratioZ2XY = getRatio_Z_XY(cal);
        IJ.log("Ratio: MP:"+ratioMicro2Px+"  ration XY2Z:"+ratioZ2XY);
        if(op==0)
        {
            getTotalContacted(typeCell, typeInf, minDistance, maxDistance, reg);
        }
        else{
            getEachCellProContacted(typeCell, typeInf, minDistance, maxDistance, reg);
        }
        

        
//        IJ.log("**************** LAYER LEVEL AFTER ************************");
//        byte stA = 2, optionA = 1;
//        computeShortestDistanceCell2Cell(stA, optionA);
//        computeAllContacts(1);
//        
        
        IJ.log("Finished");

    }
    
    private double getRatioMicro2Pixel(Calibration cal)
    {
        return (1 * 1.0f / cal.pixelWidth);
    }        
    private double getRatio_Z_XY(Calibration cal)
    {
        return (cal.pixelWidth * 1.0f / cal.pixelDepth);
    }  
    /**
     * 
     * @param st : 0: before, 1: after
     */
    private void computeAllContacts(int st)
    {
        float aa = 0, bb = 0, dd = 0, ab = 0, ad = 0, bd = 0;
        float na = 0, nb = 0, nd = 0;
        //ResultsTable rtContact = new ResultsTable();
        for (Cell C : popCells) 
        {
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
        rtContact.incrementCounter();
        //int count = 0;
        rtContact.setValue("na", st, na);
        rtContact.setValue("nb", st, nb);
        rtContact.setValue("nd", st, nd);
        rtContact.setValue("totalabd", st, sumnb);
        
        rtContact.setValue("pa", st, (float)(na / sumnb));
        rtContact.setValue("pb", st, (float)(nb / sumnb));
        rtContact.setValue("pd", st, (float)(nd / sumnb));
        
        rtContact.setValue("aa", st, aa);
        rtContact.setValue("bb", st, bb);
        rtContact.setValue("dd", st, dd);
        rtContact.setValue("ab", st, ab);
        rtContact.setValue("ad", st, ad);
        rtContact.setValue("bd", st, bd);
        rtContact.setValue("all", st, sum);
        
        rtContact.setValue("paa", st, (float)(aa / sum));
        rtContact.setValue("pbb", st, (float)(bb / sum));
        rtContact.setValue("pdd", st, (float)(dd / sum));
        rtContact.setValue("pab", st, (float)(ab / sum));
        rtContact.setValue("pad", st, (float)(ad / sum));
        rtContact.setValue("pbd", st, (float)(bd / sum));
        IJ.log("Nb : na=" + na + " nb=" + nb + " nd=" + nd + " sum=" + sumnb);
        IJ.log("Nb : pa=" + na / sumnb + " pb=" + nb / sumnb + " pd=" + nd / sumnb);
        IJ.log("Assos : aa=" + aa + " bb=" + bb + " dd=" + dd + " ab=" + ab + " ad=" + ad + " bd=" + bd + " all=" + sum);
        IJ.log("Assos : aa=" + aa / sum + " bb=" + bb / sum + " dd=" + dd / sum + " ab=" + ab / sum + " ad=" + ad / sum + " bd=" + bd / sum);
        rtContact.show("AllContact_"+st);
    }
    private ImageHandler drawCellTypes(boolean nuc) 
    {
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

    private ImageHandler drawCellNuclei() {
        ImageHandler draw = new ImageShort("NUC", img.sizeX, img.sizeY, img.sizeZ);

        for (Cell C : popCells) {
            C.nucleus.draw(draw, C.id);
        }

        return draw;
    }

    private ImageHandler drawCellRegions() 
    {
        ImageHandler draw = new ImageShort("REG", img.sizeX, img.sizeY, img.sizeZ);

        for (Cell C : popCells) {
            C.region.draw(draw, C.id);
        }

        return draw;
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
    
    /**
     * 
     * @param st 
     */
    public void computeDistanceCell2Cell(int st) 
    {
        int sizeV = popCells.size();
        ArrayList<Vertex> vertexes = new ArrayList<Vertex>(sizeV);
        ArrayList<DirectedEdgeC> edges = new ArrayList<DirectedEdgeC>();
        for(Cell C : popCells) 
        {
            Vertex vC = new Vertex(C.id - 1, "Cell" + C.id, C.type);
            vertexes.add(vC);
        }
        IJ.log("Nb of vertex is: " + vertexes.size());

        for (Cell C : popCells) 
        {
            if (C.type > 0) {
                for (Cell C1 : C.nei1) 
                {
                    if (C1.type > 0) 
                    {
                        int r1 = C.id - 1;
                        int c1 = C1.id - 1;
//                        Point3D centerDCell = C.nucleus.getCenterAsPoint();
//                        Point3D center1 = C1.nucleus.getCenterAsPoint();
                        double weight = 1.0;
//                        if(typeDistance==1)
//                        {
//                            weight = 1.0;
//                        }
//                        else{
//                            weight = centerDCell.distance(center1);
//                        }
                        DirectedEdgeC e = new DirectedEdgeC(vertexes.get(r1), vertexes.get(c1), weight);
                        //IJ.log("e: "+e.toString());
                        edges.add(e);
                    }

                }
            }
        }
        IJ.log("Nb of edges is: " + edges.size());

        //initiale the graph
        EdgeWeightedDigraphC G = new EdgeWeightedDigraphC(vertexes, edges);

        //vertex source = list all vertexes without delta cell and unlabelled cell
        //vertex destination = list delta cells
        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : vertexes) {
            if (v1.getType() == BETA) {
                sourceLst.add(v1);
            }
            if (v1.getType() == DELTA) {
                desLst.add(v1);
            }
        }
        IJ.log("Size source:" + sourceLst.size() + " desLst: " + desLst.size());

        for (Vertex vs : sourceLst) 
        {
            int s = vs.getId();
            DijkstraSPC sp = new DijkstraSPC(G, s);
            int minStep = Integer.MAX_VALUE;
            double minDistance = Double.POSITIVE_INFINITY;
            //int idxS = -1;
            //int idxD = -1, idxDS = -1;
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                if (sp.hasPathTo(d)) 
                {
                    double distanceInG = sp.distTo(d);
                    //StdOut.printf("-------------------------------------");
                    //StdOut.printf("%d to %d (%.2f)  ", s, d, distanceInG);
//                    int step = 0;
//                    for (DirectedEdgeC e : sp.pathTo(d)) {
//                        //StdOut.print(e + "   ");
//                        step++;
//                    }
                    if (d != s && minDistance > distanceInG) {
                        minDistance = distanceInG;
                        //idxD = d;
                        //idxDS = step;
                    }
//                    if (d != s && minStep > step) {
//                        minStep = step;
//                        idxS = d;
//                    }
                }
            }
            vs.setShortestDistance(minDistance);
            ArrayList<Integer> idxShortestLayer = new ArrayList<Integer>();
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                if (sp.hasPathTo(d) && sp.distTo(d)==minDistance) 
                {
                    idxShortestLayer.add(d);
                }
            }
            vs.idxShortestLayer = idxShortestLayer;
        }

        int sizeDesLst = desLst.size();
        ArrayUtil xAxis = new ArrayUtil(sizeDesLst);
        ArrayUtil nbLayer1 = new ArrayUtil(sizeDesLst);
        ArrayUtil nbLayer2 = new ArrayUtil(sizeDesLst);
        ArrayUtil nbLayer3 = new ArrayUtil(sizeDesLst);
        ArrayUtil nbLayer4 = new ArrayUtil(sizeDesLst);
        ArrayUtil nbLayer5 = new ArrayUtil(sizeDesLst);
        try {
            FileWriter ryt = null;
            if(st==1)
            {
                ryt = new FileWriter(dir + "method3/deltaInf_beforeShuffle.txt");
            }
            else{
                ryt = new FileWriter(dir + "method3/deltaInf_afterShuffle.txt");
            }
            BufferedWriter out = new BufferedWriter(ryt);

            for (int d = 0; d < desLst.size(); d++) {
                out.write("D" + desLst.get(d).getId() + " ");
                xAxis.putValue(d, d+1);
            }
            
            out.write("\n");
            for (int d = 0; d < desLst.size(); d++) 
            {
                Vertex vdx = desLst.get(d);
                int count = 0, c1=0, c2=0, c3=0, c4=0, c5=0;
                for (int s = 0; s < sourceLst.size(); s++) 
                {
                    Vertex vsx = sourceLst.get(s);
                    if (vsx.idxShortestLayer.contains(vdx.getId())) 
                    {
                        count++;
                        if(vsx.getShortestDistance()==1.0)
                        {
                            c1++;
                        }
                        if(vsx.getShortestDistance()==2.0)
                        {
                            c2++;
                        }
                        if(vsx.getShortestDistance()==3.0)
                        {
                            c3++;
                        }
                        if(vsx.getShortestDistance()==4.0)
                        {
                            c4++;
                        }
                        if(vsx.getShortestDistance()==5.0)
                        {
                            c5++;
                        }
                    }
                }
                nbLayer1.putValue(d, c1);
                nbLayer2.putValue(d, c2);
                nbLayer3.putValue(d, c3);
                nbLayer4.putValue(d, c4);
                nbLayer5.putValue(d, c5);
                out.write(count + " ");
            }
            
            out.close();
        } catch (IOException ex) {
            Logger.getLogger(Protrusion_.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.log("total: layer1"+nbLayer1.getSum()+" layer2: "+nbLayer2.getSum()+" layer3: "+nbLayer3.getSum()
                     + "layer4: "+nbLayer4.getSum() + " layer5: "+ nbLayer5.getSum());
        Plot layerPlot = PlotDraw.createPlot("Layer Level", "delta cell","nb cell influenced by delta cell", xAxis, nbLayer1, nbLayer2, nbLayer3);
        layerPlot.draw();
        PlotWindow plotLayer = layerPlot.show();

    }

    public Point3D rotationVector(Point3D center, double theta, double translation)
    {
        
        double x1 = center.getRoundX()*translation*Math.cos(Math.toRadians(theta)) - center.getRoundY()*translation * Math.sin(Math.toRadians(theta));
        double y1 = center.getRoundX()*translation*Math.sin(Math.toRadians(theta)) + center.getRoundY()*translation* Math.cos(Math.toRadians(theta));
        double z1 = center.getRoundZ()*10;
        Point3D p = new Point3D(x1, y1, z1);
        return p;
    }
    
    public ImageByte createImgProtrusionV4(double microDistance, byte typeInf, ArrayList<Cell> popTmp, int reg, int st)
    {
        IJ.log("Nb objs of seg image: "+ popTmp.size());
        String title = "Dilated_"+microDistance+"_"+reg+"_"+st;
        ImageHandler label = new ImageShort("seg_", img.sizeX, img.sizeY, img.sizeZ);
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
        if(verbose)
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
            return temp;
        } catch (Exception e) {
            exceptionPrinter.print(e, null, true);
        }
        return null;
        
//        if(verbose)
//        {
//            temp.show(title);
//        }
//        Objects3DPopulationGrid popObjs = new Objects3DPopulationGrid((ImageInt)dilateImg);
//        IJ.log("Nb objs of dilated image: "+ popObjs.getNbObjects());
//        return popObjs;
        
    }        
    public ArrayList<Object3DVoxels> createImgProtrusionV3(double microDistance, byte typeInf, ArrayList<Cell> popTmp, int reg, int st)
    {
        ArrayList<Object3DVoxels> arr = new ArrayList<Object3DVoxels>();
        String title = "Dilated_"+microDistance+"_"+reg+"_"+st;
        ImageHandler label = new ImageShort(title, img.sizeX, img.sizeY, img.sizeZ);
        for (Cell C : popTmp) 
        {
            Object3D region = C.region;
            Object3D nuc = C.nucleus;
            float rDilate = (float) (microDistance * ratioMicro2Px);
            Object3DVoxels tempObj;
            if(reg==0)  //observed nucleus region
            {
                tempObj = nuc.getDilatedObject(rDilate, rDilate, (float) (rDilate*ratioZ2XY));
            }
            else{        //observed cell region
                tempObj = region.getDilatedObject(rDilate, rDilate, (float) (rDilate*ratioZ2XY));
            }
            tempObj.setValue(region.getValue());
            tempObj.draw(label, region.getValue());
            arr.add(tempObj);
        }
        if(verbose)
        {
            label.show();
        }    
        
        return arr;
    }        
    public ArrayList<Object3DVoxels> createImgProtrusionV2(double microDistance, byte typeInf)
    {
        ObjectCreator3D objs;
        objs = new ObjectCreator3D(img.sizeX, img.sizeY, img.sizeZ);
//        ImageHandler imgProtrusion = new ImageByte("CellProtrusion"+microDistance,img.sizeX, img.sizeY, img.sizeZ);
        ArrayList<Object3DVoxels> arr = new ArrayList<Object3DVoxels>();
        int count = 0;
        String valReg = " ";
        for (Cell C : popCells) 
        {

            if (C.type == typeInf)
            {
                count++;
                Object3D region = C.region;
                Object3D nuc = C.nucleus;
                Point3D centerDCell = nuc.getCenterAsPoint();
                Vector3D center = new Vector3D(centerDCell);
                double ratioMicro2Pixel = 2.64;
                double distanceProtrusion = microDistance * ratioMicro2Pixel;
                double averDistance = nuc.getDistCenterMean();
                
                double totalDistance = averDistance + distanceProtrusion;
                objs.createEllipsoid(center, totalDistance, totalDistance, totalDistance, region.getValue());
                Object3DVoxels obj = objs.getObject3DVoxels(region.getValue());
                arr.add(obj);
                //what is the problem if the distance go beyond the boundary
                valReg += "   "+region.getValue()+ " ";
                
               
            }
        }
        
        IJ.log("distance : "+microDistance+" nb: "+count);
        IJ.log("val:  "+valReg);
        IJ.log("size of arr: "+ arr.size());
        return arr;
        
//        ImagePlus imgProtrusion = new ImagePlus("CellProtrusion"+microDistance, objs.getStack());
//        Objects3DPopulationGrid popObjs = new Objects3DPopulationGrid(arr);
//        IJ.log("step: "+microDistance+"  size: "+ popObjs.getNbObjects());
    } 
    
    
    public ImageInt createImgProtrusion(double microDistance)
    {
        ObjectCreator3D objs;
        objs = new ObjectCreator3D(img.sizeX, img.sizeY, img.sizeZ);
        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
        
        Point3D centerSCore = computeDegreeConnection();
        Vector3D centerCore = new Vector3D(centerSCore);
//        Voxel3D c = new Voxel3D(centerSCore, 255);
//        arr.add(c);
        //IJ.log("center shell core is: "+centerSCore);
        for (Cell C : popCells) 
        {

            if (C.type == DELTA)
            {
                //IJ.log("delta cell");
                Object3D region = C.region;
                Object3D nuc = C.nucleus;
//                ArrayList<Voxel3D> vox = region.getVoxels();
//                arr.addAll(vox);
                //Object3DVoxels re = (Object3DVoxels)region;
                //Object3D nuc = C.nucleus;
                Point3D centerDCell = nuc.getCenterAsPoint();
                double ratioMicro2Pixel = 2.64;
                double distanceProtrusion = microDistance * ratioMicro2Pixel;
                //IJ.log("distance is: "+distanceProtrusion);
                
                
                Vector3D direction = new Vector3D(centerDCell, centerSCore);
                double distanceCenterBorder = region.distPixelBorderUnit(centerDCell.getRoundX(), centerDCell.getRoundY(),
                                                                      centerDCell.getRoundZ(), centerCore);
                double totalDistance = distanceCenterBorder + distanceProtrusion;
                //IJ.log("total distance is: "+totalDistance);
                direction.normalize();
                Point3D target = new Point3D(totalDistance * direction.getX() + centerDCell.getX(),totalDistance * direction.getY() + centerDCell.getY(),
                                             totalDistance * direction.getZ()+ centerDCell.getZ());
                //IJ.log("target is: "+target);
                objs.createLine(centerDCell, target, region.getValue(), 10);
                //objs.createLine(centerDCell, centerSCore, region.getValue(), 1);
//                Voxel3D vC = new Voxel3D(centerDCell, region.getValue());
//                Voxel3D vT = new Voxel3D(target, region.getValue());
                //arr.add(vC);
                //arr.add(vT);
            }
        }
        objs.drawVoxels(arr);
        ImagePlus imgProtrusion = new ImagePlus("CellProtrusion"+microDistance, objs.getStack());
        //plusShape.setSlice((int) (cz2));
        imgProtrusion.setDisplayRange(0, 255);
        //imgProtrusion.show();
        ImageInt imgRes = ImageInt.wrap(imgProtrusion);
        //imgRes.show();
        return imgRes;
    } 
    
    /**
     * met qua, ec ec
     * @param typeObserved 
     */
    public void randomProtrusionV2(byte typeObserved, int disProtrusion)
    {
        
        ArrayList<Cell> lsCells = new ArrayList<Cell>();
        IJ.log("----------------BEFORE---------------------");
        int nbInf = 0;
        for (Cell C : popCells) 
        {
            if (C.type == DELTA) 
            {
                //int nbLabelled = 0;
                for (Cell Ci : C.nei1) 
                {
                    if (Ci.type == typeObserved && !lsCells.contains(Ci)) 
                    {
                        lsCells.add(Ci);
                        nbInf++;
                        //nbLabelled++;
                    }
                }
                //nbInf += nbLabelled;
            }
        }
        IJ.log("step 0: "+ nbInf + " nb beta cell influenced is : "+lsCells.size());
        
        IJ.log("----------------AFTER PROTRUSION---------------------");
        //compute nb of influenced
       
        ArrayList<Cell> lsCells2 = new ArrayList<Cell>();
        //create image protrusion
        ImageInt imgRes = createImgProtrusion(disProtrusion);
        imgRes.show();

        Object3DVoxels[] arrObjs = imgRes.getObjects3D();
        //IJ.log("step"+d+" size: "+arrObjs.length);
        //int nbCells = 0; 
        for (int i = 0; i < arrObjs.length; i++) 
        {
            Object3DVoxels voxx = arrObjs[i];
            int idDeltaC = voxx.getValue();
            Cell C = region2Cell.get(idDeltaC);
            //IJ.log("delta: "+idDeltaC);
            ArrayUtil ac = voxx.listValues(imgWat);
            ac = ac.distinctValues();
            //IJ.log("--------------------------------------");
            //IJ.log("delta cell:" + idDeltaC);
            for (int j = 0; j < ac.getSize(); j++) 
            {
                if (ac.getValueInt(j) != idDeltaC && ac.getValueInt(j) != 0) 
                {
                    int idx = ac.getValueInt(j);
                    Cell C1 = region2Cell.get(idx);
                    if ((C1 != null) && (C1 != C) && !C.nei1.contains(C1)) 
                    {
                        C.nei1.add(C1);
                        if(!C1.nei1.contains(C))
                        { 
                            C1.nei1.add(C);
                        }
                        if(C1.type==typeObserved && !lsCells2.contains(C1))
                        {
                            //nbCells++;
                            lsCells2.add(C1);
                        }
                    }
                }
            }
        }    
        IJ.log("step 1 (after protrusion) nb beta cell influenced is : "+(nbInf+lsCells2.size())); 
    }    
    public void computeCellInfluenced(byte typeObserved, byte typeInf, int minDistance, int maxDistance)
    {
        
        //compute xaxis 
        int sizeDis = maxDistance - minDistance + 1;
//        ArrayUtil xAxis = new ArrayUtil(sizeDis);
//        ArrayUtil nbCellInf = new ArrayUtil(sizeDis);
        ArrayList<Cell> lsCells = new ArrayList<Cell>();
        ResultsTable rtCellsPro = new ResultsTable();
        
//        for(int x1=0; x1<sizeDis; x1++)
//        {
//            xAxis.putValue(x1, x1+minDistance);
//        }
        
        IJ.log("----------------BEFORE---------------------");
        int nbInf = 0;
        for (Cell C : popCells) 
        {
            if (C.type == DELTA) 
            {
                int nbLabelled = 0;
                for (Cell Ci : C.nei1) 
                {
                    if (Ci.type == typeObserved && !lsCells.contains(Ci)) 
                    {
                        lsCells.add(Ci);
                        nbLabelled++;
                    }
                }
                nbInf += nbLabelled;
            }
        }
//        nbCellInf.putValue(0, lsCells.size());
//        rtCellsPro.incrementCounter();
//        rtCellsPro.setValue("MicroM", 0, 0);
//        rtCellsPro.setValue("nbCellInfluenced", 0, lsCells.size());
        IJ.log("step 0 (before propagation): "+ nbInf + "size list cell is :"+lsCells.size());
        
        IJ.log("----------------AFTER PROTRUSION---------------------");
        //compute nb of influenced
        int nb = 0;
        for(int d=minDistance; d<=maxDistance; d++)
        {
            ArrayList<Cell> lsCells2 = new ArrayList<Cell>();
            ArrayList<Object3DVoxels> arr = createImgProtrusionV2(d, typeInf);
            
            for (Object3DVoxels o : arr) 
            {
                int idDeltaC = o.getValue();
                Cell C = region2Cell.get(idDeltaC);
                ArrayUtil ac = o.listValues(imgWat);
                ac = ac.distinctValues();
                for (int j = 0; j < ac.getSize(); j++) 
                {
                    if (ac.getValueInt(j) != idDeltaC && ac.getValueInt(j) != 0) 
                    {
                        int idx = ac.getValueInt(j);
                        Cell C1 = region2Cell.get(idx);
                        if ((C1 != null) && (C1 != C) && !C.nei1.contains(C1) 
                            && C1.type==typeObserved && !lsCells2.contains(C1)) 
                        {
                            lsCells2.add(C1);
                        }
                    }
                }
            }
            IJ.log("distance:  "+d+ "   inscrease: "+ lsCells2.size());
            rtCellsPro.incrementCounter();
//            nbCellInf.putValue(nb, nbInf+lsCells2.size());
            rtCellsPro.setValue("MicroM", nb, d);
            rtCellsPro.setValue("nbCellInfluenced", nb, nbInf+lsCells2.size());
            nb++;
        }
        if(typeObserved==BETA)
        {
            rtCellsPro.show("Beta Cells Influenced");
        }    
        if(typeObserved==ALPHA)
        {
            rtCellsPro.show("Alpha Cells Influenced");
        }
//        Plot layerPlot = PlotDraw.createPlot("Level 1 Influenced","distance (micro m)","nb cells influenced by delta cell", xAxis, nbCellInf);
//        PlotWindow.noGridLines = true;
//        layerPlot.draw();
//        PlotWindow plotLayer = layerPlot.show();
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
//        ImageByte res = createImgProtrusionV4(minDistance, typeInf, popTmp, reg, st);
//        res.show();
        computeContacted(typeObserved, typeInf, minDistance, maxDistance, popTmp, reg, st);
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
                computeContacted(typeObserved, typeInf, minDistance, maxDistance, popTmp, reg, st);
            }
        }
    }        
    public void computeContacted(byte typeObserved, byte typeInf, int minDistance, int maxDistance, ArrayList<Cell> popTmp, int reg, int st)
    {
        IJ.log("nb delta cells: "+popTmp.size());
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
            ImageByte res = createImgProtrusionV4(d, typeInf, popTmp, reg, st);
//            res.show();
            Objects3DPopulationGrid popObjs = new Objects3DPopulationGrid((ImageInt)res);
//            IJ.log("Nb objs of dilated image: "+ popObjs.getNbObjects());
//            ArrayList<Object3DVoxels> arr = createImgProtrusionV3(d, typeInf, popTmp, reg, st);
//            IJ.log("Size of array: "+ arr.size());
//            IJ.log("nb objs: "+ popObjs.getNbObjects());
//            String debug = "Id Delta Cell: ";
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
//                IJ.log("ac : "+ac.toString());
            }
            
            int nbInfAfter = lsCells2.size();
            IJ.log("distance:  "+d+ "   inscrease: "+ lsCells2.size());
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
        if(st==0)
        {
            title += "_All";
        }  
        else{
            title += "_del_"+st;
        }
        rtCellsPro.show(title);
        try {
            rtCellsPro.saveAs(dir+title+".csv");
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
//        Plot layerPlot = PlotDraw.createPlot("Level 1 Influenced","distance (micro m)","nb cells influenced by delta cell", xAxis, nbCellInf);
//        PlotWindow.noGridLines = true;
//        layerPlot.draw();
//        PlotWindow plotLayer = layerPlot.show();
    }
    /**
     * 
     */
//    public void randomProtrusion() 
//    {
//        //cacheId = new HashMap<Integer, Integer>();
//        ObjectCreator3D objs;
//        objs = new ObjectCreator3D(img.sizeX, img.sizeY, img.sizeZ);
//        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
//        for (Cell C : popCells) 
//        {
//
//            if (C.type == DELTA && (C.id % 4 == 0))
//            {
//                Object3D region = C.region;
//                Object3D nuc = C.nucleus;
//                ArrayList<Voxel3D> vox = region.getVoxels();
//                arr.addAll(vox);
//                //Object3DVoxels re = (Object3DVoxels)region;
//                //Object3D nuc = C.nucleus;
//                Point3D centerDCell = nuc.getCenterAsPoint();
//                double ratioMicro2Pixel = 0.6;
//                double distanceProtrusion = 2 * ratioMicro2Pixel;
//                
//                Point3D centerSCore = computeDegreeConnection();
//                Vector3D centerCore = new Vector3D(centerSCore);
//                Vector3D direction = new Vector3D(centerDCell, centerSCore);
//                double distanceCenterBorder = nuc.distPixelBorderUnit(centerDCell.getX(), centerDCell.getY(),
//                                                                      centerDCell.getZ(), centerCore);
//                double totalDistance = distanceCenterBorder + distanceProtrusion;
//                
//                direction.normalize();
//                Point3D target = new Point3D(totalDistance * direction.getRoundX(),totalDistance * direction.getRoundY(),
//                                             totalDistance * direction.getRoundZ());
//                objs.createLine(centerDCell, target, region.getValue(), 10);
//                Voxel3D vC = new Voxel3D(centerDCell, region.getValue());
//                Voxel3D vT = new Voxel3D(centerSCore, region.getValue());
//                arr.add(vC);
//                arr.add(vT);
//            }
//        }
//        objs.drawVoxels(arr);
//        ImagePlus imgProtrusion = new ImagePlus("CellProtrusion", objs.getStack());
//        //plusShape.setSlice((int) (cz2));
//        imgProtrusion.setDisplayRange(0, 255);
//        //imgProtrusion.show();
//        ImageInt imgRes = ImageInt.wrap(imgProtrusion);
//        Object3DVoxels[] arrObjs = imgRes.getObjects3D();
//
//        int max = (int) imgWat.getMax();
//        ArrayList<Integer>[] neiAdj = new ArrayList[max + 1];
//        for (int i = 0; i < neiAdj.length; i++) {
//            neiAdj[i] = new ArrayList();
//        }
//        for (int i = 0; i < arrObjs.length; i++) 
//        {
//            Object3DVoxels voxx = arrObjs[i];
//            int idDeltaC = voxx.getValue();
//            ArrayUtil ac = voxx.listValues(imgWat);
//            ac = ac.distinctValues();
//            //IJ.log("--------------------------------------");
//            //IJ.log("delta cell:" + idDeltaC);
//            for (int j = 0; j < ac.getSize(); j++) {
//
//                if (ac.getValueInt(j) != idDeltaC && ac.getValueInt(j) != 0) 
//                {
//                    neiAdj[idDeltaC].add(ac.getValueInt(j));
//                    //IJ.log(idDeltaC + " --> " + ac.getValueInt(j));
//                }
//            }
//            //IJ.log("size " + neiAdj[idDeltaC].size());
//        }
//
//        int nbDeltaCells = 0;
//        //ArrayList<Integer> arrId = new ArrayList<>();
//        for (Cell C : popCells) {
//            if (C.type == DELTA) {
//                //arrId.add(C.id);
//                nbDeltaCells++;
//            }
//        }
//        ArrayUtil xAxis = new ArrayUtil(nbDeltaCells);
//        for(int k=0; k<nbDeltaCells; k++)
//        {
//            xAxis.putValue(k, k);
//        }
//        ArrayUtil nbCells1 = new ArrayUtil(nbDeltaCells);
//        ArrayUtil nbCells2 = new ArrayUtil(nbDeltaCells);
//        try {
//            FileWriter ryt = new FileWriter(dir + "neiDeltaContact_labelled.txt");
//            BufferedWriter out = new BufferedWriter(ryt);
//
//            //before protrusion
//            //int nbDeltaCells = 0;
//            for (Cell C : popCells) {
//                if (C.type == DELTA) {
//                    //nbDeltaCells++;
//                    out.write("D" + C.id + " ");
//                }
//            }
//            //ArrayUtil nbCells1 = new ArrayUtil(nbDeltaCells);
//            //ArrayUtil nbCells2 = new ArrayUtil(nbDeltaCells);
//            out.write("\n");
//            
//            IJ.log("----------------BEFORE---------------------");
//            int c1 = 0;
//            for (Cell C : popCells) 
//            {
//                if (C.type == DELTA) 
//                {
//                    IJ.log("*********delta: "+C.id);
//                    int nbLabelled = 0;
//                    for (Cell Ci : C.nei1) 
//                    {
//                        if (Ci.type > 0) 
//                        {
//                            nbLabelled++;
//                            IJ.log("id:"+Ci.id);
//                        }
//                    }
//                    
//                    out.write(nbLabelled + " ");
//                    nbCells1.putValue(c1, nbLabelled);
//                    c1++;
//                }
//            }
//            out.write("\n");
//            // index refers to region label
//            IJ.log("----------------AFTER---------------------");
//            IJ.log("size of nei adj is: " + neiAdj.length);
//            for (int i = 0; i < neiAdj.length; i++) 
//            {
//                if (neiAdj[i].isEmpty()) {
//                    continue;
//                }
//                Cell C = region2Cell.get(i);
//                if (C != null) 
//                {
//                    for (int j : neiAdj[i]) {
//                        Cell C1 = region2Cell.get(j);
//                        if ((C1 != null) && (C1 != C) && !C.nei1.contains(C1)) 
//                        {
//                            C.nei1.add(C1);
//                            if(!C1.nei1.contains(C))
//                            { 
//                                C1.nei1.add(C);
//                            }
//                            //IJ.log("ok");
//                        }
//                    }
//                }
//            }
//            
//            int c2 = 0;
//            for (Cell C : popCells) 
//            {
//                if (C.type == DELTA) 
//                {
//                    IJ.log("*********delta: "+C.id);
//                    int nbLabelled = 0;
//                    for (Cell Ci : C.nei1) 
//                    {
//                        if (Ci.type > 0) 
//                        {
//                            IJ.log("id:"+Ci.id);
//                            nbLabelled++;
//                        }
//                    }
//                    out.write(nbLabelled + " ");
//                    nbCells2.putValue(c2, nbLabelled);
//                    c2++;
//                    
//                }
//            }
//            
//            out.close();
//        } catch (IOException ex) {
//            Logger.getLogger(Protrusion_.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        Plot deltaPlot = PlotDraw.createPlot("Level 1","delta cell","nb cell influenced by delta cell", xAxis, nbCells1, nbCells2, null);
//        deltaPlot.draw();
//        PlotWindow plotD = deltaPlot.show();
//        imgRes.show();
//
//    }

    /**
     * 
     * @param st 1: before shuffle cells, 2: after shuffle cells
     * @param option : 1: compute shortest layer distance or 2: compute shortest path distance
     */
    private void computeShortestDistanceCell2Cell(byte st, byte option) 
    {
        String stStr = null;
        if(st==1)
        {
            stStr = "Before";
            IJ.log("Compute the distance with cells original");
        }
        else{
            stStr = "After";
            IJ.log("Compute the distance after shuffle cells type");
        }
        int sizeV = popCells.size();
        ArrayList<Vertex> vertexes = new ArrayList<>(sizeV);
        ArrayList<DirectedEdgeC> edges = new ArrayList<>();
        for (Cell C : popCells) 
        {
            Vertex vC = new Vertex(C.id - 1, "Cell" + C.id, C.type);
            vertexes.add(vC);
        }
        IJ.log("Nb of vertex is: " + vertexes.size());

        for (Cell C : popCells) 
        {
            if (C.type > 0) 
            {
                for (Cell C1 : C.nei1) 
                {
                    if (C1.type > 0) {
                        int r1 = C.id - 1;
                        int c1 = C1.id - 1;
                        double weight = 0.0;
                        if(option==1)     //shortest layer distance
                        {
                            weight = 1.0;
                        }
                        if(option==2)     //shortest path distance
                        {
                            Point3D center = C.nucleus.getCenterAsPoint();
                            Point3D center1 = C1.nucleus.getCenterAsPoint();
                            weight = center.distance(center1);
                        }    
                        DirectedEdgeC e = new DirectedEdgeC(vertexes.get(r1), vertexes.get(c1), weight);
                        edges.add(e);
                    }
                }
            }
        }
        IJ.log("Nb of edges is: " + edges.size());

        //initiale the graph
        EdgeWeightedDigraphC G = new EdgeWeightedDigraphC(vertexes, edges);

        //vertex source : list all beta cells - list vertexes 
        //vertex destination = list all delta cells
        ArrayList<Vertex> sourceLst = new ArrayList<Vertex>();
        ArrayList<Vertex> desLst = new ArrayList<Vertex>();
        for (Vertex v1 : vertexes) 
        {
            if (v1.getType() == BETA)
            {
                sourceLst.add(v1);
            }
            if (v1.getType() == DELTA) {
                desLst.add(v1);
            }
        }
        
        IJ.log("Size source:" + sourceLst.size() + " desLst: " + desLst.size());

        for (Vertex vs : sourceLst) 
        {
            int s = vs.getId();
            DijkstraSPC sp = new DijkstraSPC(G, s);
            //IJ.log("s: "+s);
            double minDistance = Double.POSITIVE_INFINITY;
            int idxD = -1;
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                if (sp.hasPathTo(d)) 
                {
                    double distanceInG = sp.distTo(d);
                    //StdOut.printf("-------------------------------------");
                    //StdOut.printf("%d to %d (%.2f)  ", s, d, distanceInG);
                    
                    if (d != s && minDistance > distanceInG) {
                        minDistance = distanceInG;
                        idxD = d;
                    }
                }

            }
            if(minDistance < Double.POSITIVE_INFINITY)
            {
                vs.setShortestDistance(minDistance);
            }    
            
            ArrayList<Integer> idxShortestLayer = new ArrayList<Integer>();
            Point3D center = popCells.get(s).nucleus.getCenterAsPoint();
            double minEuDistance = Double.POSITIVE_INFINITY;
            int idxMin = -1;
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                
                if (sp.hasPathTo(d) && sp.distTo(d)==minDistance) 
                {
                    idxShortestLayer.add(d);
                    
                    Point3D center1 = popCells.get(d).nucleus.getCenterAsPoint();
                    double dis = center.distance(center1);
                    if(minEuDistance > dis){
                        minEuDistance = dis;
                        idxMin = d;
                    }
                }
            }
            if(idxMin != -1)
            {
                vs.idDCShortestLayer = idxMin; 
            }  
            vs.idxShortestLayer = idxShortestLayer;
        }
        drawResult(st, stStr, sourceLst, desLst);

    }
    private void drawResult(byte st, String stStr, ArrayList<Vertex> sourceLst, ArrayList<Vertex> desLst)
    {
        ImageHandler drawBDCont = new ImageShort("BD Contact "+stStr, img.sizeX, img.sizeY, img.sizeZ);
        ImageHandler drawLayer = new ImageShort("LAYER"+stStr, img.sizeX, img.sizeY, img.sizeZ);
        int c1=0, c2=0, c3=0, c4=0, c5=0;
        for (int d = 0; d < desLst.size(); d++) 
        {
            Vertex vdx = desLst.get(d);
            if(st==1)
            {
                rtInf.incrementCounter();
            }
            
            //initial id delta cell
            int idDeltaBefore = 0;
            int idDeltaAfter = 0;
            if(st==1)
            {
                idDeltaAfter = 0;
                idDeltaBefore = vdx.getId();
                
            }
            if(st==2)  //after shuffle
            {
                idDeltaAfter = vdx.getId();
                idDeltaBefore = vdx.getId();
//                for(int k : cacheId.keySet())
//                {
//                    if(cacheId.get(k)==idDeltaAfter)
//                    {
//                        idDeltaBefore = k;
//                    }
//                }

            }
            rtInf.setValue("idDeltaBefore", d, idDeltaBefore+1);
            rtInf.setValue("idDeltaAfter", d, idDeltaAfter+1);


            Cell dx = popCells.get(vdx.getId());
            dx.region.draw(drawBDCont, idDeltaBefore+1);
            dx.region.draw(drawLayer, 10);
            int count = 0, count1 = 0;
            for (int s = 0; s < sourceLst.size(); s++) 
            {
                Vertex vsx = sourceLst.get(s);
                Cell sx = popCells.get(vsx.getId());
                sx.region.draw(drawLayer, (int)vsx.getShortestDistance());
                if (vsx.idxShortestLayer.contains(vdx.getId())) 
                {
                    count++;
                    if(vsx.getShortestDistance()==1.0)
                    {
                        c1++;
                    }
                    if(vsx.getShortestDistance()==2.0)
                    {
                        c2++;
                    }
                    if(vsx.getShortestDistance()==3.0)
                    {
                        c3++;
                    }
                    if(vsx.getShortestDistance()==4.0)
                    {
                        c4++;
                    }
                    if(vsx.getShortestDistance()==5.0)
                    {
                        c5++;
                    }
                    
                    sx.region.draw(drawBDCont, idDeltaBefore+1);
                    
                }
                if(vsx.idDCShortestLayer==vdx.getId())
                {
                    count1++;
                }
            }
            if(st==1)
            {
                rtInf.setValue("InfCellBefore", d, count);
                rtInf.setValue("InfCellAfter", d, 0);
            }
            if(st==2)
            {
                rtInf.setValue("InfCellAfter", d, count);
            }
        }
          rtLayer.incrementCounter();
          rtLayer.setValue("layer1", st-1, c1);
          rtLayer.setValue("layer2", st-1, c2);
          rtLayer.setValue("layer3", st-1, c3);
          rtLayer.setValue("layer4", st-1, c4);
          rtLayer.setValue("layer5", st-1, c5);
          if(st==2)
          {
              ArrayUtil nbLayerBefore = new ArrayUtil(5);
              ArrayUtil nbLayerAfter = new ArrayUtil(5);
              ArrayUtil xAxis = new ArrayUtil(5);
              for(int k=0; k<5; k++)
              {
                  xAxis.putValue(k, k+1);
              }
              nbLayerBefore.putValue(0, rtLayer.getValue("layer1", 0));
              nbLayerBefore.putValue(1, rtLayer.getValue("layer2", 0));
              nbLayerBefore.putValue(2, rtLayer.getValue("layer3", 0));
              nbLayerBefore.putValue(3, rtLayer.getValue("layer4", 0));
              nbLayerBefore.putValue(4, rtLayer.getValue("layer5", 0));
              
              nbLayerAfter.putValue(0, rtLayer.getValue("layer1", 1));
              nbLayerAfter.putValue(1, rtLayer.getValue("layer2", 1));
              nbLayerAfter.putValue(2, rtLayer.getValue("layer3", 1));
              nbLayerAfter.putValue(3, rtLayer.getValue("layer4", 1));
              nbLayerAfter.putValue(4, rtLayer.getValue("layer5", 1));
              ArrayUtil[] nbLayer = new ArrayUtil[2];
              nbLayer[0] = nbLayerBefore;
              nbLayer[1] = nbLayerAfter;
              rtLayer.show("layer"+stStr);
              Plot layerPlot = PlotDraw.createPlotLayer("layer influence", "layer st", "nb cell influenced", xAxis, nbLayer);
              layerPlot.draw();
              PlotWindow plotLayer = layerPlot.show();
              rtInf.show("InfluenceOfDeltaCell");
          }
          drawLayer.show();
          drawBDCont.show();
            
    }  

    /**
     * Compute the centerDCell of the clusters cells
     */
    private Point3D computeDegreeConnection() 
    {
        popS = new ArrayList<Cell>();
        ArrayList<Cell> sCells = new ArrayList<Cell>();
        sCells.addAll(popCells);
        for (Cell C : popCells) {
            int cLabelled = 0, cUnLabelled = 0;
            //double ratio = 0;
            if (C.type > 0) {
                for (Cell C1 : C.nei1) {
                    if ((C1.type > 0)) {
                        cLabelled++;
                    } else {
                        cUnLabelled++;
                    }
                }
            }
            //ratio = cUnLabelled / cLabelled;
            C.degreeLabelled = cLabelled;
            C.degreeUnlabelled = cUnLabelled;
            if (C.degreeLabelled < 8 && C.degreeLabelled > 0 && C.degreeUnlabelled > 0) {
                popS.add(C);
                sCells.remove(C);
                //IJ.log("id: "+ C.id + " labelled: "+C.degreeLabelled + " unlabelled: " + C.degreeUnlabelled);
            }
        }
        Collections.sort(popS);
        ImageHandler draw = new ImageShort("TYPE", img.sizeX, img.sizeY, img.sizeZ);

        //for (Cell C1 : popS) 
        for (int k = popS.size() - 1; k >= 0; k--) 
        {
            Cell C1 = popS.get(k);
//            IJ.log("Cell: " + C1.id + " l:" + C1.degreeLabelled + " unL: " + C1.degreeUnlabelled
//                    + " ratio: " + (float) C1.degreeUnlabelled / C1.degreeLabelled);

            C1.region.draw(draw, C1.degreeLabelled * 10);
        }
        double totalX=0.0, totalY=0.0, totalZ=0.0;
        for(Cell s : sCells)
        {
            totalX += s.nucleus.getCenterX();
            totalY += s.nucleus.getCenterY();
            totalZ += s.nucleus.getCenterZ();
        }
        int sizeSCells = sCells.size();
        Point3D target = new Point3D(totalX/sizeSCells, totalY/sizeSCells, totalZ/sizeSCells);
        return target;
        //draw.show("DEGREE");
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

    public class Cell implements Comparable<Cell> {

        int id;
        Object3D nucleus;
        Object3D region;
        byte type;
        int layer, degreeLabelled, degreeUnlabelled;
        int border = -1;
        // info
        int info1 = -1;
        int info2 = -1;

        ArrayList<Cell> nei1 = null;
        ArrayList<Cell> nei2 = null;

        public int compareTo(Cell v) {
            float ratio = (float) degreeUnlabelled / degreeLabelled;
            float ratioV = (float) v.degreeUnlabelled / v.degreeLabelled;
            if (ratio > ratioV) {
                return 1;
            } else if (ratio < ratioV) {
                return -1;
            } else {
                return 0;
            }
        }

        @Override
        public String toString() {
            return "(" + nucleus.getValue() + ", " + region.getValue() + ", " + type + ")";
        }
        private int[] computeNeiTypeHisto(ArrayList<Cell> nei) {
            int[] res = new int[4];
            // al types , includin unlabelled
            for (Cell C : nei) {
                if (C.type ==ALPHA || C.type==BETA || C.type==DELTA) 
                {
                    res[C.type]++;
                }
            }
            return res;
        }

//        private int[] computeNeiTypeHisto(ArrayList<Cell> nei) {
//            int[] res = new int[4];
//            // al types , includin unlabelled
//            for (Cell C : nei) {
//                //if (C.type > 0) {
//                res[C.type]++;
//                //}
//            }
//            return res;
//        }

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
