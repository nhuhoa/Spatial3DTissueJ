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
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Vector3D;
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
public class Protrusion_V3_ implements PlugIn
{
    
    private final int UNLABELLED = 0;
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
    String subdir = null;
    ResultsTable rtLayer = new ResultsTable();
    ResultsTable rtInf = new ResultsTable();
    ResultsTable rtContact = new ResultsTable();
    public boolean verbose = false, disp = false, disp1 = false, disp2=false, disp3=false;
    public double ratioMicro2Px=0, ratioZ2XY=0;
    int minDistance = 0, maxDistance = 35;
    int regVal = 0;
    public byte typeCellDef = 0;
    public byte typeInfDef = 0;
    public int op = 0, reg = 1;
    String[] typeDes = {"*NONE*","DELTA", "BETA","ALPHA"};
    String[] option = {"ALL_CELLS_PROTRUSION","EACH_CELL_PROTRUSION","ONE_CELL_PROTRUSION"};
    public String save_dir = " ";
    @Override
    public void run(String arg) 
    {
        if (Dialog())
        {
            IJ.log("Analysing ...");
            IJ.log("Reading input data included 2 results file Nuclei.zip & Regions.zip from " + save_dir);
            popRegions = new Objects3DPopulationGrid();
            popNuclei = new Objects3DPopulationGrid();
            popRegions.loadObjects(save_dir + "Regions.zip");
            popNuclei.loadObjects(save_dir + "Nuclei.zip");

            imgWat = wat.createSameDimensions();
            popRegions.draw(imgWat);
            //imgWat.show("Reg");
            imgSeg = wat.createSameDimensions();
            popNuclei.draw(imgSeg);
            //imgSeg.show("Nuc");

            // init cells 
            initCells2();
            IJ.log("====================================================================");
            IJ.log("Analysing ...");
            IJ.log("Nb cells in this tissue : " + popCells.size());
            IJ.log("Association 1 ...");
            computeAssoCells(imgWat, 0);
            IJ.log("------------------PROTRUSION-EXPANDED------------------");
//            byte typeInf = DELTA;
    //        
//            Calibration cal = dapi.getCalibration();
            Calibration cal = dapi.getImagePlus().getCalibration().copy();
            ratioMicro2Px = getRatioMicro2Pixel(cal);
            ratioZ2XY = getRatio_Z_XY(cal);
           
            IJ.log("Ratio micro to pixel: "+ratioMicro2Px);
            IJ.log("Ratio XY to Z: "+ratioZ2XY);
            if(op==0)
            {
                computeAllProtrusions(typeCellDef, typeInfDef, minDistance, maxDistance);
            }
            else if(op==1){
                computeContacts4EachCells(typeCellDef, typeInfDef, minDistance, maxDistance);
            }
            else{
                computeContacts4OneCells(typeCellDef, typeInfDef, minDistance, maxDistance, regVal);
            }
            IJ.log("Finished");
            IJ.log("====================================================================");
            IJ.selectWindow("Log");
            IJ.saveAs("Text", subdir+ "Log_protrusion_"+typeDes[typeInfDef]+"_"+option[op]+".txt");
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
                    save_dir = IJ.getDirectory("image");
                }
                
            }    
        }
        String[] typeObserved = {"BETA","ALPHA"};
        String[] typeInf = {"DELTA", "BETA","ALPHA"};
        
        String[] regionObs = {"NUCLEUS_REGION","CELL_REGION"};
        
        
        
        String unit = "microns";
        String unit1 = "pixel value";
        GenericDialogPlus gd = new GenericDialogPlus("3D Tissue Spatial Analysis - Protrusion");
//        gd.addMessage("3D Tissue Analysis Framework");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : input folder");
        gd.addMessage("Input : watershed image");
        gd.addMessage("Input : dapi filtered image");
        gd.addMessage("Output: cell-cell contacts.");
        gd.addMessage(" ");
//        gd.addStringField("Input Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("Watershed_Image :", titles, titles[0]);
        gd.addChoice("DAPI_Image :", titles, titles[1]);
        gd.addNumericField("Min_Distance: ", minDistance, 0 , 10, unit);
        gd.addNumericField("Max_Distance: ", maxDistance, 0 , 10, unit);
        gd.addNumericField("Region_Value (0-default)", regVal, 0, 10, unit1);
        gd.addChoice("Target_Cell : ", typeObserved, typeObserved[typeCellDef]);
        gd.addChoice("Protrusion_Cell : ", typeInf, typeInf[typeInfDef]);
        gd.addChoice("Range_Observe : ", option, option[op]);
        gd.addChoice("Region_Observe : ", regionObs, regionObs[reg]);
        gd.addCheckbox("Show_Protrusion_Image", verbose);
        gd.addCheckbox("Show_Distance_Map", disp1);
        gd.addCheckbox("Show_Observed_Cells", disp);
        gd.addCheckbox("Show_Layer_Distance", disp2);
        gd.addCheckbox("Show_Cells_Network", disp3);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return false;
        }
        
        int[] index = new int[2];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        minDistance = (int) gd.getNextNumber();
        maxDistance = (int) gd.getNextNumber();
        regVal = (int) gd.getNextNumber();
        typeCellDef = (byte)(gd.getNextChoiceIndex()+2);
        typeInfDef = (byte)(gd.getNextChoiceIndex()+1);
        op =  (int)(gd.getNextChoiceIndex());
        reg =  (int)(gd.getNextChoiceIndex());
        verbose = gd.getNextBoolean();
        disp1 = gd.getNextBoolean();
        disp = gd.getNextBoolean();
        disp2 = gd.getNextBoolean();
        disp3 = gd.getNextBoolean();
        if(gd.wasOKed())
        {
            IJ.log("\n------------------------------------------------------------------");
            IJ.log(" Min distance: "+minDistance
                   +" microns   &   max distance: "+maxDistance + " microns \n & type cell observed: "+typeDes[typeCellDef]
                   +" & type cell inf: "+typeDes[typeInfDef]+"\n&  range(0-all, 1-each, 2-one cell): "+op+"  & region (0-nuc, 1-reg): " + reg);
            if(op==2)
            {
                IJ.log("Input region value :"+regVal);
            }    
            IJ.log("\n------------------------------------------------------------------");
        } 
        ImagePlus plus;
        plus = WindowManager.getImage(wList[index[0]]);
        wat = ImageInt.wrap(plus);
        WindowManager.setTempCurrentImage(plus);
        ImagePlus plus1;
        plus1 = WindowManager.getImage(wList[index[1]]);
        dapi = ImageInt.wrap(plus1);
//        dir = IJ.getDirectory("image");
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        String subStr = "protrusion_"+typeDes[typeInfDef];
        File fl = new File(save_dir+subStr);
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                subdir = save_dir + subStr + "/";
                IJ.log("Create folder protrusion in the directory: "+save_dir);
            } else {
                subdir = save_dir;
            }
        }
        else{
            subdir = save_dir + subStr + "/";
        }
        return gd.wasOKed();
    }
    private void computeContacts4EachCells(byte typeObserved, byte typeInf, int minDistance, int maxDistance)
    {
        int st = 0;
        for (Cell C : popCells) 
        {
            if (C.type == typeInf) 
            {
                st++;
                IJ.log("\n----------------------"+st+"----------------------");
                IJ.log(typeDes[typeInf]+"  Cell Id: "+C.id);
                ArrayList<Cell> popTmp = new ArrayList<Cell>();
                popTmp.add(C);
                String desc = typeDes[typeInf]+"_"+C.id;
                computeContacts(typeObserved, typeInf, minDistance, maxDistance, popTmp, desc);
            }
        }
    }   
    private void computeContacts4OneCells(byte typeObserved, byte typeInf, int minDistance, int maxDistance, int regVal)
    {
        int st = 0;
        Cell C = region2Cell.get(regVal);
        if (C.type == typeInf) 
        {
            st++;
            IJ.log("\n----------------------"+st+"----------------------");
            IJ.log(typeDes[typeInf]+" Cell Id: "+C.id+"  region value: "+regVal);
            ArrayList<Cell> popTmp = new ArrayList<Cell>();
            popTmp.add(C);
            String desc = typeDes[typeInf]+"_"+C.id;
            computeContacts(typeObserved, typeInf, minDistance, maxDistance, popTmp, desc);
        }
        else{
            IJ.log("Input region value do not correspond to cell type "+typeDes[typeInf]);
            IJ.log("Please choose other value");
        }
    }  
    private void computeAllProtrusions(byte typeObserved, byte typeInf, int minDistance, int maxDistance)
    {
        ArrayList<Cell> popTmp = new ArrayList<Cell>();
        String idStr = " Id Cell Inf: ";
        for (Cell C : popCells) 
        {
            if (C.type == typeInf) 
            {
                idStr += "  "+C.id;
                popTmp.add(C);
            }
        }
        IJ.log(idStr);
        IJ.log("Nb cells inf: "+popTmp.size());
        String desc = "ALL_"+typeDes[typeInf]+"_cells";
        computeContacts(typeObserved, typeInf, minDistance, maxDistance, popTmp, desc);
    }        
    private void computeContacts(byte typeObserved, byte typeInf, int minDistance, int maxDistance, ArrayList<Cell> popInf, String desc)
    {
        ArrayList<Integer> idDelLs = new ArrayList<Integer>();
        ArrayList<Cell> lsCells = new ArrayList<Cell>();
        for (Cell C : popInf) 
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
        String title = "";
        title += desc;
        title += "_"+typeDes[typeObserved];
        title += "_Contact";
        
        //compute xaxis 
        int nbCells = popCells.size();
        int sizeDis = maxDistance - minDistance + 1;
//        ArrayUtil xAxis = new ArrayUtil(sizeDis);
//        ArrayUtil nbCellInf = new ArrayUtil(sizeDis);
        
        ResultsTable rtCellsPro = new ResultsTable();
        int count = 0;
        for (Cell C : popCells) 
        {
            if(C.type == typeObserved)
            {
                count++;
            } 
        }
        int nbInf = lsCells.size();
//        nbCellInf.putValue(0, lsCells.size());
        rtCellsPro.incrementCounter();
        rtCellsPro.setValue("MicroM", 0, 0);
        rtCellsPro.setValue("nbCells"+title, 0, nbInf);
        rtCellsPro.setValue("nbCells_"+typeDes[typeObserved], 0, count);
        rtCellsPro.setValue("percent_"+typeDes[typeObserved]+"_total_beta", 0,  nbInf * 100.0f / count);
        rtCellsPro.setValue("nbCells_ABD", 0, nbCells);
        rtCellsPro.setValue("percent_"+typeDes[typeObserved]+"_ABD", 0, nbInf * 100.0f / nbCells);
        int nb = 1;
        ImageFloat edm = getEDT(popInf);
        ImageFloat edm1 = edm.duplicate();
        edm1.setScale(dapi);
        if(disp1)
        {
            edm1.setTitle("EDM");
            edm1.show();
            IJ.selectWindow(edm1.getTitle());
            IJ.saveAs("Tiff", subdir+edm1.getTitle()+".tif");
        }    
        
        for(int d=minDistance; d<=maxDistance; d++)
        {
            IJ.log("------------------------------------------------------------------");
            IJ.log("Protrusion with distance: "+d+" micro m");
            ArrayList<Cell> popObs = new ArrayList<Cell>();
            ImageByte res = protrusion(d, typeInf, edm);
            if(verbose)
            {
                res.show();
                IJ.selectWindow(res.getTitle());
                IJ.saveAs("Tiff", subdir+res.getTitle()+".tif");
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
                        if ((C1 != null) && C1.type==typeObserved && !popObs.contains(C1)) 
                        {
                            popObs.add(C1);
                        }
                    }
                }
            }
            if(disp2)
            {
                calculHistogramContact(typeObserved, typeInf, popInf, popObs, d);
            }    
            if(disp3)
            {
                drawGraphJ(popInf, popObs, d);
            }    
            int nbInfAfter = popObs.size();
            IJ.log("  Distance:  "+d+ "   inscreased: "+ popObs.size());
            rtCellsPro.incrementCounter();
            rtCellsPro.setValue("MicroM", nb, d);
            rtCellsPro.setValue("nbCells"+title, nb, nbInfAfter);
            rtCellsPro.setValue("nbCells_"+typeDes[typeObserved], nb, count);
            rtCellsPro.setValue("percent_"+typeDes[typeObserved]+"_total_beta", nb, (nbInfAfter)* 100.0f/count);
            rtCellsPro.setValue("nbCells_ABD", nb, nbCells);
            rtCellsPro.setValue("percent_"+typeDes[typeObserved]+"_ABD", nb, (nbInfAfter)* 100.0f/nbCells);
            nb++;
        }
        try {
            rtCellsPro.saveAs(subdir+title+".csv");
            IJ.log("Saving the output table "+title+".csv in the "+subdir);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
        
    } 
    private void calculHistogramContact(byte sourceType, byte desType, ArrayList<Cell> popInf, ArrayList<Cell> popObs, int distance)
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
        for(Cell C : popInf)
        {
            if (!C.connected) 
            {
                //Vertex vC = new Vertex(C.id - 1, "Cell" + C.id, C.type);
                cacheId.put(idx, C.id);
                Vertex vC = new Vertex(idx, "Cell" + C.id, C.type);
                vertexes.add(vC);
                idx++;
            }
        }
        for(Cell C : popObs)
        {
            if (!C.connected) 
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
        IJ.log("popInf size "+popInf.size()+"  popObs size "+popObs.size());
        for(Cell C : popInf)
        {
            for(Cell Ci : popObs)
            {
                if(!C.nei1.contains(Ci))
                {
                    int keyC = (int) MapUtils.getKeyFromValue(cacheId, C.id);
                    int keyC1 = (int) MapUtils.getKeyFromValue(cacheId, Ci.id);
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
                        DirectedEdgeC e = new DirectedEdgeC(source, des, weight);  
                        DirectedEdgeC e1 = new DirectedEdgeC(des, source, weight);  
                        edges.add(e);
                        edges.add(e1);
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
        String str = distance +"_"+ typeDes[sourceType]+"_"+typeDes[desType];
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(vertexes);
        pop.setMask(vertexes, edges);
        ArrayUtil arr = pop.computeMinDistance(sourceLst, desLst);
//        spa.process(desLst, vertexes, edges, "G");
        
        ImageHandler drawLayer = new ImageShort("LAYER_DISTANCE_"+str, imgWat.sizeX, imgWat.sizeY, imgWat.sizeZ);
        
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
        drawLayer.closeImagePlus();
        IJ.log("Saving the output image "+drawLayer.getTitle() +".tif into the folder: "+subdir);
    }
    private void drawGraphJ(ArrayList<Cell> popInf, ArrayList<Cell> popObs, int distanceProtrusion)
    {
        ObjectCreator3D objs;
        objs = new ObjectCreator3D(imgWat.sizeX, imgWat.sizeY, imgWat.sizeZ);
        int rad = 6;
        for (Cell C : popCells) 
        {
            Vector3D cen = C.nucleus.getCenterAsVector();
            int rad1 = rad + (int)C.degreeLabelled/2;
            objs.createSphere(cen.getX(), cen.getY(), cen.getZ(), rad1, getCellColor(C), false);
        }
        for (Cell C : popCells) 
        {
            if (C.type != UNLABELLED) 
            {
                Object3D nuc = C.nucleus;
                Vector3D cen = nuc.getCenterAsVector();
                
                for (Cell C1 : C.nei1) 
                {
                    if (C1.type != UNLABELLED) 
                    {
                        int r1 = C.id - 1;
                        int c1 = C1.id - 1;
                        Object3D nuc1 = C1.nucleus;
                        Vector3D cen1 = nuc1.getCenterAsVector();
                        int distance = 2;
                        int val = getContactColor(C.type, C1.type);
                        objs.createLine(cen, cen1, val, distance);
                    }

                }
            }
        }
        
        for(Cell C : popInf)
        {
            Vector3D cen = C.nucleus.getCenterAsVector();
            for(Cell Ci : popObs)
            {
                if(!C.nei1.contains(Ci))
                {
                    int r1 = C.id - 1;
                    int c1 = Ci.id - 1;
                    Object3D nuc1 = Ci.nucleus;
                    Vector3D cen1 = nuc1.getCenterAsVector();
                    int distance = 2;
                    int val = getContactColor(C.type, Ci.type);
                    objs.createLine(cen, cen1, val, distance);
                }    
            }    
        }    
        
        ImagePlus graph = new ImagePlus("CellsNetwork_"+distanceProtrusion, objs.getStack());
        graph.show();
        IJ.selectWindow(graph.getTitle());
        IJ.saveAs(graph, "Tiff", subdir+"CellsNetwork.tif");
        graph.close();
        IJ.log("Saving the output image CellsNetwork.tif into the folder: "+subdir);
        
    }
    private int getCellColor(Cell C)
    {
        int val = 0;
        if (C.type != UNLABELLED) 
        {
            if(C.type==DELTA)
            {
                val = 50;
            }
            else if(C.type==BETA)
            {
                val = 100;
            }
            else
            {
                val = 150;
            }
        }    
        return val;
    }
    private int getContactColor(byte cType1, byte cType2)
    {
        int val = 0; 
        if(cType1==DELTA && cType2==DELTA)
        {
            val = 50;
        }
        else if(cType1==BETA && cType2==BETA)
        {
            val = 100;
        }
        else if(cType1==ALPHA && cType2==ALPHA)
        {
            val = 150;
        }
        else
        {
            val = 200;
        }
        return val;
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
        ImageHandler drawConnected = new ImageShort("connected_component", imgWat.sizeX, imgWat.sizeY, imgWat.sizeZ);
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
//                IJ.log("Size of max cluster is: "+lsCells.size());
            }    
            
        }
//        drawConnected.show();
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
    private ImageFloat getEDT(ArrayList<Cell> popTmp)
    {
        ImageHandler label = new ImageShort("segmented_objs_", wat.sizeX, wat.sizeY, wat.sizeZ);
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
            IJ.selectWindow(label.getTitle());
            IJ.saveAs("Tiff", subdir+label.getTitle()+".tif");
        }
//        label.setCalibration(dapi.getCalibration());
        int nbCPUs = ThreadUtil.getNbCpus();
        ImageFloat edm = null;
//        edm = EDT.run(label, 0, (float)dapi.getScaleXY(), (float)dapi.getScaleZ(), true, nbCPUs);
        edm = EDT.run(label, 0, (float)wat.getScaleXY(), (float)wat.getScaleZ(), true, nbCPUs);
        return edm;
    }        
    public ImageByte protrusion(double microDistance, byte typeInf, ImageFloat edm)
    {
        try {
            float rDilate = (float) (microDistance * ratioMicro2Px);
//            IJ.log("rDilate"+rDilate);
//            float radiusZ = (float) (rDilate*ratioZ2XY);
//            ImageFloat edm = EDT.run(label, 0, 1, rDilate / radiusZ, true, nbCPUs);
            ImageByte temp = edm.threshold(rDilate, true, false);
//            edm.flush();
//            edm = null;
            temp.setOffset(wat);
            temp.setScale(wat);
            temp.setTitle("dilated_"+microDistance);
            return temp;
        } catch (Exception e) {
            exceptionPrinter.print(e, null, true);
        }
        return null;
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
            
    private ArrayList<Integer>[] computeAsso(ImageHandler img, int BorderValue) 
    {
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
                        for (int j = 0; j < tab.size(); j++) {
                            int val = (int) tab.getValue(j);
                            if (!tab.hasOnlyValuesInt(nei[val])) {
                                for (int i = 0; i < tab.size(); i++) {
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
