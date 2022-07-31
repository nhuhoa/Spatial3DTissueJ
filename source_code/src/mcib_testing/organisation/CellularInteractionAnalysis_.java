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
//import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import mcib3d.geom.Object3D;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Vector3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;
import mcib_testing.Utils.CellSpatialAnalysis;
import stdlib.DirectedEdgeC;
import stdlib.MapUtils;
import stdlib.MathStat;
import stdlib.Queue;
import stdlib.Vertex;

/**
 *
 * @author tranhoa
 */
public class CellularInteractionAnalysis_ implements ij.plugin.PlugIn
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
//    public String dir = null, subdir = null;
    public boolean contactCal = true, neighborCal=true, hisContact=true, cellsNetwork=true;
    public boolean show = false;
    public String save_dir = "", input_dir="";
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
        ImagePlus temp = WindowManager.getImage(wList[wList.length-1]);
        if (null !=temp){
            WindowManager.setTempCurrentImage(temp);
            save_dir = IJ.getDirectory("image");
            input_dir = IJ.getDirectory("image");
        }
        GenericDialogPlus gd = new GenericDialogPlus("Cellular Interaction Analysis");
//        gd.addMessage("3D Tissue Analysis Framework");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : watershed image");
        gd.addMessage("Output: cell-cell contact.");
        gd.addMessage(" ");
        gd.addDirectoryField("Input_Dir: ", input_dir, 20);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        
//        gd.addStringField("Input Dir: ", input_dir);
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addChoice("Watershed_Image ", titles, titles[0]);
        gd.addCheckbox("cell-cell_contact_analysis", true);
        gd.addCheckbox("layer_contact", true);
        gd.addCheckbox("histogram_contact", true);
        gd.addCheckbox("cells_network_visualization", true);
//        gd.addCheckbox("nb clusters analysis", true);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        input_dir = gd.getNextString();
        save_dir = gd.getNextString();
        if (" ".equals(input_dir)) 
        {
                return;
        }
        if (" ".equals(save_dir)) 
        {
                return;
        }
        
        int[] index = new int[1];
        index[0] = gd.getNextChoiceIndex();
        contactCal  = gd.getNextBoolean();
        neighborCal = gd.getNextBoolean();
        hisContact = gd.getNextBoolean();
        cellsNetwork = gd.getNextBoolean();
        
        ImagePlus plus;
        plus = WindowManager.getImage(wList[index[0]]);
        img = ImageInt.wrap(plus);
        
        if (!input_dir.endsWith("/"))
            input_dir = input_dir + "/";
        IJ.log("Noted: input dir: where you keep list objects Regions.zip and Nuclei.zip");
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
//        WindowManager.setTempCurrentImage(plus);
//        dir = IJ.getDirectory("image");
        File fl = new File(save_dir+"raw_data_stat");
        if (!fl.exists()) 
        {
            if (fl.mkdir()) {
                save_dir = save_dir + "raw_data_stat/";
                IJ.log("Create a result folder : raw_data_stat");
            } else {
                save_dir = save_dir;
            }
        }
        else{
            save_dir = save_dir + "raw_data_stat/";
        }
   
        IJ.log("====================================================================");
        IJ.log("Analysing ...");
        IJ.log("Input watershed image : "+img.getTitle());
        IJ.log("Reading data from " + input_dir);
        IJ.log("Load list objects "+" Regions.zip"+ " and Nuclei.zip from "+input_dir);
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
//        testContact();
//        printCellContacts();
        if(contactCal)
        {
            computeAllContacts(); 
        }    
        if(neighborCal)
        {
            getAverageNbNeigbors();
        } 
        if(hisContact)
        {
            calculHistogramContact((byte)BETA, (byte)ALPHA);
        } 
        if(cellsNetwork)
        {
            drawGraphJ();
        }    
        IJ.log("Finished");
        IJ.log("====================================================================");
        IJ.selectWindow("Log");
        IJ.saveAs("Text", save_dir+ "Log_Cellular_Interaction_Analysis.txt");
        
    }
    private void testContact()
    {
        ResultsTable rtCont = new ResultsTable();
        ImageHandler draw = new ImageShort("TEST_", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        int nb = 0;
        for (Cell C : popCells) 
        {
            if(C.type==ALPHA)
            {
                C.nucleus.draw(draw, 50);
                int c = 0;
                rtCont.incrementCounter();
                rtCont.setValue("a_id", nb, C.id);
                for(Cell ne : C.nei1)
                {
                    if(ne.type==BETA)
                    {
                        ne.nucleus.draw(draw, ne.id);
                        c++;
                        rtCont.setValue("b_id_"+c, nb, ne.id);
                    }    
                }
                nb++;
                
            } 
            
        }
        rtCont.show("neighbourAlpha");
        draw.show();
        
    }        
    private void computeAllContacts()
    {
        ResultsTable rtContact = new ResultsTable();
        ResultsTable rtTheoretical = new ResultsTable();
//        float taa = 0, tbb = 0, tdd = 0, tab = 0, tad = 0, tbd = 0, tsum = 0;
        float aa = 0, bb = 0, dd = 0, ab = 0, ad = 0, bd = 0;
        float na = 0, nb = 0, nd = 0;
        int nbUnlabelled = 0;
        for (Cell C : popCells) 
        {
            if (C.type == ALPHA) 
            {
                na++;
                int[] ne = C.computeNei1TypeHisto();
                aa += ne[ALPHA-1];
                ab += ne[BETA-1];
                ad += ne[DELTA-1];
            } else if (C.type == BETA) {
                nb++;
                int[] ne = C.computeNei1TypeHisto();
                ab += ne[ALPHA-1];
                bb += ne[BETA-1];
                bd += ne[DELTA-1];
            } else if (C.type == DELTA) {
                nd++;
                int[] ne = C.computeNei1TypeHisto();
                ad += ne[ALPHA-1];
                bd += ne[BETA-1];
                dd += ne[DELTA-1];
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
        
        rtContact.setValue("total_cells", 0, popCells.size());
        rtContact.setValue("nb_unlabelled_cells", 0, nbUnlabelled);
        rtContact.setValue("Nb_alpha", 0, na);
        rtContact.setValue("Nb_beta", 0, nb);
        rtContact.setValue("Nb_delta", 0, nd);
        rtContact.setValue("total_AlphaBetaDelta", 0, sumnb);
        
        rtContact.setValue("proportion_alpha", 0, pa);
        rtContact.setValue("proportion_beta", 0, pb);
        rtContact.setValue("proportion_delta", 0, pd);
        
        rtContact.setValue("alpha_alpha", 0, aa);
        rtContact.setValue("beta_beta", 0, bb);
        rtContact.setValue("delta_delta", 0, dd);
        rtContact.setValue("alpha_beta", 0, ab);
        rtContact.setValue("alpha_delta", 0, ad);
        rtContact.setValue("beta_delta", 0, bd);
        rtContact.setValue("all_contacts", 0, sum);
        
        rtTheoretical.incrementCounter();
        rtTheoretical.setValue("total_cells", 0, popCells.size());
        rtTheoretical.setValue("nb_unlabelled_cells", 0, nbUnlabelled);
        rtTheoretical.setValue("Nb_alpha", 0, na);
        rtTheoretical.setValue("Nb_beta", 0, nb);
        rtTheoretical.setValue("Nb_delta", 0, nd);
        rtTheoretical.setValue("total_AlphaBetaDelta", 0, sumnb);
        rtTheoretical.setValue("proportion_alpha", 0, pa);
        rtTheoretical.setValue("proportion_beta", 0, pb);
        rtTheoretical.setValue("proportion_delta", 0, pd);
        rtTheoretical.setValue("p_alpha_alpha", 0, pa * pa);
        rtTheoretical.setValue("p_beta_beta", 0, pb * pb);
        rtTheoretical.setValue("p_delta_delta", 0, pd * pd);
        rtTheoretical.setValue("p_alpha_beta", 0, 2 * pa * pb);
        rtTheoretical.setValue("p_alpha_delta", 0, 2 * pa * pd);
        rtTheoretical.setValue("p_beta_delta", 0, 2 * pb * pd);
        
        rtContact.setValue("p_alpha_alpha", 0, aa / sum);
        rtContact.setValue("p_beta_beta", 0, bb / sum);
        rtContact.setValue("p_delta_delta", 0, dd / sum);
        rtContact.setValue("p_alpha_beta", 0, ab / sum);
        rtContact.setValue("p_alpha_delta", 0, ad / sum);
        rtContact.setValue("p_beta_delta", 0, bd / sum);
        rtContact.setValue("p_homotypic_contacts", 0, (aa + bb + dd) / sum);
        rtContact.setValue("p_heterotypic_contacts", 0, (ab + ad + bd) / sum);
        if(show)
        {
            rtContact.show("Observed_raw_data");
            rtTheoretical.show("Theoretical_Random_Model");
        }    
        
        try {
            rtContact.saveAs(save_dir+"AllContact.csv");
            rtTheoretical.saveAs(save_dir+"Theoretical_Random_Model.csv");
            IJ.log("Saving the output table AllContact.csv in the folder: "+save_dir);
            IJ.log("Saving the output table Theoretical_Random_Model.csv in the folder: "+save_dir);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
    }
    private void getAverageNbNeigbors()
    {
        ArrayList<Cell> popS = (ArrayList<Cell>) popCells.clone();
//        IJ.log("size of pop centroid before : "+ popS.size());
        Iterator<Cell> it = popS.iterator();
        ResultsTable rTotal = new ResultsTable();
        ResultsTable rA = new ResultsTable();
        ResultsTable rB = new ResultsTable();
        ResultsTable rD = new ResultsTable();
        int maxLabelledDegree = 0;
        Cell tmp = null;
//        int st = 0;
        int cA = 0, cB = 0, cD = 0, cT = 0;
        int nA = 0, nB = 0, nD = 0;
        int st = 0;
        rTotal.incrementCounter();
        rA.incrementCounter();
        rB.incrementCounter();
        rD.incrementCounter();
        while(it.hasNext())
        {
            Cell Ci = it.next();
            int cLabelled = 0;
            if (Ci.type != UNLABELLED) 
            {
                st++;
                if(Ci.type==ALPHA)
                {
                    nA++;
                }    
                if(Ci.type==BETA)
                {
                    nB++;
                }  
                if(Ci.type==DELTA)
                {
                    nD++;
                }  
                for (Cell C1 : Ci.nei1) 
                {
                    if ((C1.type != UNLABELLED)) 
                    {
                        cLabelled++;
                    } 
                }
                Ci.degreeLabelled = cLabelled;
            }
        }
        double a[] = new double[nA];
        double b[] = new double[nB];
        double d[] = new double[nD];
        double t[] = new double[st];
        for(Cell ci : popCells)
        {
            if(ci.type!=UNLABELLED)
            {
                t[cT] = ci.degreeLabelled;
                cT++;
                if(maxLabelledDegree < ci.degreeLabelled)
                {
                    maxLabelledDegree = ci.degreeLabelled;
                }    
            }
            if(ci.type==ALPHA)
            {
                a[cA] = ci.degreeLabelled;
                cA++;
                rA.setValue("id_"+ci.id+"_", 0, ci.degreeLabelled);
            }    
            if(ci.type==BETA)
            {
                b[cB] = ci.degreeLabelled;
                cB++;
                rB.setValue(ci.id+"_", 0, ci.degreeLabelled);
            }  
            if(ci.type==DELTA)
            {
                d[cD] = ci.degreeLabelled;
                cD++;
                rD.setValue(ci.id+"_", 0, ci.degreeLabelled);
            }  
        }
        rTotal.setValue("Total_Direct_Interaction", 0, MathStat.getSum(t)/2);
        rTotal.setValue("Max_Interaction_Cell", 0, maxLabelledDegree);
        rTotal.setValue("Average_Total", 0, MathStat.getMean(t));
        rTotal.setValue("Std_Total", 0, MathStat.getStdDev(t));
        rTotal.setValue("Average_Alpha", 0, MathStat.getMean(a));
        rTotal.setValue("Std_Alpha", 0, MathStat.getStdDev(a));
        rTotal.setValue("Average_Beta", 0, MathStat.getMean(b));
        rTotal.setValue("Std_Beta", 0, MathStat.getStdDev(b));
        rTotal.setValue("Average_Delta", 0, MathStat.getMean(d));
        rTotal.setValue("Std_Delta", 0, MathStat.getStdDev(d));
        
        
        rA.setValue("average_alpha", 0, MathStat.getMean(a));
        rA.setValue("Std_Alpha", 0, MathStat.getStdDev(a));
        rB.setValue("average_beta", 0, MathStat.getMean(b));
        rB.setValue("Std_Beta", 0, MathStat.getStdDev(b));
        rD.setValue("average_delta", 0, MathStat.getMean(d));
        rD.setValue("Std_Delta", 0, MathStat.getStdDev(d));
        if(show)
        {
            rA.show("Aver_A");
            rB.show("Aver_B");
            rD.show("Aver_D");
            rTotal.show("Aver_ABD");
        } 
        
        try {
            rTotal.saveAs(save_dir+"aver_direct_interaction_AlphaBetaDdelta.csv");
            rA.saveAs(save_dir+"aver_direct_interaction_Alpha.csv");
            rB.saveAs(save_dir+"aver_direct_interaction_Beta.csv");
            rD.saveAs(save_dir+"aver_direct_interaction_Delta.csv");
            IJ.log("Saving the output table aver_direct_interaction_AlphaBetaDdelta.csv ");
            IJ.log("Saving the output table aver_direct_interaction_Alpha ");
            IJ.log("Saving the output table aver_direct_interaction_Beta.csv ");
            IJ.log("Saving the output table aver_direct_interaction_Delta.csv ");
            IJ.log(" Into the folder: "+save_dir);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
//        rTotal.setValue("max_", 0, maxLabelledDegree);
//        rTotal.setValue("average_", 0, sumNeib/count);
//        rTotal.show("Nb Neighbors");
//        IJ.log("size of pop cell at the shell core : "+ popS.size());
//        IJ.log("max degree labelled: "+maxLabelledDegree);
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
    private void calculHistogramContact(byte sourceType, byte desType)
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
        String str = typeObserved[sourceType]+"_"+typeObserved[desType];
        int numPoints = 1000;
        int numRandomSamples = 100;
        double env = 0.05;
        CellSpatialAnalysis spa = new CellSpatialAnalysis(numPoints, numRandomSamples, env, save_dir);
        spa.process(sourceLst, desLst, vertexes, edges, str);
//        spa.process(desLst, vertexes, edges, "G");
        
//        ImageHandler drawLayer = new ImageShort("LAYER_DISTANCE_"+typeObserved[sourceType]+"_"+typeObserved[desType], img.sizeX, img.sizeY, img.sizeZ);
//        for (Vertex vd : desLst) 
//        {
//            Cell dx = popCells.get(cacheId.get(vd.getId()));
//            dx.region.draw(drawLayer, 1);
//        }
//        for (Vertex vs : sourceLst) 
//        {
//            Cell ds = popCells.get(cacheId.get(vs.getId()));
//            if(vs.getShortestDistance()!=-1)
//            {
//                ds.region.draw(drawLayer, (int)vs.getShortestDistance()+1);
//            }    
//        }
//        drawLayer.show();
//        IJ.selectWindow(drawLayer.getTitle());
//        IJ.saveAs(drawLayer.getImagePlus(), "Tiff", subdir+drawLayer.getTitle()+".tif");
//        IJ.log("Saving the output image "+drawLayer.getTitle() +" into the folder: "+subdir);
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
    private void printCellContacts()
    {
        String[] typeObs = {"del", "be","al"};
        ResultsTable rtContactInfo = new ResultsTable();
        
        int c = 0;
        for (Cell C : popCells) 
        {
            if (C.type != UNLABELLED) 
            {
                rtContactInfo.incrementCounter();
                int nbNei = 0;
                rtContactInfo.setValue("cell_id", c, C.id+"__"+typeObs[C.type-1]);
                if(C.nei1.size()>0)
                {
                    for(Cell ne : C.nei1)
                    {
                        if(ne.type != UNLABELLED)
                        {
                            nbNei++;
                            rtContactInfo.setValue("nei_"+nbNei, c, ne.id+"__"+typeObs[ne.type-1]);
                        }    
                        
                    }    
                } 
                c++;
                
            }
        }
        rtContactInfo.show("CellContactInfos");
        
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
    private void drawGraphJ()
    {
        ObjectCreator3D objs;
        objs = new ObjectCreator3D(img.sizeX, img.sizeY, img.sizeZ);
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
        
        ImagePlus graph = new ImagePlus("CellsNetwork", objs.getStack());
        graph.show();
        IJ.selectWindow(graph.getTitle());
        IJ.saveAs(graph, "Tiff", save_dir+"CellsNetwork.tif");
        IJ.log("Saving the output image CellsNetwork.tif into the folder: "+save_dir);
        
    }
    
}
