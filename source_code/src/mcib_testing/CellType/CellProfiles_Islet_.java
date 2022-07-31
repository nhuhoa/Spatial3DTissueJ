/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
//import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import mcib3d.geom.Object3D;
//import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
//import mcib3d.image3d.ImageByte;
//import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
//import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.regionGrowing.Watershed3DVoronoi;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;
//import DecimalFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Vector3D;
//import stdlib.MathStat;


/**
 *
 * @author hoatran
 */

public class CellProfiles_Islet_ implements ij.plugin.PlugIn {
    private final int UNLABELED = 1;
    private final int MARKER = 255;
    int MARKER_VAL = 2;
    double wat_radius = 3;
//    private Objects3DPopulation popRegions = null;
    HashMap<Integer, Cell> region2Cell = null;
    Objects3DPopulation popRegions = null;
    HashMap<Integer, Cell> nucleus2Cell=null;
    Objects3DPopulation popNuclei = null;
    ArrayList<Cell> popCells = null;
    String dir = null;
    private double ratioMarker = 0.2;
//    private float min = 0, maxInside = 4, maxOutside = 2;
//    private final int IN_OUTSIDE = 0, OUTSIDE_NUCLEUS = 1, INSIDE_NUCLEUS = 2, WATERSHED_REGION = 3;
    public boolean showUnlabelled = true;
    
    String save_dir = " ";
    
    //        UNLABELED = 0, ALPHA = 3, BETA = 2, DELTA = 1
    
    
    public void run(String arg) 
    {
//        String[] regionObs = {"INSIDE_OUTSIDE NUCLEUS", "OUTSIDE NUCLEUS", "INSIDE NUCLEUS","WATERSHED REGION"};
          // 2, 3, 4, 5 pixel val
//        String[] markerObs = {"MYC", "Ki67"};
//        int marker_type = 0;
        
        int[] wList = WindowManager.getIDList();
//        if (wList==null) {
//            IJ.error("No images are open.");
//            return;
//        }
        if(wList==null || wList.length<3)
        {
            IJ.error("At least segmented nuclei, and filtered markers image should be open!");
            return;
        }  
        
        String[] titles = new String[wList.length+1];
        titles[0] = "*None*";
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i+1] = imp!=null?imp.getTitle():"";
        }
//        if(wList.length==2){
//            titles = new String[wList.length+1];
//            for (int i=0; i<wList.length; i++) {
//                ImagePlus imp = WindowManager.getImage(wList[i]);
//                titles[i] = imp!=null?imp.getTitle():"";
//            }
//            titles[wList.length] = none;
//        }else{
//            titles = new String[wList.length];
//            for (int i=0; i<wList.length; i++) {
//                ImagePlus imp = WindowManager.getImage(wList[i]);
//                titles[i] = imp!=null?imp.getTitle():"";
//            }
//            titles[wList.length] = none;
//        }
        
//        boolean save = true; 
        GenericDialog gd = new GenericDialog("Cell Type Detection");
//        gd.addMessage("3D Tissue Spatial Analysis");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        "organization applied to \nthe study of the islets of Langerhans");
//        gd.addMessage("Input: nuclei segmented, watershed, labelF");
//        gd.addMessage("Output: Cell Type Image.");
        gd.addMessage("    Spatial3DTissueJ  ");
        gd.addStringField("Save Dir: ", save_dir);
        gd.addChoice("SEG NUC image :  ", titles, titles[0]);
//        gd.addChoice("Watershed: ", titles, titles[1]);

//        gd.addChoice("BIN DELTA: ", titles, titles[1]);
        gd.addChoice("FILTERED DELTA: ", titles, titles[0]);
        
//        gd.addChoice("BIN BETA: ", titles, titles[0]);
        gd.addChoice("FILTERED BETA: ", titles, titles[0]);
        
//        gd.addChoice("BIN ALPHA: ", titles, titles[0]);
        gd.addChoice("FILTERED ALPHA: ", titles, titles[0]);
        
        
//        gd.addChoice("Region Observed : ", regionObs, regionObs[reg]);
//        gd.addChoice("Marker Observed : ", markerObs, markerObs[marker_type]);
        gd.addNumericField("Percent_Marker : ", (double)ratioMarker, 1);
        gd.addNumericField("Watershed_Cell_Radius : ", wat_radius, 0);
//        gd.addNumericField("Max_Distance_Inside : ", maxInside, 1);
//        gd.addNumericField("Max_Distance_Outside : ", maxOutside, 3);
//        gd.addCheckbox("Save the output into folder", true);
//        gd.addCheckbox("Show unlabelled nucleus", true);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        
        save_dir = gd.getNextString();
        
        if (" ".equals(save_dir)) 
        {
                return;
        }

        File wdir = new File(save_dir);
        
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
//                IJ.showMessage("Working directory error");
//                return;
        }

        int[] index = new int[4];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        index[3] = gd.getNextChoiceIndex();
        
        
//        marker_type =  (int)(gd.getNextChoiceIndex());
        ratioMarker = (double) gd.getNextNumber();
        wat_radius = (double) gd.getNextNumber();
//        min = (float) gd.getNextNumber();
//        maxInside = (float) gd.getNextNumber();
//        maxOutside = (float) gd.getNextNumber();
//        save  = gd.getNextBoolean();
//        showUnlabelled  = gd.getNextBoolean();
        IJ.log("Ratio marker to assign a cell type: "+ratioMarker);
//        IJ.log("Detecting cell type: "+markerObs[marker_type]);
        ImageInt imgSeg = null;
        if(index[0]!=0)
        {
            imgSeg = ImageInt.wrap(WindowManager.getImage(wList[index[0]-1]));
            IJ.log("Seg: " + imgSeg.getTitle());
        }else{
            IJ.log("Need to input segmented image results, input error...");
        }
        ImageInt imgWat = null;
//        if(none.equals(wList[index[1]])){ //(wList.length==2) || 
//            IJ.log("Dont exist wat img. compute it later");
//        }else{
//            imgWat = ImageInt.wrap(WindowManager.getImage(wList[index[1]]));
//            IJ.log("Wat: " + imgWat.getTitle());
//        }
        
        
        
        
//        ArrayList<ImageInt> lsbin = new ArrayList<ImageInt>();
        ArrayList<ImageInt> lsraw = new ArrayList<ImageInt>();
        
        String[] markerTypes = {"DELTA", "BETA","ALPHA"};
        ArrayList<String> marker_titles = new ArrayList<String>();
        for(int i=1; i<=3; i=i+1){
            if(index[i]!=0)
            {   marker_titles.add(markerTypes[i-1]);
                IJ.log("idx: "+index[i]+" img: "+wList[index[i]-1]);
    //            IJ.log("idx2: "+index[2]+" img: "+wList[index[2]]);
    //            IJ.log("idx3: "+index[10]+" img: "+wList[index[3]]);
    //            ImageInt imgLabel = ImageInt.wrap(WindowManager.getImage(wList[index[i]]));
                ImageInt imgRawLabel = ImageInt.wrap(WindowManager.getImage(wList[index[i]-1]));
    //            imgLabel.setTitle("BIN_"+markerTypes[count]);
                IJ.log("Input filtered image: " + imgRawLabel.getTitle());
    //            lsbin.add(imgLabel);
                lsraw.add(imgRawLabel);
                
            }
            
        }   
        
//        IJ.log("Nb images bin is: " + lsbin.size());
        IJ.log("Nb images raw is: " + lsraw.size());
//        IJ.log("Input data: wat: "+imgWat.getTitle()+" seg: "+imgSeg.getTitle()+" label: "+imgLabel.getTitle());
//        IJ.log("Ratio marker: "+ ratioMarker);
//        IJ.log("Region observed: "+regionObs[reg]);
//        MARKER_VAL = marker_type + 2;
//        IJ.log("Observed marker: " + markerObs[marker_type] + "  MARKER_VAL: " + MARKER_VAL);
        IJ.log("Save dir is: " + save_dir);
//        IJ.log("Observe distance inside of nuc region:   min  :  "+min+"    max: "+maxInside);
//        IJ.log("Observe distance outside of nuc region:   min  :  "+min+"    max: "+maxOutside);
//        WindowManager.setTempCurrentImage(imgSeg.getImagePlus());
//        dir = IJ.getDirectory("image");
//        IJ.log("Initialization Tissue Analysis ...");
        analysis_v2(imgSeg, imgWat, lsraw, marker_titles, save_dir);
        
        
    }
    public void analysis_v2(ImageInt imgSeg, ImageInt imgWat,
                         ArrayList<ImageInt> lsraw, //ArrayList<ImageInt> lsbin, 
                         ArrayList<String> marker_titles, String save_dir){
        double background_threshold = 0;
        
//        ImageInt imgLabel;
//        HashMap<Integer, Cell> region2Cell=null;
        
        String image_fn = imgSeg.getTitle();
        Pattern p = Pattern.compile("_SEG");
        Matcher m = p.matcher(image_fn); 
        image_fn = m.replaceAll("");
        IJ.log(image_fn);
        IJ.log("Generate watershed cell segmentation image using radius from nucleus: " + wat_radius);
//        boolean save = false;
//        // Change the save flag using a random function, big dataset, save several data only
//        int random_int = (int)(Math.random() * 15);
//        if(random_int==8){
//            save = true;
//        }
        boolean save = true;
        if(imgWat==null){
            float radMax = (float)wat_radius;
            IJ.log("Watershed radius: "+radMax);
            if(imgSeg.getMax() > background_threshold){
                imgWat = WatershedVoronoi(imgSeg.getImagePlus(), image_fn, save_dir, radMax, false, save);
                imgWat.setTitle(image_fn+"_WAT");
                imgWat.show();
                IJ.selectWindow(image_fn+"_WAT");
                IJ.saveAs("Tiff", save_dir + image_fn+"_WAT" + ".tif");
            } else{
                IJ.log("Background image, skip this step");
            }
        }
        
 
        initCells(imgSeg, imgWat, image_fn, save_dir);
        computeAssoCells(imgWat, 0);
////        TO DO: wat, cell networks
        IJ.log("Detecting Cell Type ...");
        IJ.log(image_fn);
//        IJ.log("Nb images bin is: " + lsbin.size());
//        IJ.log("Nb filtered images is: " + lsraw.size());
        computeCellTypeV2(lsraw, marker_titles, save_dir, image_fn); 
        IJ.log("Completed");
        IJ.log("====================================================================\n");
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
    private ImageInt WatershedVoronoi(ImagePlus seedPlus, String image_fn, String save_dir,
                                  float radMax, boolean showLines, boolean save) {
        IJ.log("Computing watershed, cell region detection");
//        save_dir=save_dir+"SEG_WAT/";
        boolean showEDT = false;
//        float radMax = 0;
//        long t = System.currentTimeMillis();
        Watershed3DVoronoi watershed3DVoronoi = new Watershed3DVoronoi(ImageInt.wrap(seedPlus), radMax);
        watershed3DVoronoi.setLabelSeeds(false);
        ImageInt wat = watershed3DVoronoi.getVoronoiZones(showEDT);
        String save_fn = image_fn + "_SEG_WAT";
        wat.setTitle(save_fn);
        
        if (showLines) watershed3DVoronoi.getVoronoiLines(true).show("VoronoiLines");
        if(save)
        {
            wat.show();
            IJ.selectWindow(wat.getTitle());
            IJ.run(wat.getImagePlus(), "3-3-2 RGB", "");
            IJ.saveAs("Tiff", save_dir + save_fn + ".tif");
            IJ.log("Save the result as " + save_fn);
//            IJ.log("Save the result into the folder : " + save_dir);
            
        }  
//        IJ.log("Finished in " + (System.currentTimeMillis() - t) + " ms.");
        return wat;
    }
    private void initCells(ImageInt nucLabel, ImageInt regionLabel, String image_fn, String save_dir) 
    {
        popNuclei = new Objects3DPopulation(nucLabel);
        popRegions = new Objects3DPopulation(regionLabel, 1); // exclude value 1 used by borders

        popCells = new ArrayList<Cell>(popRegions.getNbObjects());

        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer, Cell>(popNuclei.getNbObjects());

        // get nucleus label for each region
        int c = 1;
        //int count = 0;
        for (Object3D region : popRegions.getObjectsList()) 
        {
            int nuc = (int) region.getPixModeNonZero(nucLabel);
            //IJ.log("nuc " + nuc);
            if(nuc==-1){ continue;}
            else
            {
                Cell cell = new Cell();
                cell.region = region;
                cell.nucleus = popNuclei.getObjectByValue(nuc);
                popCells.add(cell);
                cell.id = c++;
                region2Cell.put(region.getValue(), cell);
                nucleus2Cell.put(nuc, cell);
            }
            
                  
        }
        IJ.log("Number of cells is: " + popCells.size());
        
//        popNuclei.saveObjects(save_dir+image_fn+"_SEG.zip");
//        popRegions.saveObjects(save_dir+image_fn+"_WAT.zip");
        
    }
    
    
    private void initCells(ImageInt nucLabel,
                            Objects3DPopulation popRegions,
                            Objects3DPopulation popNuclei,
                            HashMap<Integer, Cell> nucleus2Cell,
                            ArrayList<Cell> popCells) 
    {
        popNuclei = new Objects3DPopulation(nucLabel);
//        popRegions = new Objects3DPopulation(regionLabel, 1); // exclude value 1 used by borders

//        popCells = new ArrayList<Cell>(popRegions.getNbObjects());
        popCells = new ArrayList<Cell>(popNuclei.getNbObjects());

//        region2Cell = new HashMap<Integer, Cell>(popRegions.getNbObjects());
        nucleus2Cell = new HashMap<Integer, Cell>(popNuclei.getNbObjects());

        // get nucleus label for each region
        int c = 1;
        //int count = 0;
        for (Object3D region : popNuclei.getObjectsList()) 
        {
            int nuc = (int) region.getPixModeNonZero(nucLabel);
            //IJ.log("nuc " + nuc);
            if(nuc==-1){ continue;}
            else
            {
                Cell cell = new Cell();
//                cell.region = region;
//                cell.nucleus = popNuclei.getObjectByValue(nuc);
                cell.nucleus = region;
                popCells.add(cell);
                cell.id = c++;
//                region2Cell.put(region.getValue(), cell);
//                nucleus2Cell.put(nuc, cell);
            }
            
                  
        }
        IJ.log("Number of cells is: " + popCells.size());
    }
    private void drawCellTypes(boolean nuc, ImageHandler imgSeg, ImageHandler imgRawLabel,
                                        String image_fn, String save_dir, boolean save) {
        
        ResultsTable nodes = new ResultsTable();
        ResultsTable edges = new ResultsTable();
        if(save){
            IJ.log("Save cell type image to file");
            String image_fn_nuc = image_fn + "_CELLTYPE_NUC";
            String image_fn_wat = image_fn + "_CELLTYPE_WAT";
//            String image_fn_intensity = image_fn + "_INTENSITY";
            ImageHandler draw_nuc = new ImageShort(image_fn_nuc, imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
            ImageHandler draw_cell = new ImageShort(image_fn_wat, imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
//            ImageHandler draw_intensity = new ImageShort(image_fn_intensity, imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
            for (Cell C : popCells) {
//                if (nuc) {
                    C.nucleus.draw(draw_nuc, C.type);
//                } else {
                    C.region.draw(draw_cell, C.type);
//                }
//                  if(C.type!=UNLABELED){
//                      C.nucleus.draw(draw_intensity, (int)C.nucleus.getPixMeanValue(imgRawLabel));
//                  }
            }
            draw_nuc.show();
            IJ.selectWindow(image_fn_nuc);
            IJ.saveAs("Tiff", save_dir + image_fn_nuc + ".tif");
            draw_cell.show();
            IJ.selectWindow(image_fn_wat);
            IJ.saveAs("Tiff", save_dir + image_fn_wat + ".tif");
//            draw_intensity.show();
//            IJ.selectWindow(image_fn_intensity);
//            IJ.saveAs("Tiff", save_dir + image_fn_intensity + ".tif");
            
        }
        
//        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
        int c = 0;
        int k = 0;
        for (Cell C : popCells) {
            String[] markerObs2 = {"UNLABELLED", "MYC", "Ki67"};//"UNLABELLED"=1, "MYC"=2, "Ki67"=3
            // First extract all cells to a csv file
            nodes.incrementCounter();
//            nodes.setValue("cell_name", c, C.nucleus.getName());
            nodes.setValue("cell_type_id", c, C.nucleus.getType()); // CD3, CD4, CD8, CD57
            nodes.setValue("cell_type", c, markerObs2[(int)(C.nucleus.getType()-1)]); // CD3, CD4, CD8, CD57
            nodes.setValue("cell_id", c, C.nucleus.getValue()); //intensity id
            nodes.setValue("x", c, Math.round(C.nucleus.getCenterX()));
            nodes.setValue("y", c, Math.round(C.nucleus.getCenterY()));
            nodes.setValue("vol", c, Math.round(C.nucleus.getVolumePixels()));
            nodes.setValue("mean_intensity", c, Math.round(C.nucleus.getPixMeanValue(imgRawLabel)));
            nodes.setValue("median_intensity", c, Math.round(C.nucleus.getPixMedianValue(imgRawLabel)));
            c = c + 1;
            
            // Second extract all edges of cells interaction to a csv file
            for (Cell C1 : C.nei1) 
            {
                String cell_interaction_type = "homo_contact";
                if(C.type!=C1.type){
                    cell_interaction_type = "hetero_contact";
                }
                edges.incrementCounter();
                edges.setValue("from", k, C.nucleus.getValue());
                edges.setValue("to", k, C1.nucleus.getValue()); 
//                edges.setValue("from_cell_name", k, C.nucleus.getName());
//                edges.setValue("to_cell_name", k, C1.nucleus.getName()); 
                edges.setValue("cell_interaction", k, cell_interaction_type);
                k = k + 1;
            }
        }
//        IJ.selectWindow(image_fn);
//        IJ.saveAs("Tiff", save_dir + image_fn + ".tif");
        
        try {
            nodes.saveAs(save_dir+image_fn+"_nodes.csv");
            edges.saveAs(save_dir+image_fn+"_edges.csv");
        } catch (IOException ex) {
            IJ.log("Have the problem of save nodes data");
        }
        
        
        
    }
    private void computeCellTypeV2(ArrayList<ImageInt> lsraw, ArrayList<String> marker_titles, //ArrayList<ImageInt> lsbin, 
                                   String save_dir, String image_fn) 
    {
        
        IJ.log("Computing cell type...");
        ResultsTable nodes = new ResultsTable();
        ResultsTable edges = new ResultsTable();
        int c = 0;
        int k = 0;
        IJ.log("Nb cells is: "+ popCells.size());
        for(int m=0; m<marker_titles.size(); m++){
                IJ.log(marker_titles.get(m));
        }  
        DecimalFormat df = new DecimalFormat("#.##");
        for (Cell C : popCells) {
            if(C.region==null) continue;
            nodes.incrementCounter();
//            IJ.log("NUC: " + C.nucleus.getValue());
            nodes.setValue("cell_id", c, C.nucleus.getValue()); //index
            nodes.setValue("x", c, Math.round(C.nucleus.getCenterX()));
            nodes.setValue("y", c, Math.round(C.nucleus.getCenterY()));
            nodes.setValue("pixvol", c, Math.round(C.nucleus.getVolumePixels()));
 
            for(int m=0; m<marker_titles.size(); m++){
//                ArrayUtil list = C.region.listValues(lsbin.get(m));
//                int nbM = list.countValueAbove(0);
//                double vol = C.nucleus.getVolumePixels();
//                double coverage = nbM /vol;
                nodes.setValue(marker_titles.get(m)+"_mean_intensity_nuc", c, Math.round(C.nucleus.getPixMeanValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_median_intensity_nuc", c, Math.round(C.nucleus.getPixMedianValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_mean_intensity_cellzone", c, Math.round(C.region.getPixMeanValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_median_intensity_cellzone", c, Math.round(C.region.getPixMedianValue(lsraw.get(m))));
//                nodes.setValue(markerTypes[m]+"_mean_coverage_nuc", c, df.format(nbM /C.nucleus.getVolumePixels()));
//                nodes.setValue(markerTypes[m]+"_mean_coverage_cellzone", c, df.format(nbM /C.region.getVolumePixels()));
            }    
            
            c = c + 1;
            
            // Second extract all edges of cells interaction to a csv file
            for (Cell C1 : C.nei1) 
            {
                edges.incrementCounter();
                edges.setValue("from", k, C.nucleus.getValue());
                edges.setValue("to", k, C1.nucleus.getValue()); 
                k = k + 1;
            }
        }
//        IJ.selectWindow(image_fn);
//        IJ.saveAs("Tiff", save_dir + image_fn + ".tif");
        
        try {
            String node_fn = save_dir+image_fn+"_cellNodes.csv";
            String edge_fn = save_dir+image_fn+"_cellsEdges.csv";
            IJ.log(node_fn);
            IJ.log(edge_fn);
            nodes.saveAs(node_fn);
            edges.saveAs(edge_fn);
        } catch (IOException ex) {
            IJ.log("Have the problem of save cell profiles info results");
        }
        
        
    }
    /**
     * Description: compute the cell type by counting all stainning inside nucleus
     * @param label
     * @param ratio
     */
//    private void computeCellType(ImageInt label, String save_dir, String image_fn) 
//    {
//        int cM=0, cUnlabelled=0;
//        ResultsTable cellM = new ResultsTable();
//        if(label==null){
//            IJ.log("Dont exist marker for this image");
//            for (Cell C : popCells) 
//            {
//                C.type = (byte)UNLABELED;
//                C.nucleus.setType(C.type);
//                C.nucleus.setName("Nuc" + C.id);
//                C.nucleus.setValue(C.id);
//            } 
//            cUnlabelled = popCells.size();
//        } else{
//            for (Cell C : popCells) 
//            {
//                Object3D nucleus = C.nucleus;
//                double vol = nucleus.getVolumeUnit();
//    //                Object3D region = C.region;
//    //                int reVal = region.getValue();
//                int reVal = nucleus.getValue();
//                
//                if(nucleus==null) continue;
//                ArrayUtil list = nucleus.listValues(label);
//
//                int nbM=0; int nbU=0;
//                for(int i=0;i<list.size();i++)
//                {
//                    if(list.getValue(i)==MARKER) nbM++;             
//                }
//                C.type = (byte)UNLABELED;
//                if(nbM /vol >= ratioMarker)
//                {
//                    C.type = (byte)MARKER_VAL;
//                } 
//
//                
//                if(C.type==UNLABELED){cUnlabelled++;}
//                if(C.type==MARKER_VAL){
//                    cM++;
//                }
//
//                nucleus.setType(C.type);
//                nucleus.setName("Nuc" + C.id);
//                nucleus.setValue(C.id);
//    //                region.setType(C.type);
//    //                region.setName("Reg" + C.id);
//    //                region.setValue(C.id);
//            }
//
//            
//        }
//        cellM.incrementCounter();
//        cellM.setValue("cell_name", 0, image_fn);
//        cellM.setValue("marker_cells", 0, cM);
//        cellM.setValue("unlabelled_cells", 0, cUnlabelled);
//        cellM.setValue("total_cells", 0, popCells.size());
//        cellM.setValue("image_width_X", 0, label.sizeX);
//        cellM.setValue("image_height_Y", 0, label.sizeY);
//        IJ.log("Nb of each cells type: " + cM + "   Unlabelled: "+cUnlabelled);
//        try {
//            cellM.saveAs(save_dir+image_fn+"_celltype.csv");
//        } catch (IOException ex) {
//            IJ.log("Have the problem of saving data");
//        }
////        ImageHandler unLabelledCell;
////        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);
////        for (Cell C : popCells) {
////            Object3D nucleus = C.nucleus;
////            if(C.type==UNLABELED)
////            {  
////               nucleus.draw(unLabelledCell, 255);
////            }
////        }
////        unLabelledCell.show();
//    }
    
    private ArrayList<Cell> getNbNeigbors(ArrayList<Cell> popCells)
    {
        IJ.log("Computing nb neighbors for each cell");
        Iterator<Cell> it = popCells.iterator();
        while(it.hasNext())
        {
            Cell Ci = it.next();
            int cLabelled = 0;
            Ci.degreeLabelled = cLabelled;
        }
        
        
       
        return popCells;
    }
    private int getCellColor(Cell C)
    {
        int val = 0;
        if (C.type != UNLABELED) 
        {
            val = MARKER_VAL;
//            if(C.type==DELTA)
//            {
//                val = 50;
//            }
//            else if(C.type==BETA)
//            {
//                val = 100;
//            }
//            else
//            {
//                val = 150;
//            }
        }    
        return val;
    }
    
    
    // String[] markerObs = {"CD3", "CD4", "CD57","CD8"};  // 2, 3, 4, 5 pixel val
    private void drawGraphJ(ImageHandler img, String image_fn, String save_dir)
    {
        ObjectCreator3D objs;
        objs = new ObjectCreator3D(img.sizeX, img.sizeY, img.sizeZ);
        int rad = 6;
        
        // Draw cells
        for (Cell C : popCells) 
        {
            Vector3D cen = C.nucleus.getCenterAsVector();
            int rad1 = rad + (int)C.degreeLabelled;
            int colorVal = 1;
            if (C.type != UNLABELED) 
            {
                colorVal = MARKER_VAL;
            }    
            objs.createSphere(cen.getX(), cen.getY(), cen.getZ(), rad1, colorVal, false);
        }
        
        // Draw connection between cells
        for (Cell C : popCells) 
        {
            Object3D nuc = C.nucleus;
                Vector3D cen = nuc.getCenterAsVector();
                
                for (Cell C1 : C.nei1) 
                {
                    Object3D nuc1 = C1.nucleus;
                    Vector3D cen1 = nuc1.getCenterAsVector();
                    int distance = 2;
                    int val = getContactColor(C.type, C1.type);
                    objs.createLine(cen, cen1, val, distance);

                }
        }
        
        ImagePlus graph = new ImagePlus("CellsNetwork", objs.getStack());
        graph.show();
        IJ.selectWindow(graph.getTitle());
        IJ.saveAs(graph, "Tiff", save_dir+image_fn+"_CellsNetwork.tif");
//        IJ.log("Saving the output image CellsNetwork.tif into the folder: "+save_dir);
        
    }
    
    private int getContactColor(byte cType1, byte cType2)
    {
        int val = 0; 
        if((cType1==UNLABELED && cType2==MARKER_VAL) || (cType1==MARKER_VAL && cType2==UNLABELED))
        {
            val = 7;  //heterogeneous contact
        }
        else if(cType1==MARKER_VAL && cType2==MARKER_VAL)
        {
            val = MARKER_VAL; //homogeneous contact
        }
        else
        {
            val = 1; // 2 UNLABELED cells, homogeneous contact
        }
        return val;
    }
}
