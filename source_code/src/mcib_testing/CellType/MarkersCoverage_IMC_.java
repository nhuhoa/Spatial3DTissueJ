/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
//import ij.io.Opener;
import ij.measure.ResultsTable;
import java.io.File;
import java.io.IOException;
//import java.awt.event.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
//import ij.plugin.FolderOpener;
//import ij.process.*;
//import ij.measure.Calibration;
//import ij.plugin.FolderOpener;
import ij.util.*;
//import ij.plugin.frame.Recorder;
import java.text.DecimalFormat;
import java.util.ArrayList;
//import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
//import java.util.List;
import mcib3d.geom.Object3D;
//import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
//import mcib3d.geom.Voxel3D;
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

public class MarkersCoverage_IMC_ implements ij.plugin.PlugIn {
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
    public String curr_dir = " ";
    private double ratioMarker = 0.2;
//    private float min = 0, maxInside = 4, maxOutside = 2;
//    private final int IN_OUTSIDE = 0, OUTSIDE_NUCLEUS = 1, INSIDE_NUCLEUS = 2, WATERSHED_REGION = 3;
    public boolean showUnlabelled = true;
//    private static String none = "*None*";
    
    public String save_dir = " ";
    
    public void run(String arg) 
    {
        
//        boolean save = true; 
//        curr_dir = "/Users/miu/Documents/workspace/projects_BCCRC/IMC/XP1487_IMCpipeline/results";
////        String suffixe = ".tif";
//        if (!curr_dir.endsWith("/"))
//            curr_dir = curr_dir + "/";
//        String binary_marker_dir = curr_dir+"BIN/";
//        String raw_marker_dir = curr_dir+"FILTERED/";
//        IJ.log("Current working directory: " + curr_dir);
        
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return;
        }
        if(wList==null || wList.length<3)
        {
            IJ.error("At least 3 imgs : segmented nuclei, a stack of filtered images, and a stack of binary markers images should be open.");
            return;
        }  

        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        boolean save = true; 
        GenericDialog gd = new GenericDialog("Cell Profiles Estimation");
        gd.addMessage(" ");
        gd.addMessage("    Spatial3DTissueJ  ");
        gd.addMessage("Input:  nuclei segmented");
        gd.addMessage("        stack of filtered images");
        gd.addMessage("        stack of binary markers images");
        gd.addMessage("Output: cell profiles");
        gd.addMessage("Output: cell network");
        gd.addMessage(" ");
        gd.addStringField("Save Dir: ", save_dir);
        gd.addChoice("NUC SEGMENTED IMAGE : ", titles, titles[0]);
        gd.addChoice("FILTERED MARKERS STACK : ", titles, titles[0]);
        gd.addChoice("BINARY MARKERS STACK :  ", titles, titles[0]);
        gd.addNumericField("Percent_Marker : ", ratioMarker, 1);
        gd.addCheckbox("Save the output into folder", true);
        
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        
        save_dir = gd.getNextString();
        if (" ".equals(save_dir)) 
        {
                return;
        }
        int[] index = new int[3];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        ratioMarker = (double) gd.getNextNumber();
        save  = gd.getNextBoolean();
        


//        ratioMarker = (double) gd.getNextNumber();
//        wat_radius = (double) gd.getNextNumber();
//        IJ.log("DEBUG is: " + binary_marker_dir);
//        IJ.log("DEBUG is: " + raw_marker_dir);
        IJ.log("Watershed radius is: " + wat_radius);
//        min = (float) gd.getNextNumber();
//        maxInside = (float) gd.getNextNumber();
//        maxOutside = (float) gd.getNextNumber();
//        save  = gd.getNextBoolean();
//        showUnlabelled  = gd.getNextBoolean();
//        IJ.log("Ratio marker to assign a cell type: "+ratioMarker);
//        IJ.log("Detecting cell type: "+markerObs[marker_type]);
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage(wList[index[0]]));
        IJ.log("Seg: " + imgSeg.getTitle());
        
        ImagePlus impFiltered = WindowManager.getImage(wList[index[1]]);
        IJ.log("Filtered markers images: " + impFiltered.getTitle());
        ImagePlus impBin = WindowManager.getImage(wList[index[2]]);
        IJ.log("Binary images: " + impBin.getTitle());
//        
//        ImageInt imgWat = null;
//        if(none.equals(wList[index[1]])){ //(wList.length==2) || 
//            IJ.log("Dont exist wat img. compute it later");
//        }else{
//            imgWat = ImageInt.wrap(WindowManager.getImage(wList[index[1]]));
//            IJ.log("Wat: " + imgWat.getTitle());
//        }
//        save_dir = curr_dir+"CELL_TYPE/";
        File fl = new File(save_dir);
        if (!fl.exists()) 
        {
            fl.mkdir();
            IJ.log("Creating save dir folder to store results");
            IJ.log(save_dir);
        }
        
//        IJ.log("Initialization Tissue Analysis ...");
        ImageInt imgWat = null;
        analysis_v2(imgSeg, imgWat, impFiltered, impBin, save_dir);
        
        
    }
    


    
    public void analysis_v2(ImageInt imgSeg, ImageInt imgWat,
                         ImagePlus impFiltered, ImagePlus impBin, 
                         String save_dir){
        ArrayList<ImageInt> lsbin = new ArrayList<ImageInt>();
        ArrayList<ImageInt> lsraw = new ArrayList<ImageInt>();
        double background_threshold = 0;
        ImageStack s = impFiltered.getStack();
        ArrayList<String> marker_titles = new ArrayList<String>();
        for (int i=1; i<s.getSize(); i++) {
            ImagePlus sp = new ImagePlus(s.getSliceLabel(i),s.getProcessor(i));
            String img_title = s.getSliceLabel(i);
            IJ.log("Input filtered images: "+img_title);
            if (img_title.endsWith(".tif")){
                img_title = img_title.substring(0, img_title.length()-4);
            }
//            directory += "/";
            marker_titles.add(img_title);
            ImageInt imgRaw = ImageInt.wrap(sp);
            lsraw.add(imgRaw);
        }
//        ArrayList<String> bin_marker_titles = new ArrayList<String>();
        ImageStack sb = impBin.getStack();
        for (int i=1; i<sb.getSize(); i++) {
            ImagePlus spb = new ImagePlus(""+sb.getSliceLabel(i),sb.getProcessor(i));
            ImageInt imgBIN = ImageInt.wrap(spb);
            lsbin.add(imgBIN);
        }
        IJ.log("DEBUG FILTER SIZE" + lsraw.size());
        IJ.log("DEBUG BIN SIZE" + lsbin.size());
//        ImageInt imgLabel;
//        HashMap<Integer, Cell> region2Cell=null;
        
        String image_fn = imgSeg.getTitle();
        Pattern p = Pattern.compile("_SEG");
        Matcher m = p.matcher(image_fn); 
        image_fn = m.replaceAll("");
        IJ.log(image_fn);
        IJ.log("Generate watershed cell segmentation image using radius from nucleus: " + wat_radius);

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
        IJ.log("Eastimating cell profiles ...");
        IJ.log(image_fn);
        IJ.log("Nb images bin is: " + lsbin.size());
        IJ.log("Nb images raw is: " + lsraw.size());
        computeCellTypeV2(lsbin, lsraw, marker_titles, save_dir, image_fn); 
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
    private void computeCellTypeV2(ArrayList<ImageInt> lsbin, ArrayList<ImageInt> lsraw, 
                                   ArrayList<String> marker_titles, 
                                   String save_dir, String image_fn) 
    {
        
        IJ.log("Computing cell type...");
        ResultsTable nodes = new ResultsTable();
        ResultsTable edges = new ResultsTable();
        int c = 0;
        int k = 0;
        IJ.log("Nb cells is: "+ popCells.size());
        
//        DecimalFormat df = new DecimalFormat("#.###");
        for (Cell C : popCells) {
            if(C.region==null) continue;
            nodes.incrementCounter();
//            IJ.log("NUC: " + C.nucleus.getValue());
            nodes.setValue("cell_id", c, C.nucleus.getValue()); //index
            nodes.setValue("x", c, Math.round(C.nucleus.getCenterX()));
            nodes.setValue("y", c, Math.round(C.nucleus.getCenterY()));
            nodes.setValue("z", c, Math.round(C.nucleus.getCenterZ()));
            nodes.setValue("nuc_area_pixel", c, Math.round(C.nucleus.getVolumePixels()));
//            nodes.setValue("nuc_area_unit", c, Math.round(C.nucleus.getVolumeUnit()));
            nodes.setValue("cell_area_pixel", c, Math.round(C.region.getVolumePixels()));
//            nodes.setValue("cell_area_unit", c, Math.round(C.region.getVolumeUnit()));
 
            for(int m=0; m<lsbin.size(); m++){
                ArrayUtil list = C.region.listValues(lsbin.get(m));

                int nbM = list.countValueAbove(0);
//                double vol = C.nucleus.getVolumePixels();
                nodes.setValue(marker_titles.get(m)+"_mean_intensity_nuc", c, Math.round(C.nucleus.getPixMeanValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_median_intensity_nuc", c, Math.round(C.nucleus.getPixMedianValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_mean_intensity_cellzone", c, Math.round(C.region.getPixMeanValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_median_intensity_cellzone", c, Math.round(C.region.getPixMedianValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_pct_coverage_nuc", c, Math.round(100*nbM /C.nucleus.getVolumePixels()));
                nodes.setValue(marker_titles.get(m)+"_pct_coverage_cellzone", c, Math.round(100*nbM /C.region.getVolumePixels()));
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
            String node_fn = save_dir+image_fn+"_cell_profiles.csv";
            String edge_fn = save_dir+image_fn+"_cells_interaction_network.csv";
            IJ.log(node_fn);
            IJ.log(edge_fn);
            nodes.saveAs(node_fn);
            edges.saveAs(edge_fn);
        } catch (IOException ex) {
            IJ.log("Have the problem of save cell profiles data");
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
    public String[] trimFileList(String[] rawlist) {
        String filter=".tif";
        if (rawlist==null)
                return null;
        int count = 0;
        for (int i=0; i< rawlist.length; i++) {
                String name = rawlist[i];
                if (name.startsWith(".")||name.equals("Thumbs.db"))
                        rawlist[i] = null;
                else
                        count++;
        }
        if (count==0) return null;
        String[] list = rawlist;
        if (count<rawlist.length) {
                list = new String[count];
                int index = 0;
                for (int i=0; i< rawlist.length; i++) {
                        if (rawlist[i]!=null && rawlist[i].contains(filter))
                            IJ.log("File name is: "+rawlist[i]);
                            list[index++] = rawlist[i];
                }
        }
        return list;
    }
}



// list all files from this path  
//    public String[] listFilesFromFolder(String directory) {
//        if (directory==null || directory.length()==0) {
//                IJ.error("No directory specified.     ");
//                return null;
//        }
//        File file = new File(directory);
//        String[] list = file.list();
//        if (list==null) {
//            String parent = file.getParent();
//            if (parent!=null) {
//                    file = new File(parent);
//                    list = file.list();
//                    IJ.log("Number of files: "+list.length);
//            }
//            if (list!=null){
//                    directory = parent;
//            }
//            else {
//                    IJ.log("Directory not found: "+directory);
//                    return null;
//            }
//        }
//        IJ.log("Directory is:  "+directory);
//        if (!(directory.endsWith("/")||directory.endsWith("\\")))
//            directory += "/";
//        //remove subdirectories from list
//        ArrayList fileList = new ArrayList();
//        for (int i=0; i<list.length; i++) {
//            File f = new File(directory+list[i]);
//            if (!f.isDirectory())
//                        fileList.add(list[i]);
//                        IJ.log(": "+list[i]);
//                
//        }
//        if (fileList.size()<list.length)
//            list = (String[])fileList.toArray(new String[fileList.size()]);
//        
//        list = trimFileList(list);
//        if (list==null)
//            return null;
//            
//        return list;
//
//    }
//    public ArrayList<ImageInt> loadImage(String[] list, String directory) {
//        ArrayList<ImageInt> ls_imgs = new ArrayList<ImageInt>();
//        try {
//            for (int i=0; i<list.length; i++) {
//                    Opener opener = new Opener();
//                    opener.setSilentMode(true);
//                    IJ.redirectErrorMessages(true);
////                    ImagePlus imp = opener.openTempImage(directory, list[i]);
//                    ImagePlus imp = opener.openTiff(directory, list[i]);
//                    IJ.redirectErrorMessages(false);
//                    if (imp!=null) {
//                            IJ.log("Image name: "+ list[i]);
//                            IJ.log("Image depth: "+imp.getBitDepth());
//                            ImageInt img = ImageInt.wrap(imp);
//                            ls_imgs.add(img);
////                            break;
//                    }
//            }
//        }  catch(OutOfMemoryError e) {
//                IJ.outOfMemory("FolderOpener");
//                return(null);
//        }  
//        return ls_imgs;
//        
//    }




        
//        ImagePlus imp1 = Opener.openUsingBioFormats(raw_marker_dir+"C01.tif");
//        ImageInt imgLabel = ImageInt.wrap(imp1);
//        ImagePlus imp2 = Opener.openUsingBioFormats(raw_marker_dir+"C02.tif");
//        ImageInt imgLabel2 = ImageInt.wrap(imp2);
//        ArrayList<ImageInt> lsraw = new ArrayList<ImageInt>();
//        lsraw.add(imgLabel);
//        lsraw.add(imgLabel2);
//        IJ.log("Debug list: " + lsraw.size());
//        ImagePlus bin_stack = FolderOpener.open(binary_marker_dir, " file=BIN");
//        IJ.run("Stack to Images", "");
//        
//        IJ.run("Image Sequence...", "");
        
//        Opener opener = new Opener();  
//        String imageFilePath = raw_marker_dir+"C03.tif";
//        ImagePlus imp = opener.openImage(imageFilePath);
//        IJ.log("Bit depth: "+imp.getBitDepth());
//        FolderOpener fo = new FolderOpener();
//        ImagePlus raw_stack = fo.openFolder(raw_marker_dir);
//        ImagePlus raw_stack = FolderOpener.open(raw_marker_dir, " file=0");
//        IJ.log("Debug list: " + raw_stack.getSlice());
//        raw_stack.
//        IJ.run("Stack to Images", "");
//        String prefix = "C0";
//        for (int i=1; i<bin_stack.getSlice(); i++) {
////            String fn=binary_marker_dir+prefix + "_BIN.tif";
////            IJ.log(fn);
////            ImagePlus imp = IJ.openImage(fn);	
//            
////            if (imp!=null) {
////                IJ.log("Loading image: C0"+i+"_BIN.tif");
////                ImageInt imgLabel = ImageInt.wrap(imp);
////                lsbin.add(imgLabel);
////            }
//            ImageInt imgLabel = ImageInt.wrap(imp);
//            lsbin.add(imgLabel);
//        }
//        String dir = "/Users/miu/Documents/workspace/projects_BCCRC/IMC/XP1487_IMCpipeline/results/FILTERED/";
//        String[] list = listFilesFromFolder(dir);
//        IJ.log("Loading raw marker images from directory: " + raw_marker_dir);
////        ArrayList<ImageInt> lsraw = loadImage(list, raw_marker_dir);
//        IJ.log("Loading binary marker images from directory: " + binary_marker_dir);
//        IJ.log("Noted that binary images have same name as raw image but end with _BIN.tif: " + binary_marker_dir);
//        String[] list_filtered_fns = new String[list.length];
//        for(int i=0; i<list.length; i++){
//            list_filtered_fns[i]=list[i].substring(1,list[i].length()-4);
//            IJ.log("Filtered fn: "+list_filtered_fns[i]);
//        }
        
//        ImagePlus imp = IJ.openImage("/Users/miu/Documents/workspace/projects_BCCRC/IMC/XP1487_IMCpipeline/results/FILTERED/C01.tif");
//        IJ.log("testing: blah"+imp.getBitDepth());
//        ArrayList<ImageInt> lsbin = loadImage(list_filtered_fns, binary_marker_dir);
//        
//        IJ.log("Nb images bin is: " + lsbin.size());
//        IJ.log("Nb images raw is: " + lsraw.size());
//        if(lsbin.size()!=lsraw.size()){
//            IJ.error("Double check the list of input markers images in raw folder and binary folder");
//        }
//        IJ.log("Input data: wat: "+imgWat.getTitle()+" seg: "+imgSeg.getTitle()+" label: "+imgLabel.getTitle());
//        IJ.log("Ratio marker: "+ ratioMarker);
//        IJ.log("Region observed: "+regionObs[reg]);
//        MARKER_VAL = marker_type + 2;
//        IJ.log("Observed marker: " + markerObs[marker_type] + "  MARKER_VAL: " + MARKER_VAL);
//        IJ.log("Save dir is: " + save_dir);
