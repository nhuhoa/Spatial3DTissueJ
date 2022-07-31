/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;

/**
 *
 * @author tranhoa
 */
public class DetectCellType_V2 implements ij.plugin.PlugIn 
{
    private final int UNLABELED = 0, ALPHA = 3, BETA = 2, DELTA = 1;
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;
    public String save_dir = "";
    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    private double ratioMarker = 0.1, rDB = 0.5, rAB = 0.9, rDA=0.9;
    private float min = 0, maxInside = 1, maxOutside = 3;
    private final int IN_OUTSIDE = 0, OUTSIDE_NUCLEUS = 1, INSIDE_NUCLEUS = 2, WATERSHED_REGION = 3;
    public boolean showUnlabelled = true;
    @Override
    public void run(String arg) 
    {
        String[] regionObs = {"INSIDE_OUTSIDE NUCLEUS", "OUTSIDE NUCLEUS", "INSIDE NUCLEUS","WATERSHED REGION"};
        int reg = 0;
        
        int[] wList = WindowManager.getIDList();
//        if (wList==null) {
//            IJ.error("No images are open.");
//            return;
//        }
        if(wList==null || wList.length<3)
        {
            IJ.error("At least 3 imgs : watershed, segmented nuclei, and SLIC labelF should be open.");
            return;
        }  

//        String[] titles = new String[wList.length];
//        for (int i=0; i<wList.length; i++) {
//            ImagePlus imp = WindowManager.getImage(wList[i]);
//            titles[i] = imp!=null?imp.getTitle():"";
//        }
        String[] titles = new String[wList.length+1];
        titles[0] = "*None*";
        
        ImagePlus temp = null;
        int count=0;
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i+1] = imp!=null?imp.getTitle():"";
            if (null !=imp){
                count++;
                if(count==1){
                    temp = imp;
                    WindowManager.setTempCurrentImage(temp);
                    save_dir = IJ.getDirectory("image");
                }
                
            }    
        }
        boolean save = true; 
        GenericDialogPlus gd = new GenericDialogPlus("3D Cell Type Detection");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("3D Tissue Spatial Analysis");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        "organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input: nuclei segmented");
        gd.addMessage("Input: watershed image");
        gd.addMessage("Input: SLIC composite label");
        gd.addMessage("Output: cell type image");
        gd.addMessage(" ");
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addChoice("SLIC_composite_img : ", titles, titles[0]);
        gd.addChoice("Watershed_CellZone_Img : ", titles, titles[1]);
        gd.addChoice("Nuclei_Segmented_Img :  ", titles, titles[2]);
//        gd.addChoice("BIN ALPHA: ", titles, titles[0]);
        gd.addChoice("FILTERED ALPHA: ", titles, titles[0]);
        //        gd.addChoice("BIN BETA: ", titles, titles[0]);
        gd.addChoice("FILTERED BETA: ", titles, titles[0]);
        gd.addChoice("FILTERED DELTA: ", titles, titles[0]);
//        gd.addNumericField("Ratio_Delta_Beta : ", rDB, 1);
//        gd.addNumericField("Ratio_Alpha_Beta : ", rAB, 1);
//        gd.addNumericField("Ratio_Delta_Alpha : ", rDA, 1);
        gd.addChoice("Region_Observed : ", regionObs, regionObs[reg]);
        gd.addNumericField("Percent_marker_coverage : ", ratioMarker, 1);
        gd.addNumericField("Min_Distance : ", min, 0);
        gd.addNumericField("Max_Distance_Inside : ", maxInside, 1);
        gd.addNumericField("Max_Distance_Outside : ", maxOutside, 3);
        gd.addCheckbox("Save output into folder", true);
        gd.addCheckbox("Show unlabelled cell type", true);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return;
        }
        
//        int[] index = new int[3];
//        index[0] = gd.getNextChoiceIndex();
//        index[1] = gd.getNextChoiceIndex();
//        index[2] = gd.getNextChoiceIndex();
        int[] index = new int[6];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        index[3] = gd.getNextChoiceIndex();
        index[4] = gd.getNextChoiceIndex();
        index[5] = gd.getNextChoiceIndex();

        reg =  (int)(gd.getNextChoiceIndex());
        ratioMarker = (double) gd.getNextNumber();
//        rDB = (double) gd.getNextNumber();
//        rAB = (double) gd.getNextNumber();
//        rDA = (double) gd.getNextNumber();
        min = (float) gd.getNextNumber();
        maxInside = (float) gd.getNextNumber();
        maxOutside = (float) gd.getNextNumber();
        save  = gd.getNextBoolean();
        showUnlabelled  = gd.getNextBoolean();
        
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        ImageInt imgLabel = ImageInt.wrap(WindowManager.getImage(wList[index[0]-1]));
        ImageInt imgWat = ImageInt.wrap(WindowManager.getImage(wList[index[1]-1]));
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage(wList[index[2]-1]));
        
        ArrayList<ImageInt> lsraw = new ArrayList<ImageInt>();
        String[] markerTypes = {"ALPHA","BETA","DELTA"};
        ArrayList<String> marker_titles = new ArrayList<String>();
        for(int i=3; i<=5; i=i+1){
            if(index[i]!=0)
            {   //marker_titles.add(markerTypes[i-3]);
                IJ.log("marker: "+markerTypes[i-3]);
                IJ.log("idx: "+index[i]+" img: "+wList[index[i]-1]);
    //            IJ.log("idx2: "+index[2]+" img: "+wList[index[2]]);
    //            IJ.log("idx3: "+index[10]+" img: "+wList[index[3]]);
    //            ImageInt imgLabel = ImageInt.wrap(WindowManager.getImage(wList[index[i]]));
                ImageInt imgRawLabel = ImageInt.wrap(WindowManager.getImage(wList[index[i]-1]));
    //            imgLabel.setTitle("BIN_"+markerTypes[count]);
                IJ.log("Input filtered image: " + imgRawLabel.getTitle());
                marker_titles.add(imgRawLabel.getTitle());
    //            lsbin.add(imgLabel);
                lsraw.add(imgRawLabel);
                
            }
            
        }   
        
//        IJ.log("Nb images bin is: " + lsbin.size());
        
        
        IJ.log("====================================================================");
        IJ.log("Input data: wat: "+imgWat.getTitle()+" seg: "+imgSeg.getTitle()+" label: "+imgLabel.getTitle());
        IJ.log("Ratio marker: "+ ratioMarker+" rDB: "+ rDB+" rDA: "+rDA+ "  rAB: "+rAB);
        IJ.log("Region observed: "+regionObs[reg]);
        IJ.log("Observe distance inside of nuc region:   min  :  "+min+"    max: "+maxInside);
        IJ.log("Observe distance outside of nuc region:   min  :  "+min+"    max: "+maxOutside);
//        WindowManager.setTempCurrentImage(imgWat.getImagePlus());
//        dir = IJ.getDirectory("image");
        IJ.log("Initialization Tissue Analysis ...");
        initCells(imgSeg, imgWat);
        
        IJ.log("Detecting Cell Type ...");
        if(reg==IN_OUTSIDE)
        {
            computeCellType_V1(imgSeg, imgWat, imgLabel, false);
        }   
        else if(reg==OUTSIDE_NUCLEUS)
        {
            computeCellType_V2(imgSeg, imgWat, imgLabel);
        }
        else if(reg==INSIDE_NUCLEUS)
        {
            computeCellType_V3(imgSeg, imgWat, imgLabel);
        }    
        else   //WATERSHED_REGION
        {
            computeCellType_V4(imgLabel);
        }  
        
        ImageHandler cellTypeWat = drawCellTypes(false, imgSeg);
        cellTypeWat.show("TYPE_WAT");
        ImageHandler cellTypeNuc = drawCellTypes(true, imgSeg);
        cellTypeNuc.show("TYPE_NUC");
        
        
        if(save)
        {
            IJ.log("Save the result as Nuclei.zip and Regions.zip");
            IJ.log("Save the result into the folder : "+save_dir);
            popNuclei.saveObjects(save_dir + "Nuclei.zip");
            popRegions.saveObjects(save_dir + "Regions.zip");
        }    
        else{
            IJ.log("Not save the result");
        }
        IJ.log("Association1 ...");
        computeAssoCells(imgWat, 0);// 0 for thomas´s cell seg; -1 for farsight
        extract_cell_profiles(lsraw, marker_titles);
        IJ.log("Finished");
        IJ.log("====================================================================");
    }
    
    private ImageHandler computeCenterPartNuc(ImageHandler imgSeg)
    {
        ImageHandler imgSeg1 = new ImageShort("seg_", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        for (Cell C : popCells) {
            Object3D nuc = C.nucleus;
            ArrayList<Voxel3D> objs = nuc.getVoxels();
            ArrayList<Voxel3D> objsTmp = new ArrayList<Voxel3D>();
            int zMax = nuc.getZmax();
            int zMin = nuc.getZmin();
            for(Voxel3D vox : objs)
            {
                if(vox.getZ() < zMax && vox.getZ() > zMin)
                {
                    objsTmp.add(vox);
                }    
            }
            Object3DVoxels newObjs = new Object3DVoxels(objsTmp);
            newObjs.draw(imgSeg1);
//            IJ.log("    Size of objs: "+objs.size()+"   of tmp: "+objsTmp.size());
        }
        return imgSeg1;
    }        
    /**
     * description: compute the staining at 10% outside of nucleus
     * @param seg 
     */
    private void computeCellType_V2(ImageInt seg, ImageInt wat, ImageInt label)
    {
        ImageHandler imgSeg1 = computeCenterPartNuc(seg);
//        imgSeg.show();
        boolean inverse = true;
        int threshold = 0;
        ImageFloat r = EDT.run(imgSeg1, threshold, inverse, Runtime.getRuntime().availableProcessors());
        if (r != null) 
        {   
            ImageFloat r2 = r.duplicate(); 
            ImageByte imgTh = r2.thresholdRangeExclusive(min, maxInside);
//            normalizeDistanceMap(r2, imgSeg.threshold(threshold, inverse, true), (float) 0.1);
            imgTh.intersectMask(wat);
            imgTh.show("distance outside_"+min+"_"+maxInside);
            ResultsTable mapR = new ResultsTable();
            int nb = 0;
            for (Cell C : popCells) 
            {
                Object3D nuc = C.nucleus;
                if(nuc==null) continue;
                double volNuc = nuc.getVolumePixels();
                Object3D region = C.region;
                int reVal = region.getValue();
                mapR.incrementCounter();
                mapR.setValue("cell", C.id, reVal);
                ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
                ArrayList<Voxel3D> objs = region.getVoxels();
                for(Voxel3D vox : objs)
                {
                    if(r2.getPixel((int)vox.getX(), (int)vox.getY(), (int)vox.getZ()) > 0)
                    {
                        arr.add(vox);
                    }    
                }
                Object3DVoxels oo = new Object3DVoxels(arr);
                ArrayUtil ac1 = oo.listValues(label); 
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<ac1.size();i++)
                {
                    if(ac1.getValue(i)==ALPHA)nbA++;
                    if(ac1.getValue(i)==BETA)nbB++;
                    if(ac1.getValue(i)==DELTA)nbD++;                
                }
                mapR.setValue("nbA", C.id, nbA);
                mapR.setValue("nbB", C.id, nbB);
                mapR.setValue("nbD", C.id, nbD);
                mapR.setValue("vOutside", C.id, arr.size());
                mapR.setValue("volReg", C.id, region.getVolumePixels());
                mapR.setValue("volNuc", C.id, volNuc);
                C.type = UNLABELED;
                int typeCell = getCellType(nbA, nbB, nbD, arr.size(), 0.1, 1, 1);
                C.type = (byte)typeCell;
                nuc.setType(C.type);
                nuc.setName("Nuc" + C.id);
                nuc.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
                if(C.type != UNLABELED)
                {
                    nb++;
                }    
                mapR.setValue("method_outside", C.id, typeCell);
            }
            IJ.log("__________Nb cells type: "+nb);
            mapR.show("Mapping Cell Type");
        }
        
        
        
    }        
    private void computeCellType_V3(ImageInt seg, ImageInt wat, ImageInt label)
    {
        ImageHandler imgSeg1 = computeCenterPartNuc(seg);
//        imgSeg.show();
        
        boolean inverse = false;
        int threshold = 0;
        ImageFloat r = EDT.run(imgSeg1, threshold, inverse, Runtime.getRuntime().availableProcessors());
        if (r != null) 
        {   
            ImageFloat r2 = r.duplicate(); 
//            normalizeDistanceMap(r2, imgSeg.threshold(threshold, inverse, true), (float) 0.4);
            ImageByte imgTh = r2.thresholdRangeExclusive(min, maxInside);
            imgTh.intersectMask(wat);
            imgTh.show("distance inside_"+min+"_"+maxInside);
            ResultsTable mapR = new ResultsTable();
            int nb = 0;
            for (Cell C : popCells) 
            {
                Object3D nuc = C.nucleus;
                if(nuc==null) continue;
                double volNuc = nuc.getVolumePixels();
                Object3D region = C.region;
                int reVal = region.getValue();
                mapR.incrementCounter();
                mapR.setValue("cell", C.id, reVal);
                ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
                ArrayList<Voxel3D> objs = nuc.getVoxels();
                for(Voxel3D vox : objs)
                {
                    if(r2.getPixel((int)vox.getX(), (int)vox.getY(), (int)vox.getZ()) > 0)
                    {
                        arr.add(vox);
                    }    
                }
                Object3DVoxels oo = new Object3DVoxels(arr);
                ArrayUtil ac1 = oo.listValues(label); 
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<ac1.size();i++)
                {
                    if(ac1.getValue(i)==ALPHA)nbA++;
                    if(ac1.getValue(i)==BETA)nbB++;
                    if(ac1.getValue(i)==DELTA)nbD++;                
                }
                mapR.setValue("nbA", C.id, nbA);
                mapR.setValue("nbB", C.id, nbB);
                mapR.setValue("nbD", C.id, nbD);
                mapR.setValue("vObserved", C.id, arr.size());
                mapR.setValue("volNuc", C.id, volNuc);
                mapR.setValue("volReg", C.id, region.getVolumePixels());
                C.type = UNLABELED;
                int typeCell = getCellType(nbA, nbB, nbD, arr.size(), 0.2, 1, 1);
                C.type = (byte)typeCell;
                nuc.setType(C.type);
                nuc.setName("Nuc" + C.id);
                nuc.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
                if(C.type != UNLABELED)
                {
                    nb++;
                }  
                mapR.setValue("method_inside", C.id, typeCell);
            }
            IJ.log("__________Nb cells type: "+nb);
            mapR.show("Mapping Cell Type");
        }
        
    }  
    /**
     * des: outside + inside
     * @param seg
     * @param wat
     * @param label
     * @param min
     * @param max 
     */
    private void computeCellType_V1(ImageInt seg, ImageInt wat, ImageInt label, boolean reObs)
    {
//        Calibration cal = seg.getCalibration();
        ImageHandler imgSeg = null;
        if(reObs)
        {
            imgSeg = computeCenterPartNuc(seg);
        }   
        else{
            imgSeg = seg.duplicate();
        }
         
        boolean inverseI = false;  //inside segmented nuclei region
        int threshold = 0;
        ImageFloat rI = EDT.run(imgSeg, threshold, inverseI, Runtime.getRuntime().availableProcessors());
        
        boolean inverseO = true;  //outside segmented nuclei region
        ImageFloat rO = EDT.run(imgSeg, threshold, inverseO, Runtime.getRuntime().availableProcessors());
        if(rI != null && rO != null)
        {
            ImageFloat rI2 = rI.duplicate(); 
            ImageByte imgThI = rI2.thresholdRangeExclusive(min, maxInside);
            imgThI.intersectMask(wat);
//            imgThI.show("distance inside_"+min+"_"+maxInside);
            ImageFloat rO2 = rO.duplicate(); 
            ImageByte imgThO = rO2.thresholdRangeExclusive(min, maxOutside);
//            normalizeDistanceMap(r2, imgSeg.threshold(threshold, inverse, true), (float) 0.1);
            imgThO.intersectMask(wat);
//            imgThO.show("distance outside_"+min+"_"+maxOutside);
            
            ResultsTable mapR = new ResultsTable();
            int nb = 0;
            int nA=0, nB=0, nD = 0;
            for (Cell C : popCells) 
            {
                Object3D nuc = C.nucleus;
                if(nuc==null) continue;
                double volNuc = nuc.getVolumeUnit();
                Object3D region = C.region;
                int reVal = region.getValue();
                mapR.incrementCounter();
                mapR.setValue("cell", C.id, reVal);
                ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
                ArrayList<Voxel3D> objs = region.getVoxels();
                for(Voxel3D vox : objs)
                {
                    if(imgThO.getPixel((int)vox.getX(), (int)vox.getY(), (int)vox.getZ()) > 0
                    || imgThI.getPixel((int)vox.getX(), (int)vox.getY(), (int)vox.getZ()) > 0)
                    {
                        arr.add(vox);
                    }    
                }
                Object3DVoxels oo = new Object3DVoxels(arr);
                ArrayUtil ac1 = oo.listValues(label); 
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<ac1.size();i++)
                {
                    if(ac1.getValue(i)==ALPHA)nbA++;
                    if(ac1.getValue(i)==BETA)nbB++;
                    if(ac1.getValue(i)==DELTA)nbD++;                
                }
                mapR.setValue("nbA", C.id, nbA);
                mapR.setValue("nbB", C.id, nbB);
                mapR.setValue("nbD", C.id, nbD);
                mapR.setValue("vObserved", C.id, arr.size());
                mapR.setValue("volNuc", C.id, volNuc);
                mapR.setValue("volReg", C.id, region.getVolumePixels());
                //C.type = UNLABELED;
                int typeCell = UNLABELED;
//                typeCell = getCellType(nbA, nbB, nbD, arr.size(), ratioMarker, rDB, rAB);
//                typeCell = getCellTypeV2(nbA, nbB, nbD, arr.size(), ratioMarker);
                typeCell = getCellTypeV3(nbA, nbB, nbD, arr.size(), ratioMarker);
                C.type = (byte)typeCell;
                
                 
                nuc.setType(C.type);
                nuc.setName("Nuc" + C.id);
                nuc.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
                if(C.type != UNLABELED)
                {
                    nb++;
                    if(C.type==BETA)
                    {nB++;}    
                    else if(C.type==ALPHA)
                    {nA++;}   
                    else
                    {nD++;}   
                }  
                mapR.setValue("method_in_outside", C.id, typeCell);
            }
            if(showUnlabelled)
            {
                ImageHandler unLabelledCell;
                unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);
                for (Cell C : popCells) {
                    Object3D nucleus = C.nucleus;
                    if(C.type==UNLABELED)
                    {  
                       nucleus.draw(unLabelledCell, 255);
                    }
                }
                unLabelledCell.show();
            }    
            IJ.log("Nb cells type : " + nb);
            IJ.log("Nb alpha : " + nA);
            IJ.log("Nb beta : " + nB);
            IJ.log("Nb delta : " + nD);
            mapR.show("Mapping Cell Type");
        }
        
            
    }
    /**
     * Description: compute the cell type by counting all stainning inside nucleus
     * @param label
     * @param ratio
     * @param rDB
     * @param rAB 
     */
    private void computeCellType_V4(ImageInt label) 
    {
        int cA=0, cB=0, cD=0, cUnlabelled=0, cColocAB=0, cColocAD=0, cColocBD=0;
        ResultsTable cellM = new ResultsTable();

            for (Cell C : popCells) 
            {
                Object3D nucleus = C.nucleus;
                double vol = nucleus.getVolumeUnit();
                Object3D region = C.region;
                int reVal = region.getValue();
                cellM.incrementCounter();
                cellM.setValue("cell", C.id, reVal);
                cellM.setValue("id", C.id, C.id);
                if(nucleus==null) continue;
                ArrayUtil list = nucleus.listValues(label);
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<list.size();i++) //list.getSize()
                {
                    if(list.getValue(i)==ALPHA)nbA++;
                    if(list.getValue(i)==BETA)nbB++;
                    if(list.getValue(i)==DELTA)nbD++;                
                }
                C.type = UNLABELED;
                int typeCell = getCellType(nbA, nbB, nbD, vol, ratioMarker, rDB, rAB);
                C.type = (byte)typeCell;
                cellM.setValue("method1", C.id, typeCell);
                if(typeCell==UNLABELED){cUnlabelled++;}
                if(typeCell==ALPHA){
                    cA++;
                    IJ.log("  ");
//                    IJ.log("al: "+100 * (nbA/vol));
                }
                if(typeCell==BETA){
                    cB++;
                    IJ.log("  ");
//                    IJ.log("be: "+100 * (nbB/vol));
                }
                if(typeCell==DELTA){
                    cD++;
                    IJ.log("  ");
//                    IJ.log("del: "+100 * (nbD/vol));
                }
                nucleus.setType(C.type);
                nucleus.setName("Nuc" + C.id);
                nucleus.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
            }

        IJ.log("Nb of each cells type: A: "+cA+"   B: "+cB+"   D:"+cD+"   Unlabelled: "+cUnlabelled);
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);
        for (Cell C : popCells) {
            Object3D nucleus = C.nucleus;
            if(C.type==UNLABELED)
            {  
               nucleus.draw(unLabelledCell, 255);
            }
        }
        unLabelledCell.show();
    }

    private void extract_cell_profiles(ArrayList<ImageInt> lsraw, 
                                   ArrayList<String> marker_titles){
        IJ.log("Extracting cell profiles...");
        ResultsTable nodes = new ResultsTable();
        ResultsTable edges = new ResultsTable();
        int c = 0;
        int k = 0;
//      int UNLABELED = 0, ALPHA = 3, BETA = 2, DELTA = 1;
        String[] label_desc = {"UNLABELED", "DELTA", "BETA","ALPHA"};
        IJ.log("Noted: predefined cell type value is: UNLABELED = 0, ALPHA = 3, BETA = 2, DELTA = 1");
        IJ.log("Nb cells is: "+ popCells.size());
        IJ.log("Nb images raw is: " + lsraw.size());
        IJ.log("Nb marker is: " + marker_titles.size());
////        DecimalFormat df = new DecimalFormat("#.###");
        for (Cell C : popCells) {
//            IJ.log("DEBUG cell type"+C.type);
//            IJ.log("DEBUG cell type"+C.nucleus.getValue());
//            IJ.log("DEBUG cell type"+C.nucleus.getCenterX());
//            IJ.log("DEBUG cell type"+C.nucleus.getVolumePixels());
//            IJ.log("DEBUG cell type"+C.region.getVolumePixels());
            
            if(C.region==null) continue;
            nodes.incrementCounter();
//            IJ.log("NUC: " + C.nucleus.getValue());
            nodes.setValue("cell_id", c, C.nucleus.getValue()); //index
            nodes.setValue("cell_id1", c, C.id); //index
            nodes.setValue("x", c, Math.round(C.nucleus.getCenterX()));
            nodes.setValue("y", c, Math.round(C.nucleus.getCenterY()));
            nodes.setValue("z", c, Math.round(C.nucleus.getCenterZ()));
//            nodes.setValue("cell_type", c, label_desc[C.type]);
            nodes.setValue("cell_type", c, C.type);
            nodes.setValue("nuc_area_pixel", c, Math.round(C.nucleus.getVolumePixels()));
//            nodes.setValue("nuc_area_unit", c, Math.round(C.nucleus.getVolumeUnit()));
            nodes.setValue("cellzone_area_pixel", c, Math.round(C.region.getVolumePixels()));
//            nodes.setValue("cell_area_unit", c, Math.round(C.region.getVolumeUnit()));
 
            for(int m=0; m<lsraw.size(); m++){
                
//                double vol = C.nucleus.getVolumePixels();
                nodes.setValue(marker_titles.get(m)+"_mean_intensity_nuc", c, Math.round(C.nucleus.getPixMeanValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_median_intensity_nuc", c, Math.round(C.nucleus.getPixMedianValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_mean_intensity_cellzone", c, Math.round(C.region.getPixMeanValue(lsraw.get(m))));
                nodes.setValue(marker_titles.get(m)+"_median_intensity_cellzone", c, Math.round(C.region.getPixMedianValue(lsraw.get(m))));
                
            }    
            
            c = c + 1;
            
            // Second extract all edges of cells interaction to a csv file
            for (Cell C1 : C.nei1) 
            {
                edges.incrementCounter();
                edges.setValue("from_cell_id", k, C.nucleus.getValue());
                edges.setValue("to_cell_id", k, C1.nucleus.getValue()); 
                k = k + 1;
            }
        }
//        IJ.selectWindow(image_fn);
//        IJ.saveAs("Tiff", save_dir + image_fn + ".tif");
        
        try {
            String node_fn = save_dir+"cell_profiles.csv";
            String edge_fn = save_dir+"cell2cell_interaction_network.csv";
            IJ.log(node_fn);
            IJ.log(edge_fn);
            nodes.saveAs(node_fn);
            edges.saveAs(edge_fn);
        } catch (IOException ex) {
            IJ.log("Have the problem of save cell profiles data");
        }
        
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

    
    private void computeAssoCells(ImageHandler img, int borderValue) 
    {
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
    private ArrayList<Integer>[] computeAsso(ImageHandler img, int BorderValue) 
    {
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
    private void initCells(ImageInt nucLabel, ImageInt regionLabel) 
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
        IJ.log("Number of cells is: "+popCells.size());
    }
    public void normalizeDistanceMap(ImageFloat distanceMap, ImageInt mask, float ratioUse) 
    {
        int count = 0;
        VoxS[] idx = new VoxS[mask.countMaskVolume()];
        double volume = idx.length;
        for (int z = 0; z < distanceMap.sizeZ; z++) {
            for (int y = 0; y < distanceMap.sizeY; y++) {
                for (int x = 0; x < distanceMap.sizeX; x++) {
                    if (mask.getPixelInt(x, y, z) != 0) {
                        idx[count] = new VoxS(distanceMap.getPixel(x, y, z), x, y, z);
                        count++;
                    }
                }    
            }
        }
        Arrays.sort(idx);

        for (int i = 0; i < idx.length - 1; i++) {
            // gestion des repetitions
            if (idx[i + 1].distance == idx[i].distance) {
                int j = i + 1;
                while (j < (idx.length - 1) && idx[i].distance == idx[j].distance) {
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
        for (VoxS idx1 : idx) 
        {
            if((float) (idx1.index / volume) <= ratioUse)
            {
                distanceMap.setPixel(idx1.x, idx1.y, idx1.z, (float) (idx1.index / volume) * (float)1000.0);
            }
            else{
                distanceMap.setPixel(idx1.x, idx1.y, idx1.z, 0);
            }
            
        }
    }
    
    //for human or in case alpha, delta is very low density
    public int getCellTypeV2(int vA, int vB, int vD, double vol, double rMarker)
    {
        if(vD/vol >= 0.1)
        {
            return DELTA;
        }
        else if(vA/vol >= rMarker)
        {
            return ALPHA;
        } 
        else if(vB > vA && vB > vD && vB/vol >= rMarker)
        {
            return BETA;
        }
        else
        {
            return UNLABELED;
        }  
        
    }
    //for monkey
    public int getCellTypeV3(int vA, int vB, int vD, double vol, double rMarker)
    {
//        double rDB = 0.6, rDA = 0.9;
        if(vD >= vA * rDA && vD >= vB*rDB && vD/vol >= rMarker)
        {
            return DELTA;
        }
//        if(vD/vol >= rMarker)
//        {
//            return DELTA;
//        }
        else if(vA >= vB * rAB && vA > vD && vA/vol >= rMarker)
        {
            return ALPHA;
        } 
        else if(vB > vA && vB > vD && vB/vol >= rMarker)
        {
            return BETA;
        }
        else
        {
            return UNLABELED;
        }  
        
    }
    public int getCellType(int vA, int vB, int vD, double vol, double rMarker, double rDB, double rAB)
    {
        if(vA >= vB && vA > vD && vA/vol >= rMarker)
        {
            return ALPHA;
        } 
        else if(vD >= vA && vD >= vB && vD/vol >= rMarker)
        {
            return DELTA;
        }
        else if(vB > vA && vB > vD && vB/vol >= rMarker)
        {
            if(vD >= vB * rDB && vD>vA && vD/vol >= rMarker)  //0.5
            {
                return DELTA;
            }
            if(vA>=vB * rAB && vA>vD && vA/vol >= rMarker)
            { 
                return ALPHA;
            }
            else{
                return BETA;
            }
        }
        else
        {
            return UNLABELED;
        }  
        
    }
}
