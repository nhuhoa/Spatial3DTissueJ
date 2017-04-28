/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */
public class TypeSLIC_ implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;

    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;
    int thresh = 0;
    int thres1 = 5;
    double errRate = 0.1;
    //ImageHandler[] signals;
    
    ImageHandler labelDAPI, labelAlpha, labelBeta, labelDelta;
    ImageInt imgSeg;
    ArrayList<Cell> popCells = null;
    
    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    ImageInt imgWat = null;
    @Override
    public void run(String arg) {
//        labelDAPI = ImageHandler.wrap(WindowManager.getImage("label-re-C4-dapi.tif"));
//        labelAlpha = ImageHandler.wrap(WindowManager.getImage("label-re-C1-alpha.tif"));
//        labelBeta = ImageHandler.wrap(WindowManager.getImage("label-re-C2-beta.tif"));
//        labelDelta = ImageHandler.wrap(WindowManager.getImage("label-re-C1-delta.tif"));
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        
        WindowManager.setTempCurrentImage(plus);
        String dir = IJ.getDirectory("image");
        imgWat = ImageInt.wrap(plus);
        imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
        
      

        IJ.log("Initialization Tissue Analysis ...");
        initCells(imgSeg, imgWat);
        IJ.log("test distance map ...");
        
        testDomain();
        //this.computeTypeCellsNucleusSLIC();
                
        //computeTypeCellsLayer(10);
        /*popNuclei.saveObjects(dir + "Nuclei.zip");
        popRegions.saveObjects(dir + "Regions.zip");
        IJ.log("save the result");
        Manager3DType_ manager = new Manager3DType_();
        manager.directory=dir;
        //ArrayList<Object3D> sortNuc = popNuclei.getObjectsList();
        //Collections.sort(sortNuc, new CompareNucleus());
        //Objects3DPopulation popNucleiSort=new Objects3DPopulation();
        //popNucleiSort.addObjects(sortNuc);
        manager.addObjects3DPopulation(popNuclei);
        manager.save();
        */
        //drawCellTypes(false).show("TYPE");
        
        IJ.log("Finished");
    }  
    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE", labelDAPI.sizeX, labelDAPI.sizeY, labelDAPI.sizeZ);

        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.type);
            } else {
                C.region.draw(draw, C.type);
            }
        }

        return draw;
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
        //IJ.log("number of error: " + count);
    }

    public ImageFloat getDomain(ImageHandler img)
    {
        boolean inverse = true;
        int threshold = 1;
        ImageFloat r = EDT.run(img, threshold, inverse, Runtime.getRuntime().availableProcessors());
        if (r != null) 
        {   
            ImageFloat r2 = r.duplicate(); 
            normalizeDistanceMap(r2, img.threshold(threshold, inverse, true), (float) 0.08);
            return r2;
        }
        return null;

    }  
    public void testDomain()
    {
        ImageFloat imgBorder = getDomain(imgSeg);
        for (Cell C : popCells)
        {
            Object3D region = C.region;
            region.setLabelImage(imgWat);
            ArrayList<Voxel3D> cont = region.getContours();
            Object3DVoxels vos = new Object3DVoxels();
            vos.addVoxels(cont);
            vos.draw(imgBorder, 255);
            region.setLabelImage(null);
        } 
        for(Cell C : popCells)
        {
            Object3D nuc = C.nucleus;
            
            nuc.draw(imgBorder, 200);
        }
//        for(Cell C : popCells)
//        {
//            Object3D region = C.region;
//            Object3DVoxels intersectR = new Object3DVoxels();
//            intersectR.addVoxelsIntersection(region, region);
//        }
        ImageFloat imgBorder1 = imgBorder.duplicate();
        imgBorder1.show("result contour");
        imgBorder.intersectMask(imgWat);
        imgBorder.show("result final");
                
//        ImageHandler imgReg = new ImageShort("region", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
//        Cell C = popCells.get(5);
//        Object3D region = C.region;
//        region.draw(imgReg, 1);
//        imgBorder.intersectMask((ImageInt)imgReg);
//        imgBorder.show("result final");
    }
    private void computeTypeCellsNucleusSLIC() 
    {
        ImageFloat imgBorder = getDomain(imgSeg);
        double[] histype = {0, 0, 0, 0};
        for (Cell C : popCells) 
        {
            Object3D nucleus = C.nucleus;
            Object3D region = C.region;
            
            if(nucleus==null) continue;
            
            ArrayUtil ac1 = nucleus.listValues(labelDelta); 
            int nbTotalC1 = ac1.countValueAbove(thres1);
            ArrayUtil ac2 = nucleus.listValues(labelBeta);            
            int nbTotalC2 = ac2.countValueAbove(thres1);
            ArrayUtil ac3 = nucleus.listValues(labelAlpha);            
            int nbTotalC3 = ac3.countValueAbove(thres1);
            //ArrayUtil ac4 = nucleus.listValues(labelDAPI);            
            //int nbTotalC4 = ac4.countValueAbove(thres1);
            //int nbTotal = nbTotalC1 + nbTotalC2 + nbTotalC3 + nbTotalC4;
            
            int val1 = (int)nucleus.getPixMeanValue(labelDelta);
            int val2 = (int)nucleus.getPixMeanValue(labelBeta);
            int val3 = (int)nucleus.getPixMeanValue(labelAlpha);
            
            
            /*ArrayList<regionCell> arrList = new ArrayList<regionCell>();
            regionCell r1 = new regionCell(val1, nbTotalC1, 1);
            regionCell r2 = new regionCell(val2, nbTotalC2, 2);
            regionCell r3 = new regionCell(val3, nbTotalC3, 3);
            arrList.add(r1);
            arrList.add(r2);
            arrList.add(r3);
            double maxValMedian = 0, maxVolume = 0; 
            int indexM = 0;
            int indexV = 0;
            for(int k=0; k<arrList.size(); k++)
            {
                if(maxVolume < arrList.get(k).getVolume())
                {
                    maxVolume = arrList.get(k).getVolume();
                    indexV = arrList.get(k).getIndex1();
                } 
                if(maxValMedian < arrList.get(k).getMValue())
                {
                    maxValMedian = arrList.get(k).getMValue();
                    indexM = arrList.get(k).getIndex1();
                }
            }    
            if(indexV==indexM && indexV >0)
            */
            
            if(nbTotalC1 > nbTotalC2 && nbTotalC1 > nbTotalC3)
            {
                C.type = DELTA;
            }
            else if(nbTotalC2 > nbTotalC1 && nbTotalC2 > nbTotalC3)
            {
                C.type = BETA;
            }
            else if(nbTotalC3 > nbTotalC1 && nbTotalC3 > nbTotalC2)
            {
                C.type = ALPHA;
            }
            else if(nbTotalC1 == nbTotalC2 && nbTotalC1 > nbTotalC3)
            {
                if(val1 > val2 && val1 > val3){C.type = DELTA;}
                if(val2 > val1 && val2 > val3){C.type = BETA;}
            }
            else if(nbTotalC1 == nbTotalC3 && nbTotalC1 > nbTotalC2)
            {
                if(val1 > val2 && val1 > val3){C.type = DELTA;}
                if(val3 > val1 && val3 > val2){C.type = ALPHA;}
            }
            else if(nbTotalC2 == nbTotalC3 && nbTotalC2 > nbTotalC1)
            {
                if(val2 > val1 && val2 > val3){C.type = BETA;}
                if(val3 > val1 && val3 > val2){C.type = ALPHA;}
            }
            else{
                IJ.log("-----------------------------------------------------------");
                ImageHandler imgReg = new ImageShort("region", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
                region.draw(imgReg, 1);
                imgBorder.intersectMask((ImageInt)imgReg);
                int nbA=0, nbB=0, nbD=0;
                for (int z = 0; z < imgReg.sizeZ; z++) {
                    for (int y = 0; y < imgReg.sizeY; y++) {
                        for (int x = 0; x < imgReg.sizeX; x++) {
                            if(imgReg.getPixel(x, y, z)>0) 
                            {
                                if(labelDelta.getPixel(x, y, z)>thres1) 
                                {
                                    nbD++;
                                }
                                if(labelBeta.getPixel(x, y, z)>thres1) 
                                {
                                    nbB++;
                                }
                                if(labelAlpha.getPixel(x, y, z)>thres1) 
                                {
                                    nbA++;
                                }

                            }

                        }    
                    }
                }
                IJ.log("nb Domain 8%: A: " + nbA + " B: " + nbB + " D: " +nbD);
                
                IJ.log("Volume nuc: del: " + nbTotalC1 + " beta: " + nbTotalC2 + " al: " + nbTotalC3);
                IJ.log("Pix mean nuc is: del: " + val1 + " beta:" + val2 + " al: " + val3);
                ArrayUtil rac1 = region.listValues(labelDelta);            
                int rnbTotalC1 = rac1.countValueAbove(thres1);
                ArrayUtil rac2 = region.listValues(labelBeta);            
                int rnbTotalC2 = rac2.countValueAbove(thres1);
                ArrayUtil rac3 = region.listValues(labelAlpha);            
                int rnbTotalC3 = rac3.countValueAbove(thres1);
                int rval1 = (int)region.getPixMeanValue(labelDelta);
                int rval2 = (int)region.getPixMeanValue(labelBeta);
                int rval3 = (int)region.getPixMeanValue(labelAlpha);
                IJ.log("Volume reg: del: " + rnbTotalC1 + " beta: " + rnbTotalC2 + " al: " + rnbTotalC3);
                IJ.log("Pix mean reg is: del: " + rval1 + " beta:" + rval2 + " al: " + rval3);
                IJ.log("-----------------------------------------------------------");
            }
            
            histype[C.type]++;
            nucleus.setType(C.type);
            nucleus.setName("Nuc" + C.id);
            nucleus.setValue(C.id);
            region.setType(C.type);
            region.setName("Reg" + C.id);
            region.setValue(C.id);
        }
        IJ.log("Number of cell type is : " + histype[UNLABELED] + " " + histype[ALPHA] + " " + histype[BETA] + " " + histype[DELTA]);
        
        
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", labelDAPI.sizeX, labelDAPI.sizeY, labelDAPI.sizeZ);

        for (Cell C : popCells) {
            if(C.type==UNLABELED){
               Object3D nucleus = C.nucleus;
               nucleus.draw(unLabelledCell, 255);
            }     
        }
        unLabelledCell.show("Unlabelled");
        unLabelledCell.setTitle("unlabelled");
        //unLabelledCell.save(dir);


       /*double[] y = {0, 1, 2, 3};
        PlotWindow.noGridLines = false; // draw grid lines
        Plot plot = new Plot("Plot","X Axis","Y Axis", y, histype);
        plot.setLimits(0, 5, 5, 900);
        plot.setLineWidth(1);
        plot.changeFont(new Font("Helvetica", Font.PLAIN, 16));
        plot.setColor(Color.blue);
        plot.show();
       */        
        IJ.log("NB : unlabelled : " + histype[UNLABELED] + " alpha: " + histype[ALPHA] + " beta: " + histype[BETA] 
                + " delta: " + histype[DELTA]);
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
    private void computeTypeCellsNucleusSLICV2() 
    {
        double[] histype = {0, 0, 0, 0};
        for (Cell C : popCells) 
        {
            Object3D nucleus = C.nucleus;
            Object3D region = C.region;
            // get list of regions inside nucleus
            //IJ.log("Nuc="+nucleus+" label="+label);
            if(region==null) continue;
            
            ArrayUtil ac1 = region.listValues(labelDelta);            
            int nbTotalC1 = ac1.countValueAbove(thres1);
            ArrayUtil ac2 = region.listValues(labelBeta);            
            int nbTotalC2 = ac2.countValueAbove(thres1);
            ArrayUtil ac3 = region.listValues(labelAlpha);            
            int nbTotalC3 = ac3.countValueAbove(thres1);
            //ArrayUtil ac4 = region.listValues(labelDAPI);            
            //int nbTotalC4 = ac4.countValueAbove(thres1);
            //int nbTotal = nbTotalC1 + nbTotalC2 + nbTotalC3 + nbTotalC4;
            
            int nbA=0; int nbB=0; int nbD=0;
            if(nbTotalC1 > nbTotalC2 && nbTotalC1 > nbTotalC3)
            {
                C.type = DELTA;
            }
            if(nbTotalC2 > nbTotalC1 && nbTotalC2 > nbTotalC3)
            {
                C.type = BETA;
            }
            if(nbTotalC3 > nbTotalC1 && nbTotalC3 > nbTotalC2)
            {
                C.type = ALPHA;
            }
            histype[C.type]++;
            nucleus.setType(C.type);
            nucleus.setName("Nuc" + C.id);
            nucleus.setValue(C.id);
            region.setType(C.type);
            region.setName("Reg" + C.id);
            region.setValue(C.id);
        }
        IJ.log("Number of cell type is : " + histype[UNLABELED] + " " + histype[ALPHA] + " " + histype[BETA] + " " + histype[DELTA]);
        
        
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", labelDAPI.sizeX, labelDAPI.sizeY, labelDAPI.sizeZ);

        for (Cell C : popCells) {
            if(C.type==UNLABELED){
               Object3D nucleus = C.nucleus;
               nucleus.draw(unLabelledCell, 255);
            }     
        }
        unLabelledCell.show("Unlabelled");
        unLabelledCell.setTitle("unlabelled");
        //unLabelledCell.save(dir);


       /*double[] y = {0, 1, 2, 3};
        PlotWindow.noGridLines = false; // draw grid lines
        Plot plot = new Plot("Plot","X Axis","Y Axis", y, histype);
        plot.setLimits(0, 5, 5, 900);
        plot.setLineWidth(1);
        plot.changeFont(new Font("Helvetica", Font.PLAIN, 16));
        plot.setColor(Color.blue);
        plot.show();
       */        
        IJ.log("NB : unlabelled : " + histype[UNLABELED] + " alpha: " + histype[ALPHA] + " beta: " + histype[BETA] 
                + " delta: " + histype[DELTA]);
    }
    
    private class Cell {

        int id;
        Object3D nucleus;
        Object3D region;
        byte type;
        int layer;

        ArrayList<Cell> nei1 = null;
        ArrayList<Cell> nei2 = null;

        @Override
        public String toString() {
            return "(" + nucleus.getValue() + ", " + region.getValue() + ", " + type + ")";
        }

        private int[] computeNeiType(ArrayList<Cell> nei) {
            int[] res = new int[3];
            for (Cell C : nei) {
                if (C.type > 0) {
                    res[C.type - 1]++;
                }
            }
            return res;
        }

        public int[] computeNei1Type() {
            return computeNeiType(nei1);
        }

        public int[] computeNei2Type() {
            return computeNeiType(nei2);
        }

        public int[] computeNeiRangeType(double dist) {
            return computeNeiType(getCellsRange(dist));
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
    }

}
