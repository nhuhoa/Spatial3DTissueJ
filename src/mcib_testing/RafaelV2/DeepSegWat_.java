/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import java.util.TreeMap;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.utils.ThreadUtil;

/**
 *
 * @author thomasb
 */
public class DeepSegWat_ implements PlugIn {

    double mean = 0;
    int nbcpu = ThreadUtil.getNbCpus();
    @Override
    public void run(String string) 
    {
        ImagePlus plus = WindowManager.getImage("dapi-seg-wat.tif");
        ImageInt imgWat = ImageInt.wrap(plus);
        Objects3DPopulation popRegions = new Objects3DPopulation(imgWat, 1);
        ImageHandler w = imgWat.createSameDimensions();
        popRegions.draw(w);
        ImageInt wat = (ImageInt)w;
        TreeMap<Integer, int[]> bounds = wat.getBounds(false);
        ImageInt[] wats = wat.crop3D(bounds);
        
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
        TreeMap<Integer, int[]> boundsSeg = imgSeg.getBounds(false);
//        ImageInt[] segs = imgSeg.crop3D(boundsSeg);
        ImageInt[] segs1 = imgSeg.crop3D(bounds);
        
        ImagePlus plus1 = WindowManager.getImage("C4-dapi.tif");
        ImageInt raw = ImageInt.wrap(plus1);
        ImageInt[] raws = raw.crop3D(bounds);
//        ImageByte[] thEx = new ImageByte[raws.length]; 
//        ImageByte[] thIn = new ImageByte[raws.length]; 
//        AutoThresholder thr = new AutoThresholder();
        mean = raw.getMean();
        
//        ArrayList<ImageByte> thExErrorS = new ArrayList<ImageByte>();
//        ArrayList<ImageByte> thExErrorB = new ArrayList<ImageByte>();
//        int countErrS = 0, countErrB = 0;
        IJ.log("size raw: "+raws.length + "size wat: "+wats.length + "seg size: "+segs1.length);
        
        ImageByte[] th = new ImageByte[raws.length];
        Integer[] labels = new Integer[raws.length];
        labels = (Integer[]) (bounds.keySet().toArray(labels));
        
        for (int i = 0; i < raws.length; i++) {
            int lab = labels[i];
            if (lab > 0) {
                th[i] = processRegion(lab, wats[i], raws[i], segs1[i], true);
            } else {
                th[i] = new ImageByte("th", wats[i].sizeX, wats[i].sizeY, wats[i].sizeZ);
            }
        }
        
        ImageShort res = ImageShort.merge3DBinary(th, wat.sizeX, wat.sizeY, wat.sizeZ);
        res.show("Deep_Seg");
        
//        ImageByte mask0 = wats[0].thresholdRangeInclusive((1), (1));
//        //mask0.erode(1, 2);
//        int[] hist0=raws[0].getHistogram(null);            
//        double thre0=thr.getThreshold(AutoThresholder.Method.Otsu, hist0);
//        //IJ.log("i "+thre);
//        thEx[0] = raws[0].thresholdAboveExclusive((float) thre0);
//        thIn[0] = raws[0].thresholdAboveInclusive((float) thre0);
//        thEx[0].intersectMask(mask0);
//        thIn[0].intersectMask(mask0);
//        
//        
//        
//        
//        for (int i = 1; i < raws.length-1; i++) 
//        {            
//            if(i==5)
//            {
//                wats[i].show("mask");
//                raws[i].show("raw before threhold");
//            } 
//            ImageByte mask = wats[i].thresholdRangeInclusive((i+1), (i+1));
//            
//            //mask.erode(2, 2);
//            int[] hist=raws[i].getHistogram(null);            
//            double thre=thr.getThreshold(AutoThresholder.Method.Otsu, hist);
////            IJ.log("i "+thre);
//            thEx[i] = raws[i].thresholdAboveExclusive((float) thre);
//            thIn[i] = raws[i].thresholdAboveInclusive((float) thre);
//            if(i==5)
//            {
//                mask.show("mask");
//                thEx[i].show("ex");
//                thIn[i].show("in");
//                segs1[i].show("seg");
//            } 
//            thEx[i].intersectMask(mask);
//            thIn[i].intersectMask(mask);
//            
////            IJ.log("size raw: "+raws[i].sizeX+" "+raws[i].sizeY + " "+raws[i].sizeZ);
////            IJ.log("size thEx: "+thEx[i].sizeX+" "+thEx[i].sizeY + " "+thEx[i].sizeZ);
////            IJ.log("size Seg1: "+segs1[i].sizeX+" "+segs1[i].sizeY + " "+segs1[i].sizeZ);
////            IJ.log("size Seg "+segs[i].sizeX+" "+segs[i].sizeY + " "+segs[i].sizeZ);
//            
////            int volOtsu = 0;
////            for (int z = 0; z < thEx[i].sizeZ; z++) {
////                for (int y = 0; y < thEx[i].sizeY; y++) {
////                    for (int x = 0; x < thEx[i].sizeX; x++) {
////                        if(thEx[i].getPixel(x, y, z)!=0)
////                        {
////                            volOtsu++;
////                        }
////                    }
////                }
////            }    
//            
////            Object3DVoxels[] objs1 =  thEx[i].getObjects3D();
////            for(int k=0; k < objs1.length; k++)
////            {
////                volOtsu += objs1[k].getVolumeUnit();
////            }
////            Object3DVoxels[] objs =  segs[i].getObjects3D();
////            for(int k=0; k < objs.length; k++)
////            {
////                volSeg += objs[k].getVolumeUnit();
////            }
////            float ratio = ((float)volSeg/(float)volOtsu);
////            if(ratio < 0.9)
////            {
////                countErrS++;
////                thExErrorS.add(thEx[i]);
////            }    
////            if(ratio > 4.0)
////            {
////                countErrB++;
////                thExErrorB.add(thEx[i]);
////            }  
//            
//            
//        }
        
        
        //}
//        IJ.log(" Nb err: smaller than otsu :" + countErrS + " bigger than otsu: "+countErrB);
//        int countS = 0;
//        ImageByte[] thExErrorS1 = new ImageByte[thExErrorS.size()]; 
//        for(ImageByte im : thExErrorS)
//        {
//            thExErrorS1[countS] = im;
//            countS++;
//            
//        }
//        
//        int countB = 0;
//        ImageByte[] thExErrorB1 = new ImageByte[thExErrorB.size()]; 
//        for(ImageByte im : thExErrorB)
//        {
//            thExErrorB1[countB] = im;
//            countB++;
//            
//        }
//        ImageShort resEx = ImageShort.merge3DBinary(thEx, wat.sizeX, wat.sizeY, wat.sizeZ);
//        resEx.show("Otsu Exclusive");
//        ImageShort resIn = ImageShort.merge3DBinary(thIn, wat.sizeX, wat.sizeY, wat.sizeZ);
//        resIn.show("Otsu Inclusive");
//        ImageShort resErrorS = ImageShort.merge3DBinary(thExErrorS1, wat.sizeX, wat.sizeY, wat.sizeZ);
//        resErrorS.show("errorS");
//        ImageShort resErrorB = ImageShort.merge3DBinary(thExErrorB1, wat.sizeX, wat.sizeY, wat.sizeZ);
//        resErrorB.show("errorB");
        
    }
    private ImageByte processRegion(int i, ImageInt wat, ImageInt raw, ImageInt seg, boolean debug) 
    {
        ImageByte th;
        AutoThresholder thr = new AutoThresholder();
        ImageByte mask = seg.thresholdRangeInclusive(i, i);
        ImageInt ma = null;
        //int[] hist = raw.getHistogram(ma, 65536, 0, 65535);
        double minI = raw.getMin();
        double maxI = raw.getMax();
        double sbin = (maxI - minI) / 256;
        int[] hist = raw.getHistogram(ma, 256, minI, maxI);
        double thre = thr.getThreshold(AutoThresholder.Method.Otsu, hist);
        double thre16 = minI + thre/2 * sbin;
        //
        IJ.log("-----------------------------------------------------------");
        if (debug) {
            IJ.log("Threshold for region " + i + " : " + thre + " " + thre16);
//            if (i == 10) {
//                raw.duplicate().show("raw_" + i);
//                mask.duplicate().show("mask_" + i);
//                wat.duplicate().show("wat_" + i);
//                seg.duplicate().show("seg_" + i);
//            }
        }
        if (thre16 < mean) {
            thre16 = 65535;
        }
        th = raw.thresholdAboveExclusive((float) thre16);
//        th = raw.thresholdAboveExclusive((float) thre);
        th.intersectMask(mask);
        ImageByte res = (ImageByte)FastFilters3D.filterIntImage(th, 7, 5, 5, 2, nbcpu, true);
//        if (i == 10) 
//        {
//            th.show("threshold img_"+i);
//            res.show("threshold_closing_"+i);
//        }
//        Object3DVoxels O = new Object3DVoxels(th);
//        if (!O.isConnex()) {
//            ArrayList<Object3DVoxels> conn = O.getConnexComponents();
//            //IJ.log(O + " is not connex " + conn.size());
//            int max = 0;
//            int maxi = 0;
//            for (int o = 0; o < conn.size(); o++) {
//                if (conn.get(o).getVolumePixels() > max) {
//                    max = conn.get(o).getVolumePixels();
//                    maxi = o;
//                }
//            }
//            for (int o = 0; o < conn.size(); o++) {
//                if (o != maxi) {
//                    conn.get(o).draw(th, 0);
//                }
//            }
//        }
        
//        int volOtsu = 0;
//        for (int z = 0; z < th.sizeZ; z++) {
//            for (int y = 0; y < th.sizeY; y++) {
//                for (int x = 0; x < th.sizeX; x++) {
//                    if(th.getPixel(x, y, z)!=0)
//                    {
//                        volOtsu++;
//                    }
//                }
//            }
//        } 
//        IJ.log("After intersect: "+volOtsu);
//        double volSeg = 0;
//        for (int z = 0; z < seg.sizeZ; z++) {
//            for (int y = 0; y < seg.sizeY; y++) {
//                for (int x = 0; x < seg.sizeX; x++) {
//                    if(seg.getPixel(x, y, z)!=0)
//                    {
//                        volSeg++;
//                    }
//                }
//            }
//        }
//        int volIntersect = 0;
//        for (int z = 0; z < wat.sizeZ; z++) {
//            for (int y = 0; y < wat.sizeY; y++) {
//                for (int x = 0; x < wat.sizeX; x++) {
//                    if(th.getPixel(x, y, z)!=0 && seg.getPixel(x, y, z)!=0)
//                    {
//                        volIntersect++;
//                    }
//                }
//            }
//        }    
//
//        float ratioDSC = (float)(2 * volIntersect)/(float)(volSeg + volOtsu);
//
//        IJ.log("  Region: "+ i +" threshold: "+thre+ "   vol Otsu: "+volOtsu + "   vol Seg: "+volSeg
//              +"  ratio DSC: "+ ratioDSC + "  vol intersection: " + volIntersect);
        return th;
        
        // need to remove small region
        
    }

}
