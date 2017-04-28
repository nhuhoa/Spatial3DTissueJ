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
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.image3d.regionGrowing.AllRegionsAssociation;
import mcib3d.image3d.regionGrowing.AssociationRegion;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */
public class Validate_Nuclei_ implements ij.plugin.PlugIn 
{

    private Objects3DPopulation pop, popNuclei = null;
    ImageInt imgLabel;
    String dir = null;
    ImageHandler img1, img2, img3, img4;
    int thresh = 0;
    int thres1 = 15;
    double errRate = 0.1;
    private final int UNLABELED = 0;
    private final int ALPHA = 1;
    private final int BETA = 2;
    private final int DELTA = 3;
    private final int DAPI = 4;
    
    public void run(String arg) 
    {
//        img1 = ImageHandler.wrap(WindowManager.getImage("C1-delta.tif"));
//        img2 = ImageHandler.wrap(WindowManager.getImage("C2-beta.tif"));
//        img3 = ImageHandler.wrap(WindowManager.getImage("C1-alpha.tif"));
        //img4 = ImageHandler.wrap(WindowManager.getImage("C4-dapi.tif"));
        ImagePlus plus1 = WindowManager.getImage("C1-delta.tif");
        img1 = ImageHandler.wrap(plus1);
        ImagePlus plus2 = WindowManager.getImage("C2-beta.tif");
        img2 = ImageHandler.wrap(plus2);
        ImagePlus plus3 = WindowManager.getImage("C1-alpha.tif");
        img3 = ImageHandler.wrap(plus3);
        ImagePlus plus4 = WindowManager.getImage("C4-dapi.tif");
        img4 = ImageHandler.wrap(plus4);
        //dir = IJ.getDirectory("image");
        imgLabel = ImageInt.wrap(WindowManager.getImage("label.tif"));
        
        IJ.log("Beginning, normalise the value of pixels");
        //validateCellTypeV3(imgLabel);
        normalizeValuePixels(img1, true);
        normalizeValuePixels(img2, true);
        normalizeValuePixels(img3, true);
        normalizeValuePixels(img4, true);
        img1.show("after");
        img2.show("after");
        img3.show("after");
        img4.show("after");
        IJ.log("Start create SLIC image");
        validateCellTypeV4(imgLabel);
        IJ.log("Finished");
    }
   
    private void test()
    {
        img1.multiplyByValue((float) (1000.0 / img1.getMax()));
        img2.multiplyByValue((float) (1000.0 / img2.getMax()));
        img3.multiplyByValue((float) (1000.0 / img3.getMax()));
        img4.multiplyByValue((float) (1000.0 / img4.getMax()));
        
        img1.resetStats(null);
        img2.resetStats(null);
        img3.resetStats(null);
        img4.resetStats(null);
        
        IJ.log("img1 max: " + img1.getMax());
        IJ.log("img2 max: " + img2.getMax());
        IJ.log("img3 max: " + img3.getMax());
        IJ.log("img4 max: " + img4.getMax());
    }
    public void normalizeValuePixels(ImageHandler imgOrig, boolean excludeZeros) 
    {
        // int count = 0;
        ArrayList<VoxS> idxList = new ArrayList<VoxS>();
        //VoxEVF[] idx = new VoxEVF[mask.countMaskVolume()];
        //double volume = idx.length;
        double minDist = Double.NEGATIVE_INFINITY;
        if (excludeZeros) 
        {
            minDist = 0;
        }
        
        
        for (int x = 0; x < imgOrig.sizeX; x++) {
            for (int y = 0; y < imgOrig.sizeY; y++) {
                for (int z = 0; z < imgOrig.sizeZ; z++) {
                    if (imgOrig.getPixel(x, y, z) > minDist) {
                        idxList.add(new VoxS(imgOrig.getPixel(x, y, z), x, y, z));
                        //idx[count] = new VoxEVF(distanceMap.pixels[z][xy], xy, z);
                        //count++;
                    } else {
                        imgOrig.setPixel(x, y, z, 1.0f);
                    }
                }    
            }
        }
        if (idxList.isEmpty()) {
            return;
        }
        VoxS[] idx = new VoxS[idxList.size()];
        idx = (VoxS[]) idxList.toArray(idx);
        double volume = idx.length;
        Arrays.sort(idx);
        double maxV = idx[idx.length - 1].distance;
        IJ.log("Max value is : " + maxV);
        for (int i = 0; i < idx.length; i++) 
        {
            imgOrig.setPixel(idx[i].x, idx[i].y, idx[i].z, (float) (1000.0 * idx[i].distance / maxV));
        }
    }
    private void validateCellTypeV3(ImageInt regionLabel)
    {
        pop = new Objects3DPopulation(regionLabel);
        ImageHandler label1 = new ImageShort("label1", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label2 = new ImageShort("label2", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label3 = new ImageShort("label3", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label4 = new ImageShort("label4", img1.sizeX, img1.sizeY, img1.sizeZ);

        ImageHandler labelF = new ImageShort("labelF", img1.sizeX, img1.sizeY, img1.sizeZ);

//        img1.multiplyByValue((float) (1000 / img1.getMax()));
//        img2.multiplyByValue((float) (1000 / img2.getMax()));
//        img3.multiplyByValue((float) (1000 / img3.getMax()));
//        img4.multiplyByValue((float) (1000 / img4.getMax()));
        
        int countE=0, countR = 0;
        for (int i = 0; i < pop.getNbObjects(); i++) 
        {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            double valC1 = O.getPixMedianValue(img1);
            double valC2 = O.getPixMedianValue(img2);
            double valC3 = O.getPixMedianValue(img3);
            double valC4 = O.getPixMedianValue(img4);
            
            ArrayUtil ac1 = O.listValues(img1);            
            int nbTotalC1 = ac1.countValueAbove(thres1);
            ArrayUtil ac2 = O.listValues(img2);            
            int nbTotalC2 = ac2.countValueAbove(thres1);
            ArrayUtil ac3 = O.listValues(img3);            
            int nbTotalC3 = ac3.countValueAbove(thres1);
            ArrayUtil ac4 = O.listValues(img4);            
            int nbTotalC4 = ac4.countValueAbove(thres1);
            int nbTotal = nbTotalC1 + nbTotalC2 + nbTotalC3 + nbTotalC4;
            
            ArrayList<regionCell> arrList = new ArrayList<regionCell>();
            regionCell r1 = new regionCell(valC1, nbTotalC1, 1);
            regionCell r2 = new regionCell(valC2, nbTotalC2, 2);
            regionCell r3 = new regionCell(valC3, nbTotalC3, 3);
            regionCell r4 = new regionCell(valC4, nbTotalC4, 4);
            arrList.add(r1);
            arrList.add(r2);
            arrList.add(r3);
            arrList.add(r4);
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
                /*if(maxValMedian < arrList.get(k).getMedian())
                {
                    maxValMedian = arrList.get(k).getMedian();
                    indexM = k;
                } */
            }    
            
            //int max1 = arrayList.get(arrayList.size() - 1);
            //int max2 = arrayList.get(arrayList.size() - 2);
            //IJ.log("nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
              //      "max: " + max1 + " " + max2);
//            O.draw(label1, (int) valC1);
//            O.draw(label2, (int) valC2);
//            O.draw(label3, (int) valC3);
//            O.draw(label4, (int) valC4);
                              
            if ((valC1 > thresh) && (valC1 > valC2) && (valC1 > valC3) && (valC1 > valC4)) 
            {
                //O.draw(labelF, 1);
                //O.draw(label1, 255);
                if((nbTotalC1 == maxVolume) || (nbTotalC1 < maxVolume && ((maxVolume-nbTotalC1) < (errRate * maxVolume))) )
                {
                    O.draw(labelF, 1);
                    O.draw(label1, 255);
                }    
                else{
                    if(arrList.get(indexV-1).getMValue()>0)
                    {
                        countR++;
                        O.draw(labelF, indexV);
                        if(indexV==2){O.draw(label2, 255); IJ.log("C1 --> C2");}
                        if(indexV==3){O.draw(label3, 255); IJ.log("C1 --> C3");}
                        if(indexV==4){O.draw(label4, 255); IJ.log("C1 --> C4");}
                    }
                    else{
                        O.draw(labelF, 5);
                        countE++;
                        IJ.log("C1, nbTotal: " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
                        "Max: " + maxVolume + " idx: "+ indexV + "  ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
                    }
                    
                }
//                if(nbTotalC1 < maxVolume && ((maxVolume-nbTotalC1) > (errRate * maxVolume)))
//                {
//                    countE++;
//                    IJ.log("C1, nbTotal: " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
//                    "   Max: " + maxVolume + " idx: "+ indexV + "  ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
//                }
            }
            if ((valC2 > thresh) && (valC2 > valC1) && (valC2 > valC3) && (valC2 > valC4)) {
                //O.draw(labelF, 2);
                //O.draw(label2, 255);
                if((nbTotalC2 == maxVolume) || (nbTotalC2 < maxVolume && ((maxVolume-nbTotalC2) < (errRate * maxVolume))) )
                {
                    O.draw(labelF, 2);
                    O.draw(label2, 255);
                }    
                else{
                    if(arrList.get(indexV-1).getMValue()>0)
                    {
                        countR++;
                        O.draw(labelF, indexV);
                        if(indexV==1){O.draw(label1, 255); IJ.log("C2 --> C1");}
                        if(indexV==3){O.draw(label3, 255); IJ.log("C2 --> C3");}
                        if(indexV==4){O.draw(label4, 255); IJ.log("C2 --> C4");}
                    }
                    else
                    {
                        countE++;
                        O.draw(labelF, 5);
                        IJ.log("C2, nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
                        " Max: " + maxVolume + " idx: "+ indexV + "  ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
                    }
                    
                }
//                if(nbTotalC2 < maxVolume && ((maxVolume-nbTotalC2) > errRate * maxVolume))
//                {
//                    countE++;
//                    IJ.log("C2, nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
//                    " Max: " + maxVolume + " idx: "+ indexV + "  ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
//                }
            }
            if ((valC3 > thresh) && (valC3 > valC1) && (valC3 > valC2) && (valC3 > valC4)) {
//                O.draw(labelF, 3);
//                O.draw(label3, 255);
//                if(nbTotalC3 < maxVolume && ((maxVolume-nbTotalC3) > errRate * maxVolume))
//                {
//                    countE++;
//                    IJ.log("C3, nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
//                    "   Max: " + maxVolume + " idx: "+ indexV + "   ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
//                }
                if((nbTotalC3 == maxVolume) || (nbTotalC3 < maxVolume && ((maxVolume-nbTotalC3) < (errRate * maxVolume))) )
                {
                    O.draw(labelF, 3);
                    O.draw(label3, 255);
                }    
                else{
                    if(arrList.get(indexV-1).getMValue()>0)
                    {
                        countR++;
                        O.draw(labelF, indexV);
                        if(indexV==1){O.draw(label1, 255); IJ.log("C3 --> C1");}
                        if(indexV==2){O.draw(label2, 255); IJ.log("C3 --> C2");}
                        if(indexV==4){O.draw(label4, 255); IJ.log("C3 --> C4");}
                    }
                    else
                    {
                        countE++;
                        O.draw(labelF, 5);
                        IJ.log("C3, nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
                        " Max: " + maxVolume + " idx: "+ indexV + "  ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
                    }
                    
                }
            }
            if ((valC4 > thresh) && (valC4 > valC1) && (valC4 > valC2) && (valC4 > valC3)) {
//                O.draw(labelF, 4);
//                O.draw(label4, 255);
//                if(nbTotalC4 < maxVolume && ((maxVolume-nbTotalC4) > errRate * maxVolume))
//                {
//                    countE++;
//                    IJ.log("C4, nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
//                    "  Max: " + maxVolume + " idx: "+ indexV + "   ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
//                }
                if((nbTotalC4 == maxVolume) || (nbTotalC4 < maxVolume && ((maxVolume-nbTotalC4) < (errRate * maxVolume))) )
                {
                    O.draw(labelF, 4);
                    O.draw(label4, 255);
                }    
                else{
                    if(arrList.get(indexV-1).getMValue()>0)
                    {
                        countR++;
                        O.draw(labelF, indexV);
                        if(indexV==1){O.draw(label1, 255); IJ.log("C4 --> C1");}
                        if(indexV==3){O.draw(label3, 255); IJ.log("C4 --> C3");}
                        if(indexV==2){O.draw(label2, 255); IJ.log("C4 --> C2");}
                    }
                    else
                    {
                        countE++;
                        O.draw(labelF, 5);
                        IJ.log("C4, nbTotal " + nbTotalC1 + " " + nbTotalC2 + " " + nbTotalC3 + " " + nbTotalC4 + 
                        " Max: " + maxVolume + " idx: "+ indexV + "  ValC:  " + valC1+"  "+valC2+"  "+valC3+"  "+valC4);
                    }
                    
                }
            }
            
//            
//            if ((valC1 > thresh) && (nbTotalC1 > nbTotalC2) && (nbTotalC1 > nbTotalC3) && (nbTotalC1 > nbTotalC4)) {
//                O.draw(labelF, 1);
//                O.draw(label1, 255);                
//            }
//            if ((valC2 > thresh) && (nbTotalC2 > nbTotalC1) && (nbTotalC2 > nbTotalC3) && (nbTotalC2 > nbTotalC4)) {
//                O.draw(labelF, 2);
//                O.draw(label2, 255);
//            }
//            if ((valC3 > thresh) && (nbTotalC3 > nbTotalC1) && (nbTotalC3 > nbTotalC2) && (nbTotalC3 > nbTotalC4)) {
//                O.draw(labelF, 3);
//                O.draw(label3, 255);
//            }
//            if ((valC4 > thresh) && (nbTotalC4 > nbTotalC1) && (nbTotalC4 > nbTotalC2) && (nbTotalC4 > nbTotalC3)) {
//                O.draw(labelF, 4);
//                O.draw(label4, 255);
//            }
//            
            
            
        }
        IJ.log("Count error: " + countE + " resolve: " + countR);
        label1.show();
        label2.show();
        label3.show();
        label4.show();
        labelF.show();
        
    
    
    }
    
    private void validateCellType(ImageInt regionLabel)
    {
        pop = new Objects3DPopulation(regionLabel);
        ImageHandler label1 = new ImageShort("label1", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label2 = new ImageShort("label2", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label3 = new ImageShort("label3", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label4 = new ImageShort("label4", img1.sizeX, img1.sizeY, img1.sizeZ);

        ImageHandler labelF = new ImageShort("labelF", img1.sizeX, img1.sizeY, img1.sizeZ);

        img1.multiplyByValue((float) (1000 / img1.getMax()));
        img2.multiplyByValue((float) (1000 / img2.getMax()));
        img3.multiplyByValue((float) (1000 / img3.getMax()));
        img4.multiplyByValue((float) (1000 / img4.getMax()));

        for (int i = 0; i < pop.getNbObjects(); i++) {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            double valC1 = O.getPixMedianValue(img1);
            double valC2 = O.getPixMedianValue(img2);
            double valC3 = O.getPixMedianValue(img3);
            double valC4 = O.getPixMedianValue(img4);

            if ((valC1 > thresh) && (valC1 > valC2) && (valC1 > valC3) && (valC1 > valC4)) {
                O.draw(labelF, 1);
                O.draw(label1, 255);   
                O.setType(ALPHA);
            }
            if ((valC2 > thresh) && (valC2 > valC1) && (valC2 > valC3) && (valC2 > valC4)) {
                O.draw(labelF, 2);
                O.draw(label2, 255);
                O.setType(BETA);
            }
            if ((valC3 > thresh) && (valC3 > valC1) && (valC3 > valC2) && (valC3 > valC4)) {
                O.draw(labelF, 3);
                O.draw(label3, 255);
                O.setType(DELTA);
            }
            if ((valC4 > thresh) && (valC4 > valC1) && (valC4 > valC2) && (valC4 > valC3)) {
                O.draw(labelF, 4);
                O.draw(label4, 255);
                O.setType(DAPI);
            }
            
        }
        
        AllRegionsAssociation a = new AllRegionsAssociation();
        int count = 0;
        for (int i = 0; i < (pop.getNbObjects()/500); i++) 
        {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            
            if(O.getType()==DAPI)
            {
                count++;
                ArrayList<Voxel3D> listVoxels = O.getVoxels();
                //IJ.log("count : " + count + "size: " + listVoxels.size());
                for(Voxel3D v : listVoxels)
                {
                    ArrayUtil tab = regionLabel.getNeighborhoodCross3D((int)v.getX(), (int)v.getY(), (int)v.getZ());
                    int val = (int) regionLabel.getPixel((int)v.getX(), (int)v.getY(), (int)v.getZ());
                    tab = tab.distinctValues();
                    //IJ.log("size of tab: "+ tab.getSize());
                    for (int j = 0; j < tab.getSize(); j++) {
                            int val2 = tab.getValueInt(j);
                            if (val2 != val) {
                                
                                AssociationRegion asso = new AssociationRegion();
                                asso.addRegion(val);
                                asso.addRegion(val2);
                                boolean f = a.addAssoRegion(asso);
                                if(f) {IJ.log("ok");}
                                
                            }
                    }
                
                }
            }  
        }
        
        IJ.log("Number of nuclei is: " + count);
        ArrayList<AssociationRegion> listAsso = a.getListAssociation();
        IJ.log("Number of association region is : "+ listAsso.size());
        
        
        for(AssociationRegion as : listAsso)
        {   
            String resStr = as.toString();
            IJ.log("ass is: "+ resStr);
            
            int begin = as.getFirst();
            Object3D obj = pop.getObjectByValue(begin);
            IJ.log("type of first region is : "+obj.getType());
            for (int i = 1; i < as.size(); i++) {
                int nei = as.getRegion(i);
                Object3D next = pop.getObjectByValue(nei);
                IJ.log("type of neig is: " + next.getType());
            }
        }
                
                
//        label1.show();
//        label2.show();
//        label3.show();
//        label4.show();
//        labelF.show();
        
    }
    
    private void validateCellTypeV4(ImageInt regionLabel)
    {
        pop = new Objects3DPopulation(regionLabel);
        ImageHandler label1 = new ImageShort("label1-delta", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label2 = new ImageShort("label2-beta", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label3 = new ImageShort("label3-alpha", img1.sizeX, img1.sizeY, img1.sizeZ);
        ImageHandler label4 = new ImageShort("label4-dapi", img1.sizeX, img1.sizeY, img1.sizeZ);

        ImageHandler labelF = new ImageShort("labelF", img1.sizeX, img1.sizeY, img1.sizeZ);

        for (int i = 0; i < pop.getNbObjects(); i++) {
            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
            double valC1 = O.getPixMedianValue(img1);
            double valC2 = O.getPixMedianValue(img2);
            double valC3 = O.getPixMedianValue(img3);
            double valC4 = O.getPixMedianValue(img4);

            if ((valC1 > thresh) && (valC1 > valC2) && (valC1 > valC3) && (valC1 > valC4)) {
                O.draw(labelF, 1);
                O.draw(label1, 255);   
                O.setType(ALPHA);
            }
            if ((valC2 > thresh) && (valC2 > valC1) && (valC2 > valC3) && (valC2 > valC4)) {
                O.draw(labelF, 2);
                O.draw(label2, 255);
                O.setType(BETA);
            }
            if ((valC3 > thresh) && (valC3 > valC1) && (valC3 > valC2) && (valC3 > valC4)) {
                O.draw(labelF, 3);
                O.draw(label3, 255);
                O.setType(DELTA);
            }
            if ((valC4 > thresh) && (valC4 > valC1) && (valC4 > valC2) && (valC4 > valC3)) {
                O.draw(labelF, 4);
                O.draw(label4, 255);
                O.setType(DAPI);
            }
            
        } 
                
        label1.show();
        label2.show();
        label3.show();
        label4.show();
        labelF.show();
        
    }
    
}
