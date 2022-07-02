/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Evaluation;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.util.ArrayList;
import java.util.Iterator;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
/**
 *
 * @author tranhoa
 */
public class CompareDetection_ implements PlugIn
{
    public double distMin = 1, distMax = 10;
    @Override
    public void run(String string) 
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
        if (wList.length < 3) {
            IJ.error("Open at least 3 images.");
            return;
        }
//        String none = "*None*";
//        titles[wList.length] = none;

        GenericDialog gd = new GenericDialog("Evaluation Segmentation");
        gd.addChoice("Method_1_Result:", titles, titles[0]);
        gd.addChoice("Method_2_Result:", titles, titles[1]);
        gd.addChoice("Original_Image:", titles, titles[2]);
        gd.addNumericField("Min distance: ", distMin, 1);
        gd.addNumericField("Max distance: ", distMax, 10);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        int[] index = new int[3];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        distMin = (int) gd.getNextNumber();
        distMax = (int) gd.getNextNumber();
        ImageHandler imgJaza = ImageHandler.wrap(WindowManager.getImage(wList[index[0]]));
        ImageHandler imgH = ImageHandler.wrap(WindowManager.getImage(wList[index[1]]));
        ImageHandler imgOri = ImageHandler.wrap(WindowManager.getImage(wList[index[2]]));
        IJ.log("Initialing ...");
        IJ.log("----------------------------------------------------------");
        ImageHandler imgOrigin = imgOri.duplicate();
        compareDetection(imgJaza, imgH, imgOrigin);
        IJ.log("----------------------------------------------------------");
        IJ.log("Finished ...");
    }
    
    private void compareDetection(ImageHandler imgJaza, ImageHandler imgH, ImageHandler imgOrigin)
    {
        IJ.log("Min distance: "+ distMin + "    Max distance: "+ distMax);
        IJ.log("----------------------------------------------------------");
        
        Objects3DPopulation popJ = new Objects3DPopulation((ImageInt)imgJaza);
        Objects3DPopulation popH = new Objects3DPopulation((ImageInt)imgH);
        int sizeJ = popJ.getNbObjects();
        int sizeH = popH.getNbObjects();
        IJ.log("Nb detected objs in "+imgJaza.getTitle()+" is: "+sizeJ+" in "+imgH.getTitle()+" is: "+sizeH);
        IJ.log("----------------------------------------------------------");
//        Objects3DPopulation popJ1 = new Objects3DPopulation((ImageInt)imgJaza);
//        Objects3DPopulation popH1 = new Objects3DPopulation((ImageInt)imgH);
        int count = 0;
        //int sizePop = popH.getNbObjects();
        ImageHandler labelJ = new ImageShort("label_correspond_"+imgJaza.getTitle(), imgJaza.sizeX, imgJaza.sizeY, imgJaza.sizeZ);
        ImageHandler labelJErr = new ImageShort("label_err_"+imgJaza.getTitle(), imgJaza.sizeX, imgJaza.sizeY, imgJaza.sizeZ);
        ImageHandler labelH = new ImageShort("label_correspond_"+imgH.getTitle(), imgJaza.sizeX, imgJaza.sizeY, imgJaza.sizeZ);
        ImageHandler labelHErr = new ImageShort("label_err_"+imgH.getTitle(), imgJaza.sizeX, imgJaza.sizeY, imgJaza.sizeZ);
        
        ArrayUtil arrH = new ArrayUtil(popH.getNbObjects());
        for (int k = 0; k < popH.getNbObjects(); k++) 
        {
            arrH.putValue(k, popH.getObject(k).getValue());
        }
        float dilatedRad = 2;
        float dilatedRadErr = 4;
        float dilatedRadErr1 = 6;
        for (int k = 0; k < popJ.getNbObjects(); k++) 
        {
            Object3D O = popJ.getObject(k);
            Object3D neiObj;
            neiObj = closestCenterV2(O, popH, distMin, distMax);
            boolean flag = false;
            if(neiObj!=null)
            {
                for(int j=0;j<arrH.getSize();j++)
                {
                    if(arrH.getValue(j)==neiObj.getValue())
                    {
                        flag = true;
                        arrH.putValue(j, 0);
                        neiObj.draw(labelH,neiObj.getValue());
                        count++;
                        IJ.log("("+imgJaza.getTitle()+")      obj:"+O.getValue()
                                +"   ----  "+" obj: "+neiObj.getValue()+ "      ("+imgH.getTitle()+")");
                        Point3D C1 = neiObj.getCenterAsPoint();
                        ArrayList<Voxel3D> arr1 = new ArrayList<Voxel3D>();
                        Voxel3D vox1 = new Voxel3D(C1, 5000);
                        arr1.add(vox1);
                        Object3DVoxels obj1 = new Object3DVoxels(arr1);
                        Object3D obj11 = obj1.getDilatedObject(dilatedRad, dilatedRad, dilatedRad);
                        obj11.draw(imgOrigin, 5000);
                    }    
                }    
                
                Point3D C2 = O.getCenterAsPoint();
                ArrayList<Voxel3D> arr2 = new ArrayList<Voxel3D>();
                Voxel3D vox2 = new Voxel3D(C2, 6000);
                arr2.add(vox2);
                Object3DVoxels obj2 = new Object3DVoxels(arr2);
                Object3D obj22 = obj2.getDilatedObject(dilatedRad, dilatedRad, dilatedRad);
                obj22.draw(imgOrigin, 6000);
                O.draw(labelJ, O.getValue());
            }
            if(!flag)
            {
                O.draw(labelJErr,O.getValue());
                Point3D C3 = O.getCenterAsPoint();
                ArrayList<Voxel3D> arr3 = new ArrayList<Voxel3D>();
                Voxel3D vox3 = new Voxel3D(C3, 15000);
                arr3.add(vox3);
                Object3DVoxels obj3 = new Object3DVoxels(arr3);
                Object3D obj33 = obj3.getDilatedObject(dilatedRadErr, dilatedRadErr, dilatedRad);
                obj33.draw(imgOrigin, 15000);
                IJ.log("----------------------------------------------------------");
                IJ.log("Not correspond in ("+ imgJaza.getTitle()+") obj val : " + O.getValue());
                
            }
            
        }
        IJ.log("----------------------------------------------------------");
        int count1 = 0;
        for(int j=0;j<arrH.getSize();j++)
        {
            if(arrH.getValue(j)!=0)
            {
                count1++;
                Object3D objH = popH.getObjectByValue((int)arrH.getValue(j));
                objH.draw(labelHErr, objH.getValue());
                Point3D C4 = objH.getCenterAsPoint();
                ArrayList<Voxel3D> arr4 = new ArrayList<Voxel3D>();
                Voxel3D vox4 = new Voxel3D(C4, 16000);
                arr4.add(vox4);
                Object3DVoxels obj4 = new Object3DVoxels(arr4);
                Object3D obj44 = obj4.getDilatedObject(dilatedRadErr, dilatedRadErr, dilatedRad);
                obj44.draw(imgOrigin, 16000);
                IJ.log("----------------------------------------------------------");
                IJ.log("Not correspond in ("+ imgH.getTitle()+")   obj val : " + objH.getValue());
            }    
        } 
        
        IJ.log("-------------------SUMMARY---------------------------------------");
        IJ.log("Nb objs have problem in "+ imgJaza.getTitle()+" is: "+(sizeJ-count));
        IJ.log("Nb objs have problem in "+ imgH.getTitle()+" is: "+(sizeH-count)+" count: "+count1);
        
        labelJ.show();
        labelH.show();
        labelJErr.show();
        labelHErr.show();
        imgOrigin.show("RemarkedSegmentation");
        IJ.log("----------------------------------------------------------");
        IJ.log("Nb objs correspondant in 2 images : "+count);
        
    }
    public Object3D closestCenterV2(Object3D O, Objects3DPopulation pop, double distMin, double distMax) 
    {
        ArrayList<Object3D> objects = pop.getObjectsList();
        Point3D P = O.getCenterAsPoint();
        Object3D res = null;
        Object3D tmp;
        double d;

        for (Iterator e = objects.iterator(); e.hasNext();) {
            tmp = (Object3D) e.next();
            d = tmp.distPixelCenter(P.getX(), P.getY(), P.getZ());
            if ((d < distMax) && (d > distMin)) {
                distMax = d;
                res = tmp;
            }
        }

        return res;
    }
   
}
