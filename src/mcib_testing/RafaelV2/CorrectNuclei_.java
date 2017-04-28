/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;


import mcib_testing.Rafael.*;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.ResultsTable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DPoint;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */

public class CorrectNuclei_ implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;

    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    
    ImageInt imgLabel;
    String dir = null;

    @Override
    public void run(String arg) 
    {
        //ImagePlus plus = WindowManager.getImage("dapi-seg-wat.tif");
        
        //ImageInt imgWat = ImageInt.wrap(plus);
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
        //ImageInt imgNucleiOrigin = ImageInt.wrap(WindowManager.getImage("C4-dapi.tif"));
        
        //ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg-SLIC-25000.tif"));
        dir = IJ.getDirectory("image");
        imgLabel = ImageInt.wrap(WindowManager.getImage("label-re-C4-dapi.tif"));
    
        IJ.log("Correct the result of segmentation");
        this.correctNucleusSLIC(imgSeg, imgLabel);
        //this.computeTypeCellsNucleusSLIC(imgLabel);
        IJ.log("Finished");

    }


    private void correctNucleusSLIC(ImageInt nucLabel, ImageInt label)
    {
        IJ.log("Correct the result of segmentation using image label SLIC, result of SLIC");
        //int[] histype = {0};
        int rm = 0; 
        int ac = 0;
        popNuclei = new Objects3DPopulation(nucLabel);
        IJ.log("Number of cells detected at the first time : " + popNuclei.getNbObjects());
        
        //for(Object3D nuc : popNuclei.getObjectsList())
        for (int i = 0; i < popNuclei.getNbObjects(); i++)
        {
            Object3D nuc = popNuclei.getObject(i);
            int nbColoc = 0; 
            ArrayUtil list = nuc.listValues(label);
            int nbTotal = list.getSize();
            for(int k=0; k<nbTotal; k++)
            {
                if(list.getValue(k)!=0)
                {
                    nbColoc++;
                }
                
            }
            if(nbColoc >= nbTotal * 0.75)
            {
                //IJ.log("Percentage of exact segment: " + nb_dapi +" " + nbTotal);
                ac++;
                
            }
            else
            {
                //IJ.log("Error of segmentation, not exact a nucleus");
                rm++;
                popNuclei.removeObject(i);
           
            }
        }  
        IJ.log("Number of error of segmentation is: " + rm);
        IJ.log("Number of cells detected after the correction : " + popNuclei.getNbObjects());
        ImageHandler nucLabelDraw;
        nucLabelDraw = new ImageShort("DAPI Correct", nucLabel.sizeX, nucLabel.sizeY, nucLabel.sizeZ);

        for (Object3D nuc : popNuclei.getObjectsList()) {
            //C.nucleus.draw(draw, C.id);
            if(nuc==null) continue; 
            nuc.draw(nucLabelDraw, nuc.getValue());
        }
        nucLabelDraw.setTitle(nucLabel.getTitle()+"-correct");
        nucLabelDraw.show();
        
    }  
    
    
}
