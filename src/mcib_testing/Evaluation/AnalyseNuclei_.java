/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Evaluation;


import ij.IJ;
import ij.WindowManager;
import java.util.ArrayList;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib_testing.Utils.Cell;

/**
 *
 * @author ttnhoa
 */
public class AnalyseNuclei_ implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
    private final int DAPI = 4;
    private final int COLOCAB = 23;
    private final int COLOCAD = 13;
    private final int COLOCBD = 12;
   
    
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    
    //ImageInt imgLabel;
    String dir = null;
    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;
    ImageInt imgSeg = null;
    @Override
    public void run(String arg) 
    {
        imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
  
        IJ.log("Initialization Tissue Analysis ...");
        //ImageHandler im = removeObjs(imgSeg);
        //removeObjs(imgSeg);
        separateObjects(imgSeg);
        IJ.log("Finished");

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
                if(popNuclei.getObjectByValue(nuc)!=null)
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
            
                  
        }
        //IJ.log("number of error: " + count);
        IJ.log("Number of cells is: "+popCells.size());
    }
    
    
    private void separateObjects(ImageInt nucLabel)
    {
        //ImageInt cloneNuc = nucLabel.duplicate();
        Objects3DPopulation popNuclei1 = new Objects3DPopulation(nucLabel);
        
        for(Object3D nuc : popNuclei1.getObjectsList())
        {
            // separating objects
            //set value 0 for contour
            //float radDilate = (float)0.5;
            ArrayList<Voxel3D> cont = nuc.getContours();
            Object3DVoxels contObj = new Object3DVoxels(cont);
            //Object3DVoxels contDilated = contObj.getDilatedObject(radDilate, radDilate, radDilate);
            //ArrayList<Voxel3D> lsDilated = contDilated.getVoxels(); 
            contObj.draw(nucLabel, 0);
//            for(Voxel3D vox : cont)
//            {
//                Point3D p = vox.getPosition();
//                cloneNuc.setPixel(p, 0);
//            }
            
            
        }
        
        nucLabel.show();
    }
    
    private ImageHandler removeObjs(ImageInt nucLabel)
    {
        ImageHandler label = new ImageShort("label", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        ImageHandler labelBorder = new ImageShort("label border", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        ImageHandler labelSmall = new ImageShort("label small", imgSeg.sizeX, imgSeg.sizeY, imgSeg.sizeZ);
        Objects3DPopulation popNuc = new Objects3DPopulation(nucLabel);
        Objects3DPopulation popNucV2 = new Objects3DPopulation(nucLabel);
        int c = 0, cV = 0, cUse=0;
        IJ.log("Nb obj before:"+ popNucV2.getNbObjects());
        for(Object3D nuc : popNuc.getObjectsList())
        {
            if(nuc.getVolumePixels()<1000)
            {
                cV++;
                nuc.draw(labelSmall, nuc.getValue());
            }
            else if((nuc.getVolumePixels()<2000) && (nuc.getXmin()<5 || nuc.getXmax()<5 || nuc.getXmax()>(imgSeg.sizeX-5) 
                        || nuc.getYmin()<5 || nuc.getYmax()<5 || nuc.getYmax()>(imgSeg.sizeY-5)))
            {
                    c++;
                    nuc.draw(labelBorder, nuc.getValue());
            }
            else{
                cUse++;
                nuc.draw(label, nuc.getValue());
            }
        }
        IJ.log("Nb of small objs: "+cV);
        IJ.log("Nb of border objs: "+c);
        IJ.log("Nb of objs used: "+cUse);
//        label.show();
//        labelBorder.show();
//        labelSmall.show();
        return label;
        
    }

   
}

