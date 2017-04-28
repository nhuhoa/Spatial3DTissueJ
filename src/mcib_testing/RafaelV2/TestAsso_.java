/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

import ij.IJ;
import ij.WindowManager;
import java.util.ArrayList;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.regionGrowing.AllRegionsAssociation;
import mcib3d.image3d.regionGrowing.AssociationRegion;

/**
 *
 * @author tranhoa
 */
public class TestAsso_ implements ij.plugin.PlugIn 
{
    
    private Objects3DPopulation pop, popNuclei = null;
    ImageInt imgLabel;
    
    public void run(String arg) 
    {
        IJ.log("Start test all region");
        imgLabel = ImageInt.wrap(WindowManager.getImage("label.tif"));
        testAssos(imgLabel);
        IJ.log("Finish");
    }
    private void testAssos(ImageInt regionLabel)
    {
        //pop = new Objects3DPopulation(regionLabel);
        AllRegionsAssociation asso = new AllRegionsAssociation();
        asso.computeAllRegionsAssociationPairs(regionLabel, -1);

        IJ.log("asso " + asso.getNbLabels() + " " + asso.getListAssociation().size() + " " + asso.getAssociation(asso.getListAssociation().size() / 2));

        int max = 0;
        int r1 = -1;
        int r2 = -1;
        ArrayList<AssociationRegion> list = asso.getListAssociation();
        for (AssociationRegion A : list) {
            int i1 = A.getRegion(0);
            int i2 = A.getRegion(1);
            IJ.log("---------------------------------------------------");
            IJ.log("index 2 regions: " + i1 +" " + i2);
            Object3DVoxels O1 = (Object3DVoxels) pop.getObjectByValue(i1);
            Object3DVoxels O2 = (Object3DVoxels) pop.getObjectByValue(i2);
            if((O1==null)||(O2==null)) continue;
            int co = O1.edgeContact(O2, 1);
            if (co > max) {
                max = co;
                r1 = i1;
                r2 = i2;
                IJ.log("Max contact " + max + " " + r1 + " " + r2);
            }
            IJ.log("---------------------------------------------------");
        }
    }
    
}
