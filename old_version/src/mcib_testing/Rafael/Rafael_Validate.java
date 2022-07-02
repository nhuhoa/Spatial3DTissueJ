/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Rafael;

import ij.IJ;
import ij.ImagePlus;
import java.util.ArrayList;
import java.util.Collections;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;

/**
 *
 * @author thomasb
 */
public class Rafael_Validate implements ij.plugin.PlugIn {

    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    ImageHandler img;

    @Override
    public void run(String arg) {
        ImagePlus plus = IJ.getImage();
        //String directory = plus.getFileInfo().directory;
        img = ImageHandler.wrap(plus);
        String directory = IJ.getDirectory("image");

        IJ.log("Reading data from " + directory);
        //popRegions = new Objects3DPopulation();
        popNuclei = new Objects3DPopulation();
        //popRegions.loadObjects(directory + "Regions.zip");
        popNuclei.loadObjects(directory + "SEG/Nuclei.zip");

        // read file and update from user
        Manager3DType manager = new Manager3DType();
        manager.directory = directory+"SEG/";
        ArrayList<Object3D> sortNuc = popNuclei.getObjectsList();
        Collections.sort(sortNuc, new CompareNucleus());
        Objects3DPopulation popNucleiSort = new Objects3DPopulation();
        popNucleiSort.addObjects(sortNuc);
        manager.addObjects3DPopulation(popNucleiSort);
        manager.load();

        IJ.log("Finished");

    }
}
