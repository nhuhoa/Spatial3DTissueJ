/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Segmentation;



import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.regionGrowing.Watershed3DVoronoi;

/**
 * Description of the Class
 *
 * @author thomas
 * @created 16 avril 2006
 */
public class WatershedSegmentation_ implements PlugIn {
    // seeds
    ImagePlus seedPlus = null;
    float radMax = 0;
    boolean showLines = false;
    boolean showEDT = false;
    boolean save = true;
    String dir = null;
    @Override
    public void run(String arg) {
        if (Dialog())
            IJ.log("====================================================================");
            WatershedVoronoi();
            IJ.log("====================================================================");
    }

    private void WatershedVoronoi() {
        IJ.log("");
        long t = System.currentTimeMillis();
        Watershed3DVoronoi watershed3DVoronoi = new Watershed3DVoronoi(ImageInt.wrap(seedPlus), radMax);
        watershed3DVoronoi.setLabelSeeds(false);
        ImageInt wat = watershed3DVoronoi.getVoronoiZones(showEDT);
        wat.show("dapi-seg-wat");
        if (showLines) watershed3DVoronoi.getVoronoiLines(true).show("VoronoiLines");
        if(save)
        {
            wat.setTitle("dapi-seg-wat");
            IJ.selectWindow(wat.getTitle());
            IJ.saveAs("Tiff", dir+"dapi-seg-wat.tif");
            IJ.log("Save the result as dapi-seg-wat.tif");
            IJ.log("Save the result into the folder : "+dir);
            
        }  
        IJ.log("Finished in " + (System.currentTimeMillis() - t) + " ms.");
    }

    private boolean Dialog() 
    {
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return false;
        }
        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i] = imp!=null?imp.getTitle():"";
        }
        ImagePlus img = WindowManager.getImage(wList[0]);
        Calibration calibration = img.getCalibration();
        String unit = "pixel";
        if (calibration != null) unit = calibration.getUnits();
        GenericDialog gd = new GenericDialog("3D Cell Zone");
        gd.addMessage("3D Tissue Spatial Analysis");
        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : nuclei segmented image");
        gd.addMessage("Output: watershed - cell zone image.");
        gd.addMessage(" ");
        gd.addChoice("Nuclei Segmented Image : ", titles, titles[0]);
        gd.addNumericField("Radius_Max (0 for no max)", radMax, 2, 10, unit);
        gd.addCheckbox("Show_EDT", showEDT);
        gd.addCheckbox("Show_Lines", showLines);
        gd.addCheckbox("Save the output into folder", true);
        gd.showDialog();
        int idxImg = gd.getNextChoiceIndex();
        radMax = (float) gd.getNextNumber();
        if (radMax == 0) radMax = Float.MAX_VALUE;
        showEDT = gd.getNextBoolean();
        showLines = gd.getNextBoolean();
        save = gd.getNextBoolean();
        seedPlus = WindowManager.getImage(wList[idxImg]);
        WindowManager.setTempCurrentImage(seedPlus);
        dir = IJ.getDirectory("image");
        return gd.wasOKed();
    }
}