/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Segmentation;



import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import java.io.File;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.regionGrowing.Watershed3DVoronoi;

/**
 * Description of the Class
 *
 * @author based on script from thomas boudier, with modification from Hoa Tran
 * @modified 31 July 2022
 */
public class WatershedSegmentation_ implements PlugIn {
    // seeds
    ImagePlus seedPlus = null;
    float radMax = 3;
    boolean showLines = false;
    boolean showEDT = false;
    boolean save = true;
//    String dir = null;
    public String save_dir = "";
    @Override
    public void run(String arg) {
        if (Dialog()){
            IJ.log("====================================================================");
            WatershedVoronoi();
            IJ.log("====================================================================");
        }
        
    }

    private void WatershedVoronoi() {
        IJ.log("");
        long t = System.currentTimeMillis();
        Watershed3DVoronoi watershed3DVoronoi = new Watershed3DVoronoi(ImageInt.wrap(seedPlus), radMax);
        watershed3DVoronoi.setLabelSeeds(false);
        ImageInt wat = watershed3DVoronoi.getVoronoiZones(showEDT);
        
        
        if(wat!=null){
            if (showLines) watershed3DVoronoi.getVoronoiLines(true).show("VoronoiLines");
            String wat_fn = "dapi-seg-wat";
            wat.show(wat_fn);
            if(save)
            {
                wat.setTitle(wat_fn);
                IJ.selectWindow(wat_fn);
                IJ.saveAs("Tiff", save_dir+"dapi-seg-wat.tif");
                IJ.log("Save the result as dapi-seg-wat.tif");
                IJ.log("Save the result into the folder : "+save_dir);

            }  
        }else{
            IJ.log("Issue with cell zone estimation, please double check input segmented image");
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
        ImagePlus temp = WindowManager.getImage(wList[wList.length-1]);
        if (null !=temp){
            WindowManager.setTempCurrentImage(temp);
            save_dir = IJ.getDirectory("image");
        }
        
        ImagePlus img = WindowManager.getImage(wList[0]);
        Calibration calibration = img.getCalibration();
        String unit = "pixel";
        if (calibration != null) unit = calibration.getUnits();
        GenericDialogPlus gd = new GenericDialogPlus("3D Cell Zone Estimation");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("3D Tissue Spatial Analysis");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        " organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : nuclei segmented image");
        gd.addMessage("Output: watershed - cell zone image.");
        gd.addMessage(" \n");
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("Nuclei Segmented Image : ", titles, titles[0]);
        gd.addNumericField("Radius_Max (0 for no max)", radMax, 2, 10, unit);
        gd.addCheckbox("Show_EDT", showEDT);
        gd.addCheckbox("Show_Lines", showLines);
        gd.addCheckbox("Save the output into folder", true);
        gd.showDialog();
        
        if (gd.wasCanceled()) return false;
        
        
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return false;
        }
        int idxImg = gd.getNextChoiceIndex();
        radMax = (float) gd.getNextNumber();
        if (radMax == 0) radMax = Float.MAX_VALUE;
        showEDT = gd.getNextBoolean();
        showLines = gd.getNextBoolean();
        save = gd.getNextBoolean();
        seedPlus = WindowManager.getImage(wList[idxImg]);
//        WindowManager.setTempCurrentImage(seedPlus);
//        dir = IJ.getDirectory("image");
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        return gd.wasOKed();
    }
}