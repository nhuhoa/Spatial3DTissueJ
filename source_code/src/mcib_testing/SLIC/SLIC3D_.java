/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.SLIC;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
//import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import java.io.File;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;


/**
 *
 * @author thomasb
 */
public class SLIC3D_ implements PlugIn 
{    
    private int min = 80;
    private int max = 300;
    private int stepX = 16;
    private int stepY = 16;
    private int stepZ = 8;
    private int m2 = 10000;
    private int ite = 10;
    public String dir = null;
    public String save_dir = "";
    
    @Override
    public void run(String string) {
        
//        int nbima = WindowManager.getImageCount();
        
        
//        for (int i = 0; i < nbima; i++) {
//            imgs[i] = ImageHandler.wrap(WindowManager.getImage(i + 1));
//        }
        
        
        int[] wList = WindowManager.getIDList();
        if (wList==null) {
            IJ.error("No images are open.");
            return;
        }
        
        String[] titles = new String[wList.length+1];
        titles[0] = "*None*";
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            titles[i+1] = imp!=null?imp.getTitle():"";
        }
        ImagePlus temp = WindowManager.getImage(wList[wList.length-1]);
        if (null !=temp){
            WindowManager.setTempCurrentImage(temp);
            save_dir = IJ.getDirectory("image");
        }else{
            save_dir = "";
        }
        
        GenericDialogPlus gd = new GenericDialogPlus("3D SLIC Clustering");
        gd.addMessage("    Spatial3DTissueJ  ");
//        gd.addMessage("See and quote reference:\n A novel toolbox to investigate tissue\nspatial" +
//        "organization applied to \nthe study of the islets of Langerhans");
        gd.addMessage("Input : markers staining channels");
        gd.addMessage("Output: SLIC clustering label");
        gd.addMessage(" \n");
//        gd.addStringField("Save Dir: ", save_dir);
        gd.addDirectoryField("Save_Dir: ", save_dir, 20);
        gd.addChoice("Marker_1 ", titles, titles[0]);
        gd.addChoice("Marker_2 ", titles, titles[0]);
        gd.addChoice("Marker_3 ", titles, titles[0]);
        gd.addChoice("Marker_4 ", titles, titles[0]);
        gd.addNumericField("min_size", min, 0, 10, "");
        gd.addNumericField("max_size", max, 0, 10, "");
//        gd.addNumericField("Ratio M", m2, 0);
        gd.addNumericField("nb_iterations", ite, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        save_dir = gd.getNextString();
        if ("".equals(save_dir)) 
        {
                return;
        }
        
        int[] index = new int[4];
        index[0] = gd.getNextChoiceIndex();
        index[1] = gd.getNextChoiceIndex();
        index[2] = gd.getNextChoiceIndex();
        index[3] = gd.getNextChoiceIndex();
        min = (int) gd.getNextNumber();
        max = (int) gd.getNextNumber();
//        m2 = (int) gd.getNextNumber();
        ite = (int) gd.getNextNumber();
        
        if (!save_dir.endsWith("/"))
            save_dir = save_dir + "/";
        File wdir = new File(save_dir);
        if (!wdir.exists()) { //!wdir.isDirectory() ||  || !wdir.canRead()
            wdir.mkdirs();
        }
        

        IJ.log("SLIC Multi-SuperVoxels Clustering ...");
        IJ.log("Input Images: ");
        int nb = 0;
        for(int k=0; k<index.length; k++)
        {
            if(index[k]!=0)
            {
                nb++;
            }    
        } 
        ImageHandler[] imgs = new ImageHandler[nb];
        IJ.log("Nb input images: " + nb);
        int count = 0;
        for(int k=0; k<index.length; k++)
        {
            if(index[k]!=0)
            {
                imgs[count] = ImageHandler.wrap(WindowManager.getImage(wList[index[k]-1]));
                IJ.log("Img "+(count+1)+"  : "+imgs[count].getTitle());
                count++;
            }    
        }  
        if (imgs.length==0) {
            IJ.error("Need to choose at least one input image");
            return;
        }
//        WindowManager.setTempCurrentImage(imgs[0].getImagePlus());
//        dir = IJ.getDirectory("image");
        
        
        ImageInt slic = getSLICClustering(imgs);
        slic.show();
//        IJ.selectWindow(slic.getTitle());
        IJ.saveAs("Tiff", save_dir+slic.getTitle()+".tif");
        IJ.log("Save SLIC clustering result as "+slic.getTitle()+".tif ");
        IJ.log("into the folder: "+save_dir);
        
//        ImageHandler img = imgs[0];
//        ImageHandler[] labelsF = new ImageHandler[imgs.length];
//        for (int c = 0; c < imgs.length; c++) {
//            labelsF[c] = new ImageShort("labelF-" + imgs[c].getTitle(), img.sizeX, img.sizeY, img.sizeZ);
//        }
//        
//        Objects3DPopulation pop = new Objects3DPopulation(slic);
//        for (int i = 0; i < pop.getNbObjects(); i++) {
//            Object3DVoxels O = (Object3DVoxels) pop.getObject(i);
//            for (int c = 0; c < imgs.length; c++) {
//                O.draw(labelsF[c], (int) O.getPixMedianValue(imgs[c]));
//            }
//        }
//        for (int c = 0; c < imgs.length; c++) 
//        {
//            labelsF[c].show();
//        }

//        // TEST POP
//        SLICPopulation popSLIC = new SLICPopulation(slic);
//        popSLIC.getHistogram(img).getPlot().show();
//        
//        ImageInt maxlo = new ImageShort("maxLocal", img.sizeX, img.sizeY, img.sizeZ);
//        for (int o = 0; o < popSLIC.getNbObjects(); o++) {
//            IJ.showStatus("Processing " + o);
//            Object3D O = popSLIC.getObject(o);
//            double V = O.getPixMeanValue(img);
//            ArrayList<Integer> ne = popSLIC.getNeighbors(o);
//            boolean maxlocal = true;
//            for (int i : ne) {
//                if (popSLIC.getObjectByValue(i).getPixMeanValue(img) > V) {
//                    maxlocal = false;
//                    break;
//                }
//            }
//            if (maxlocal) 
//            {
//                O.draw(maxlo, 255);
//                O.setType(1);
//            }
//        }
//        
//        for (int o = 0; o < popSLIC.getNbObjects(); o++) {
//            Object3D O = popSLIC.getObject(o);
//            if (O.getType() == 0) {
//                continue;
//            }
//            ArrayList<Integer> ne = popSLIC.getNeighbors(o);
//            for (int i : ne) {
//                popSLIC.getObjectByValue(i).draw(maxlo, 128);
//            }            
//        }
//        
//        maxlo.show();
    }
    private ImageInt getSLICClustering(ImageHandler[] imgs)
    {
//        Calibration cal = imgs[0].getCalibration();
        Calibration cal = imgs[0].getImagePlus().getCalibration().copy();
        int nbcube = (int) (max / (cal.pixelWidth * cal.pixelHeight * cal.pixelDepth));
        double ratio = cal.pixelDepth / cal.pixelWidth;
        SLIC3D slicP = new SLIC3D(min, max, nbcube, ratio, m2);
        ImageInt slic = slicP.getSLICImage(imgs, ite);
        slic.setTitle("combined_markers_SLIC");
        return slic;
    }        
    
}
