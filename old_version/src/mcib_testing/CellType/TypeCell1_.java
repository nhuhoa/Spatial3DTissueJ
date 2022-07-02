/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import java.util.ArrayList;
import java.util.HashMap;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;
import mcib_testing.Utils.Cell;

/**
 *
 * @author ttnhoa
 */
public class TypeCell1_ implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 3;
    private final int BETA = 2;
    private final int DELTA = 1;
//    private final int COLOCAB = 23;
//    private final int COLOCAD = 13;
//    private final int COLOCBD = 12;
   
    
    private Objects3DPopulation popRegions = null;
    private Objects3DPopulation popNuclei = null;

    //ImageHandler[] signals;
    
    ImageInt imgLabel;
    String dir = null;
    ArrayList<Cell> popCells = null;

    HashMap<Integer, Cell> region2Cell;
    HashMap<Integer, Cell> nucleus2Cell;

    @Override
    public void run(String arg) 
    {
        imgLabel = ImageInt.wrap(WindowManager.getImage("labelF.tif"));
        ImagePlus plus;
        plus = WindowManager.getImage("dapi-seg-wat.tif");
        
        WindowManager.setTempCurrentImage(plus);
        dir = IJ.getDirectory("image");
        //IJ.log("dir is: "+dir);
        double ratio = 0.1;
        double rDB = 0.8;
        double rAB = 0.85;
        ImageInt imgWat = ImageInt.wrap(plus);
        ImageInt imgSeg = ImageInt.wrap(WindowManager.getImage("dapi-seg.tif"));
        GenericDialog gd = new GenericDialog(" COMPUTE CELL TYPE ");
        gd.addNumericField(" Ratio ", ratio, 0);
//        gd.addNumericField(" Ratio DB ", rDB, 0);
//        gd.addNumericField(" Ratio AB ", rAB, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        ratio = (double) gd.getNextNumber();
//        rDB = (double) gd.getNextNumber();
//        rAB = (double) gd.getNextNumber();
        IJ.log("Initialization Tissue Analysis ...");
        initCells(imgSeg, imgWat);
        
        
        IJ.log("Type ...");
        IJ.log("Vol percent: "+ratio + "  rDB: "+rDB +"  rAB: "+rAB);
//        this.computeTypeCellsNucleus(imgLabel, ratio, rDB, rAB);
        this.computeTypeCellsNucleus(imgLabel, ratio, 1, 1);
        
        IJ.log("Save the result");
        popNuclei.saveObjects(dir + "Nuclei.zip");
        popRegions.saveObjects(dir + "Regions.zip");
        
        
        // draw 
//        drawCellTypes(false).show("TYPE_WAT");
//        drawCellTypes(true).show("TYPE NUC");
        
//        IJ.log("Association1 ...");
//        computeAssoCells(imgWat, 0);
//        computeAllContacts("DIRECT CONTACT");
        IJ.log("Finished");

    }
//    private void computeAllContacts(String name)
//    {
//        float aa = 0, bb = 0, dd = 0, ab = 0, ad = 0, bd = 0;
//        float na = 0, nb = 0, nd = 0;
//        ResultsTable rtContact = new ResultsTable();
//        for (Cell C : popCells) 
//        {
//            if (C.type == ALPHA) {
//                na++;
//                int[] ne = C.computeNei1TypeHisto();
//                aa += ne[ALPHA-1];
//                ab += ne[BETA-1];
//                ad += ne[DELTA-1];
//            } else if (C.type == BETA) {
//                nb++;
//                int[] ne = C.computeNei1TypeHisto();
//                ab += ne[ALPHA-1];
//                bb += ne[BETA-1];
//                bd += ne[DELTA-1];
//            } else if (C.type == DELTA) {
//                nd++;
//                int[] ne = C.computeNei1TypeHisto();
//                ad += ne[ALPHA-1];
//                bd += ne[BETA-1];
//                dd += ne[DELTA-1];
//            }
//        }
//        aa /= 2;
//        bb /= 2;
//        dd /= 2;
//        ab /= 2;
//        ad /= 2;
//        bd /= 2;
//        float sum = aa + bb + dd + ab + ad + bd;
//        float sumnb = na + nb + nd;
//        rtContact.incrementCounter();
//        int count = 0;
//        rtContact.setValue("unlabelled", count, (popCells.size()-sumnb));
//        rtContact.setValue("na", count, na);
//        rtContact.setValue("nb", count, nb);
//        rtContact.setValue("nd", count, nd);
//        rtContact.setValue("totalabd", count, sumnb);
//        
//        rtContact.setValue("pa", count, na / sumnb);
//        rtContact.setValue("pb", count, nb / sumnb);
//        rtContact.setValue("pd", count, nd / sumnb);
//        
//        rtContact.setValue("aa", count, aa);
//        rtContact.setValue("bb", count, bb);
//        rtContact.setValue("dd", count, dd);
//        rtContact.setValue("ab", count, ab);
//        rtContact.setValue("ad", count, ad);
//        rtContact.setValue("bd", count, bd);
//        rtContact.setValue("all", count, sum);
//        
//        rtContact.setValue("paa", count, aa / sum);
//        rtContact.setValue("pbb", count, bb / sum);
//        rtContact.setValue("pdd", count, dd / sum);
//        rtContact.setValue("pab", count, ab / sum);
//        rtContact.setValue("pad", count, ad / sum);
//        rtContact.setValue("pbd", count, bd / sum);
//        IJ.log("Nb : na=" + na + " nb=" + nb + " nd=" + nd + " sum=" + sumnb);
//        IJ.log("Nb : pa=" + na / sumnb + " pb=" + nb / sumnb + " pd=" + nd / sumnb);
//        IJ.log("Assos : aa=" + aa + " bb=" + bb + " dd=" + dd + " ab=" + ab + " ad=" + ad + " bd=" + bd + " all=" + sum);
//        IJ.log("Assos : aa=" + aa / sum + " bb=" + bb / sum + " dd=" + dd / sum + " ab=" + ab / sum + " ad=" + ad / sum + " bd=" + bd / sum);
//        rtContact.show("AllContact_"+name);
//    } 
    private void computeAllContacts(String name)
    {
        float aa = 0, bb = 0, dd = 0, ab = 0, ad = 0, bd = 0;
        float na = 0, nb = 0, nd = 0;
        ResultsTable rtContact = new ResultsTable();
        ResultsTable rtH = new ResultsTable();
        for (Cell C : popCells) 
        {
            if (C.type == ALPHA) {
                na++;
                int[] ne = C.computeNei1TypeHisto();
                aa += ne[ALPHA-1];
                ab += ne[BETA-1];
                ad += ne[DELTA-1];
            } else if (C.type == BETA) {
                nb++;
                int[] ne = C.computeNei1TypeHisto();
                ab += ne[ALPHA-1];
                bb += ne[BETA-1];
                bd += ne[DELTA-1];
            } else if (C.type == DELTA) {
                nd++;
                int[] ne = C.computeNei1TypeHisto();
                ad += ne[ALPHA-1];
                bd += ne[BETA-1];
                dd += ne[DELTA-1];
            }
        }
        aa /= 2;
        bb /= 2;
        dd /= 2;
        ab /= 2;
        ad /= 2;
        bd /= 2;
        float sum = aa + bb + dd + ab + ad + bd;
        float sumnb = na + nb + nd;
        rtContact.incrementCounter();
        
        int count = 0;
        rtContact.setValue("na", count, na);
        rtContact.setValue("nb", count, nb);
        rtContact.setValue("nd", count, nd);
        rtContact.setValue("totalabd", count, sumnb);
        
        rtContact.setValue("pa", count, na / sumnb);
        rtContact.setValue("pb", count, nb / sumnb);
        rtContact.setValue("pd", count, nd / sumnb);
        
        rtContact.setValue("aa", count, aa);
        rtContact.setValue("bb", count, bb);
        rtContact.setValue("dd", count, dd);
        rtContact.setValue("ab", count, ab);
        rtContact.setValue("ad", count, ad);
        rtContact.setValue("bd", count, bd);
        rtContact.setValue("all", count, sum);
        
        rtH.incrementCounter();
        rtH.setValue("aa", count, aa);
        rtH.setValue("bb", count, bb);
        rtH.setValue("dd", count, dd);
        rtH.setValue("ab", count, ab);
        rtH.setValue("ad", count, ad);
        rtH.setValue("bd", count, bd);
        rtH.setValue("all", count, sum);
        
        rtH.setValue("homotypic_alpha", count, (aa/(aa+ab+ad)));
        rtH.setValue("heterotypic_alpha", count, ((ab+ad)/(aa+ab+ad)));
        
        rtH.setValue("homotypic_beta", count, (bb/(bb+ab+bd)));
        rtH.setValue("heterotypic_beta", count, ((ab+bd)/(bb+ab+bd)));
        
        rtH.setValue("homotypic_delta", count, (dd/(dd+ad+bd)));
        rtH.setValue("heterotypic_delta", count, ((bd+ad)/(dd+ad+bd)));
        
        rtContact.setValue("paa", count, aa / sum);
        rtContact.setValue("pbb", count, bb / sum);
        rtContact.setValue("pdd", count, dd / sum);
        rtContact.setValue("pab", count, ab / sum);
        rtContact.setValue("pad", count, ad / sum);
        rtContact.setValue("pbd", count, bd / sum);
        
        IJ.log("Nb : na=" + na + " nb=" + nb + " nd=" + nd + " sum=" + sumnb);
        IJ.log("Nb : pa=" + na / sumnb + " pb=" + nb / sumnb + " pd=" + nd / sumnb);
        IJ.log("Assos : aa=" + aa + " bb=" + bb + " dd=" + dd + " ab=" + ab + " ad=" + ad + " bd=" + bd + " all=" + sum);
        IJ.log("Assos : aa=" + aa / sum + " bb=" + bb / sum + " dd=" + dd / sum + " ab=" + ab / sum + " ad=" + ad / sum + " bd=" + bd / sum);
        rtContact.show("AllContact_"+name);
        rtH.show("Homo-Hetero Contact");
    }
    private void computeAssoCells(ImageHandler img, int borderValue) {

        ArrayList<Integer>[] nei = computeAsso(img, borderValue);

        // index refers to region label
        for (int i = 0; i < nei.length; i++) {
            if (nei[i].isEmpty()) {
                continue;
            }
            Cell C = region2Cell.get(i);
            if (C != null) {
                C.nei1 = new ArrayList<Cell>();
                for (int j : nei[i]) {
                    Cell C1 = region2Cell.get(j);
                    if ((C1 != null) && (C1 != C)) {
                        C.nei1.add(C1);
                    }
                }
            }
        }
    }
    private ArrayList<Integer>[] computeAsso(ImageHandler img, int BorderValue) {
        //ImagePlus imp = WindowManager.getCurrentImage();
        //ImageHandler img = ImageHandler.wrap(imp);
        int max = (int) img.getMax();

        ArrayList<Integer>[] nei = new ArrayList[max + 1];
        for (int i = 0; i < nei.length; i++) {
            nei[i] = new ArrayList();
        }

        for (int x = 0; x < img.sizeX; x++) {
            for (int y = 0; y < img.sizeY; y++) {
                for (int z = 0; z < img.sizeZ; z++) {
                    if (((BorderValue >= 0) && (img.getPixel(x, y, z) == BorderValue)) || (BorderValue < 0)) {
                        ArrayUtil tab = img.getNeighborhood3x3x3(x, y, z);
                        tab = tab.distinctValues();
                        for (int j = 0; j < tab.getSize(); j++) {
                            int val = (int) tab.getValue(j);
                            if (!tab.hasOnlyValuesInt(nei[val])) {
                                for (int i = 0; i < tab.getSize(); i++) {
                                    if (!nei[val].contains((int) tab.getValue(i))) {
                                        nei[val].add((int) tab.getValue(i));
//                                        if (val == 9) {
//                                            IJ.log(" " + tab);
//                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return nei;
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
                Cell cell = new Cell();
                cell.region = region;
                cell.nucleus = popNuclei.getObjectByValue(nuc);
                popCells.add(cell);
                cell.id = c++;
                region2Cell.put(region.getValue(), cell);
                nucleus2Cell.put(nuc, cell);
            }
            
                  
        }
        //IJ.log("number of error: " + count);
        IJ.log("Number of cells is: "+popCells.size());
    }

    public int getCellType2(int vA, int vB, int vD, double vol, double ratio, double rDB, double rAB)
    {
        //double ratioDel = 0.25;
        //double ratioDel = 0.2;
        if(vA >= vB && vA > vD && vA/vol >= ratio)
        {
            return ALPHA;
        } 
        else if(vB > vA && vB > vD && vB/vol >= ratio)
        {
            if(vD >= vB * rDB && vD>vA && vD/vol >= ratio)  //0.5
            {
                return DELTA;
            }
            if(vA>=vB * rAB && vA>vD && vA/vol >= ratio)
            { 
                return ALPHA;
            }
            else{
                return BETA;
            }
        }
        else if(vD >= vA && vD >= vB && vD/vol >= ratio)
        {
            return DELTA;
        }
        
        else
        {
            return UNLABELED;
        }  
        
    }
//    public int getCellType(int vA, int vB, int vD)
//    {
//        if(vA > vB && vA > vD)
//        {
//            if(vD >= vA * 0.5)
//            {
//                return DELTA;
//            }
//            else if(vB >= vA*0.95)
//            {
//                return COLOCAB;
//            }
//            else if(vD >= vA*0.7)
//            {
//                return COLOCAD;
//            }
//            else{
//                return ALPHA;
//            }
//        }
//        else if(vB > vA && vB > vD)
//        {
//            if(vD >= vB * 0.45)
//            {
//                return DELTA;
//            }
//            else if(vA >= vB*0.8)
//            {
//                return ALPHA;
//            }
//            else if(vA >= vB*0.6)
//            {
//                return COLOCAB;
//            }
//            else if(vD >= vB*0.35)
//            {
//                return COLOCBD;
//            }
//            else{
//                return BETA;
//            }
//        }
//        else if(vD > vA && vD > vB)
//        {
//            return DELTA;
//        }
//        else if(vA==vB && vA>vD)
//        { 
//            return COLOCAB;
//        }
//        else if((vA==vD && vA>vB))
//        {
//            return COLOCAD;
//        }
//        else if(vD==vB && vD>vA)
//        {
//            return COLOCBD;
//        } 
//        else
//        {
//            return UNLABELED;
//        }  
//        
//    }      
    private void computeTypeCellsNucleus(ImageInt label, double ratio, double rDB, double rAB) 
    {
        int cA=0, cB=0, cD=0, cUnlabelled=0, cColocAB=0, cColocAD=0, cColocBD=0;
        ResultsTable cellM = new ResultsTable();

            for (Cell C : popCells) 
            {
                Object3D nucleus = C.nucleus;
                double vol = nucleus.getVolumeUnit();
                Object3D region = C.region;
                int reVal = region.getValue();
                cellM.incrementCounter();
                cellM.setValue("cell", C.id, reVal);
                cellM.setValue("id", C.id, C.id);
                if(nucleus==null) continue;
                ArrayUtil list = nucleus.listValues(label);
                int nbA=0; int nbB=0; int nbD=0; int nbU=0;
                for(int i=0;i<list.getSize();i++)
                {
                    if(list.getValue(i)==ALPHA)nbA++;
                    if(list.getValue(i)==BETA)nbB++;
                    if(list.getValue(i)==DELTA)nbD++;                
                }
                C.type = UNLABELED;
                //double ratio = 0.2;
//                if(nbD>0 && (nbD/vol)>0.15)
//                {
//                    IJ.log("");
//                    IJ.log("Cell val: "+nucleus.getValue()+" NbD: "+nbD+" nbA: "+nbA+" nbB: "+nbB+" percent: "+(nbD/vol));
//                    if(nbB>0)
//                    {
//                        IJ.log("------percent delta/beta : "+(nbD/nbB));
//                    }
//                }
                int typeCell = getCellType2(nbA, nbB, nbD, vol, ratio, rDB, rAB);
                C.type = (byte)typeCell;
                cellM.setValue("method1", C.id, typeCell);
                if(typeCell==UNLABELED){cUnlabelled++;}
                if(typeCell==ALPHA){
                    cA++;
//                    IJ.log("  ");
//                    IJ.log("al: "+100 * (nbA/vol));
                }
                if(typeCell==BETA){
                    cB++;
//                    IJ.log("  ");
//                    IJ.log("be: "+100 * (nbB/vol));
                }
                if(typeCell==DELTA){
                    cD++;
//                    IJ.log("  ");
//                    IJ.log("del: "+100 * (nbD/vol));
                }
//                if(typeCell==COLOCAB){cColocAB++;}
//                if(typeCell==COLOCAD){cColocAD++;}
//                if(typeCell==COLOCBD){cColocBD++;}
                String content = nbA +"  "+nbB +"  "+nbD +"  "+typeCell +"\n";
                //out.write(content);
                nucleus.setType(C.type);
                nucleus.setName("Nuc" + C.id);
                nucleus.setValue(C.id);
                region.setType(C.type);
                region.setName("Reg" + C.id);
                region.setValue(C.id);
            }

        //cellM.show("Method1");
        
        IJ.log("Nb of each cells type: A: "+cA+"   B: "+cB+"   D:"+cD+"   Unlabelled: "+cUnlabelled);
        ImageHandler unLabelledCell;
        unLabelledCell = new ImageShort("Unlabelled", label.sizeX, label.sizeY, label.sizeZ);
//        ImageHandler typeImg;
//        typeImg = new ImageShort("TYPE", label.sizeX, label.sizeY, label.sizeZ);
        for (Cell C : popCells) {
            Object3D nucleus = C.nucleus;
            if(C.type==UNLABELED)
            {  
               nucleus.draw(unLabelledCell, 255);
            }
//            else
//            {
//                nucleus.draw(typeImg, C.type);
//            }
        }
        unLabelledCell.show();
//        typeImg.show();
    }

    
    private ImageHandler drawCellTypes(boolean nuc) {
        ImageHandler draw = new ImageShort("TYPE_", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);
//        ImageHandler drawDelta = new ImageShort("DELTA", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);
//        ImageHandler drawDeltaCenter = new ImageShort("DELTA CENTER", imgLabel.sizeX, imgLabel.sizeY, imgLabel.sizeZ);
        ArrayList<Voxel3D> arr = new ArrayList<Voxel3D>();
        for (Cell C : popCells) {
            if (nuc) {
                C.nucleus.draw(draw, C.type);
            } else {
//                if(C.type==DELTA)
//                {
//                    C.region.draw(drawDelta, C.type);
//                    Point3D center = C.nucleus.getCenterAsPoint();
//                    Voxel3D v = new Voxel3D(center, 255);
//                    arr.add(v);
//                }
                C.region.draw(draw, C.type);
            }
        }
//        Object3DVoxels obj = new Object3DVoxels(arr);
//        obj.draw(drawDeltaCenter);
//        drawDelta.show();
//        drawDeltaCenter.show();
        return draw;
    }

}
