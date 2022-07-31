/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.SLIC;

import ij.IJ;
import ij.measure.Calibration;
import java.util.ArrayList;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.ObjectCreator3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageShort;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author thomasb
 */
public class SLIC3D {

    int stepX, stepY, stepZ;
    int nbreg;
    ImageHandler[] imgs;
    int min, max, nbcube;
    double ratio;
    int m2, ite;

    public SLIC3D(int min, int max, int nbcube, double ratio, int m2) {
        this.min = min;
        this.max = max;
        this.ratio = ratio;
        this.nbcube = nbcube;
        this.m2 = m2;
    }

    public ImageInt getSLICImage(ImageHandler[] imgs, int nbit) {
        ImageHandler initSLICImg = initSLIC(imgs);
        initSLICImg.show();
        return processSLIC(imgs, initSLICImg, nbit);
    }

    private ImageHandler initSLIC(ImageHandler[] imgs) {

//        IJ.log("DEBUG: nbcube: " +nbcube);
        stepX = (int) Math.floor(Math.pow(nbcube * ratio, 0.333));
        // 2D 
        if (imgs[0].sizeZ == 1) {
            stepX = (int) Math.floor(Math.pow(nbcube * ratio, 0.5));
        }
        stepY = stepX;
        stepZ = (int) (stepX / ratio);

        // case 2D
        if (imgs[0].sizeZ == 1) {
            stepZ = 1;
        }
//        IJ.log("raw image sizeX: " + imgs[0].sizeX + " sizeY: " + imgs[0].sizeY + " sizeZ: " + imgs[0].sizeZ);    
        IJ.log("grid stepX: " + stepX + " stepY: " + stepY + " stepZ: " + stepZ);
//        ObjectCreator3D create = new ObjectCreator3D(imgs[0].sizeX, imgs[0].sizeY, imgs[0].sizeZ);
//        ObjectCreator3D_V2 create = new ObjectCreator3D_V2(imgs[0].sizeX, imgs[0].sizeY, imgs[0].sizeZ);
        ImageHandler img = new ImageShort("creator", imgs[0].sizeX, imgs[0].sizeY, imgs[0].sizeZ);
        nbreg = 1;
        for (int z = stepZ / 2; z <= imgs[0].sizeZ - stepZ / 2; z += stepZ) {
            for (int x = stepX / 2; x <= imgs[0].sizeX - stepX / 2; x += stepX) {
                for (int y = stepY / 2; y <= imgs[0].sizeY - stepY / 2; y += stepY) {
//                    create.createBrick(x, y, z, stepX / 2, stepY / 2, stepZ / 2, nbreg++);
                    
//                    create.createBrick(x, y, z, stepX / 2, stepY / 2, stepZ / 2, (float)nbreg);
                    img = createBrickV2(img, x, y, z, stepX / 2, stepY / 2, stepZ / 2, nbreg++);
                }
            }
        }
        IJ.log("Number of bricks: "+ nbreg);
        // create.getImageHandler().show();
        Objects3DPopulation test = new Objects3DPopulation((ImageInt) img);
        IJ.log("Nb regions: " + test.getNbObjects());
//        return create.getImageHandler();
        return img;
    }

    private ImageInt processSLIC(ImageHandler[] imgs, ImageHandler create, int ite) 
    {
        ImageHandler img = imgs[0];
        Objects3DPopulation pop = new Objects3DPopulation((ImageInt) create);
        ImageInt label = null;
        ImageFloat dist;
        //ImageLabeller labeller = new ImageLabeller();
        long t0 = System.currentTimeMillis();
        for (int it = 0; it < ite; it++) {
            IJ.log("SLIC " + it + " / " + (ite + 1) + ", " + pop.getNbObjects()+" regions");

            label = new ImageShort("label", img.sizeX, img.sizeY, img.sizeZ);
            label.fill(0);

            dist = new ImageFloat("dist", img.sizeX, img.sizeY, img.sizeZ);
            dist.fill(Float.MAX_VALUE);

            float S2 = (img.sizeXYZ / nbreg);
            double[] valsC = new double[imgs.length];

            for (int i = 0; i < pop.getNbObjects(); i++) {
                Object3D O = pop.getObject(i);
                Point3D C = O.getCenterAsPoint();
                for (int c = 0; c < imgs.length; c++) {
                    valsC[c] = O.getPixMedianValue(imgs[c]);
                }
                for (int z = C.getRoundZ() - stepZ; z < C.getRoundZ() + stepZ; z++) {
                    for (int x = C.getRoundX() - stepX; x < C.getRoundX() + stepX; x++) {
                        for (int y = C.getRoundY() - stepY; y < C.getRoundY() + stepY; y++) {
                            if (img.contains(x, y, z)) {
                                double dc2 = 0;
                                for (int c = 0; c < imgs.length; c++) {
                                    dc2 += (valsC[c] - imgs[c].getPixel(x, y, z)) * (valsC[c] - imgs[c].getPixel(x, y, z));
                                }
                                double ds2 = C.distanceSquare(new Point3D(x, y, z));
                                double D2 = dc2 + (ds2 / S2) * m2;
                                if (D2 < dist.getPixel(x, y, z)) {
                                    dist.setPixel(x, y, z, (float) D2);
                                    label.setPixel(x, y, z, i);
                                }
                            }
                        }
                    }
                }
            }
            // relabel, test connexity
            pop = new Objects3DPopulation();
//            pop.addImage(label, -1, new Calibration());
            pop.addImage(label, -1);
            int nb = pop.getNbObjects();
            for (int i = pop.getNbObjects() - 1; i >= 0; i--) 
            {
                ArrayList<Object3DVoxels> con = ((Object3DVoxels) (pop.getObject(i))).getConnexComponents();
                if ((con.size() > 1) || (pop.getObject(i).getValue() == 0)) 
                {
                    //IJ.log(pop.getObject(i) + " is not connex " + pop.getObject(i).getVolumePixels() + " " + con.size());
                    for (Object3DVoxels OO : con) 
                    {
                        OO.setValue(nb++);
                        pop.addObject(OO);
                        OO.draw(label);
                    }
                    pop.removeObject(i);
                }
            }

            pop.updateNamesAndValues();

            // CHECK VOLUMES
            ArrayList<Object3D> toRemove = new ArrayList();
            for (int i = pop.getNbObjects() - 1; i >= 0; i--) {
                Object3D OO = pop.getObject(i);
                if (OO.getVolumeUnit() < min) {
                    int ad = findAdjacentLableQuick(OO, label);
                    int ni = 100;
                    while ((ad == -1) && (ni > 0)) {
                        ad = findAdjacentLableQuick(OO, label);
                        ni--;
                    }
                    if (ad >= 0) {
                        //IJ.log("small " + OO + " " + OO.getVolumePixels() + " " + ad);
                        OO.draw(label, ad);
                        mergeRegions(pop, i, ad);
                        toRemove.add(OO);
                    }
                }
            }

            for (Object3D OO : toRemove) {
                pop.removeObject(OO);
            }
        }

        return label;
    }
    /**
     * Create a brick in Z direction with centre and radius
     *
     * @param centerx centre X
     * @param centery centre Y
     * @param centerz centre Z
     * @param rx      radius X
     * @param ry      radius Y
     * @param rz      radius Z
     * @param value   the value to draw
     */
    public ImageHandler createBrickV2(ImageHandler img, int centerx, int centery, int centerz, double rx, double ry, double rz, float value) {
        int startx = (int) (centerx - rx);
        if (startx < 0) {
            startx = 0;
        }
        int starty = (int) (centery - ry);
        if (starty < 0) {
            starty = 0;
        }
        int startz = (int) (centerz - rz);
        if (startz < 0) {
            startz = 0;
        }
        int endx = (int) (centerx + rx);
        if (endx >= img.sizeX) {
            endx = img.sizeX - 1;
        }
        int endy = (int) (centery + ry);
        if (endy >= img.sizeY) {
            endy = img.sizeY - 1;
        }
        int endz = (int) (centerz + rz);
        if (endz >= img.sizeZ) {
            endz = img.sizeZ - 1;
        }
        for (int k = startz; k <= endz; k++) {
            for (int j = starty; j <= endy; j++) {
                for (int i = startx; i <= endx; i++) {
                    img.setPixel(i, j, k, value);
                }
            }
        }
        return img;
    }

    private int findAdjacentLableQuick(Object3D obj, ImageHandler label) 
    {

        ArrayList<Voxel3D> cont = obj.getContours();
        if (cont.isEmpty()) {
//            IJ.log("Empty cont " + obj);
            return -1;
        }
        Voxel3D V = cont.get((int) (Math.random() * cont.size()));
        ArrayUtil tab = label.getNeighborhoodCross3D(V.getRoundX(), V.getRoundY(), V.getRoundZ());

        int p = 0;
//        while ((p < tab.getSize()) && (tab.getValue(p) == obj.getValue())) {
        while ((p < tab.size()) && (tab.getValue(p) == obj.getValue())) {
            p++;
        }

        int res = -1;
//        if (p < tab.getSize()) {
        if (p < tab.size()) {
            res = tab.getValueInt(p);
        }

        return res;
    }

    private void mergeRegions(Objects3DPopulation pop, int a, int b) {
        Object3D A = pop.getObject(a);
        Object3DVoxels B = (Object3DVoxels) pop.getObjectByValue(b);
        if(B!=null && A!=null)
        {
            B.addVoxels(A.getVoxels());
        }    
        
    }
}
