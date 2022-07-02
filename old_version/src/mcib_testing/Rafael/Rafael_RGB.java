/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Rafael;

import ij.IJ;
import ij.WindowManager;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;

/**
 *
 * @author thomasb
 */
public class Rafael_RGB implements ij.plugin.PlugIn {

    private final int UNLABELED = 0;
    private final int ALPHA = 1;
    private final int BETA = 2;
    private final int DELTA = 3;

    ImageHandler[] signals;

    @Override
    public void run(String arg) {
        ImageInt imgAlpha = ImageInt.wrap(WindowManager.getImage("alpha.tif"));
        ImageInt imgBeta = ImageInt.wrap(WindowManager.getImage("beta.tif"));
        ImageInt imgDelta = ImageInt.wrap(WindowManager.getImage("delta.tif"));
        ImageInt imgDAPI=ImageInt.wrap(WindowManager.getImage("DAPI-seg.tif"));
        ImageByte binNuc=imgDAPI.thresholdAboveInclusive(1);
        
        // FILTERS
        
        imgAlpha.intersectMask(binNuc);
        imgBeta.intersectMask(binNuc);
        imgDelta.intersectMask(binNuc);
        
        

        int[] TA = new int[8];
        int[] TB = new int[8];
        int[] TC = new int[8];
        TA[0] = (int) imgAlpha.getMax();
        TB[0] = (int) imgBeta.getMax();
        TC[0] = (int) imgDelta.getMax();
        int[] errors;
//        errors = computeN0(imgAlpha, imgBeta, imgDelta, TA[0], TB[0], TC[0]);
//        for (int i = 0; i < errors.length; i++) {
//            IJ.log("" + errors[i]);
//        }

        int error = 0;
        int step = 20;
        while (error != -1) {
            error = -1;
            IJ.log("");
            TA[1] = TA[0] - step;
            TB[1] = TB[0];
            TC[1] = TC[0];
            TA[2] = TA[0];
            TB[2] = TB[0] - step;
            TC[2] = TC[0];
            TA[3] = TA[0];
            TB[3] = TB[0];
            TC[3] = TC[0] - step;
            TA[4] = TA[0] - step;
            TB[4] = TB[0] - step;
            TC[4] = TC[0];
            TA[5] = TA[0] - step;
            TB[5] = TB[0];
            TC[5] = TC[0] - step;
            TA[6] = TA[0];
            TB[6] = TB[0] - step;
            TC[6] = TC[0] - step;
            TA[7] = TA[0] - step;
            TB[7] = TB[0] - step;
            TC[7] = TC[0] - step;
            int maxCount = 0;
            int maxId = -1;
            int maxcoloc = 0;
            int maxpc = 0;
            for (int i = 1; i < 8; i++) {
                errors = computeN0(imgAlpha, imgBeta, imgDelta, TA[i], TB[i], TC[i]);
                int countCell = errors[1] + errors[2] + errors[3];
                int coloc = errors[4];
                int pc;
                if (coloc == 0) {
                    pc = Integer.MAX_VALUE;
                } else {
                    pc = countCell / coloc;
                }
                if ((countCell >= maxCount) && (coloc<1000)) {
                    maxCount = countCell;
                    maxId = i;
                    maxcoloc = coloc;
                    maxpc = pc;
                }
            }
            if (maxId > 0) {
                TA[0] = TA[maxId];
                TB[0] = TB[maxId];
                TC[0] = TC[maxId];
                error = 0;
            } else {
                if (step > 1) {
                    step = step / 2;
                    error = 0;
                }
            }
            IJ.log("" + maxCount + " : " + TA[0] + " " + TB[0] + " " + TC[0] + " " + maxcoloc + " " + maxpc);
        }
        IJ.log("Finished");
    }

    private boolean checkTT(ImageInt A, ImageInt B, ImageInt C, int TA, int TB, int TC) {
        int w = A.sizeXYZ;
        for (int i = 0; i < w; i++) {
            int a = A.getPixelInt(i);
            int b = B.getPixelInt(i);
            int c = C.getPixelInt(i);
            if ((a >= TA) && (b >= TB)) {
                return false;
            }
            if ((b >= TB) && (c >= TC)) {
                return false;
            }
            if ((a >= TA) && (c >= TC)) {
                return false;
            }
        }
        return true;
    }

    private int[] computeN0(ImageInt A, ImageInt B, ImageInt C, int TA, int TB, int TC) {
        int[] NN = new int[5];
        int w = A.sizeXYZ;
        for (int i = 0; i < w; i++) {
            int a = A.getPixelInt(i);
            int b = B.getPixelInt(i);
            int c = C.getPixelInt(i);
            // noise
            if ((a < TA) && (b < TB) && (c < TC)) {
                NN[0]++;
            }
            // A, B,C
            if ((a >= TA) && (b < TB) && (c < TC)) {
                NN[1]++;
            }
            if ((a < TA) && (b >= TB) && (c < TC)) {
                NN[2]++;
            }
            if ((a < TA) && (b < TB) && (c >= TC)) {
                NN[3]++;
            }
            // coloc
            if (((a >= TA) && (b >= TB)) || ((b >= TB) && (c >= TC)) || ((a >= TA) && (c >= TC))) {
                NN[4]++;
            }

        }
        return NN;
    }

}
