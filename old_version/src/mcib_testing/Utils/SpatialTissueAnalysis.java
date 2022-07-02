/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Utils;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DLabel;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;
/**
 *
 * @author tranhoa
 */
public class SpatialTissueAnalysis {
    private final int numEvaluationPoints;
    private final int numRandomSamples;
    private final double distHardCore;
    //int imaspots = 0;
    //int imamask = 1;
    //String[] names;
    //Point3D[] evaluationPoints;
    private final double env;
    // JEAN DB
    private double sdi_F = Double.NaN;
    private double sdi_G = Double.NaN;
    private double sdi_H = Double.NaN;
    private String desc = "";
    //private double sdi_K;
    ImageHandler randomPop;

    int nbBins = 1000;

    public SpatialTissueAnalysis(int numPoints, int numRandomSamples, double distHardCore, double env) {
        this.numEvaluationPoints = numPoints;
        this.numRandomSamples = numRandomSamples;
        this.distHardCore = distHardCore;
        this.env = env;
    }

//    public void process(ImagePlus plusSpots, ImagePlus plusMask, String functions, boolean verbose, boolean show, boolean save) {
//        process(ImageInt.wrap(plusSpots), ImageInt.wrap(plusMask), functions, verbose, show, save);
//    }
//
//    public void processAll(ImageHandler plusSpots, ImageHandler plusMask, boolean verbose, boolean show, boolean save) {
//        process(plusSpots, plusMask, "FGH", verbose, show, save);
//    }

//    public void processAll(ImagePlus plusSpots, ImagePlus plusMask, boolean verbose, boolean show, boolean save) {
//        process(plusSpots, plusMask, "FGH", verbose, show, save);
//    }
    /**
     * ttnhoa
     * Process the data in the existing list of elements (grid)
     * @param plusSpots
     * @param plusMask
     * @param lsCenterCells : list of object center 
     * @param functions
     * @param verbose
     * @param show
     * @param save 
     */
    public void process(ImagePlus plusSpots, ImagePlus plusMask, ArrayList<Point3D> lsCenterCells, String functions, boolean verbose, boolean show, boolean save, String dir) {
        processInGrid(ImageInt.wrap(plusSpots), ImageInt.wrap(plusMask), lsCenterCells, functions, verbose, show, save, dir);
    }
    
    //ttnhoa
    public void process(Objects3DPopulationGrid pop, ImagePlus plusMask, ArrayList<Point3D> lsCenterCells, String functions, boolean verbose, boolean show, boolean save, String directory, String desc) 
    {
        this.desc = desc;
        processInGrid(pop, ImageInt.wrap(plusMask), lsCenterCells, functions, verbose, show, save, directory);
    }
    /**
     * ttnhoa
     * @description: process elements in an existing list of elements (grid), instead of all space
     * @param plusSpots
     * @param plusMask
     * @param lsCenterCells
     * @param functions
     * @param verbose
     * @param show
     * @param save 
     */
    public void processInGrid(ImageHandler plusSpots, ImageHandler plusMask, ArrayList<Point3D> lsCenterCells, String functions, boolean verbose, boolean show, boolean save, String dir) 
    {
        Calibration calibration = plusSpots.getCalibration();
        if (calibration == null) 
        {
            IJ.log("Image not calibrated");
            calibration = new Calibration();
            calibration.setUnit("pix");
            calibration.pixelWidth = 1;
            calibration.pixelHeight = 1;
            calibration.pixelDepth = 1;
        }

        ImageInt inImage = (ImageInt) plusSpots;
        ImageInt segImage;
        if (inImage.isBinary(0)) {
            if (verbose) {
                IJ.log("Segmenting image...");
            }
            inImage = inImage.threshold(0, false, true);
            ImageLabeller labels = new ImageLabeller(false);
            segImage = labels.getLabels(inImage);
            if (verbose) {
                segImage.show("Labelled Image");
            }
        } else {
            segImage = (ImageInt) inImage.duplicate();
        }
        segImage.setCalibration(calibration);

        int nbspots;

        Objects3DPopulationGrid pop = new Objects3DPopulationGrid();
        ImageInt maskHandler = (ImageInt) plusMask;
        Object3D mask = new Object3DLabel(maskHandler, (int) maskHandler.getMax());
        mask.setCalibration(calibration);
        pop.setMask(mask);
        pop.addImage(segImage, calibration);
        pop.setCalibration(calibration);

        // random sample
        Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
        poprandom.setCalibration(calibration);
        poprandom.setMask(mask);
        poprandom.createRandomPopulationInGrid(lsCenterCells, pop.getNbObjects(), distHardCore);
        //poprandom.createRandomPopulation(pop.getNbObjects(), distHardCore);
        randomPop = maskHandler.createSameDimensions();
        randomPop.setCalibration(calibration);
        poprandom.draw(randomPop);
        //randomPop.show("random");

        if ((plusMask.getCalibration() == null) || (!plusMask.getCalibration().scaled())) {
            if (verbose) {
                IJ.log("mask not calibrated, calibrating ...");
            }
            plusMask.setCalibration(calibration);
            plusMask.getImagePlus().updateAndRepaintWindow();
        }

        nbspots = pop.getNbObjects();
        if (verbose) {
            IJ.log("Computing spatial statistics, please wait ...");
        }

        if (verbose) {
            IJ.log("Nb Spot=" + nbspots);
            IJ.log("Volume mask=" + mask.getVolumeUnit());
            IJ.log("Density=" + (nbspots / mask.getVolumeUnit()));
        }

        if (functions.contains("F")) {
            //processF(pop, mask, verbose, show, save);
            processFFunctionInGrid(pop, mask, lsCenterCells, verbose, show, save, dir);
        }

        if (functions.contains("G")) {
            processGFunctionInGrid(pop, mask, lsCenterCells, verbose, show, save, dir);
        }

        if (functions.contains("H")) {
            //processH(pop, mask, verbose, show, save);
            processHFunctionInGrid(pop, mask, lsCenterCells, verbose, show, save);
            
        }
    }
    private Point3D[] createEvaluationPoints(int numPoints, Objects3DPopulationGrid population) {
        Point3D[] evaluationPoints = new Point3D[numPoints];
        for (int i = 0; i < numPoints; ++i) {
            evaluationPoints[i] = population.getRandomPointInMask();
        }

        return evaluationPoints;
    }
    public double getSdi_F() {
        return sdi_F;
    }

    public double getSdi_G() {
        return sdi_G;
    }

    public double getSdi_H() {
        return sdi_H;
    }

    public ImageHandler getRandomSample() {
        return randomPop;
    }
    /**
     * ttnhoa
     * @description Random elements in an existing list of elements (grid), instead of all space 
     * @param pop
     * @param mask
     * @param lsCenterCells
     * @param verbose
     * @param show
     * @param save 
     */
    private void processFFunctionInGrid(Objects3DPopulationGrid pop, Object3D mask, ArrayList<Point3D> lsCenterCells, boolean verbose, boolean show, boolean save, String dir) 
    {
        Calibration calibration = mask.getCalibration();
        int nbspots = pop.getNbObjects();
        IJ.log("*******F_function_ Nb Spot=" + nbspots);
        ArrayUtil distances;
        // F
        Point3D[] evaluationPoints;
        ArrayUtil observedDistancesF;
        ArrayUtil[] sampleDistancesF;
        ArrayUtil xEvalsF;
        ArrayUtil observedCDF;
        ArrayUtil averageCDF;

//        evaluationPoints = createEvaluationPoints(numEvaluationPoints, pop);
//        evaluationPoints = createEvaluationPointsInList(lsCenterCells);
        evaluationPoints = createEvaluationPointsInList(numEvaluationPoints, pop, lsCenterCells);
        observedDistancesF = pop.computeDistancesV2(evaluationPoints);
        observedDistancesF.sort();
        observedCDF = CDFTools.cdf(observedDistancesF);
        xEvalsF = new ArrayUtil(numRandomSamples * numEvaluationPoints);
        sampleDistancesF = new ArrayUtil[numRandomSamples];

        ResultsTable rtDistance = new ResultsTable();
        for (int i = 0; i < observedDistancesF.getSize(); i++) 
        {
            rtDistance.incrementCounter();
            rtDistance.setValue("observed_X", i, observedDistancesF.getValue(i));
            rtDistance.setValue("observed_Y", i, observedCDF.getValue(i));
        }
        for (int i = 0; i < numRandomSamples; i++) {
            if (verbose) {
                IJ.showStatus("Random population F " + i);
            }
            Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
            poprandom.setCalibration(calibration);
            poprandom.setMask(mask);
            poprandom.createRandomPopulationInGrid(lsCenterCells, nbspots, distHardCore);
            //poprandom.createRandomPopulation(nbspots, distHardCore);
            poprandom.createKDTreeCenters();
            distances = poprandom.computeDistancesV2(evaluationPoints);
            distances.sort();
            sampleDistancesF[i] = distances;//           
            xEvalsF.insertValues(i * numEvaluationPoints, distances);
        }
        xEvalsF.sort();
        averageCDF = CDFTools.cdfAverage(sampleDistancesF, xEvalsF);

        ResultsTable rtSamples = new ResultsTable();
        for (int i = 0; i < numRandomSamples; i++) 
        {
            if (verbose) {
                IJ.showStatus("Random population F " + i);
            }
            Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
            poprandom.setCalibration(calibration);
            poprandom.setMask(mask);
            poprandom.createRandomPopulationInGrid(lsCenterCells, nbspots, distHardCore);
            //poprandom.createRandomPopulation(nbspots, distHardCore);
            poprandom.createKDTreeCenters();
            distances = poprandom.computeDistancesV2(evaluationPoints);
            distances.sort();
            sampleDistancesF[i] = distances;
            if(i==0)
            {
                for (int k = 0; k < sampleDistancesF[0].getSize(); k++) 
                {
                    rtSamples.incrementCounter();
                    rtSamples.setValue("sF_"+0, k, sampleDistancesF[0].getValue(k));
                }
            }
            else
            {
                for (int k = 0; k < sampleDistancesF[i].getSize(); k++) 
                {
                    rtSamples.setValue("sF_"+i, k, sampleDistancesF[i].getValue(k));
                }
            }
            
        }
        

        // plot
        Plot plotF = null;
        if (show || save) {
            plotF = createPlot(xEvalsF, sampleDistancesF, observedDistancesF, observedCDF, averageCDF, "Euclidean-Distance-F");
            
//            ResultsTable rtF = new ResultsTable();
//            for(int k=0; k<observedDistancesF.getSize(); k++)
//            {
//                rtF.incrementCounter();
//                rtF.setValue("observedX", k, observedDistancesF.getValue(k));
//                rtF.setValue("observedY", k, observedCDF.getValue(k));
//                rtF.setValue("xaxis", k, observedDistancesF.getValue(k));
//                rtF.setValue("xaxis", k, observedCDF.getValue(k));
//                rtF.setValue("xaxis", k, averageCDF.getValue(k));
//            }
        }

        if (show) {
            PlotWindow.noGridLines = true;
            plotF.draw();
//            rtDistance.show("observedF_Eu"+desc);
//            rtSamples.show("samplesF_Eu"+desc);
//            PlotWindow plotW = plotF.show();
        }
        if (save) {
            PlotWindow plotW = plotF.show();
            if (plotW != null) {
                try {
                    plotW.getResultsTable().saveAs(dir + "StatsPlot-F-Eu"+desc+".csv");
                    IJ.selectWindow(plotW.getTitle());
                    IJ.saveAs("PNG", dir+"F_function_Euclidean_Distance.png");
                    plotW.close();
                } catch (IOException ex) {
                    Logger.getLogger(SpatialTissueAnalysis.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
            try {
                rtDistance.saveAs(dir+"observedF_Eu"+desc+".csv");
                rtSamples.saveAs(dir+"samplesF_Eu"+desc+".csv");
            } catch (IOException ex) {
                IJ.log("Have the problem of saving data");
            }
        }
        
        IJ.log("--- F Function ---");
        IJ.log("");
        sdi_F = CDFTools.SDI(observedDistancesF, sampleDistancesF, averageCDF, xEvalsF);
        IJ.log("SDI F=" + sdi_F);
        IJ.log("");
        
    }
    
    /**
     * ttnhoa
     * @description Random elements in an existing list of elements (grid), instead of all space 
     * @param pop
     * @param mask
     * @param verbose
     * @param show
     * @param save 
     */
    private void processGFunctionInGrid(Objects3DPopulationGrid pop, Object3D mask, ArrayList<Point3D> lsCenterCells, boolean verbose, boolean show, boolean save, String dir) 
    {
        Calibration calibration = mask.getCalibration();
        int nbspots = pop.getNbObjects();
        IJ.log("*******  G_function_ Nb Spot=" + nbspots);
        ArrayUtil distances;
        
        // G
        ArrayUtil observedDistancesG;
        ArrayUtil observedCDG;
        observedDistancesG = pop.distancesAllClosestCenter();
        observedDistancesG.sort();
        observedCDG = CDFTools.cdf(observedDistancesG);
        ResultsTable rtDistance = new ResultsTable();
        for (int i = 0; i < observedDistancesG.getSize(); i++) 
        {
            rtDistance.incrementCounter();
            rtDistance.setValue("observed_X", i, observedDistancesG.getValue(i));
            rtDistance.setValue("observed_Y", i, observedCDG.getValue(i));
        }
        //G 
        ArrayUtil xEvalsG;
        ArrayUtil[] sampleDistancesG;
        ArrayUtil averageCDG;

        xEvalsG = new ArrayUtil(numRandomSamples * nbspots);
        sampleDistancesG = new ArrayUtil[numRandomSamples];
        
        
        for (int i = 0; i < numRandomSamples; i++) {
            if (verbose) {
                IJ.showStatus("Random population G " + i);
            }
            Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
            poprandom.setCalibration(calibration);
            poprandom.setMask(mask);
            //poprandom.createRandomPopulation(nbspots, distHardCore);
            poprandom.createRandomPopulationInGrid(lsCenterCells, nbspots, distHardCore);
            //poprandom.createKDTreeCenters();
            distances = poprandom.distancesAllClosestCenter();
            distances.sort();
            sampleDistancesG[i] = distances;
            xEvalsG.insertValues(i * nbspots, distances);
        }
        xEvalsG.sort();
        averageCDG = CDFTools.cdfAverage(sampleDistancesG, xEvalsG);

        // G function
        ResultsTable rtSamples = new ResultsTable();
        for (int i = 0; i < numRandomSamples; i++) {
            if (verbose) {
                IJ.showStatus("Random population G " + i);
            }
            Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
            poprandom.setCalibration(calibration);
            poprandom.setMask(mask);
            poprandom.createRandomPopulationInGrid(lsCenterCells, nbspots, distHardCore);
            //poprandom.createRandomPopulation(nbspots, distHardCore);
            poprandom.createKDTreeCenters();
            distances = poprandom.distancesAllClosestCenter();
            distances.sort();
            sampleDistancesG[i] = distances;//    
            if(i==0)
            {
                for (int k = 0; k < sampleDistancesG[0].getSize(); k++) 
                {
                    rtSamples.incrementCounter();
                    rtSamples.setValue("sG_"+0, k, sampleDistancesG[0].getValue(k));
                }
            }
            else
            {
                for (int k = 0; k < sampleDistancesG[i].getSize(); k++) 
                {
                    rtSamples.setValue("sG_"+i, k, sampleDistancesG[i].getValue(k));
                }
            }
        }

        // plot
        Plot plotG = null;
        if (show || save) {
//            plotG = createPlot(xEvalsG, sampleDistancesG, observedDistancesG, observedCDG, averageCDG, "G");
            plotG = createPlot(xEvalsG, sampleDistancesG, observedDistancesG, observedCDG, averageCDG, "Euclidean-Distance-G");
        }

        if (show) {
            PlotWindow.noGridLines = true;
            plotG.draw();
//            rtDistance.show("observedG_Eu"+desc);
//            rtSamples.show("samplesG_Eu"+desc);
//            PlotWindow plotW = plotG.show();
        }
        if (save) {
            PlotWindow plotW = plotG.show();
            if (plotW != null) {
                try {
                    plotW.getResultsTable().saveAs(dir + "StatsPlot-G-Eu"+desc+".csv");
                    IJ.selectWindow(plotW.getTitle());
                    IJ.saveAs("PNG", dir+"G_function_Euclidean_Distance.png");
                    plotW.close();
                } catch (IOException ex) {
                    Logger.getLogger(SpatialTissueAnalysis.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
            try {
                rtDistance.saveAs(dir+"observedG_Eu"+desc+".csv");
                rtSamples.saveAs(dir+"samplesG_Eu"+desc+".csv");
            } catch (IOException ex) {
                IJ.log("Have the problem of saving data");
            }
        }
        
        IJ.log("--- G Function ---");
        sdi_G = CDFTools.SDI(observedDistancesG, sampleDistancesG, averageCDG, xEvalsG);
        IJ.log("SDI G = " + sdi_G);
        IJ.log("------------------");
        
    }
    /**
     * ttnhoa
     * @description Random elements in an existing list of elements (grid), instead of all space 
     * @param pop
     * @param mask
     * @param lsCenterCells
     * @param verbose
     * @param show
     * @param save 
     */
    private void processHFunctionInGrid(Objects3DPopulationGrid pop, Object3D mask, ArrayList<Point3D> lsCenterCells, boolean verbose, boolean show, boolean save) {
        Calibration calibration = mask.getCalibration();
        int nbspots = pop.getNbObjects();
        ArrayUtil distances;

        // H
        ArrayUtil observedDistancesH;
        ArrayUtil observedCDH;

        observedDistancesH = pop.distancesAllCenter();
        observedDistancesH.sort();
        observedCDH = CDFTools.cdf(observedDistancesH);

        // H
        ArrayUtil xEvalsH;
        ArrayUtil[] sampleDistancesH;
        ArrayUtil averageCDH;

        xEvalsH = new ArrayUtil(numRandomSamples * nbspots);
        sampleDistancesH = new ArrayUtil[numRandomSamples];
        for (int i = 0; i < numRandomSamples; i++) {
            if (verbose) {
                IJ.showStatus("Random population H " + i);
            }
            Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
            poprandom.setCalibration(calibration);
            poprandom.setMask(mask);
            poprandom.createRandomPopulationInGrid(lsCenterCells, nbspots, distHardCore);
            //poprandom.createRandomPopulation(nbspots, distHardCore);
            distances = poprandom.distancesAllCenter();
            distances.sort();
            sampleDistancesH[i] = distances;
            xEvalsH.insertValues(i * nbspots, distances);
        }
        xEvalsH.sort();
        averageCDH = CDFTools.cdfAverage(sampleDistancesH, xEvalsH);

        // H function
        for (int i = 0; i < numRandomSamples; i++) {
            if (verbose) {
                IJ.showStatus("Random population H " + i);
            }
            Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
            poprandom.setCalibration(calibration);
            poprandom.setMask(mask);
            poprandom.createRandomPopulationInGrid(lsCenterCells, nbspots, distHardCore);
            //poprandom.createRandomPopulation(nbspots, distHardCore);
            distances = poprandom.distancesAllCenter();
            distances.sort();
            sampleDistancesH[i] = distances;
        }

        // plot
        Plot plotH = null;
        if (show || save) {
            plotH = createPlot(xEvalsH, sampleDistancesH, observedDistancesH, observedCDH, averageCDH, "H");
        }

        if (show) {
            plotH.draw();
        }
        if (save) {
            PlotWindow plotW = plotH.show();
            if (plotW != null) {
                try {
                    plotW.getResultsTable().saveAs(IJ.getDirectory("home") + "StatsPlot-H.csv");
                    IJ.selectWindow(plotW.getTitle());
                    IJ.saveAs("PNG", IJ.getDirectory("home") +"H_function_EuclideanDistance.png");
                    plotW.close();
                } catch (IOException ex) {
                    Logger.getLogger(SpatialTissueAnalysis.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

        }

        IJ.log("--- H Function ---");
        sdi_H = CDFTools.SDI(observedDistancesH, sampleDistancesH, averageCDH, xEvalsH);
        IJ.log("SDI H=" + sdi_H);
    }

    private Plot createPlot(ArrayUtil xEvals, ArrayUtil[] sampleDistances, ArrayUtil observedDistances, ArrayUtil observedCD, ArrayUtil averageCD, String function) {
        double plotmaxX = observedDistances.getMaximum();
        double plotmaxY = observedCD.getMaximum();

        // low env      
        double max = xEvals.getMaximum();
        ArrayUtil xEvals0 = new ArrayUtil(nbBins);
        for (int i = 0; i < nbBins; i++) {
            xEvals0.addValue(i, ((double) i) * max / ((double) nbBins));
        }
        // get the values
        ArrayUtil samplespc5 = CDFTools.cdfPercentage(sampleDistances, xEvals0, env / 2.0);
        ArrayUtil samplespc95 = CDFTools.cdfPercentage(sampleDistances, xEvals0, 1.0 - env / 2.0);
        // get the limits        
        if (xEvals0.getMaximum() > plotmaxX) {
            plotmaxX = xEvals0.getMaximum();
        }
        if (samplespc5.getMaximum() > plotmaxY) {
            plotmaxY = samplespc5.getMaximum();
        }
        if (samplespc95.getMaximum() > plotmaxY) {
            plotmaxY = samplespc95.getMaximum();
        }
        if (xEvals.getMaximum() > plotmaxX) {
            plotmaxX = xEvals.getMaximum();
        }
        if (averageCD.getMaximum() > plotmaxY) {
            plotmaxY = averageCD.getMaximum();
        }
        if (observedCD.getMaximum() > plotmaxY) {
            plotmaxY = observedCD.getMaximum();
        }
        if (observedDistances.getMaximum() > plotmaxX) {
            plotmaxX = observedDistances.getMaximum();
        }
        // create the plot
        Plot plot = new Plot(function + "-function", "distance", "cumulated frequency");
        plot.setLimits(0, plotmaxX, 0, plotmaxY);

        // enveloppe  for e.g 10 % at 5 and 95 %
        plot.setColor(Color.green);
        plot.addPoints(xEvals0.getArray(), samplespc5.getArray(), Plot.LINE);

        // high envxEvals.getMaximum
        plot.setColor(Color.black);
        plot.addPoints(xEvals0.getArray(), samplespc95.getArray(), Plot.LINE);

        // average
        plot.setColor(Color.red);
        plot.addPoints(xEvals.getArray(), averageCD.getArray(), Plot.LINE);

        // observed
        plot.setColor(Color.blue);
        plot.addPoints(observedDistances.getArray(), observedCD.getArray(), Plot.LINE);

        return plot;
    }
    /**
     * ttnhoa
     * @description: process elements in an existing list of elements (grid), instead of all space
     * @param plusSpots
     * @param plusMask
     * @param lsCenterCells
     * @param functions
     * @param verbose
     * @param show
     * @param save 
     */
    public void processInGrid(Objects3DPopulationGrid pop, ImageHandler plusMask, ArrayList<Point3D> lsCenterCells, String functions, boolean verbose, boolean show, boolean save, String dir) 
    {
        Calibration calibration = plusMask.getCalibration();
        if (calibration == null) 
        {
            IJ.log("Image not calibrated");
            calibration = new Calibration();
            calibration.setUnit("pix");
            calibration.pixelWidth = 1;
            calibration.pixelHeight = 1;
            calibration.pixelDepth = 1;
        }

//        ImageInt inImage = (ImageInt) plusSpots;
//        ImageInt segImage;
//        if (inImage.isBinary(0)) {
//            if (verbose) {
//                IJ.log("Segmenting image...");
//            }
//            inImage = inImage.threshold(0, false, true);
//            ImageLabeller labels = new ImageLabeller(false);
//            segImage = labels.getLabels(inImage);
//            if (verbose) {
//                segImage.show("Labelled Image");
//            }
//        } else {
//            segImage = (ImageInt) inImage.duplicate();
//        }
//        segImage.setCalibration(calibration);
//
//        int nbspots;
//
//        Objects3DPopulation pop = new Objects3DPopulation();
        ImageInt maskHandler = (ImageInt) plusMask;
        Object3D mask = new Object3DLabel(maskHandler, (int) maskHandler.getMax());
        mask.setCalibration(calibration);
//        pop.setMask(mask);
//        pop.addImage(segImage, calibration);
//        pop.setCalibration(calibration);

        // random sample
        Objects3DPopulationGrid poprandom = new Objects3DPopulationGrid();
        poprandom.setCalibration(calibration);
        poprandom.setMask(mask);
        poprandom.createRandomPopulationInGrid(lsCenterCells, pop.getNbObjects(), distHardCore);
        //poprandom.createRandomPopulation(pop.getNbObjects(), distHardCore);
        randomPop = maskHandler.createSameDimensions();
        randomPop.setCalibration(calibration);
        poprandom.draw(randomPop);
        //randomPop.show("random");

        if ((plusMask.getCalibration() == null) || (!plusMask.getCalibration().scaled())) {
//            if (verbose) {
//                IJ.log("mask not calibrated, calibrating ...");
//            }
            plusMask.setCalibration(calibration);
            plusMask.getImagePlus().updateAndRepaintWindow();
        }

        int nbspots = pop.getNbObjects();
        if (verbose) {
            IJ.log("Computing spatial statistics, please wait ...");
        }
        IJ.log("*******Nb Spot=" + nbspots);
//        if (verbose) {
//            IJ.log("Nb Spot=" + nbspots);
//            IJ.log("Volume mask=" + mask.getVolumeUnit());
//            IJ.log("Density=" + (nbspots / mask.getVolumeUnit()));
//        }

        if (functions.contains("F")) {
            //processF(pop, mask, verbose, show, save);
            processFFunctionInGrid(pop, mask, lsCenterCells, verbose, show, save, dir);
        }

        if (functions.contains("G")) {
            processGFunctionInGrid(pop, mask, lsCenterCells, verbose, show, save, dir);
        }

        if (functions.contains("H")) {
            //processH(pop, mask, verbose, show, save);
            processHFunctionInGrid(pop, mask, lsCenterCells, verbose, show, save);
            
        }
    }

    private Point3D[] createEvaluationPointsInList(int numPoints, Objects3DPopulationGrid population, ArrayList<Point3D> lsCenterCells) {
        Point3D[] evaluationPoints = new Point3D[numPoints];
        for (int i = 0; i < numPoints; ++i) {
            evaluationPoints[i] = population.getRandomPointInList(lsCenterCells);
        }

        return evaluationPoints;
    }
    private Point3D[] createEvaluationPointsInList(ArrayList<Point3D> lsCenterCells) {
        int numPoints = lsCenterCells.size();
        Point3D[] evaluationPoints = new Point3D[numPoints];
        for (int i = 0; i < numPoints; ++i) {
            evaluationPoints[i] = lsCenterCells.get(i);
        }
        return evaluationPoints;
    }

}
