/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Utils;

import stdlib.Vertex;
import stdlib.DirectedEdgeC;
import ij.IJ;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.ResultsTable;
import java.io.IOException;
import java.util.ArrayList;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;
import stdlib.PlotDraw;

/**
 *
 * @author tranhoa
 */
public class CellSpatialAnalysis 
{
    
    private final int numEvaluationPoints;
    private final int numRandomSamples;
    int nbBins = 1000;
    private double sdi_F = Double.NaN;
    private double sdi_G = Double.NaN;
    private double sdi_H = Double.NaN;
    private final double env;
    private String directory = null;
    private boolean show = false;
    public CellSpatialAnalysis(int numPoints, int numRandomSamples, double env, String directory) {
        this.numEvaluationPoints = numPoints;
        this.numRandomSamples = numRandomSamples;
        this.env = env;
        this.directory = directory;
    }
    public CellSpatialAnalysis() {
        this.numEvaluationPoints = 0;
        this.numRandomSamples = 0;
        this.env = 0.05;
    }
    
    public void process(ArrayList<Vertex> sourceLst, ArrayList<Vertex> desLst, ArrayList<Vertex> lsTotalVerx, ArrayList<DirectedEdgeC> lsTotalEdges, String str)
    {
        ResultsTable rtContact = new ResultsTable();
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(lsTotalVerx);
        pop.setMask(lsTotalVerx, lsTotalEdges);
        ArrayUtil arr = pop.computeMinDistance(sourceLst, desLst);
        rtContact.incrementCounter();
        int distanceLength = (int)arr.getMaximum()+1;
        int[] res = new int[distanceLength];
        for (int i = 0; i < arr.getSize(); i++) 
        {
            res[arr.getValueInt(i)]++;
        }
        for(int j=1; j< res.length; j++)
        {
            if(res[j]!=0)
            {
                rtContact.setValue("d_"+j, 0, res[j]);
            }
            else{
                rtContact.setValue("d_"+j, 0, 0);
            }
        }
//        rtContact.show("Histogram_"+str);
        try {
            rtContact.saveAs(directory+"Histogram_Contact_"+str+".csv");
            IJ.log("Saving the output table "+"Histogram_Contact_"+str+".csv"+" in the folder: "+directory);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
        
    }
    
    public void process(ArrayList<Vertex> lsObservedVerx, Cells3DPopulation popTotal, String functions)
    {
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(lsObservedVerx);
        pop.setMask(popTotal.getVertexes(), popTotal.getEdges());
        if (functions.contains("F")) {
            process_F_Function(pop, popTotal);
        }
        if (functions.contains("G")) {
            process_G_Function(pop, popTotal);
        }
//        ArrayUtil arr = new ArrayUtil(2);
//        arr.addValue(0, sdi_F);
//        arr.addValue(1, sdi_G);
//        return arr;
    }
    public void process(ArrayList<Vertex> lsObservedVerx, ArrayList<Vertex> lsTotalVerx, ArrayList<DirectedEdgeC> lsTotalEdges, String functions)
    {
        Cells3DPopulation pop = new Cells3DPopulation();
        pop.setVertexes(lsObservedVerx);
        pop.setMask(lsTotalVerx, lsTotalEdges);
        //int nbspots = pop.getNbObjects();
        // random sample
//        Cells3DPopulation popRandom = new Cells3DPopulation();
//        popRandom.setMask(lsTotalVerx, lsTotalEdges);
//        popRandom.getRandomCellInPopulation(nbspots, lsTotalVerx);
        
        Cells3DPopulation popTotal = new Cells3DPopulation(lsTotalVerx, lsTotalEdges);
        if (functions.contains("F")) {
            //processF(pop, mask, verbose, show, save);
            process_F_Function(pop, popTotal);
        }
        if (functions.contains("G")) {
            //processF(pop, mask, verbose, show, save);
            process_G_Function(pop, popTotal);
        }
    }
    private void process_G_Function(Cells3DPopulation pop, Cells3DPopulation popTotal)
    {
        int nbspots = pop.getNbObjects();
//        ResultsTable rtContact = new ResultsTable();
//        ResultsTable rtDistance = new ResultsTable();
        ArrayUtil distances;
        ArrayList<Vertex> lsTotalVerx = popTotal.getVertexes();
        ArrayList<DirectedEdgeC> lsTotalEdges = popTotal.getEdges();
        // G
        ArrayUtil observedDistancesG;
        ArrayUtil observedCDG;
        observedDistancesG = pop.distancesAllClosestCenter();
        observedDistancesG.sort();
        observedCDG = CDFTools.cdf(observedDistancesG);
        
        //G 
        ArrayUtil xEvalsG;
        ArrayUtil[] sampleDistancesG;
        ArrayUtil averageCDG;
        xEvalsG = new ArrayUtil(numRandomSamples * nbspots);
        ResultsTable rtDistance = new ResultsTable();
        for (int i = 0; i < observedDistancesG.getSize(); i++) 
        {
            rtDistance.incrementCounter();
            rtDistance.setValue("observed_Cell_X", i, observedDistancesG.getValue(i));
            rtDistance.setValue("observed_Cell_Y", i, observedCDG.getValue(i));
        }
//        rtContact.incrementCounter();
//        rtDistance.incrementCounter();
//        int distanceLength = 20;
//        int[] res = new int[distanceLength];
//        for (int i = 0; i < observedDistancesG.getSize(); i++) 
//        {
//            rtDistance.setValue((i+1)+"_", 0, observedDistancesG.getValueInt(i));
//            res[observedDistancesG.getValueInt(i)]++;
//        }
//        for(int j=0; j< res.length; j++)
//        {
//            if(res[j]!=0)
//            {
//                rtContact.setValue("d_"+j, 0, res[j]);
//            }
//            else{
//                rtContact.setValue("d_"+j, 0, 0);
//            }
//        }
        sampleDistancesG = new ArrayUtil[numRandomSamples];
        // G function
        for (int i = 0; i < numRandomSamples; i++) 
        {
            IJ.showStatus("Random population G " + i);
//            IJ.log("Random population G " + i);
            Cells3DPopulation popRandom = new Cells3DPopulation();
            popRandom.setMask(lsTotalVerx, lsTotalEdges);
            popRandom.getRandomCellInPopulation(nbspots, lsTotalVerx);
            distances = popRandom.distancesAllClosestCenter();
            distances.sort();
            sampleDistancesG[i] = distances;
            xEvalsG.insertValues(i * nbspots, distances);
        }
        xEvalsG.sort();
        averageCDG = CDFTools.cdfAverage(sampleDistancesG, xEvalsG);
        ResultsTable rtSamples = new ResultsTable();
        for (int i = 0; i < numRandomSamples; i++) 
        {
            IJ.showStatus("Random population G " + i);
//            IJ.log("Random population G " + i);
            Cells3DPopulation popRandom = new Cells3DPopulation();
            popRandom.setMask(lsTotalVerx, lsTotalEdges);
            popRandom.getRandomCellInPopulation(nbspots, lsTotalVerx);
            distances = popRandom.distancesAllClosestCenter();
            distances.sort();
            sampleDistancesG[i] = distances;
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
//            rtContact.incrementCounter();
//            int[] resRandom = new int[distanceLength];
//            for (int k = 0; k < sampleDistancesG[i].getSize(); k++) 
//            {
//                resRandom[sampleDistancesG[i].getValueInt(k)]++;
//            }
//            for(int j1=0; j1< distanceLength; j1++)
//            {
//                if(resRandom[j1]!=0)
//                {
//                    rtContact.setValue("d_"+j1, (i+1), resRandom[j1]);
//                }
//                else{
//                    rtContact.setValue("d_"+j1, (i+1), 0);
//                }
//            }    
        }
        Plot plotG = null;
        plotG = PlotDraw.createPlot(xEvalsG, sampleDistancesG, observedDistancesG, observedCDG, averageCDG, "Cell-Distance-G",nbBins,env);
//        plotG.draw();
//        plotG.show();
        PlotWindow.noGridLines = true;
        plotG.draw();
        PlotWindow plotW = plotG.show();
        IJ.selectWindow(plotW.getTitle());
        IJ.saveAs("PNG", directory+"G_function_Cell_Distance.png");
        plotW.close();
        if(show)
        {
            rtDistance.show("observedG_CellDistance");
            rtSamples.show("samplesG_CellDistance");  
        } 
        
        try {
            rtDistance.saveAs(directory+"observedG_CellDistance.csv");
            rtSamples.saveAs(directory+"samplesG_CellDistance.csv");
            IJ.log("Saving the output table in the folder: "+directory);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
        IJ.log("--- G Function Cell Distance---");
        sdi_G = CDFTools.SDI(observedDistancesG, sampleDistancesG, averageCDG, xEvalsG);
        IJ.log("SDI G=" + sdi_G);
    }
    private void process_F_Function(Cells3DPopulation pop, Cells3DPopulation popTotal)
    {
        int nbspots = pop.getNbObjects();
//        ResultsTable rtContact = new ResultsTable();
//        rtContact.incrementCounter();
        ArrayList<Vertex> lsTotalVerx = popTotal.getVertexes();
        ArrayList<DirectedEdgeC> lsTotalEdges = popTotal.getEdges();
        ArrayUtil distances;
        
        Vertex[] evaluationPoints;
        ArrayUtil observedDistancesF;
        ArrayUtil[] sampleDistancesF;
        ArrayUtil observedCDF;
        ArrayUtil xEvalsF;
        ArrayUtil averageCDF;
        evaluationPoints = createEvaluationPoints(numEvaluationPoints, popTotal);
        observedDistancesF = pop.computeDistances(evaluationPoints);
        observedDistancesF.sort();
        observedCDF = CDFTools.cdf(observedDistancesF);
        
        int distanceLength = 20;
//        int[] res = new int[distanceLength];
        
        ResultsTable rtDistance = new ResultsTable();
        for (int i = 0; i < observedDistancesF.getSize(); i++) 
        {
            rtDistance.incrementCounter();
            rtDistance.setValue("observed_Cell_X", i, observedDistancesF.getValue(i));
            rtDistance.setValue("observed_Cell_Y", i, observedCDF.getValue(i));
        }
//        for(int j=0; j< res.length; j++)
//        {
//            if(res[j]!=0)
//            {
//                rtContact.setValue("d_"+j, 0, res[j]);
//            }
//            else{
//                rtContact.setValue("d_"+j, 0, 0);
//            }
//            
//        }    
//        
        xEvalsF = new ArrayUtil(numRandomSamples * numEvaluationPoints);
        sampleDistancesF = new ArrayUtil[numRandomSamples];
        for (int i = 0; i < numRandomSamples; i++) 
        {
            IJ.showStatus("Random population F " + i);
//            IJ.log("Random population F " + i);
            Cells3DPopulation popRandom = new Cells3DPopulation();
            popRandom.setMask(lsTotalVerx, lsTotalEdges);
            popRandom.getRandomCellInPopulation(nbspots, lsTotalVerx);
            distances = popRandom.computeDistances(evaluationPoints);
            distances.sort();
            sampleDistancesF[i] = distances;
            xEvalsF.insertValues(i * numEvaluationPoints, distances);
        }
        xEvalsF.sort();
        averageCDF = CDFTools.cdfAverage(sampleDistancesF, xEvalsF);
        
        ResultsTable rtSamples = new ResultsTable();
        for (int i = 0; i < numRandomSamples; i++) 
        {
            IJ.showStatus("Random population F " + i);
//            IJ.log("Random population F " + i);
            Cells3DPopulation popRandom = new Cells3DPopulation();
            popRandom.setMask(lsTotalVerx, lsTotalEdges);
            popRandom.getRandomCellInPopulation(nbspots, lsTotalVerx);
            distances = popRandom.computeDistances(evaluationPoints);
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
//            rtContact.incrementCounter();
//            int[] resRandom = new int[distanceLength];
//            for (int k = 0; k < sampleDistancesF[i].getSize(); k++) 
//            {
//                resRandom[sampleDistancesF[i].getValueInt(k)]++;
//            }
//            for(int j1=0; j1< distanceLength; j1++)
//            {
//                if(resRandom[j1]!=0)
//                {
//                    rtContact.setValue("d_"+j1, (i+1), resRandom[j1]);
//                }
//                else{
//                    rtContact.setValue("d_"+j1, (i+1), 0);
//                }
//            }    
        }


        Plot plotF = null;
        plotF = PlotDraw.createPlot(xEvalsF, sampleDistancesF, observedDistancesF, observedCDF, averageCDF, "Cell-Distance-F",nbBins,env);
        PlotWindow.noGridLines = true;
        plotF.draw();
        PlotWindow plotW = plotF.show();
        IJ.selectWindow(plotW.getTitle());
        IJ.saveAs("PNG", directory+"F_function_Cell_Distance.png");
        plotW.close();
        if(show)
        {
            rtDistance.show("observedF_CellDistance");
            rtSamples.show("samplesF_CellDistance");
        }    
        
        try {
            rtDistance.saveAs(directory+"observedF_CellDistance.csv");
            rtSamples.saveAs(directory+"samplesF_CellDistance.csv");
            IJ.log("Saving the output table in the folder: "+directory);
        } catch (IOException ex) {
            IJ.log("Have the problem of saving data");
        }
        IJ.log("--- F Function_ Cell Distance ---");
        IJ.log("");
        sdi_F = CDFTools.SDI(observedDistancesF, sampleDistancesF, averageCDF, xEvalsF);
        IJ.log("SDI F=" + sdi_F);
    }
    
    private Vertex[] createEvaluationPoints(int numPoints, Cells3DPopulation popTotal) {
        Vertex[] evaluationPoints = new Vertex[numPoints];
        //Cells3DPopulation evaluatePop = new Cells3DPopulation();
        if(popTotal.getNbObjects()==0)
        {
            return null;
        }    
        for (int i = 0; i < numPoints; ++i) 
        {
            evaluationPoints[i] = popTotal.getRandomCell();
        }
        return evaluationPoints;
    }
//    private Vertex[] createEvaluationPointsV2(int numPoints, Cells3DPopulation popTotal) 
//    {
//        Vertex[] evaluationPoints = new Vertex[numPoints];
//        //Cells3DPopulation evaluatePop = new Cells3DPopulation();
//        if(popTotal.getNbObjects()==0)
//        {
//            return null;
//        }    
//        for (int i = 0; i < numPoints; ++i) 
//        {
//            evaluationPoints[i] = popTotal.getObject(i);
//        }
//        return evaluationPoints;
//    }
//    public void test(ArrayList<Vertex> vertexes, ArrayList<DirectedEdgeC> edges)
//    {
//        pop = new Cells3DPopulation(vertexes, edges);
//        
//        IJ.log("Nb obj: " + pop.getNbObjects());
//        
//        Vertex s = pop.getRandomCell();
//        Vertex d = pop.getRandomCell();
//        IJ.log("source : "+(s.getId()+1)+" des: "+(d.getId()+1));
//        double distance = pop.computeDistance(s, d);
//        IJ.log(" distance is: "+distance);
//        
//    } 
//    public void test2(ArrayList<Vertex> sourceLst, ArrayList<Vertex> desLst)
//    {
//        ArrayUtil arr = pop.computeMinDistance(sourceLst, desLst);
//        IJ.log("layer distance is: "+arr.toString());
//        double closestE = pop.getClosest(sourceLst.get(1));
//        IJ.log("closest is: "+closestE);
//    }
//    public void test3()
//    {
//        Cells3DPopulation popRandom = new Cells3DPopulation();
//        popRandom.getRandomCellInPopulation(4, pop.getVertexes());
//        IJ.log("pop random: "+popRandom.toString());
//    }        
    
}
