/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stdlib;

import ij.gui.Plot;
import java.awt.Color;
import mcib3d.utils.ArrayUtil;
import mcib3d.utils.CDFTools;

/**
 *
 * @author tranhoa
 */
public class PlotDraw 
{
    public static Plot createPlot(String function, String xStr, String yStr, ArrayUtil xAxis, ArrayUtil nbCells1, ArrayUtil nbCells2, ArrayUtil nbCells3)
    {
        double plotmaxX = xAxis.getMaximum() + 1;
        double plotmaxY = 0.0;
        double max1 = nbCells1.getMaximum();
        double max2 = nbCells2.getMaximum();
        if(max1>max2)
        {
            plotmaxY = max1 + 2;
        }
        else{
            plotmaxY = max2 + 2;
        }
        Plot plot = new Plot(function, xStr, yStr);
        plot.setLimits(0, plotmaxX, 0, plotmaxY);
        plot.setColor(Color.red);
        plot.addPoints(xAxis.getArray(), nbCells1.getArray(), Plot.LINE);

        plot.setColor(Color.blue);
        plot.addPoints(xAxis.getArray(), nbCells2.getArray(), Plot.LINE);

        if(nbCells3!=null)
        {
            plot.setColor(Color.green);
            plot.addPoints(xAxis.getArray(), nbCells3.getArray(), Plot.LINE);
        }
        return plot;

    }
    public static Plot createPlot(ArrayUtil xEvals, ArrayUtil[] sampleDistances, ArrayUtil observedDistances, ArrayUtil observedCD, ArrayUtil averageCD, String function, int nbBins, double env) {
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
        plot.setColor(Color.green);
        plot.addPoints(xEvals0.getArray(), samplespc95.getArray(), Plot.LINE);

        // average
        plot.setColor(Color.red);
        plot.addPoints(xEvals.getArray(), averageCD.getArray(), Plot.LINE);

        // observed
        plot.setColor(Color.blue);
        plot.addPoints(observedDistances.getArray(), observedCD.getArray(), Plot.LINE);

        return plot;
    }
    public static Plot createPlot(String function, String xStr, String yStr, ArrayUtil xAxis, ArrayUtil nbCells1)
    {
        double plotmaxX = xAxis.getMaximum();
        double plotmaxY = 0.0;
        double max1 = nbCells1.getMaximum();
        plotmaxY = max1 + 2;
        Plot plot = new Plot(function, xStr, yStr);
        plot.setLimits(0, plotmaxX, 0, plotmaxY);
        plot.setColor(Color.red);
        plot.setAxes(true, true, true, true, false, false, (int)plotmaxX/5, 0);
        plot.addPoints(xAxis.getArray(), nbCells1.getArray(), Plot.LINE);
        return plot;

    } 
    public static Plot createPlotLayer(String function, String xStr, String yStr, ArrayUtil xAxis, ArrayUtil[] nbCellsArr)
    {
        double plotmaxX = xAxis.getMaximum() + 1;
        double plotmaxY = 0.0;
        double sizeArr = nbCellsArr.length;
        for(int i=0; i < sizeArr; i++)
        {
            if(plotmaxY < nbCellsArr[i].getMaximum())
            {
                plotmaxY = nbCellsArr[i].getMaximum();
            }
        }
        plotmaxY = plotmaxY + 2; 
        Plot plot = new Plot(function, xStr, yStr);
        plot.setLimits(0, plotmaxX, 0, plotmaxY);
        plot.setColor(Color.red);
        plot.addPoints(xAxis.getArray(), nbCellsArr[0].getArray(), Plot.LINE);

        plot.setColor(Color.blue);
        plot.addPoints(xAxis.getArray(), nbCellsArr[1].getArray(), Plot.LINE);

        return plot;

    } 
}
