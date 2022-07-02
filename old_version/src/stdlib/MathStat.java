/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stdlib;

import java.util.Arrays;

/**
 *
 * @author tranhoa
 */
public class MathStat {

//    public static double[] data;
//    public static int size;   

//    public MathStat(double[] data) 
//    {
//        this.data = data;
//        size = data.length;
//    }   

    public static double getSum(double[] data)
    {
        int size = data.length;
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum;
    }
    public static double getMean(double[] data)
    {
        int size = data.length;
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return sum/size;
    }

    public static double getVariance(double[] data)
    {
        int size = data.length;
        double mean = getMean(data);
        double temp = 0;
        for(double a :data)
            temp += (a-mean)*(a-mean);
        return temp/size;
    }

    public static double getStdDev(double[] data)
    {
        return Math.sqrt(getVariance(data));
    }

    public static double median(double[] data) 
    {
       Arrays.sort(data);

       if (data.length % 2 == 0) 
       {
          return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
       } 
       else 
       {
          return data[data.length / 2];
       }
    }
}

