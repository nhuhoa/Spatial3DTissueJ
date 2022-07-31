/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;

import ij.IJ;
import ij.measure.ResultsTable;

/**
 *
 * @author tranhoa
 */
public class AnalyseResult_ implements ij.plugin.PlugIn
{
    public void run(String arg) 
    {
        IJ.log("Analysing result......");
        fusionTable();
        IJ.log("Finished");
        
    }
    public void fusionTable()
    {
        ResultsTable total = new ResultsTable();
        ResultsTable r1 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/final_test/SEG/method1/method1.csv");
        ResultsTable r2 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/final_test/SEG/method2/method2.csv");
        ResultsTable r3 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/final_test/SEG/method3/method3.csv");
        ResultsTable r4 = ResultsTable.open2("/home/tranhoa/Raphael/data_test/testV2/final_test/SEG/method4/method4.csv");
        int cA=0, cB=0, cD=0, cAB=0, cAD=0, cBD=0;
        for(int i=0; i<r1.size(); i++)
        {
            double val1 = r1.getValue("cell", i);
            double val2 = r2.getValue("cell", i);
            double val3 = r3.getValue("cell", i);
            double val4 = r3.getValue("cell", i);
            IJ.log("val1: "+val1 +" val2: "+val2+" val3:"+val3 + " val4: "+val4);
            if(val1==val2 && val1==val3 && val1==val4)
            {
                double type1 = r1.getValue("method1", i);
                double type2 = r2.getValue("method2", i);
                double type3 = r3.getValue("method3", i);
                double type4 = r4.getValue("method4", i);
                IJ.log("t1: "+type1 +" t2: "+type2+" t3:"+type3 + " t4: "+type4);
                double typeC = computeCellType(type1, type2, type3, type4);
                if(typeC==1){cD++;}
                if(typeC==2){cB++;}
                if(typeC==3){cA++;}
                if(typeC==12){cBD++;}
                if(typeC==13){cAD++;}
                if(typeC==23){cAB++;}
                
                total.incrementCounter();
                total.setValue("cell", i, val1);
                total.setValue("method1", i, type1);
                total.setValue("method2", i, type2);
                total.setValue("method3", i, type3);
                total.setValue("method4", i, type4);
                total.setValue("type", i, typeC);
            }
        }
        IJ.log("nb of cell type: A: "+cA+" B: "+cB+" D: "+cD+" AB: "+cAB+" AD: "+cAD+" BD: "+cBD);
        total.show("Result");
    }
    public double computeCellType(double type1, double type2, double type3, double type4)
    {
        if(     (type3==type1 && type3==type2 && type3==type4) ||  
                (type3==type1 && type3==type2 && type3!=type4) ||
                (type3==type4 && type3==type2 && type3!=type1) ||
                (type3==type4 && type3==type1 && type3!=type2) ||
                (type3==type4) || (type3==type1) || (type3==type2) )        
        {
            return type3;
        }
        else if(type1==type4 && type1==type2){
            return type1;
        }
        else{
            return 0;
        }
    
    }        
    
}
