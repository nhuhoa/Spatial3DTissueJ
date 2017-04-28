/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.RafaelV2;

/**
 *
 * @author tranhoa
 */

public class regionCell 
    {
        double regM; 
        double regVolume;
        int index1;
        String label;
        public regionCell(double regM, double regVolume, int index)
        {        
            this.regM = regM;
            this.regVolume = regVolume;
            this.index1 = index;
        }
        public double getMValue(){ return regM;}
        public void setMValue(double val)
        {
            this.regM = val;
        }
        public double getVolume(){ return regVolume;}
        public void setValVolume(double val)
        {
            this.regVolume = val;
        }
        public int getIndex1(){ return index1;}
        public void setIndex1(int index)
        {
            this.index1 = index;
        }
        public void setLabel(String label){
            this.label = label;
        }
        public String getLabel(){ return label;}
        
    }