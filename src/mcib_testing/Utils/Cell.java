/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Utils;

import java.util.ArrayList;
import mcib3d.geom.Object3D;

/**
 *
 * @author tranhoa
 */
public class Cell 
{
        public int id;
        public Object3D nucleus;
        public Object3D region;
        public byte type;
        public int layer, degreeLabelled, degreeUnlabelled;
        public int layer_distance = 0;
        public boolean marked = false;        //for clustering only
        public boolean connected = false;
        public boolean seen = false;          //for erode and dilation
        public int currLevel = -1;
        public boolean lastnode = false;
        public ArrayList<Cell> nei1 = null;
        public ArrayList<Cell> nei2 = null;

             
        @Override
        public String toString() {
            return "(" + nucleus.getValue() + ", " + region.getValue() + ", " + type + ")";
        }

        public int[] computeNeiType(ArrayList<Cell> nei) {
            int[] res = new int[3];
            for (Cell C : nei) {
                if (C.type > 0) {
                    res[C.type - 1]++;
                }
            }
            return res;
        }
        public int[] computeNei1TypeHisto() {
            return computeNeiType(nei1);
        }
        public boolean hasContact(int type, ArrayList<Cell> nei) {
            for (Cell N : nei) {
                if (N.type == type) {
                    return true;
                }
            }
            return false;
        }
        
        public boolean hasContact1(int type) {
            return hasContact(type, nei1);
        }

        public boolean hasContact2(int type) {
            return hasContact(type, nei2);
        }
//        public void setDistance(int layerDistance)
//        {
//            layer_distance = layerDistance;
//        }            
//        public int getDistance()
//        {
//            return layer_distance;
//        }        

}
