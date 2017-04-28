/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stdlib;

import java.util.ArrayList;

/**
 *
 * @author tranhoa
 */

public class Vertex implements Comparable<Vertex>{
    public final int id;
    public double index;
    public final String name;
    public byte typeCell;
    public double disShortestLayer;
    public double disShortestPath;
    public int idDCShortestLayer = -1;
    public ArrayList<Integer> idxShortestLayer = null;
    public int idDCShortestPath = -1;
  
    public Vertex(int id, String name, byte typeC) 
    {
        this.id = id;
        this.name = name;
        this.typeCell = typeC;
        this.disShortestLayer = -1;
        this.disShortestPath = -1;
    }
    public int getId() {
        return id;
    }
    public String getName() {
        return name;
    }
    public void setType(byte newType){
        this.typeCell = newType;
    }
    public byte getType(){
        return typeCell;
    }
    public void setDistanceLayer(double distanceLayer)
    {
        this.disShortestLayer = distanceLayer;
    }
    public double getDistanceLayer()
    {
        return this.disShortestLayer;
    }
    public void setShortestDistance(double distanceEu)
    {
        this.disShortestPath = distanceEu;
    }
    public double getShortestDistance()
    {
        return this.disShortestPath;
    }
    
  
  @Override
    public boolean equals(Object obj) {
        if (this == obj)
          return true;
        if (obj == null)
          return false;
        if (getClass() != obj.getClass())
          return false;
        Vertex other = (Vertex) obj;
        if (id < 0) {
          if (other.id < 0)
            return false;
        } else if (id!=other.id)
          return false;
        return true;
    }

    public int compareTo(Vertex v) {
        if (this.disShortestPath > v.getShortestDistance()) {
            return 1;
        } else if (this.disShortestPath < v.getShortestDistance()) {
            return -1;
        } else {
            return 0;
        }
    }
    @Override
    public String toString() {
          return name;
    }
  
} 