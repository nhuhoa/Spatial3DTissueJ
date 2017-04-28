/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stdlib;

import mcib3d.geom.Point3D;

/**
 *
 * @author tranhoa
 */
public class CentroidC implements Comparable<CentroidC>
{
    private int id;
    public double index;
    private byte typeCell;
    private float shortestDis = -1;
    private Point3D centroid = null;
    public CentroidC(int id, float dis) 
    {
        this.id = id;
        this.shortestDis = dis;
    }
    public CentroidC(int id, byte typeC, Point3D pos) 
    {
        this.id = id;
        this.typeCell = typeC;
        this.centroid = pos; 
    }
    public int getId() {
        return id;
    }
    public byte getType(){
        return typeCell;
    }
    public void setShortestDistance(float distanceEu)
    {
        this.shortestDis = distanceEu;
    }
    public float getShortestDistance()
    {
        return this.shortestDis;
    }
    public void setCentroid(Point3D cent)
    {
        this.centroid = cent;
    }
    public Point3D getCentroid()
    {
        return this.centroid;
    }
    public int compareTo(CentroidC v) {
        if (shortestDis > v.getShortestDistance()) {
            return 1;
        } else if (shortestDis < v.getShortestDistance()) {
            return -1;
        } else {
            return 0;
        }
    }
}
