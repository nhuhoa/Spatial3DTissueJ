/*
 * objects3d population in graph
 * cell distance instead of euclidean distance
 * 
 */
package mcib_testing.Utils;

import stdlib.Vertex;
import stdlib.EdgeWeightedDigraphC;
import stdlib.DirectedEdgeC;
import stdlib.DijkstraSPC;
import ij.IJ;
import ij.measure.Calibration;
import java.util.ArrayList;
import static javafx.application.Platform.exit;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageInt;
import mcib3d.utils.ArrayUtil;

/**
 *
 * @author tranhoa
 */
public class Cells3DPopulation 
{
    ArrayList<Vertex> vertexes;
    ArrayList<DirectedEdgeC> edges;
    EdgeWeightedDigraphC G, Gmask;
    
    public Cells3DPopulation()
    {
        vertexes = new ArrayList<Vertex>();
        edges = new ArrayList<DirectedEdgeC>();
    }
    public Cells3DPopulation(ArrayList<Vertex> vs, ArrayList<DirectedEdgeC> ed)
    {
        this.vertexes = vs;
        this.edges = ed;
        G = new EdgeWeightedDigraphC(vertexes, edges);
    }
    public void initGraph()
    {
        G = new EdgeWeightedDigraphC(vertexes, edges);
    }
    public void setMask(ArrayList<Vertex> vTotal, ArrayList<DirectedEdgeC> edTotal)
    {
        Gmask = new EdgeWeightedDigraphC(vTotal, edTotal);
    }        
    public ArrayList<Vertex> getVertexes()
    {
        return vertexes;
    }
    public ArrayList<DirectedEdgeC> getEdges()
    {
        return edges;
    }
    public void setVertexes(ArrayList<Vertex> lsVertexes)
    {
        this.vertexes = lsVertexes;
    }   
    
    public void addVertexes(Vertex v)
    {
        vertexes.add(v);
    }  
    public int getNbObjects() {
        return vertexes.size();
    }
    public Vertex getObject(int i) {
        return vertexes.get(i);
    }
    /**
     * min distance between pair of vertex
     * @param source
     * @param des 
     */
    public double computeDistance(Vertex source, Vertex des)
    {
        int s = source.getId();
        DijkstraSPC sp = new DijkstraSPC(Gmask, s);
        //IJ.log("s: "+s);
        double minDistance = Double.POSITIVE_INFINITY;
        //int idxD = -1;
        int d = des.getId();
        if(d==s)
        {
            return 0;
        }    
        if (sp.hasPathTo(d)) 
        {
            double distanceInG = sp.distTo(d);
            if (minDistance > distanceInG) {
                minDistance = distanceInG;
                //idxD = d;
            }
        }
        if(minDistance < Double.POSITIVE_INFINITY)
        {
            return minDistance;
        }    
        else{
            return -1;
        }
        
    }        
    
    
    /**
     * min distance from each element in the source list 
     * to all element in the destination list
     * @param sourceLst
     * @param desLst 
     */
    public int computeLayerDistance(ArrayList<Vertex> sourceLst, ArrayList<Vertex> desLst)
    {
        //ArrayList<Integer> idxShortestLayer = new ArrayList<Integer>();
        ArrayUtil arr = new ArrayUtil(sourceLst.size());
        double maxDistanceLayer = -1;
        for (Vertex vs : sourceLst) 
        {
            int s = vs.getId();
            double minDis = Double.POSITIVE_INFINITY;
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                double dis = computeDistance(vs, vd);
                if(dis!=-1 && dis < minDis)
                {
                    minDis = dis;
                }    
            }
            if(minDis<Double.POSITIVE_INFINITY)
            {
                vs.setDistanceLayer(minDis);
                if(maxDistanceLayer < minDis)
                {
                    maxDistanceLayer = minDis;
                }    
            }
        }
        return (int)maxDistanceLayer;
    }
    
    public void computeShortestDistance(ArrayList<Vertex> sourceLst, ArrayList<Vertex> desLst)
    {
        //ArrayList<Integer> idxShortestLayer = new ArrayList<Integer>();
        ArrayUtil arr = new ArrayUtil(sourceLst.size());
        for (Vertex vs : sourceLst) 
        {
            vs.setShortestDistance(-1);
            int s = vs.getId();
            double minDis = Double.POSITIVE_INFINITY;
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                double dis = computeDistance(vs, vd);
                if(dis!=-1 && dis < minDis)
                {
                    minDis = dis;
                }    
            }
            if(minDis<Double.POSITIVE_INFINITY)
            {
                vs.setShortestDistance(minDis);
            }
        }
    }
    /**
     * min distance from each element in the source list 
     * to all element in the destination list
     * @param sourceLst
     * @param desLst 
     */
    public ArrayUtil computeMinDistance(ArrayList<Vertex> sourceLst, ArrayList<Vertex> desLst)
    {
        //ArrayList<Integer> idxShortestLayer = new ArrayList<Integer>();
        ArrayUtil arr = new ArrayUtil(sourceLst.size());
        int count = 0;
        for (Vertex vs : sourceLst) 
        {
            int s = vs.getId();
            double minDis = Double.POSITIVE_INFINITY;
            for (Vertex vd : desLst) 
            {
                int d = vd.getId();
                double dis = computeDistance(vs, vd);
                if(dis!=-1 && dis < minDis)
                {
                    minDis = dis;
                }    
            }
            if(minDis<Double.POSITIVE_INFINITY)
            {
                vs.setShortestDistance(minDis);
                arr.putValue(count, minDis);
            }
            else{
                arr.putValue(count, 0);
                IJ.log("Not connected component");
            }
            count++;
        }
        arr.sort();
        return arr;
    }
    
    /**
     * 
     * @param source
     * @return 
     */
    public double getClosest(Vertex source, boolean in)
    {
        int s = source.getId();
        DijkstraSPC sp = null;
        if(in)  //cell only, if true
        {
            sp = new DijkstraSPC(G, s);
        }
        else{     //cell + unlabelled for domain
            if(Gmask==null)
            {
                IJ.log("Forget to set gmask");
                return -1;
            }    
            sp = new DijkstraSPC(Gmask, s);
        }
        
        //IJ.log("s: "+s);
        double minDistance = Double.POSITIVE_INFINITY;
        //int idxD = -1;
        Vertex tmp = null;
        ArrayList<Vertex> lsVertexes = this.getVertexes();
        for(Vertex v : lsVertexes)
        {
            int vId = v.getId();
            if(vId != s && sp.hasPathTo(vId))
            {
                double distanceInG = sp.distTo(vId);
                if (minDistance > distanceInG) 
                {
                    minDistance = distanceInG;
                    //idxD = vId;
                    tmp = v;
                }
                
            }
        }
        if(minDistance < Double.POSITIVE_INFINITY)
        {
            return minDistance;
        }    
        else{
            IJ.log("Not connected component, don't find out the closest element");
            return -1; 
        }
        //return tmp;
    }  
    
    public ArrayUtil distancesAllClosestCenter() 
    {
        int nb = this.getNbObjects();
        ArrayUtil tab = new ArrayUtil(nb);
        for (int i = 0; i < nb; i++) 
        {
            double dis = this.getClosest(this.getObject(i), false);
            if(dis!=-1)
            {
                tab.putValue(i, dis);
            }
            else{
                IJ.log("Not connected component, can not compute the shortest distance between 2 cells");
                tab.putValue(i, 0);
            }
        }
        return tab;
    }
    public void addImage(ImageInt seg, int threshold) 
    {
        seg.resetStats(null);
        int min = (int) seg.getMinAboveValue(threshold);
        int max = (int) seg.getMax();
        if (max == 0) {
            IJ.log("No objects found");
            return;
        }
        //IJ.log("mm "+min+" "+max);
        // iterate in image  and constructs objects
        //calibration = cali;
        ArrayList<Voxel3D>[] objectstmp = new ArrayList[max - min + 1];
        for (int i = 0; i < max - min + 1; i++) {
            objectstmp[i] = new ArrayList<Voxel3D>();
        }
        int pix;
        int sz = seg.sizeZ;
        int sy = seg.sizeY;
        int sx = seg.sizeX;

        for (int k = 0; k < sz; k++) {
            for (int j = 0; j < sy; j++) {
                for (int i = 0; i < sx; i++) {
                    pix = seg.getPixelInt(i, j, k);
                    if (pix > threshold) {
                        objectstmp[pix - min].add(new Voxel3D(i, j, k, pix));
                    }
                }
            }
        }
        // ARRAYLIST   
        for (int i = 0; i < max - min + 1; i++) {
            if (!objectstmp[i].isEmpty()) {
                //Object3DVoxels ob = new Object3DVoxels(objectstmp[i]);
                //ob.setLabelImage(null);// the image can be closed anytime
                //ob.setCalibration(cali);
                //ob.setName("Obj" + (i + 1));
                Vertex v = new Vertex(i, "Cell " +(i + 1), (byte)0); // 0 : UNLABELLED
                addVertexes(v);
            }
        }
    }
    
    /**
     * for reference point
     * @return 
     */
    public Vertex getRandomCell() 
    {
        if(vertexes.isEmpty())
        {
            return null;
        }    
        int ranInit = (int) (Math.random() * vertexes.size());
        Vertex vr = vertexes.get(ranInit);
        return vr;
    }
    
    /**
     * for reference point
     * @return 
     */
    public Vertex getRandomCell(ArrayList<Vertex> lsVertexes) 
    {
        if(lsVertexes.isEmpty())
        {
            return null;
        }    
        int ranInit = (int) (Math.random() * lsVertexes.size());
        Vertex vr = lsVertexes.get(ranInit);
        return vr;
    }
    
    public void getRandomCellInPopulation(int nb, ArrayList<Vertex> lsVertexes) 
    {
        if(lsVertexes.isEmpty())
        {
            exit();
        }   
        for (int i = 0; i < nb; i++) 
        {
            boolean flag = false; 
            Vertex vr;
            while(!flag)
            {
                int ranInit = (int) (Math.random() * lsVertexes.size());
                vr = lsVertexes.get(ranInit);
                ArrayList<Vertex> lsTmp = getVertexes();
                if(!lsTmp.contains(vr))
                {
                    addVertexes(vr);
                    flag = true;
                } 
            }    
            
        }
        
    }
    
    public ArrayUtil computeDistances(Vertex[] evaluationPoints) {
        final int numPoints = evaluationPoints.length;
        ArrayUtil arr;
        Vertex P;
        
        //Vertex cl;
        arr = new ArrayUtil(numPoints);
        for (int i = 0; i < numPoints; i++) {
            P = evaluationPoints[i];
            double minDistance = getClosest(P, false);
            arr.putValue(i, minDistance);
        }

        return arr;
    }
    
    public String toString()
    {
        String str = "";
        for(Vertex v : vertexes)
        {
            str += v.getId()+"  ";
        }
        return str;
    }        

}
