package stdlib;


import stdlib.DirectedEdgeC;
import ij.IJ;
import java.util.ArrayList;
import stdlib.Bag;
import stdlib.Stack;

/*************************************************************************
 *  Compilation:  javac EdgeWeightedDigraph.java
 *  Execution:    java EdgeWeightedDigraph V E
 *  Dependencies: Bag.java DirectedEdge.java
 *
 *  An edge-weighted digraph, implemented using adjacency lists.
 *
 *************************************************************************/

/**
 *  The <tt>EdgeWeightedDigraph</tt> class represents a edge-weighted
 *  digraph of vertices named 0 through <em>V</em> - 1, where each
 *  directed edge is of type {@link DirectedEdgeC} and has a real-valued weight.
 *  It supports the following two primary operations: add a directed edge
 *  to the digraph and iterate over all of edges incident from a given vertex.
 *  It also provides
 *  methods for returning the number of vertices <em>V</em> and the number
 *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which 
 *  is a vertex-indexed array of @link{Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the edges incident from a given vertex, which takes
 *  time proportional to the number of such edges.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/44sp">Section 4.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class EdgeWeightedDigraphC {
    public ArrayList<Vertex> vertexes;
    public ArrayList<DirectedEdgeC> edges;
    private int V;
    private int E;
    private Bag<DirectedEdgeC>[] adj;
    
    /**
     * Initializes an empty edge-weighted digraph with <tt>V</tt> vertices and 0 edges.
     * param V the number of vertices
     * @throws java.lang.IllegalArgumentException if <tt>V</tt> < 0
     */
    public EdgeWeightedDigraphC(int V) {
        if (V < 0) throw new IllegalArgumentException("Number of vertices in a Digraph must be nonnegative");
        this.V = V;
        this.E = 0;
        adj = (Bag<DirectedEdgeC>[]) new Bag[V];
        for (int v = 0; v < V; v++)
            adj[v] = new Bag<DirectedEdgeC>();
    }

    /**
     * Initializes a random edge-weighted digraph with <tt>V</tt> vertices and <em>E</em> edges.
     * param V the number of vertices
     * param E the number of edges
     * @throws java.lang.IllegalArgumentException if <tt>V</tt> < 0
     * @throws java.lang.IllegalArgumentException if <tt>E</tt> < 0
     */
//    public EdgeWeightedDigraphC(int V, int E) {
//        this(V);
//        if (E < 0) throw new IllegalArgumentException("Number of edges in a Digraph must be nonnegative");
//        for (int i = 0; i < E; i++) {
//            int v = (int) (Math.random() * V);
//            int w = (int) (Math.random() * V);
//            double weight = Math.round(100 * Math.random()) / 100.0;
//            Vertex f = new Vertex(v, "cell_"+v);
//            Vertex t = new Vertex(w, "cell_"+w);
//            DirectedEdgeC e = new DirectedEdgeC(f, t, weight);
//            addEdge(e);
//        }
//    }

    /**  
     * Initializes an edge-weighted digraph from an input stream.
     * The format is the number of vertices <em>V</em>,
     * followed by the number of edges <em>E</em>,
     * followed by <em>E</em> pairs of vertices and edge weights,
     * with each entry separated by whitespace.
     * @param vertexes
     * @param in the input stream
     * @throws java.lang.IndexOutOfBoundsException if the endpoints of any edge are not in prescribed range
     * @throws java.lang.IllegalArgumentException if the number of vertices or edges is negative
     */
//    public EdgeWeightedDigraph(In in) {
//        this(in.readInt());
//        int E = in.readInt();
//        if (E < 0) throw new IllegalArgumentException("Number of edges must be nonnegative");
//        for (int i = 0; i < E; i++) {
//            int v = in.readInt();
//            int w = in.readInt();
//            if (v < 0 || v >= V) throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V-1));
//            if (w < 0 || w >= V) throw new IndexOutOfBoundsException("vertex " + w + " is not between 0 and " + (V-1));
//            double weight = in.readDouble();
//            addEdge(new DirectedEdge(v, w, weight));
//        }
//    }
    
//    public EdgeWeightedDigraphC(In in) 
//    {
//        //this(in.readInt());
//        setV(in.readInt());
//        System.out.println("V:"+this.V());
//        while (in.hasNextLine() && !in.isEmpty()) 
//        {
//            int v = in.readInt();
//            int w = in.readInt();
//            double weight = in.readDouble();
//            //System.out.println("v: "+v+" w: " + w + " weight: "+weight);
//            Vertex f = new Vertex(v, "cell_"+v);
//            Vertex t = new Vertex(w, "cell_"+w);
//            addEdge(new DirectedEdgeC(f, t, weight));
//        }
//    }
    
    /**
     * 
     * @param vertexes
     * @param edges 
     */
    public EdgeWeightedDigraphC(ArrayList<Vertex> vertexes, ArrayList<DirectedEdgeC> edges)
    {
        this(vertexes.size());
        this.vertexes = vertexes;
        this.edges = edges;
//        IJ.log("V value is: "+this.V());
//        IJ.log("Nb edges: "+ edges.size());
        //this.E = 0;
        for(DirectedEdgeC ei : edges)
        {
            addEdge(ei);
        }
        
    }

    /**
     * Initializes a new edge-weighted digraph that is a deep copy of <tt>G</tt>.
     * @param G the edge-weighted graph to copy
     */
    public EdgeWeightedDigraphC(EdgeWeightedDigraphC G) {
        this(G.V());
        this.E = G.E();
        for (int v = 0; v < G.V(); v++) {
            // reverse so that adjacency list is in same order as original
            Stack<DirectedEdgeC> reverse = new Stack<DirectedEdgeC>();
            for (DirectedEdgeC e : G.adj[v]) {
                reverse.push(e);
            }
            for (DirectedEdgeC e : reverse) {
                adj[v].add(e);
            }
        }
    }

    /**
     * Returns the number of vertices in the edge-weighted digraph.
     * @return the number of vertices in the edge-weighted digraph
     */
    public int V() {
        return V;
    }
    public void setV(int V) {
        this.V = V;
    }

    /**
     * Returns the number of edges in the edge-weighted digraph.
     * @return the number of edges in the edge-weighted digraph
     */
    public int E() {
        return E;
    }
    public void setE(int E) {
        this.E = E;
    }

    // throw an IndexOutOfBoundsException unless 0 <= v < V
    private void validateVertex(int v) {
        if (v < 0 || v >= V)
            throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V-1));
    }

    /**
     * Adds the directed edge <tt>e</tt> to the edge-weighted digraph.
     * @param e the edge
     * @throws java.lang.IndexOutOfBoundsException unless endpoints of edge are between 0 and V-1
     */
    private void addEdge(DirectedEdgeC e) {
        int v = e.from().getId();
        int w = e.to().getId();
        validateVertex(v);
        validateVertex(w);
        adj[v].add(e);
        E++;
    }


    /**
     * Returns the directed edges incident from vertex <tt>v</tt>.
     * @return the directed edges incident from vertex <tt>v</tt> as an Iterable
     * @param v the vertex
     * @throws java.lang.IndexOutOfBoundsException unless 0 <= v < V
     */
    public Iterable<DirectedEdgeC> adj(int v) {
        validateVertex(v);
        return adj[v];
    }

    /**
     * Returns the number of directed edges incident from vertex <tt>v</tt>.
     * This is known as the <em>outdegree</em> of vertex <tt>v</tt>.
     * @return the outdegree of vertex <tt>v</tt>
     * @param v the vertex
     * @throws java.lang.IndexOutOfBoundsException unless 0 <= v < V
     */
    public int outdegree(int v) {
        validateVertex(v);
        return adj[v].size();
    }

    /**
     * Returns all directed edges in the edge-weighted digraph.
     * To iterate over the edges in the edge-weighted graph, use foreach notation:
     * <tt>for (DirectedEdge e : G.edges())</tt>.
     * @return all edges in the edge-weighted graph as an Iterable.
     */
    public Iterable<DirectedEdgeC> edges() {
        Bag<DirectedEdgeC> list = new Bag<DirectedEdgeC>();
        for (int v = 0; v < V; v++) {
            for (DirectedEdgeC e : adj(v)) {
                list.add(e);
            }
        }
        return list;
    } 

    /**
     * Returns a string representation of the edge-weighted digraph.
     * This method takes time proportional to <em>E</em> + <em>V</em>.
     * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
     *   followed by the <em>V</em> adjacency lists of edges
     */
    public String toString() {
        String NEWLINE = System.getProperty("line.separator");
        StringBuilder s = new StringBuilder();
        s.append(V + " " + E + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (DirectedEdgeC e : adj[v]) {
                s.append(e + "  ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }

    /**
     * Unit tests the <tt>EdgeWeightedDigraph</tt> data type.
     */
//    public static void main(String[] args) {
//        //In in = new In(args[0]);
//        In in = new In("/home/ttnhoa/NetBeansProjects/MyGraph/data/mediumEWG_2.txt");
//        EdgeWeightedDigraphC G = new EdgeWeightedDigraphC(in);
//        StdOut.println(G);
//    }

}
