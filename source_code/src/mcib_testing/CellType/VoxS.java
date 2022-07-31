/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.CellType;

/**
 *
 * @author tranhoa
 */

public class VoxS implements Comparable<VoxS> {

    public double distance;
    public double index;
    public int x, y, z;

    public VoxS(double distance, int x, int y, int z) {
        this.distance = distance;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public VoxS(double distance, double index, int x, int y, int z) {
        this.distance = distance;
        this.index = index;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    @Override
    public int compareTo(VoxS v) {
        if (distance > v.distance) {
            return 1;
        } else if (distance < v.distance) {
            return -1;
        } else {
            return 0;
        }
    }
}