/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing;

import ij.IJ;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

/**
 *
 **
 * /**
 * Copyright (C) 2008- 2012 Thomas Boudier and others
 *
 *
 *
 * This file is part of mcib3d
 *
 * mcib3d is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author thomas
 */
public class AssociationMitosis {

    public double[][] matrix;
    int[][] order;
    int[] association;

    public void initAll(int nb) {
        int a = nb;
        matrix = new double[a][a];
        association = new int[a];
        order = new int[a][a];
        for (int i = 0; i < a; i++) {
            matrix[i] = new double[a];
            order[i] = new int[a];
            association[i] = -1;
            for (int j = 0; j < a; j++) {
                matrix[i][j] = -1;
                order[i][j] = -1;
            }
        }
    }

    public void initAsso() {
        int a = association.length;
        for (int i = 0; i < a; i++) {
            association[i] = -1;
        }
    }

    public void setAssociation(int i, int j, double val) {
        matrix[i][j] = val;
    }

    public void setAssociation(int[] asso) {
        association = asso;
    }

    public double[][] duplicateMatrix(double[][] mat) {
        int a = mat.length;
        double[][] copy = new double[a][a];
        for (int i = 0; i < a; i++) {
            System.arraycopy(mat[i], 0, copy[i], 0, a);
        }
        return copy;
    }

    public int[] duplicateAssociation(int[] asso) {
        int a = asso.length;
        int[] copy = new int[a];
        System.arraycopy(asso, 0, copy, 0, a);
        return copy;
    }

    public void deleteNonBijectif() {
        int a = matrix.length;
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < a; j++) {
                if ((matrix[i][j] < 0) || (matrix[j][i] < 0)) {
                    matrix[i][j] = -1;
                    matrix[j][i] = -1;
                }
            }
        }
    }

    public void computeOrders() {
        ArrayList list = new ArrayList();
        int a = matrix.length;
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < a; j++) {
                if (matrix[i][j] >= 0) {
                    list.add(matrix[i][j]);
                }
            }
        }
        Collections.sort(list);
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < a; j++) {
                if (matrix[i][j] >= 0) {
                    order[i][j] = list.indexOf(matrix[i][j]);
                }
            }
        }
    }

    private ArrayList<Integer> availableAssociation(double[][] mat, int n) {
        ArrayList<Integer> asso = new ArrayList();
        for (int i = 0; i < mat.length; i++) {
            if (mat[n][i] >= 0) {
                asso.add(i);
            }
        }
        return asso;
    }

    private int randomAssociation(double[][] mat, int n, int i0) {
        ArrayList<Integer> assos = this.availableAssociation(mat, n);
        int nb = assos.size();
        if (nb <= 1) {
            if (nb == 0) {
                return -1;
            } else {
                return assos.get(0);
            }
        } else {
            Random r = new Random();
            int id = r.nextInt(nb);
            while (assos.get(id) == i0) {
                id = r.nextInt(nb);
            }
            return assos.get(id);
        }
    }

    private void keepUnique(double[][] mat, int i0, int j0) {
        int a = mat.length;
        for (int i = 0; i < a; i++) {
            if (i != i0) {
                mat[i][j0] = -1;
            }
        }
        for (int j = 0; j < a; j++) {
            if (j != j0) {
                mat[i0][j] = -1;
            }
        }
    }

    public void randomAssociation(double[][] mat) {
        int a = mat.length;
        // shuffle array
        ArrayList<Integer> rand = new ArrayList(a);
        for (int i = 0; i < a; i++) {
            rand.add(i);
        }
        Collections.shuffle(rand);
        for (int i = 0; i < a; i++) {
            int ii = rand.get(i);
            ArrayList<Integer> list = this.availableAssociation(mat, ii);
            if (list.size() > 1) {
                //IJ.log((ii+1) + " has " + list.size() + " associations");
            }
            int j = this.randomAssociation(mat, ii, -1);
            if (j >= 0) {
                association[ii] = j;
                association[j] = ii;
                this.keepUnique(mat, ii, j);
                this.keepUnique(mat, j, ii);
            }

        }
    }

    public int[] getMinCostRandomAssociation(int nbIte) {
        int[] res = null;
        double Cmin = Double.MAX_VALUE;
        double[][] B;
        for (int i = 0; i < nbIte; i++) {
            B = duplicateMatrix(matrix);
            initAsso();
            randomAssociation(B);
            double cost = getCostAssociation();
            if (cost < Cmin) {
                //System.out.println("");
                //System.out.println("cost=" + cost);
                //printAssociation(mito.association, 1);
                Cmin = cost;
                res = duplicateAssociation(association);
            }
        }
        return res;
    }

    public String printAssociation(int[] asso, int offset) {
        System.out.println();
        String out = "";
        for (int i = 0; i < asso.length; i++) {
            if (asso[i] >= 0) {
                if (i < asso[i]) {
                    out = out.concat(" " + (i + offset) + "-" + (asso[i] + offset));
                }
            }
        }
        return out;
    }
    
    public int[] getLastAssociation(){
        return association;
    }

    public double getCostAssociation() {
        double sum = 0;
        for (int i = 0; i < association.length; i++) {
            int id = association[i];
            if (id >= 0) {
                sum += matrix[i][id];
            }
        }
        return sum;
    }

    double getMinCost(int[] asso, double[][] mat) {
        double coutMin = Double.MAX_VALUE;
        int nb = asso.length;
        int i0 = 0;
        while ((asso[i0] >= 0) && (i0 < nb)) {
            i0++;
        }

        if (i0 == nb) {
            return coutMin;
        }
        ArrayList<Integer> list = availableAssociation(mat, i0);
        if (list.isEmpty()) {
            return coutMin;
        }
        for (int i = 0; i < list.size(); i++) {
            int j0 = list.get(i);
            IJ.log("" + i0 + "/" + j0);
            int[] asso1 = this.duplicateAssociation(asso);
            asso1[i0] = j0;
            asso1[j0] = i0;
            double cout = getMinCost(asso1, mat);
            if (cout < coutMin) {
                coutMin = cout;
                asso = asso1;
            }
            mat = duplicateMatrix(matrix);
        }

        printAssociation(asso, 1);
        return getCostAssociation();
    }
}
