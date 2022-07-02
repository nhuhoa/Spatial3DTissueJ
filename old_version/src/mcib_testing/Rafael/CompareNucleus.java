/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Rafael;

import java.util.Comparator;
import mcib3d.geom.Object3D;

/**
 *
 * @author thomasb
 */
public class CompareNucleus implements Comparator<Object3D> {

    @Override
    public int compare(Object3D t1, Object3D t2) {
        // if contains * (coloc) then first
        if ((t1.getComment().contains("*")) && (!t2.getComment().contains("*"))) {
            return 1;
        } else if ((!t1.getComment().contains("*")) && (t2.getComment().contains("*"))) {
            return -1;
        }
        // if both contains * or both does not contain *, sort on type
        if ((t1.getType() > 0) && (t2.getType() == 0)) {
            return 1;
        } else if ((t1.getType() == 0) && (t2.getType() > 0)) {
            return -1;
        }
        // if both no type, sort on number
        if ((t1.getType() > 0) && (t2.getType() > 0)) {
            if (t1.getType() > t2.getType()) {
                return 1;
            } else if (t1.getType() < t2.getType()) {
                return -1;
            }
        }
        // sort on number
        int nb1 = Integer.parseInt(t1.getName().substring(3));
        int nb2 = Integer.parseInt(t2.getName().substring(3));

        return (int) (Math.signum(nb1 - nb2));
    }

}
