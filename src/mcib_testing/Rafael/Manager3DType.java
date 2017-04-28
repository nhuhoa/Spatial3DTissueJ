/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mcib_testing.Rafael;

import ij.IJ;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import mcib3d.geom.Object3D;
import mcib_plugins.tools.RoiManager3D_2;

/**
 *
 * @author thomasb
 */
public class Manager3DType extends RoiManager3D_2 {
    
    String directory = "";
    
    public Manager3DType() {
        super();
        list.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent ke) {
                String key = "" + ke.getKeyChar();
                IJ.log("Key Pressed " + key);
                // sort nuclei
                if (key.equalsIgnoreCase("Z")) {
                    save();
                    return;
                }
                if ((!key.equalsIgnoreCase("0")) && (!key.equalsIgnoreCase("1")) && (!key.equalsIgnoreCase("2")) && (!key.equalsIgnoreCase("3"))) {
                    return;
                }
                
                int[] indices = list.getSelectedIndices();
                for (int i : indices) {
                    objects3D.getObject(i).setType(Integer.parseInt(key));
                    objects3D.getObject(i).setComment("");
                    if (objects3D.getObject(i).getType() > 0) {
                        updateName(i);
                    } else {
                        updateName(i);
                    }
                }
                
                list.updateUI();
            }
        });
    }
    
    public void save() 
    {
        IJ.log("Saving list");
        try {
            BufferedWriter bf = new BufferedWriter(new FileWriter(directory + "UpdatedRoi3D.list"));
            for (int i = 0; i < model.getSize(); i++) {
                bf.write(model.getElementAt(i).toString());
                bf.write("\n");
            }
            bf.close();
        } catch (IOException e) {
            IJ.log("Pb saving list " + e);
        }
    }
    
    public void load() {
        IJ.log("Loading list");
        String data;
        String[] spli;
        try {
            BufferedReader br = new BufferedReader(new FileReader(directory + "UpdatedRoi3D.list"));
            data = br.readLine();
            while (data != null) {
                spli = data.split(" ");
                String name = spli[0];
                // type
                int type = 0;
                if (spli.length > 1) {
                    String typep = spli[spli.length - 1];
                    System.out.println("load " + data + " " + name + " " + typep);
                    type = Integer.parseInt(typep.substring(1, 2));
                } else {
                    data = br.readLine();
                    continue;
                }
                Object3D obj = objects3D.getObjectByName(name);
                int ty = obj.getType();
                if (type != ty) {
                    IJ.log("Updated type for  " + name);
                    int id = objects3D.getIndexFromName(name);
                    objects3D.getObject(id).setType(type);
                    updateName(id);
                }
                // comment
                String comment = "";
                if (spli.length == 3) {
                    comment = spli[1];
                }
                String co = objects3D.getObjectByName(name).getComment();
                if (!comment.equals(co)) {
                    IJ.log("Updated comment for  " + name);
                    int id = objects3D.getIndexFromName(name);
                    objects3D.getObject(id).setComment(comment);
                    updateName(id);
                }
                data = br.readLine();
            }
        } catch (IOException e) {
            IJ.log("Pb loading list " + e);
        }
    }
}
