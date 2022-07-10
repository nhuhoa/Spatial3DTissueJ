//======================================================================================
//// Installation
//// See at https://github.com/nhuhoa/Spatial3DTissueJ/tree/master/ImageJ_plugins_jar/list_plugin_readme.txt
//// First, you can copy these 3 folders into Fiji/plugins/ or ImageJ/plugins/ and quick check if there is any duplicate plugin, ex: 3DViewerXXXvXXX//// .jar, Fiji_PluginXXX.jar. You can keep only one version for each plugin. 
//// To update plugins, you can delele the folder 3D_suite/ and replace by the most updated folder from my github, most updated plugin: Spatial3DTissueJ_v22_windows.jar

//======================================================================================
//// Preparation
//// Create a folder for each islet, ex: H1536_islet1/ as this one and put composite image of a given islet into this folder. 
//// Copy the color mapping file into this folder H1536_islet1/: cell_type_colormap.lut
When you download github files, you have this file in the github folder: https://github.com/nhuhoa/Spatial3DTissueJ/tree/master/macro/macro_windows/cell_type_colormap.lut
//// Noted: cell_type_colormap.lut is just a color map for visualization, if macro has an issue with this file, you can just remove the line 
//// And change dir variable to your directory path

//======================================================================================
//// How to run a macro here
//// I divide macro into many steps, you can run one by one step here
//// And then when you don't see any error, you can try to run entire pipeline
//// I tested macro in Windows systems and for tissue H1536_islet1
