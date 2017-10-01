# MARA

MARA ("Multiple Artifact Rejection Algorithm") is an open-source [EEGLAB](http://sccn.ucsd.edu/eeglab/) plug-in which automatizes the process of hand-labeling independent components for artifact 
rejection. 
      
The core of MARA is a supervised machine learning algorithm that learns from expert ratings of 1290 components by extracting six features from the spatial, the spectral 
and the temporal domain. Features were optimized to solve the binary classification problem "reject vs. accept". Thus, MARA is not limited to a specific type of artifact, 
and should be able to handle eye artifacts, muscular artifacts and loose electrodes equally well. 

For detailed instructions, see MARATutorial.pdf

#### Requirements

In addition to the requirements of EEGLAB, MARA needs the Matlab Statistics Toolbox, the Optimization Toolbox and the Signal Processing Toolbox.  

#### Installation 

Copy this folder into your EEGLAB plugins directory. Starting EEGLab will automatically recognize the plugin and MARA should appear in the menu Tools. 

#### Contact

If you have questions or suggestions about the toolbox, please contact: Irene Dowding, 
irenedowding16 at gmail.com
 

