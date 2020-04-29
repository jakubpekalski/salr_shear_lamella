Scripts for pattern recogniton used for the article "Shear induced ordering in systems with competing interactions, a
machine learning based study".

*CNN_activations.ipynb* -  structure recogniton with convolutional neural networks (uses Keras)

*ML_shear_T042.ipynb*   - structure recognition with dense neural networks (uses tensorflow 1)

*salr.py*         - [HOOMD-Blue](http://glotzerlab.engin.umich.edu/hoomd-blue/) script for temperature annealing in npt ensemble

*compute_.py*     - [OVITO](http://ovito.org) scripts for computing structural properties descrbied in the paper

*in.salr_shear*   - [LAMMPS](http://lammps.sandia.gov) script for imposing shear via [SLLOD] (https://lammps.sandia.gov/doc/fix_nvt_sllod.html) equations of motion. Note: a lammps extension: fix npt/sllod is required
