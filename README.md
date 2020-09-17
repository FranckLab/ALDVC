# Augmented Lagrangian Digital Volume Correlation (ALDVC)
Augmented Lagrangian Digital Volume Correlation, or we call ALDVC, is to measure full-field volumetric displacements and strains based on a hybrid local-global approach.
ALDVC is a fast, parallel-computing hybrid DVC algorithm, which combines advantages of local subset method (fast computation speed, and parallel computing) and finite-element-based global method (guarantee global kinematic compatibility and decrease noise).  

## Advantages of ALDVC algorithm
* [1] It’s a fast algorithm using distributed parallel computing.  
* [2]	Global kinematic compatibility is added as a global constraint in the form of augmented Lagrangian, and solved using Alternating Direction Method of Multipliers scheme.
* [3]	Both displacement fields and affine deformation gradients are correlated at the same time.
* [4]	No need of much manual experience about choosing displacement smoothing filters.
* [5]	Being able to compute image sequence with multiple image frames, which is especially quite useful for measuring very large deformations.

## Prerequisites & Installation
ALDVC MATLAB code was tested on MATLAB versions later than R2018a. Both single thread and parallel computing features are included in ALDVC code. Please download and unzip the code to the MATLAB working path. Then, execute the mail file main_ALDVC.m.

## ALDVC example dataset
Example images download link: 
https://uwmadison.box.com/s/kr6pjje9yfi9k9va6cbk1ayvhvbs5nvq <p>
Example results download link: 
https://uwmadison.box.com/s/wdhmysbehidobwd6oh8tn120ufjq1iov
 
## ****** ATTENTION ******  
% The "x,y,z" or "1-,2-,3-" coordinates in the ALDVC code always correspond to the 1st, 2nd and 3rd indices of Matlab workspace variable. For example, p_meas(:,1) and p_meas(:,2) are the x- & y-coordinates of scattered points.  
 
% This is a little different from some MATLAB image processing functions. 
% For example, if a 3D image has size MxNxL, in this code, we always have the image size_x=M, size_y=N, size_z=L. If you use some Matlab computer vision/image post-processing function, for example, 'imagesc3D', or 'imshow3D', or 'surf', it will reads size_x=N, size_y=M, size_z=L. 
 
% Please pay attention to this.  

## Citation
If used please cite
```bibtex
@article{Yang2020aldvc,
  title={Augmented Lagrangian Digital Volume Correlation (ALDVC)},
  author={Yang, J. and Hazlett, L. and Landauer, A. K. and Franck, C.},
  journal={Experimental Mechanics},
  volume={},
  number={},
  pages={},
  year={2020},
  Url={https://doi.org/10.1007/s11340-020-00607-3}
}
```
 
[1] Yang, J. Hazlett, L., Landauer, A. Franck, C. Augmented Lagrangian Digital Volume Correlation. Exp.Mech, 2020.  
* Full text can be requested at: 
Exp Mech Website: https://link.springer.com/article/10.1007/s11340-020-00607-3
ResearchGate: https://www.researchgate.net/publication/343676441_Augmented_Lagrangian_Digital_Volume_Correlation_ALDVC

## Contact and support
Email: aldicdvc@gmail.com;  -or- Jin Yang, jyang526@wisc.edu; -or- Prof. Christian Franck, cfranck@wisc.edu
Welcome to give us ratings and make comments at: [![View Augmented Lagrangian Digital Volume Correlation (ALDVC) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/77019-augmented-lagrangian-digital-volume-correlation-aldvc)


##

 
<p align="center">
  <img width="538" height="301" src="https://github.com/FranckLab/ALDVC/blob/master/aldvc_logo.png">
</p>



