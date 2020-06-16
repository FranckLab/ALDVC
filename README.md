# ALDVC
Augmented Lagrangian Digital Volume Correlation  - volumetric displacement and strain measurement based on a hybrid local-global approach
ALDVC is a fast, parallel-computing hybrid DVC algorithm, which combines advantages of local subset method (fast computation speed, and parallel computing) and finite-element-based global method (guarantee global kinematic compatibility and decrease noise).  

## Advantages of ALDVC algorithm
* [1] It’s a fast algorithm using distributed parallel computing.  
* [2]	Global kinematic compatibility is added as a global constraint in the form of augmented Lagrangian, and solved using Alternating Direction Method of Multipliers scheme.
* [3]	Both displacement fields and affine deformation gradients are correlated at the same time.
* [4]	No need of much manual experience about choosing displacement smoothing filters.
* [5]	Being able to compute image sequence with multiple image frames, which is especially quite useful for measuring very large deformations.

## Prerequisites & Installation
ALDVC MATLAB code was tested on MATLAB versions later than R2018a. Both single thread and parallel computing features are included in ALDVC code. Please download and unzip the code to the MATLAB working path. Then, execute the mail file main_ALDVC.m.

## DVC example dataset
### Example images download link: 
https://uwmadison.box.com/s/kr6pjje9yfi9k9va6cbk1ayvhvbs5nvq
### Example results download link: 
https://uwmadison.box.com/s/wdhmysbehidobwd6oh8tn120ufjq1iov


## Citation
[1] Yang, J. Hazlett, L., Landauer, A. Franck, C. Augmented Lagrangian Digital Volume Correlation. Exp.Mech, 2020.  
* Full text can be requested at: https://www.researchgate.net/publication/342182706_Augmented_Lagrangian_Digital_Volume_Correlation

## Contact and support
Email: aldicdvc@gmail.com;  -or- Jin Yang, jyang526@wisc.edu; -or- Prof. Christian Franck, cfranck@wisc.edu
Welcome to give us ratings and make comments at: [![View Augmented Lagrangian Digital Volume Correlation (ALDVC) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/77019-augmented-lagrangian-digital-volume-correlation-aldvc)


##

 
<p align="center">
  <img width="538" height="301" src="https://github.com/FranckLab/ALDVC/blob/master/aldvc_logo.png">
</p>



