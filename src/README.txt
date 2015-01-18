Our folder structure is the same as the structure of the Piotr Dollar Toolbox.
Simply copy the folder 'toolbox' into the toolbox. WARNING: toolboxCompile.m does override the original, but it is not strictly needed. toolboxRegCompile.m only compiles the needed .cpp files.


To run the code two scripts are provided:

RegressionRandomForest.m: This runs the training and evaluation of the standard Regression Forest.

camera_relocalization.m: This runs the code for the camera relocalization task. Training and evaluating might take some time.

The Data for the camera_relocalization.m script can be found at (we used heads, as it is the smallest dataset):

http://research.microsoft.com/en-us/projects/7-scenes/