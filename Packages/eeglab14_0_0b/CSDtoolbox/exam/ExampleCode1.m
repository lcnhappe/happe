% ====================
% CSD Toolbox Tutorial
% ====================
%
% ExampleCode1.m
% --------------
% MatLab code to supplement the CSD Toolbox online tutorial
% (URL: "http://psychophysiology.cpmc.columbia.edu/CSDtoolbox.htm")
%
% Updated: $Date: 2009/05/20 19:34:00 $ $Author: jk $
%
% Changed: $Date: 2010/07/20 12:00:00 $ $Author: jk $
%        - removed closing MapMontage figure
%        - use new figures for ERP and CSD line plots

% ------------ Step 1 -----------------------------------------------------
% understand the spherical coordinate system of the CSD toolbox
% ------------ Step 2 -----------------------------------------------------
E = textread('E31.asc','%s');                                      
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E);  
MapMontage(M);
% ------------ Step 3 -----------------------------------------------------
tic
[G,H] = GetGH(M);
toc
% ------------ Step 4 -----------------------------------------------------
D = textread('NR_C66_trr.dat');
D = D';
% ------------ Step 5 -----------------------------------------------------
tic
X = CSD (D, G, H);
toc
X = X';
figure;           % create new figure
plot(X);
plot(X(:,[14 24]));
D = D';
figure;           % create new figure 
plot(D);
plot(D(:,[14 24]));
% ------------ Step 6 -----------------------------------------------------
WriteMatrix2Text(X,'CSD_C66_trr.dat');