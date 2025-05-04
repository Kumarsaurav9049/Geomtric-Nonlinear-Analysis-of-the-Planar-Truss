% Generealized Matalab Code for the Analysis of the 2D Truss.
% Saurav Kumar Nidhi, M.Tech structural Engineering. (2311CE60)
% Dr. Vaibhav Singhal. Associate professor
% Department of Civil and Enviornmental Engineering
% Indian Institute of Technology, Patna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
close all
format short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary('output or results of this code stored in this file.txt')
disp(['*********************************************************************' ...
    '**********************************************************************']);
disp(['*******************************************Fully Generealized Matalab Code for' ...
    ' the Analysis of the Planar Truss' ...
    '**************************************']);
disp(['**Saurav Kumar Nidhi, M.Tech structural Engineering. (2311CE60)**______' ...
    '_______________________**Dr. Vaibhav Singhal. Associate Professor **'])
disp(['**Department of civil and enviornmental Engineering**________________' ...
    '_____________________________**Indian Institute of Technology, Patna**'])
disp(['*********************************************************************' ...
    '**********************************************************************']);
model.type = input(['Enter the type of the analysis which you want to analyse (linear or nonlinear):-' ...
    '(enter in lower case letter)\n'],'s');

switch model.type
    case 'linear'
        linear_analysis_of_the_2d_truss()    % Function file for the linear analysis...
                                               % (other function file are in this code)
    case 'nonlinear'
        nonlinear_analysis_of_the_2d_truss()   % Function file for the nonlinear analysis... 
                                               % in continuation with the linear analysis...
                                               % (Also conatins the many more funcrion files)  
                                               % Establish coordiantes and global dof's of nodes
                                               % not_used_nonlinear_function_file()
                                               % earleir used as the
                                               % check function function
                                               % file
                                               
end
