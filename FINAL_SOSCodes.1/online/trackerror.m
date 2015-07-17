% -------------------------------------------------------------------------
% trackerror.m 
% 
% -------------------------------------------------------------------------
% Description:
%
% Error message mainly when exiting a timer. Current information is saved. 
% -------------------------------------------------------------------------

function [] = trackerror()

global dirout finalbbox impro imraw

% stop(vidobj)               
% delete(vidobj)

% Saving data before quitting 
save finalbbox finalbbox
save impro impro
save imraw imraw

% Giving some information
display('There was an error and SOS had to stop');
display('Possible problems: animal left field of view, too high tracking speed, incorrect threholding');
display('Remember to close "Video Preview" window before starting a new experiment');
cd(dirout)
close all

end
