% -------------------------------------------------------------------------
% trackpreview.m 
% 
% -------------------------------------------------------------------------
% Description:
%
% Before tracking starts, trackpreview allows the user to tune the
% thresholding given the particular illumination conditions in the arena.
% -------------------------------------------------------------------------

function [] = trackpreview()

global vidobj  thresholdbackground  % nextbboxref goodbackground boxsize previewtimer

% Getting a snapshot and trying to find the animal by changing the threshold
mypreviewframe=getsnapshot(vidobj); 
itest=im2bw(mypreviewframe,thresholdbackground);
itest=1-itest;
itest=imclearborder(itest);
[LO,num]=bwlabel(itest,8);
Areas=regionprops(LO,'area');
area_val=[Areas.Area];
maxarea=max(area_val);
idxBig= find(maxarea == area_val);
it2=ismember(LO,idxBig);
% it2=1-it2;
% it2=imclearborder(it2);
[r c]=find(it2); 
maxc=max(c)+5;
minc=min(c)-5;
maxr=max(r)+5;
minr=min(r)-5;
% Prevent negative values at the edges
if minr<0
    minr=0;
end
if minc<0
    minc=0;
end
bbox=[minc minr maxc-minc maxr-minr];
%nextbboxref(1,:)=[bbox(1)-50  bbox(2)-50  bbox(3)+100  bbox(4)+100];

% Raw and processed images in a little box
smallImage0=imcrop(mypreviewframe,bbox);
%smallImage=imcrop(it2,bbox);
smallImage=imcrop(itest,bbox);

% % Display the result
 close all
 figure,
 subplot(2,2,1), imshow(mypreviewframe);
 subplot(2,2,2), imshow(itest);
 subplot(2,2,3), imshow(smallImage0);
 subplot(2,2,4), imshow(smallImage);


% Ask whether to change the threshold: if Y (yes), iterate another preview
entry=input('Do you want to change threshold ? [Y]/[N]  ','s');
if entry=='Y'
thresholdbackground=input('Introduce new threshold value (between 0 and 1): ');
%display('Type "return" to check the new preview');
%keyboard;
previewtimer=timer('TimerFcn','trackpreview','TaskstoExecute',1);
start(previewtimer);
end
       
end
