% -------------------------------------------------------------------------
% track.m 
% 
% -------------------------------------------------------------------------
% Description:
%
% Main code for online animal tracking. It is launched as a function with two 
% arguments: period time in seconds between frames and total number of 
% frames to track. It calls other functions such as trackgo.m, 
% trackerror.m, etc. Its output is a time sequence of raw and black-and-white
% cropped animal images, together with bounding box coordinates in the 
% arena system of reference.
% -------------------------------------------------------------------------

function [] = track(period,nframe)

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------

% Global variables
global  i vidobj threshold thresholdbackground dirout dirdest impro imraw finalbbox x y l w nextbboxref goodbackground firstbg iscrossing boxsize

close all
display('SOS online tracking')


% Directories: set the right paths
dirout='/Users/joshsalvi/Documents/Lab/Lab/Zebrafish/Video Tracking/';
%
% Create a new folder directory to save the data 
date=datestr(clock,'mmm.dd,yyyy HH.MM.SS');
dirdest = fullfile('/Users/joshsalvi/Documents/Lab/Lab/Zebrafish/Video Tracking/',date,'onlineData');
mkdir(dirdest);
cd(dirdest);

% -------------------------------------------------------------------------
% Talk to camera
vidobj=videoinput('dcam',1,'Y8_1024X768');
%
src=getselectedsource(vidobj);
set (src, 'ShutterMode', 'manual'); 
set (src, 'Shutter', 230);

% -------------------------------------------------------------------------
% View video preview of the field of view before tracking starts
preview(vidobj)%,'visible=off');  
pause(2)
%
display('Write "return" when the animal is in the arena');
keyboard;
%
% stoppreview
% closepreview
% stoppreview(vidobj)

% -------------------------------------------------------------------------
% Initializing 
impro={};
imraw={};
nextbboxref=zeros(nframe,4);


% -------------------------------------------------------------------------
% Call the preview function to iteratively optimize the threshold 
% -------------------------------------------------------------------------

firstbg= getsnapshot(vidobj); 
goodbackground=firstbg; % that will be used for subtraction afterwards
% Initial threshold value (that we can manually tune before tracking)
thresholdbackground=0.4; 

% This part of code can be skipped when we know which the threshold to use
display('Adjusting the threshold before tracking')
%
tp=timer('TimerFcn','trackpreview','TaskstoExecute',1,'ErrorFcn','trackerror');
start(tp);
wait(tp);
% display('This is the new threshold:')
% display(num2str(thresholdbackground))

% Image processing to find the cropping box around the animal
im2=im2bw(firstbg,thresholdbackground);
im2=1-im2; % If animal is white in a black background, comment this line
im2=imclearborder(im2); % Delete arena parts at the field of view border
% Find the biggest object and keep it as the only object
[LO,num]=bwlabel(im2,8);
Areas=regionprops(LO,'area');
area_val=[Areas.Area];
maxarea=max(area_val);
idxBig= find(maxarea == area_val);
it2=ismember(LO,idxBig);
% cropping box
[r c]=find(it2); 
% close all
% figure,
% imshow(it2)
% display('Biggest object detected')
% pause
boxsize=50; 
maxc=max(c)+boxsize;
minc=min(c)-boxsize;
maxr=max(r)+boxsize;
minr=min(r)-boxsize;
if minr<0
    minr=0;
end
if minc<0
    minc=0;
end
% Important data
x=minc;
y=minr;
l=maxc-minc;
w=maxr-minr;
bbox=[minc minr maxc-minc maxr-minr];


% ----------------------------------------------------------------------
% Non-trivial background reconstruction 
% ----------------------------------------------------------------------
% It implies calling trackwait.m and trackbg.m functions
% The wait timer stops when the maggot is out of the initial bounding box.
% Then is moves to the trackbg function in order to reconstruct the new bg
  
iscrossing=0;
waittimer=timer('TimerFcn','trackwait','TaskstoExecute',1,'StopFcn','trackbg');
start(waittimer);
wait(waittimer);

% Finding where the animal after background reconstruction (since it moved)
mypreviewframe=getsnapshot(vidobj); 
removebgframe=(goodbackground)-(mypreviewframe); % invert when animal is white
threshold=0.15;
itest=im2bw(removebgframe,threshold);
[LO,num]=bwlabel(itest,8);
Areas=regionprops(LO,'area');
area_val=[Areas.Area];
maxarea=max(area_val);
idxBig= find(maxarea == area_val);
it2=ismember(LO,idxBig);
[r c]=find(it2);
maxc=max(c)+5;
minc=min(c)-5;
maxr=max(r)+5;
minr=min(r)-5;
% This is the updated bbox for trackbg to know where to look
bbox=[minc minr maxc-minc maxr-minr];

% ----------------------------------------------------------------------
% Writing down the first bounding box used to speed up trackgo.m
nextbboxref(1,:)=[bbox(1)-boxsize  bbox(2)-boxsize  bbox(3)+2*boxsize  bbox(4)+2*boxsize];

% Saving the first frame just before tracking starts
% Necessary for representation and calibration purposes offline. 
myframe=getsnapshot(vidobj); 
imwrite(myframe,'firstframe.BMP');
%iii=imcrop(myframe,nextbboxref(1,:));
%figure, imshow(iii)
%pause(2)

% ----------------------------------------------------------------------
% Animal tracking starts
% ----------------------------------------------------------------------

display('Tracking starts now')
% Main timer object 
tic
i=0;
t = timer('TimerFcn','trackgo','period',period,'ExecutionMode','fixedRate','BusyMode','error','TaskstoExecute',nframe,'ErrorFcn','trackerror');
start(t);
% wait(t);
% Stopping the experiment before the timer dictates it
display('Write "return" to stop tracking & save data')
keyboard;
stop(t)
display('End of tracking');

% ----------------------------------------------------------------------
% Save cropped animal images and bounding box positions as output data
save finalbbox finalbbox
save impro impro
save imraw imraw

display('End of the experiment');

% We can make sure the animal was properly tracked (or just comment this part)
cd(dirdest);
load imraw.mat
load impro.mat
close all
figure,
for ii=1:length(imraw)
close all
subplot(1,2,1); imshow(imraw{ii})
subplot(1,2,2); imshow(impro{ii})
pause(0.25)
end
%

cd(dirout);
display('Remember to close "Video Preview" window before starting a new experiment');

end
