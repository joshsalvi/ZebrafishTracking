% -------------------------------------------------------------------
% locib.m code
%
% Description:
% This code reads a sequence of frames extracted from a behavioral movie,
% removes the background to detect the animal with proper thresholding,
% and then extracts detailed postural information such as its contour,
% its head and tail points, etc. The key is to use the curvature of the
% contour (assuming that the animal shape is properly resolved and smooth)
% to detect key points along the animal's body. Then, in further offline 
% analysis software, one can build kinematic variables (angles, speeds)
% to start establishing quantitative correlates.
%
% Input:
%
% The input data is a sequence of raw frames of the animal behaving in the
% arena captured previously. 
%
% Output:
% (a) The output of the code is a sequence of processed jpeg's : basically the
% original raw frame with the animal's head and tail positions detected
% and, little boxes illustrating the image processing. 
% (b) More importantly, the output of the code is the head, tail, centroid and
% midpoint x-y positions over time. That is the central data one needs to
% progress with the behavioral analysis. Skeletons, areas, etc. could be
% saved as well.
%
% -------------------------------------------------------------------

% -------------------------------------------------------------------

% Start with nothing (an empty workspace)
clear
close all

% Folders (where to get the input data and where to write the output)
dataDir0='/Users/agomez/Desktop/animalTrack/';
dirOut2='/Users/agomez/Desktop/animalTrack/mouse/framesPro/';
mkdir(dirOut2)
dataFolder='mouse/';
dataDir=strcat(dataDir0,dataFolder,'frames/');
outDir=dataDir;
filesArray=dir(strcat(dataDir,'frameNumber*.bmp'));
maxFile=length(filesArray);
cd(dataDir)

% load the background image for proper subtraction
% (there are several possibilities and tricks)
ibg=imread('0bg.bmp');
ibg=ibg(:,:,1);

% Initialize some stuff
fish={};
contour1={};
contour2={};
zcm=zeros(3,2); 
fileCount=0;
xhtime=[];
yhtime=[];
xctime=[];
yctime=[];
xttime=[];
yttime=[];
xmtime=[];
ymtime=[];

% flags:
flag=1; % set to 1 if you want to plot the jpeg sequence
% set to any other value otherwise

% -------------------------------------------------------------------
% The image processing party starts here!

% Temporal resolution (frames per second)
fps=30; 

% In case you want to start processing the n-th frame and end before.
iniFile=1;
finalFile=maxFile; % the last frame in my list

% In case you want to skip a few frames (coarse graining in time)
cgTime=1;

for file=iniFile:cgTime:finalFile
    
        % A useful counter
        fileCount=fileCount+1;

        % read each frame from the folder
        filename=filesArray(file).name;
        i0=imread(filename);
        i0=i0(:,:,1);
       
        % if there is a bit of the image that bothers you, you can cut it
        % (remember to do so in the background image too)
        %i0=i0(:,82:634); 
        
        % Original frame (will be used later to plot stuff)
        i00=i0;
        
        % Background subtraction (a simple way to detect the animal)
        i1=ibg-i0;
        
        % Renaming it:
        i0=i1;
      
        % Thresholding: the grey images is converted into black-and-white
        i2=im2bw(i1,0.15); 
        % The value has to be adjusted depending on the illumination, etc.
        % (Just play a little with it until it works)
        % I tuned it so that mouse-in-a-cage movie processing would work OK
        
        % Sometimes it is necessart to invert the image (black to white)        
        %i2=1-i2;
       
        % Renaming it:
        i4=i2;
        
        % Sometimes it is useful to play with these Matlab image processing
        % operations (type "help bwareaopen" for a proper explanation or/and ask me)
        % i4=bwareaopen(i2,10);
        % se = strel('diamond',1);
        % i5= imdilate(i4,se);
        % i5=imfill(i5,'holes');
        
        % Renaming it:
        i5=i4;
        
        % Sometimes the boder of the image is in the way and we can just
        % kill it
        %i5=imclearborder(1-i5);
        
        % Renaming it 
        i4=i5;
        % (note: all these redundant variables are not very optimal in terms of speed)
        

%         % When something fails, it is useful to stop the code and display what we
%         % got. This is a simple example to see how does the black-white detected
%         % animal looks like before we proceed to find its morphological properties.
%                          close all
%                          figure, imshow(i5)
%                          pause
%         % Advice: play with the image processing above and check how i5
%         % changes

        
% -------------------------------------------------------------------
% Having a nice animal shape, let's find out more about it

% Matlab gives you for free the properties of the object
i=i4;
[L,numobj] = bwlabel(i,8);
regprops=regionprops(L,'Centroid','Area');%,'Orientation','MajorAxisLength','MinorAxisLength','Perimeter');

% Find the areas of all the remaining objects in the processed image
bigobject1=find([regprops.Area]==max([regprops.Area]));
Arees0=[regprops.Area];

% I could use this as a condition to process objects of only a certain size
if Arees0(bigobject1)>0
     
% Another trick to avoid tracking dirts, curtains or artifacts: just focus
% on the largest object (biggest area), assuming it is the animal
bigobject=bigobject1(1); % in case there is more than one
% this is reminiscent from an earlier version of the code which allowed to
% track more than one animal simultaneously

% The first and most basic coordinate: its centroid position
centr=regprops(bigobject).Centroid;
xc=centr(1);
yc=centr(2);

% Now just focus on such big object and disregard the rest
[rMain, cMain]=find(bwlabel(i)==bigobject);
iL=zeros(size(i,1),size(i,2));
for k=1:length(rMain)
iL(rMain(k),cMain(k))=1;
end

% I could make sure the animal has no little holes (to avoid problems)
iL=imfill(iL,'holes'); 
% I could repeat image processing operations again if I want (smoothing, etc.)

% Next, we find its perimeter (not ordered from heat to tail)
iperim=bwperim(iL);
%
% I find those points that belong to the perimeter and write down their
% coordinates
centerIndp=find(iperim>0);
[py,px]=ind2sub(size(iperim),centerIndp);
%
% And finally, the ordered contour by walking thru the perimeter
contour=bwtraceboundary(iperim, [py(1) px(1)], 'N');
% (N=North seems important to avoid diagonals that may be dangerours)
% When that does not work, simply try other directions


% -------------------------------------------------------------------
% Now, calculate the curvature along the contour

% The X and Y coordinates of the contour are ordered to run
% the function that calculates curvature along perimeter
outlineCoords=[];
outlineCoords(:,1)=contour(:,1);
outlineCoords(:,2)=contour(:,2);
% Heuristic distance (in pxls!) between points to calculate the angles 
distCurv=10;

% Next, we calculate the angle between such points along the contour
dataLength=length(outlineCoords);
curVec=nan(1,dataLength);
xOutline=(outlineCoords(:,1));
yOutline=(outlineCoords(:,2));
%
for i=1:dataLength
index1=i;
index2=index1+distCurv;
index3=index1-distCurv;
%
length1=length(outlineCoords);
if index2<=length1
outlineCoords2 = outlineCoords(index2,:);
elseif index2>length1
indexMove = length1-index1;
indexNo = distCurv-indexMove;
outlineCoords2 = outlineCoords(indexNo,:);
end
%
if index3>=1
outlineCoords3 = outlineCoords(index3,:);
elseif index3<1
indexNo = length1+index3;
outlineCoords3 = outlineCoords(indexNo,:);
end
%
% this is what I need, for each point along
p1=outlineCoords(i,:);
p2=outlineCoords2;
p3=outlineCoords3;
angle1 = atan2((p2(1,1)-p1(1,1)),(p2(1,2) -p1(1,2)))-atan2((p3(1,1)-p1(1,1)),(p3(1,2) -p1(1,2)));             
curVec(i)=angle1;
%
end
% We can arrange the curvature a bit to plot it nicely
curVec=unwrap(curVec)-mean(unwrap(curVec));


% -------------------------------------------------------------------
% Two-endpoint detection (head and tail, in principle)

% Find the maximum positive value: should be the head or tail
idxMax0=find(curVec==max(curVec));
% Mid torsions of the animal have negative curvature and, thus, they are
% not picked (however, this part it easily customizable for other purposes)
%
% In case there is more than one maximum, just pick any
idxMax=idxMax0(1); 

% Now we find the other positive maximum, which shold be the tail
% We define an clone of the curVec function
curVec2=curVec;
% The idea is to set to zero everything around the first maximum, so that
% we can then easily find the second maximum
% (a smooth poly fit could be another option to find both peaks)
% Other improvements are feasible here.
curVec2(idxMax0)=[0];
length1=length(outlineCoords);
curVec2(length1)=[0];
% Now, we have to be careful with periodic conditions while running thru
% the curvature of the perimeter
index1=idxMax;
distCurv=length(curVec)/3; % one third of the length works
length1=length(outlineCoords);
% some points forwards
index2=index1+distCurv;
if index2<=length1
curVec2(index1:index2)=[0];
elseif index2>length1
indexMove=length1-index1;
indexNo=distCurv-indexMove;
curVec2(index1:length1)=[0];
curVec2(1:indexNo)=[0];           
end
% some points backwards
index3 = index1-distCurv;
if index3>=1
curVec2(index3:index1)=[0];
elseif index3<1
indexNo=length1+index3;
curVec2(1:index1)=[0];
curVec2(indexNo:length1)=[0];
end
%
% Finally, we find the second maximum
idxMax2=find(curVec2==max(curVec2));
% and pick any, in case there are two or more values together
idxMax2=idxMax2(1);   


% -------------------------------------------------------------------
% Building the skeleton from the endpoints (instead of the opposite) and
% then calculating the mid point (as a complementary position to the centroid)

% Now we are going to run "at different speeds" thru the perimeter,
% going from head to tail. By calculating the mid point we should draw a
% pretty cool skeleton
%
% Note 1: this option should be more efficient than "thinnig" the whole
% image. We will also avoid spurs, which can be painful.
% Note 2: in this way, even if the skel is shit, we can obtain the head 
% and tail positions.
%
zhead=[outlineCoords(idxMax,1) outlineCoords(idxMax,2)];
contourN=[];
contourS=[];

% Now it is important to be sure we "walk well":
% starting from the head, we take one direction to walk on the perimeter
try
contourN=bwtraceboundary(iperim, zhead, 'W');
catch, end
try
contourN=bwtraceboundary(iperim, zhead, 'E');
catch, end
try
contourN=bwtraceboundary(iperim, zhead, 'S');
catch, end
try
contourN=bwtraceboundary(iperim, zhead, 'N');
catch, end

% By reversing the direction of walking, we find the reversed way of going
% from head to tail:
contourS(1,:)=contourN(1,:);
for kk=2:length(contourN); 
contourS(kk,:)=contourN(length(contourN)-kk+2,:);
end
%
ztail=[outlineCoords(idxMax2,1) outlineCoords(idxMax2,2)];
% We now find how much I have to walk from each contour to get to the tail
lengthN=find(contourN(:,1)==ztail(1) & contourN(:,2)==ztail(2));
lengthS=find(contourS(:,1)==ztail(1) & contourS(:,2)==ztail(2));
vecS=[];
vecN=[];
zmid=[];
% if contourN happens to have the longest path to the tail, I will
% make a virtual contourS as large as contourN
if lengthN>lengthS
vecS=zeros(lengthN,2);
for kk=1:lengthN
vecS(kk,:)=contourS(round(lengthS*kk/lengthN+1),:);
vecN(kk,:)=contourN(kk,:);
%The skeleton is writen down
zmid(kk,:)=[(vecN(kk,2)+vecS(kk,2))/2 (vecN(kk,1)+vecS(kk,1))/2];
end
% I do the converse if contourN is shorter than contourS
else
vecN=zeros(lengthS,2);
for kk=1:lengthS
vecN(kk,:)=contourN(round(lengthN*kk/lengthS+1),:);
vecS(kk,:)=contourS(kk,:);
% Finally, the skeleton!
zmid(kk,:)=[(vecN(kk,2)+vecS(kk,2))/2 (vecN(kk,1)+vecS(kk,1))/2];
end
end

% And from the skeleton, we find the mid point coordinate
xm=zmid(round(length(zmid)/2),1);
ym=zmid(round(length(zmid)/2),2);


% -------------------------------------------------------------------
% Making sure the head and tail identities are correct 
% (based on the previous head-tail endpoint detection)

% Current guesses:
xt=ztail(2);
yt=ztail(1);
xh=zhead(2);
yh=zhead(1);

% We use a distance rule. If the animal has a particular geometry that
% univocally charaterizes its head or tail, we could use that instead.
 
% Here we (hopefully) take the max peak of the curvature from the first
% frame as a correct detection of the head.
% Otherwise, we could set it manually for the user to click on the 1st frame
if file==iniFile
xt=ztail(2);
yt=ztail(1);
xh=zhead(2);
yh=zhead(1);    
else
% and then we systematically apply "the proximity rule"
xt=ztail(2);
yt=ztail(1);
xh=zhead(2);
yh=zhead(1);
%
distToTn=sqrt((xt-oxt)^2 +  (yt-oyt)^2);
distToHn=sqrt((xh-oxt)^2 +  (yh-oyt)^2);
% find the distance and swap head and tail if necessary:           
if distToTn>distToHn
xh=ztail(2);
yh=ztail(1);
xt=zhead(2);
yt=zhead(1);
end
end
% And save current values for subsequent comparison in next frame
oxh=xh;
oyh=yh;
oxt=xt;
oyt=yt;


% -------------------------------------------------------------------
% Writing down the loci trajectories (and other stuff I may need)

xhtime=[xhtime xh];
yhtime=[yhtime yh];
xctime=[xctime xc];
yctime=[yctime yc];
xttime=[xttime xt];
yttime=[yttime yt];
xmtime=[xmtime xm];
ymtime=[ymtime ym];


% -------------------------------------------------------------------
% Plotting "biutiful" things:

if flag==1
    
close all
figure, 

subplot(3,4,[1:3 5:7 9:11]), 
imshow(i00)
axis equal 
xc=centr(1);
yc=centr(2);
hold on
plot(xh,yh,'*m')
plot(xc,yc,'*w')
plot(xt,yt,'*g')

% Some box coordinates
dz2=65;
% The biggest
centr=regprops(bigobject1).Centroid;
xc=centr(1);
yc=centr(2);

subplot(3,4,4), 
imshow(i00)
axis equal 
axis([xc-dz2 xc+dz2 yc-dz2 yc+dz2])
axis off

subplot(3,4,8), 
imshow(i5)
axis equal 
axis([xc-dz2 xc+dz2 yc-dz2 yc+dz2])
axis off

subplot(3,4,12), 
imshow(i0)
axis equal 
xc=centr(1);
yc=centr(2);
hold on
plot(contour(:,2), contour(:,1), 'linewidth', 1, 'color', 'b');
plot(xh,yh,'*m')
plot(xt,yt,'*g')
plot(xc,yc,'.w')
plot(xm,ym,'.y')
[b,a]=find(iL==1);
axis([xc-dz2 xc+dz2 yc-dz2 yc+dz2])
axis off

% And save the image for every processed frame
outfile00=strcat(dirOut2,'proFrame',num2str(fileCount),'.jpeg');
print('-f1','-djpeg',outfile00);  

display('Paused: press any key to continue')
pause

end % end of the object/area constraints   
end % end of the flag for plotting
  
% -------------------------------------------------------------------
% -------------------------------------------------------------------

end % end of the main loop


% -------------------------------------------------------------------
% Saving the positions of certain loci of interest (head, centroid and tail coordinates)
save xhtime.mat xhtime
save yhtime.mat yhtime
save xctime.mat xctime
save yctime.mat yctime
save xttime.mat xttime
save yttime.mat yttime
save xmtime.mat xmtime
save ymtime.mat ymtime
% It could all be easily saved within a unique Matlab array 
