% -------------------------------------------------------------------------
% loci.m 
%
% -------------------------------------------------------------------------
% Description:
%
% This code inputs grayscale & black-and-white animal cropped shape images
% (together with the bounding box coordinates of the cropped images) 
% from the online tracking experiments to output: the animal boundary, 
% skeleton and head, tail, centroid and midpoint positions in the arena
% (whole field of view) coordinate system. 
% The background image (firstframe) is used to callibrate the space
% and also for further sensory gradient mapping. 
% Other measures such as animal area or skeleton lengths are extracted 
% and saved as well.
% The contour and skeleton positions of the animal are saved too.
% At the end of the processing, an image of the arena and the time 
% trajectory of several points of interest is printed on the raw image.
% Finally, an animation of the whole process is displayed to identify 
% and correct possible mistakes.
% -------------------------------------------------------------------------

clear
close all

% ------------------------------------------------------------------------
% Directories
outDir='/Volumes/Promise Pegasus/Manual Backup/Lab/Videos/Zebrafish/High Speed/';
%dataDir=strcat(outDir,'onlineData/');
%cd(dataDir);
cd(outDir)
files = dir('*movieonly.mat');
fn = files(1).name;

% ------------------------------------------------------------------------
% Time resolution
fs=7; % frames per second
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Enabling manual correction mode:
correctFlag=1; % Stopping at problematic frames
%correctFlag=0; % Not going through problematic frames
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Spatial Callibration and Landmarks
% ------------------------------------------------------------------------
% Load arena image
load(fn);
bg = vidobj(1).cdata;
%bg=imread('firstframe.BMP');

% Manual click to determine spatial scale and stimulus reference positions:

figure, imshow(bg)
title('Click first at Center of the odor droplet, then 3 wells West, and finally 3 wells North')
% In the sample data, the center of the droplet was at the South-East well
% from the middle little dark mark of the 96-well plate
set(gca,'xtick',[],'ytick',[])
[X0,Y0]=ginput(1);
[XW,YW]=ginput(1);
[XN,YN]=ginput(1);
close all
figure, imshow(bg)
hold on
plot(X0,Y0,'*b')
plot(XW,YW,'*r')
plot(XN,YN,'*g')
pause(2)
display('If landmarks are incorrect, restart the program.')
display('If landmarks are correct, confirm placement and continue by typing "return".')
input('');
close all

D1=((X0-XW)^2+(Y0-YW)^2).^0.5;
D2=((X0-XN)^2+(Y0-YN)^2).^0.5;
Dm=(D1+D2)/2;
dw=3; % 3 wells distance
lw=9; % size of a well in mm
D1well=Dm/dw; 
D1mm=D1well/lw; % pixels/mm scale

% Saving temporal and spatial callibration
cd(outDir);
reference=[fs D1mm X0 Y0 XW YW XN YN];
save reference reference
cd(dataDir);
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Loading input files
% ------------------------------------------------------------------------
load impro.mat
load imraw.mat
load finalbbox

% ------------------------------------------------------------------------
% Initialization of variables and parameters
% ------------------------------------------------------------------------
% Initial and final snapshots to process
iniFile=1;
maxFile=size(impro,2);
finalFile=maxFile;
diffFiles=1+finalFile-iniFile;
% Initial arbitrary reference for head-tail classification
oxrh=0;
oyrh=0;
oxrt=0;
oyrt=0;
% Flags
UturnFlag=0;
headTailFlag=0;
countNot=0;
% Handy avatar variables
coordinates=zeros(diffFiles,8);
properties=zeros(diffFiles,7);
contour={};
skel={};
imskel={};
impoints={};
imcontour={};
bboxRef=zeros(diffFiles,4);
centAbs=zeros(diffFiles,2);
wskel=zeros(600,800);


% ------------------------------------------------------------------------
% Loop over black-and-white larval images to extract points of interest
% ------------------------------------------------------------------------
for file=iniFile:finalFile
    
   % Display fraction of total files processed
   fraction=(file-iniFile)/(finalFile-iniFile)
      
   % Try-catch loop to process all frames skipping skel spurs (see below)
   % try  
                
        % Read cropped images array 
        smallImage=impro{file};
        smallImage0=imraw{file};

        % Read the bounding box position for each image
        bbx=finalbbox(file,1);
        bby=finalbbox(file,2);
        
        % Skeleton reconstruction via extra image processing
        % Smoothing 
        h=fspecial('gaussian', 5, 1);
        w=imfilter(smallImage,h,'replicate');
        w=imfill(w,'holes');
        % Clear border trick
        w=imclearborder(w);
        % Delete small objects
        w=bwareaopen(w,10);
        % Skeleton
        % Second-the-last skel for head-tail swap flag display (see below)
        wskelOld=wskel;
        % Actual animal skeleton by simply thinning to infinity
        wskel=bwmorph(w,'thin',Inf);
        imskel{file}=wskel;
        
        % Morphological Information
        statsC=regionprops(double(w),'Centroid','Area','Perimeter','Orientation','MajorAxisLength','MinorAxisLength');
        % Geometric center of mass
        centr=statsC.Centroid;
        xrc=centr(1);
        yrc=centr(2);
        % Useful scalars
        are=statsC.Area;
        peri=statsC.Perimeter;
        ori=statsC.Orientation;
        maxis=statsC.MajorAxisLength;
        minxis=statsC.MinorAxisLength;
        % Handy proxy for skeleton length (no sqrt(2) considerations here)
        idxSkel=find(wskel==1);
        lengthSkel1=length(idxSkel);
        % Writing the data
        properties(file,:)=[0 lengthSkel1 are peri ori maxis minxis];
        % the first column will be filled later
        
        % Skeleton x-y positions
        skeleton=wskel;
        skeleton_indices=find(skeleton>0);
        [skeleton_y,skeleton_x]=ind2sub(size(skeleton),skeleton_indices);
        % Skeleton end points
        connect=zeros(1,length(skeleton_indices));
        for point=1:length(skeleton_indices)
            i=skeleton_x(point);
            j=skeleton_y(point);
            connect(point)=length(find(skeleton_x>=i-1 & skeleton_x<=i+1 & skeleton_y>=j-1 & skeleton_y<=j+1));
        end
        endPoints=find(connect==2);
        
        % Only continue with skeleton ordering if there are 2 endpoints
        % plus extra flag if the skeleton is too short (U-turn shape)
        Lthreshold=10; % in pixels
        aspectRatio=maxis/minxis;
        if (length(endPoints)==2 & aspectRatio>1.5) %& lengthSkel1>Lthreshold
            
            % Ordering the skeleton
            backbone_indices=endPoints(1);
            countwhile=0;
            while (length(backbone_indices)<length(skeleton_indices) & countwhile<1000);
                countwhile=countwhile+1;
                i=skeleton_x(backbone_indices(end));
                j=skeleton_y(backbone_indices(end));
                next_index=setdiff(find(skeleton_x>=i-1 & skeleton_x<=i+1 & skeleton_y>=j-1 & skeleton_y<=j+1),backbone_indices);
                backbone_indices=cat(2,backbone_indices,next_index);
            end
            backbone=[skeleton_y(backbone_indices),skeleton_x(backbone_indices)];
            
            % Proper skeleton length
            y=backbone(:,1); y=reshape(y,1,length(y));
            x=backbone(:,2); x=reshape(x,1,length(x));
            points=length(backbone);
            s=zeros(1,points);
            for i=2:points
                s(i)=s(i-1)+ sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
            end
            properties(file,1)=s(end);
          
            % Skeleton midpoint (as an alternative measure to centroid)
            [val,idx]=min(abs(s-s(end)/2));
            skeletonMid=backbone(idx,:);
            xrm=skeletonMid(2);
            yrm=skeletonMid(1);
          
            % Head-tail classification
            % Initial guess (50% chance) in relative coordinate system
            xrh=skeleton_x(endPoints(1));
            yrh=skeleton_y(endPoints(1));
            xrt=skeleton_x(endPoints(2));
            yrt=skeleton_y(endPoints(2));
            
            % Stopping to manually annotate head-tail at the initial frame
            % Adding problematinc frame swaps (last frame had problem)
            if ((file==iniFile) || (headTailFlag==0))  
            
            beep
            display('Head-tail manual annotation:')
            display('First, click close to head. Second, close to tail')
            
            display(['Number of endpoints = ' num2str(length(endPoints))])
            display(['Aspect ratio = ' num2str(aspectRatio)])
            
            % This part requires experimenting with the characteristics of the model organism.
            % If, after processing, all head-tail trajectory is swapped, 
            % one simply needs to clik on the other endpoint during the 
            % 1st frame manual annotation to get the classification right.
            %
            pause(1)
            %keyboard
            close all
            % Show current image frame to classify
            figure,
            imagesc(smallImage0)
            axis equal
            colormap('gray')
            title(['File number:',int2str(file)],'Color','b')
            colormap('gray')
            axis equal
            set(gca,'xtick',[],'ytick',[])
            [xTail,yTail]=ginput(1);
            oxrh=xTail;
            oyrh=yTail;  
            [xTail,yTail]=ginput(1);
            oxrt=xTail;
            oyrt=yTail;  
            close all
            %
            end

            % Distance rule to decide new head and tail
            distToTn=sqrt((xrt-oxrt)^2 + (yrt-oyrt)^2);
            distToHn=sqrt((xrh-oxrt)^2 + (yrh-oyrt)^2);
            endPointHead=endPoints(1);
            if distToTn>distToHn
                xrh=skeleton_x(endPoints(2));
                yrh=skeleton_y(endPoints(2));
                xrt=skeleton_x(endPoints(1));
                yrt=skeleton_y(endPoints(1));
                endPointHead=endPoints(2);
            end
            
            % We write down for next frame comparison
            % (notice relative spatial coordinates, so we assume fast
            % tracking)
            oxrh=xrh;
            oyrh=yrh;
            oxrt=xrt;
            oyrt=yrt;
            headTailFlag=1;
                        
            % Write head, tail centroid and midpoint points in the relative
            % system of reference
            impoints{file}=[xrc yrc xrh yrh xrt yrt xrm yrm];

            % Write centroid, head, tail and midpoint x-y positions in the
            % arena system or refence
            coordinates(file,:)=[xrc+bbx yrc+bby xrh+bbx yrh+bby xrt+bbx  yrt+bby xrm+bbx yrm+bby];
       
            % Write contour for potential further use
            wperim=bwperim(smallImage);
            imcontour{file}=wperim;
            %
            centerIndp=find(wperim>0);
            [py,px]=ind2sub(size(wperim),centerIndp);
            contour{file}.x=[px' + bbx];
            contour{file}.y=[py' + bby];
          
            % Reorder skeleton again but now from head-to-tail  
            backbone_indices=endPointHead;  
            countwhile=0;
            while (length(backbone_indices)<length(skeleton_indices) & countwhile<1000);
                countwhile=countwhile+1;
                i=skeleton_x(backbone_indices(end));
                j=skeleton_y(backbone_indices(end));
                next_index=setdiff(find(skeleton_x>=i-1 & skeleton_x<=i+1 & skeleton_y>=j-1 & skeleton_y<=j+1),backbone_indices);
                backbone_indices=cat(2,backbone_indices,next_index);
            end
           %backbone=[skeleton_y(backbone_indices),skeleton_x(backbone_indices)];
           backboneAbs=[skeleton_x(backbone_indices)+bbx,skeleton_y(backbone_indices)+bby];
           skel{file}.x=backboneAbs(:,1);
           skel{file}.y=backboneAbs(:,2);
         
        % If the skeleton had spurs or was too small we simply skip it and
        % set the headtail flag to false
        
        else
            
        headTailFlag=0; 
        if correctFlag==0; % not going thru problematic frames
            headTailFlag=1;
        end            
        
        %
        skel{file}.x=0;
        skel{file}.y=0;

        countNot=countNot+1;
        
        end
    
   % End of the try-catch  
   % catch %, ME
   % headTailFlag=0;
   % end 

end

% -------------------------------------------------------------------------
% Save files to disk
cd(outDir);
save coordinates coordinates
save properties properties
save skel skel
save contour contour


% Fraction of skipped head-tail skeletons:
framesNotPro=countNot;
allFramesPro=finalFile-iniFile;
percentatge=framesNotPro/allFramesPro*100;
display(['Skipped frames: ' num2str(percentatge) '%  (' num2str(framesNotPro) ' out of ' num2str(allFramesPro) ')'])


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Create multi-point animation with image processing display
close all
figure,
tspan=7; % to regulate the speed of the animation
for i=1:tspan:length(coordinates)% or maxFile

subplot(4,4,1)
imshow(imraw{i})
subplot(4,4,2)
imshow(impro{i})
axis off
subplot(4,4,3)
imshow(imcontour{i})
subplot(4,4,4)
imshow(imskel{i})
%
subplot(4,4,[5 16])
imagesc(bg)
colormap(gray)
title(['Frame:',int2str(i)],'Color','b')
axis equal
axis off
hold on
%
plot(skel{i}.x,skel{i}.y,'-k','lineWidth',3)
%
delayT=45; % in seconds
Di=i-delayT*fs;
if i<delayT*fs
Di=1;
end
plot(coordinates(Di:i,5),coordinates(Di:i,6),'.r','markersize',2)
plot(coordinates(Di:i,7),coordinates(Di:i,8),'.y','markersize',2)
plot(coordinates(Di:i,3),coordinates(Di:i,4),'.g','markersize',2)
%
pause(0.025)
hold off
end

close all

% -------------------------------------------------------------------------
% Create a multi-point trajectory on the raw arena image
subplot(1,1,1), 
imagesc(bg)
colormap(gray)
axis equal
axis off
hold on
plot(coordinates(:,1),coordinates(:,2),'+b','markersize',1,'linewidth',1) % centroid
plot(coordinates(:,5),coordinates(:,6),'+r','markersize',1,'linewidth',1) % tail
plot(coordinates(:,7),coordinates(:,8),'+y','markersize',1,'linewidth',1) % midpoint
plot(coordinates(:,3),coordinates(:,4),'+g','markersize',1,'linewidth',1) % head
outfile00=strcat(outDir,'trajectory.jpeg');
print('-f1','-djpeg',outfile00);



