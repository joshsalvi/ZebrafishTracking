% -------------------------------------------------------------------------
% motion.m 
%
% -------------------------------------------------------------------------
% Description:
%
% It reads, for every animal, the output tracking data and builds 
% trajectories (in physical length units), animal posture angles (bending, 
% orienation, etc.), and other variables such as speeds. 
% It also segments continuous kinematic trajectories into discrete
% behavioral modes such as "runs", "turns" and "casts".
% The whole motor information for all animals is now contained in a single
% matlab array.
% -------------------------------------------------------------------------

clear
close all

% Folders
dataDir='/Volumes/Promise Pegasus/Manual Backup/Lab/Videos/Zebrafish/High Speed/';
dirStruct=dir(strcat(dataDir,'all*'));
nDir=length(dirStruct);

% Initialization
larvaIdx=0;
motorData={};

for dirIdx=1:nDir
    
    dirName=dirStruct(dirIdx).name;
    % Positions, spatial and temporal callibration, and properties
    runStruct1=dir(strcat(dataDir,dirName,'/','coordinates*'));
    runStruct2=dir(strcat(dataDir,dirName,'/','reference*'));
    runStruct3=dir(strcat(dataDir,dirName,'/','properties*'));
    
    maxRun=length(runStruct1)
    for run=1:maxRun
        
        larvaIdx=larvaIdx+1;
        maxDist=zeros(larvaIdx,maxRun);

% ------------------------------------------------------------------------
% Continuous kinematic variables        
        
        % reading every file for every experiment
        fname=strcat(dataDir,dirName,'/',runStruct1(run).name);
        data_array=importdata(fname);
        fname2=strcat(dataDir,dirName,'/',runStruct2(run).name);
        data_array2=importdata(fname2);
        fname3=strcat(dataDir,dirName,'/',runStruct3(run).name);
        data_array3=importdata(fname3);
        
        % Error-free tracking frames
        goodFrames=find(data_array(:,1)~=0);
        %
        maxFrame=length(data_array(:,1));   
        badFrames=find(data_array(:,1)==0);
        
        % Space and time calibration
        fs=data_array2(1,1); 
        times=goodFrames/fs;
        scale=1/data_array2(1,2); %mm/pixel converson scale
        sourceXY=scale*[data_array2(1,3), data_array2(1,4)];
        % 
        % Saving data:
        motorData{larvaIdx}.fs=fs;
        motorData{larvaIdx}.scale=scale;
        motorData{larvaIdx}.sourceXY=sourceXY;
        motorData{larvaIdx}.goodFrames=goodFrames;
        
        % Morphological features
        motorData{larvaIdx}.lengthSkel=scale*data_array3(:,1);
        motorData{larvaIdx}.area=scale*scale*data_array3(:,3);
        motorData{larvaIdx}.perimeter=scale*data_array3(:,4);
        motorData{larvaIdx}.orientation0=data_array3(:,5);
        motorData{larvaIdx}.maxis=data_array3(:,6);
        motorData{larvaIdx}.minxis=data_array3(:,7);
        
%         % Raw (unfiltered) positions
%         motorData{larvaIdx}.cmXYRaw=scale*[data_array(:,1) data_array(:,2)];
%         motorData{larvaIdx}.headXYRaw=scale*[data_array(:,3) data_array(:,4)];
%         motorData{larvaIdx}.tailXYRaw=scale*[data_array(:,5) data_array(:,6)];
%         motorData{larvaIdx}.midXYRaw=scale*[data_array(:,7) data_array(:,8)];
         
        % Polynomial filter smoothing of positions
        xFilt=zeros(4,maxFrame);
        yFilt=zeros(4,maxFrame);
        nl=fs;
        nr=fs;
        pOrder=3;
        thetaFilt=zeros(1,maxFrame);
        for var=1:4
            dataIdx=var*2-1;
            for frame=nl+1:maxFrame-(nr+1)
                window=frame-nl:frame+nr;
                %the good frames window
                winIdx=setdiff(window,badFrames);
                if (length(find(winIdx>frame+1))>3 && length(find(winIdx<frame))>3)
                    px=polyfit(1:length(winIdx),squeeze(data_array(winIdx,dataIdx))',pOrder);
                    xFilt(var,frame) = polyval(px,nl+1);
                    py=polyfit(1:length(winIdx),squeeze(data_array(winIdx,dataIdx+1))',pOrder);
                    yFilt(var,frame)=polyval(py,nl+1);
                    if var==3
                      % using the filter to construct the centroid velocity direction
                       xdot=polyval(polyder(px),nl+1);
                       ydot=polyval(polyder(py),nl+1);
                       thetaFilt(frame)=atan2(ydot,xdot);
                    end
                end
            end
        end
        motorData{larvaIdx}.cmXY=scale*[xFilt(1,:)' yFilt(1,:)'];
        motorData{larvaIdx}.headXY=scale*[xFilt(2,:)',yFilt(2,:)'];
        motorData{larvaIdx}.tailXY=scale*[xFilt(3,:)',yFilt(3,:)'];
        motorData{larvaIdx}.midXY=scale*[xFilt(4,:)',yFilt(4,:)'];
        

        % Speed (in mm/s)
        % centroid 
        speed=zeros(1,maxFrame);
        for i=2:maxFrame-1;
            if length(intersect(i-1:i+1,goodFrames))==3
                dx2=(motorData{larvaIdx}.cmXY(i+1,1)-motorData{larvaIdx}.cmXY(i-1,1))^2/4;
                dy2=(motorData{larvaIdx}.cmXY(i+1,2)-motorData{larvaIdx}.cmXY(i-1,2))^2/4;     
                speed(i)=fs*sqrt(dx2+dy2);
            end
        end
        motorData{larvaIdx}.cmSpeed=speed;
        % head
        speed=zeros(1,maxFrame);
        for i=2:maxFrame-1;
            if length(intersect(i-1:i+1,goodFrames))==3
                dx2=(motorData{larvaIdx}.headXY(i+1,1)-motorData{larvaIdx}.headXY(i-1,1))^2/4;
                dy2=(motorData{larvaIdx}.headXY(i+1,2)-motorData{larvaIdx}.headXY(i-1,2))^2/4;     
                speed(i)=fs*sqrt(dx2+dy2);
            end
        end
        motorData{larvaIdx}.headSpeed=speed;
        % tail
        speed=zeros(1,maxFrame);
        for i=2:maxFrame-1;
            if length(intersect(i-1:i+1,goodFrames))==3
                dx2=(motorData{larvaIdx}.tailXY(i+1,1)-motorData{larvaIdx}.tailXY(i-1,1))^2/4;
                dy2=(motorData{larvaIdx}.tailXY(i+1,2)-motorData{larvaIdx}.tailXY(i-1,2))^2/4;     
                speed(i)=fs*sqrt(dx2+dy2);
            end
        end
        motorData{larvaIdx}.tailSpeed=speed;
        % midpoint
        speed=zeros(1,maxFrame);
        for i=2:maxFrame-1;
            if length(intersect(i-1:i+1,goodFrames))==3
                dx2=(motorData{larvaIdx}.midXY(i+1,1)-motorData{larvaIdx}.midXY(i-1,1))^2/4;
                dy2=(motorData{larvaIdx}.midXY(i+1,2)-motorData{larvaIdx}.midXY(i-1,2))^2/4;     
                speed(i)=fs*sqrt(dx2+dy2);
            end
        end
        motorData{larvaIdx}.midSpeed=speed;
        
        
        % Animal orientation
        % The tail angle is good proxy for heading, which is defined at every
        % time point and does not rely on
        % It is expressed in the lab arbitrary system or coordinates.
        dx=motorData{larvaIdx}.midXY(:,1)-motorData{larvaIdx}.tailXY(:,1);
        dy=motorData{larvaIdx}.midXY(:,2)-motorData{larvaIdx}.tailXY(:,2);
        theta=zeros(1,maxFrame);
        theta(goodFrames)=atan2(dy(goodFrames),dx(goodFrames));
        motorData{larvaIdx}.bodyTheta=theta;     
        
        
        % Rate of change in orienation
        bodyOmega=zeros(1,maxFrame);
        goodOmegaFrames=goodFrames;
        for i=2:maxFrame-1;
            %only define the derivative when there are three good frames in a row
            if length(intersect(i-1:i+1,goodFrames))~=3
                goodOmegaFrames=setdiff(goodOmegaFrames,i);
            else
                data=unwrap(motorData{larvaIdx}.bodyTheta(i-1:i+1));
                %the symmetric difference, in rad/s
                bodyOmega(i)=fs*(data(3)-data(1))/2;
            end
        end
        motorData{larvaIdx}.bodyOmega=bodyOmega;
       
      
        % Relative head bending angle.
        dx=motorData{larvaIdx}.headXY(:,1)-motorData{larvaIdx}.midXY(:,1);
        dy=motorData{larvaIdx}.headXY(:,2)-motorData{larvaIdx}.midXY(:,2);
        theta=zeros(1,maxFrame);
        theta(goodFrames)=atan2(dy(goodFrames),dx(goodFrames));
        aat=unwrap(theta)-unwrap(motorData{larvaIdx}.bodyTheta);
        motorData{larvaIdx}.headTheta=atan2(sin(aat),cos(aat));
        
        
        % Rate of change in bending
        omega=zeros(1,size(data_array,1));
        for i=2:size(data_array,1)-1;
            %only define the derivative when there are three good frames in a row
            if length(intersect(i-1:i+1,goodFrames))==3
                data=unwrap(motorData{larvaIdx}.headTheta(i-1:i+1));
                %the symmetric difference, in rad/s
                omega(i)=fs*(data(3)-data(1))/2;
            end
        end
        motorData{larvaIdx}.headOmega=omega;
                
  
% ------------------------------------------------------------------------
% Discrete behavioral modes        

% -----
% CASTS
% Definition: point intervals that correspond to large head angles
%
thresholdHead=37*pi/180; % 37 degrees body bending threshold 
aat=motorData{larvaIdx}.headTheta;
aat2=atan2(sin(aat),cos(aat));
idxBend=find(abs((aat2))>thresholdHead);  
[maxtab, mintab] = peakdet(aat2(goodFrames),0.2,goodFrames);  
relExtrem=[];
try
relExtrem=union(mintab(:,1),maxtab(:,1));
catch, end
idxBend0=intersect(idxBend,relExtrem);
%
% Frame time-stamp indexes of occurring head casts: 
motorData{larvaIdx}.idxCastEvent=idxBend0;


% -----
% TURNS
% Definition: point intervals that correspond to large body angle velocities
%
bOt=0.2; % lower reorientation speed threhold (in rad/s) to classify for turn event
idxSpin=find(abs(motorData{larvaIdx}.bodyOmega)>bOt);
idxSpin0=idxSpin;
% To cut those very large: that are errors
bOthigh=2.5; % upper threshold to prevent reorientation velocity artifacts
idxSpin=find(abs(motorData{larvaIdx}.bodyOmega)>bOt & abs(motorData{larvaIdx}.bodyOmega)<bOthigh);
% Finding the beginning and end each event
idxSpinEnd=[];
idxSpinStart=[];
zerone=zeros(1,length(motorData{larvaIdx}.bodyOmega));
zerone(idxSpin)=1;
countStart=0;
countEnd=0;
for k=1:length(idxSpin)
if (zerone(idxSpin(k))==1 & zerone(idxSpin(k)-1)==0)
countStart=countStart+1;
idxSpinStart(countStart)=idxSpin(k);
end
if (zerone(idxSpin(k))==1 & zerone(idxSpin(k)+1)==0)
countEnd=countEnd+1;
idxSpinEnd(countEnd)=idxSpin(k);
end
end
% Erase those idxSpinStart and idxSpinEnd comprising a short time interval
SpinDuration=idxSpinEnd-idxSpinStart;
lowBoundDuration=1*fs; %1 second threshold
for k=1:length(idxSpinStart)
if SpinDuration(k)<lowBoundDuration
    idxSpinStart(k)=[0];
    idxSpinEnd(k)=[0];
end
end
zerStart=find(idxSpinStart==0);
idxSpinStart(zerStart)=[];
zerEnd=find(idxSpinEnd==0);
idxSpinEnd(zerEnd)=[];
idxSpin=[];
for k=1:length(idxSpinStart)
idxSpin=[idxSpin idxSpinStart(k):idxSpinEnd(k)];
end
%
% Frame time-stamp indexes of occurring turns (start and finish):
motorData{larvaIdx}.idxTurnStart=idxSpinStart;
motorData{larvaIdx}.idxTurnEnd=idxSpinEnd;

% ------------------------------------------------------------------------
                

    end % End of the loop over runs
end % End of the loop over directories


% ------------------------------------------------------------------------
% Save all data and new variables, for every time point and experimental
% run, in an array for later use
%
%outDir=strcat(dataDir,'analysis/');
outDir=strcat(dataDir);
save(strcat(outDir,'motorData.mat'),'motorData')


