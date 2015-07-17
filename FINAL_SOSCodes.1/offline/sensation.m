% -------------------------------------------------------------------------
% sensation.m 
%
% -------------------------------------------------------------------------
% Description:
%
% It reads the experimentally reconstructed sensory landscape and maps it
% to the motor data in order to obtain: the animal local orientation 
% with respect to the stimulus gradient field direction, 
% the sensory experience at several points along the animal (potentially
% relevant when related to the corresponding sensory organs) and its 
% temporal derivatives (normalized and un-normalized).
% -------------------------------------------------------------------------

close all
clear

% Directories for sensory 
dirOut='/Users/agomez/Desktop/';
%
dataDir='/Users/agomez/Desktop/supplementary/datademo/';
cd(dataDir)


% ------------------------------------------------------------------
% Sensory landscape
% ------------------------------------------------------------------
% Spatial distribution of the stimulus intensity
load landscape.mat
data=XYsensoryMap;
% Odor concentration units: micro Molar
% Mesh spacing: 0.1 milimeters (convenient to be similar than the tracking resolution)
scaleMesh=1/10; % mm/point 
% Vertical and Horizontal coordinates
dataT=data';
Ly=size(dataT,1);
for i=1:Ly
dataTnew(i,:)=dataT(Ly-i+1,:);
end
yrealNew=(1:size(dataTnew,1))*scaleMesh;
xrealNew=(1:size(dataTnew,2))*scaleMesh;

% Center reference position (in mm)
% (in this case it should coincide with odor droplet center)
stimMax=max(max(dataTnew));
[yMN,xMN]=find(dataTnew==stimMax); % in grid numbering
xC=xMN*scaleMesh;
yC=yMN*scaleMesh;

% Local gradient intensity and direction
[gx,gy]=gradient(dataTnew);
% Direction:
angles=atan2(gy,gx);
% Gradient intensity:
mgrad=[];
for kk1=1:size(gy,1);
for kk2=1:size(gy,2);
mgrad(kk1,kk2)=sqrt(gx(kk1,kk2)*gx(kk1,kk2) + gy(kk1,kk2)*gy(kk1,kk2))/scaleMesh; %uM/mm units
end
end


% ------------------------------------------------------------------------
% Loading input files
% ------------------------------------------------------------------------
load motorData.mat
% Sequence of single-animal experiments
sequence=1:length(motorData);

% ------------------------------------------------------------------------
% Initial empty variables
sensoryData={};

% ------------------------------------------------------------------------
% Main loop over every experiment to map sensory gradient to motor trajectory
% ------------------------------------------------------------------------
for run=sequence;
    
% Display current run process
run

% ------------------------------------------------------------------------
% Time resolution
fs=motorData{1}.fs; 
% GoodFrames for derivaties (see below)
goodFrames=motorData{run}.goodFrames; 
% Source position
xsource=motorData{run}.sourceXY(1);
ysource=motorData{run}.sourceXY(2);
% Trajectory in behavioral arena coordinates
xtraj=motorData{run}.cmXY(:,1);
ytraj=motorData{run}.cmXY(:,2);
xtrajH=motorData{run}.headXY(:,1);
ytrajH=motorData{run}.headXY(:,2);
xtrajT=motorData{run}.tailXY(:,1);
ytrajT=motorData{run}.tailXY(:,2);
xtrajM=motorData{run}.midXY(:,1);
ytrajM=motorData{run}.midXY(:,2);

% ------------------------------------------------------------------------
% Align droplet centers (sensory-motor landscape mapping) via translation
deltaX=xC-xsource;
deltaY=yC-ysource;
%
xsource=xsource+deltaX;
ysource=ysource+deltaY;
xtraj=xtraj+deltaX;
ytraj=ytraj+deltaY;
xtrajH=xtrajH+deltaX;
ytrajH=ytrajH+deltaY;
xtrajT=xtrajT+deltaX;
ytrajT=ytrajT+deltaY;
xtrajM=xtrajM+deltaX;
ytrajM=ytrajM+deltaY;

% ------------------------------------------------------------------------
% Initialize grid and sensory variables
xGrid=zeros(1,length(xtraj));
yGrid=zeros(1,length(xtraj));
xGridH=zeros(1,length(xtraj));
yGridH=zeros(1,length(xtraj));
xGridT=zeros(1,length(xtraj));
yGridT=zeros(1,length(xtraj));
xGridM=zeros(1,length(xtraj));
yGridM=zeros(1,length(xtraj));
CexpCM=nan(1,length(xtraj));
CexpH=nan(1,length(xtraj));
CexpT=nan(1,length(xtraj));
CexpM=nan(1,length(xtraj));
gradCM=nan(1,length(xtraj));
gradH=nan(1,length(xtraj));
gradT=nan(1,length(xtraj));
anglesM=nan(1,length(xtraj));

    for i=1:length(xtraj) 
    % Mapping of trajectory positions to the mesh grid
    xGridCM(i)=round(xtraj(i)/scaleMesh);
    yGridCM(i)=round(ytraj(i)/scaleMesh);
    xGridH(i)=round(xtrajH(i)/scaleMesh);
    yGridH(i)=round(ytrajH(i)/scaleMesh);
    xGridT(i)=round(xtrajT(i)/scaleMesh);
    yGridT(i)=round(ytrajT(i)/scaleMesh);
    xGridM(i)=round(xtrajM(i)/scaleMesh);
    yGridM(i)=round(ytrajM(i)/scaleMesh);
        try 
        % If the kinematic position lays beyond the sensory space,
        % its sensory associated value is ignored now, and it can be 
        % reasigned later by extending the gradient reconstruction.
        CexpCM(i)=dataTnew(yGridCM(i),xGridCM(i)); 
        CexpH(i)=dataTnew(yGridH(i),xGridH(i)); 
        CexpT(i)=dataTnew(yGridT(i),xGridT(i));
        CexpM(i)=dataTnew(yGridM(i),xGridM(i));
        anglesM(i)=angles(yGridM(i),xGridM(i)); 
        gradH(i)=mgrad(yGridH(i),xGridH(i)); 
        gradCM(i)=mgrad(yGridCM(i),xGridCM(i));
        gradT(i)=mgrad(yGridT(i),xGridT(i));
        gradM(i)=mgrad(yGridM(i),xGridM(i));
        catch, end
    end

    % ---------------------------------------------------------------------
    % Calculating sensory experience time derivatives
    % Initialize: 
    Y3hReal=zeros(1,length(xtraj));
    Y3cmReal=zeros(1,length(xtraj));
    Y3tReal=zeros(1,length(xtraj));
    Y3mReal=zeros(1,length(xtraj));
    Y4hReal=zeros(1,length(xtraj));
    Y4cmReal=zeros(1,length(xtraj));
    Y4tReal=zeros(1,length(xtraj));
    Y4mReal=zeros(1,length(xtraj));
    for k=2:length(xtraj)-1;  
    % Only define the derivative when there are three good frames
    if length(intersect(k-1:k+1,goodFrames))==3
       % stimulus derivative
       Y3hReal(k)=(fs*(CexpH(k+1)-CexpH(k-1))/2);
       Y3cmReal(k)=(fs*(CexpCM(k+1)-CexpCM(k-1))/2);
       Y3tReal(k)=(fs*(CexpT(k+1)-CexpT(k-1))/2);
       Y3mReal(k)=(fs*(CexpM(k+1)-CexpM(k-1))/2);
       % normalized stimulus derivative
       Y4hReal(k)=(fs*(CexpH(k+1)-CexpH(k-1))/2)/CexpH(k);
       Y4cmReal(k)=(fs*(CexpCM(k+1)-CexpCM(k-1))/2)/CexpCM(k); 
       Y4tReal(k)=(fs*(CexpT(k+1)-CexpT(k-1))/2)/CexpT(k); 
       Y4mReal(k)=(fs*(CexpM(k+1)-CexpM(k-1))/2)/CexpM(k); 
    end
    end

% ------------------------------------------------------------------------
% Write down the mapped sensory trajectories in time
%
% Stimulus intensity
sensoryData{run}.headSens=CexpH;
sensoryData{run}.cmSens=CexpCM;
sensoryData{run}.tailSens=CexpT;
sensoryData{run}.midSens=CexpM;
% Stimulus derivative
sensoryData{run}.headSensDot=Y3hReal;
sensoryData{run}.cmSensDot=Y3cmReal;
sensoryData{run}.tailSensDot=Y3tReal;
sensoryData{run}.midSensDot=Y3mReal;
% Stimulus normalized derivative (log stimulus)
sensoryData{run}.headSensDotNorm=Y4hReal;
sensoryData{run}.cmSensDotNorm=Y4cmReal;
sensoryData{run}.tailSensDotNorm=Y4tReal;
sensoryData{run}.midSensDotNorm=Y4mReal;
% Stimulus gradient
sensoryData{run}.headGrad=gradH;
sensoryData{run}.cmGrad=gradCM;
sensoryData{run}.tailGrad=gradT;
sensoryData{run}.midGrad=gradM;
% Local stimulus direction
sensoryData{run}.stimDir=anglesM;
% Stimulus direction relative to animal orientation 
diff=(sensoryData{run}.stimDir-motorData{run}.bodyTheta);
bearing=atan2(sin(diff),cos(diff)); 
sensoryData{run}.bearing=bearing; % notice that radians are its natural units


end

% -------------------------------------------------------------------
% End of the main loop and saving sensory mapped data
save sensoryData.mat sensoryData

% -------------------------------------------------------------------
% One could add, analogously to loci.m, a live display of the sensory-motor 
% variable dynamics for visualiation purposes.


