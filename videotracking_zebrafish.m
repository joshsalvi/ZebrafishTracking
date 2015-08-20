%% load data
clear;close all;
cd '/Volumes/Promise Pegasus/Manual Backup/Lab/Videos/Zebrafish/High Speed/20150728/';
files = dir('*.avi');

for m = 1:length(files)
    
    disp(['File #' num2str(m) '/' num2str(length(files))]);
    
    warning off
    
    comments{m} = [files(m).name(1:end-4)];
    vid{m}=aviread(files(m).name);
    
end

%% choose a threshold
close all
thresh = 8;
newcropyn = 1;
for movnum = 1:length(files)
    subplot_tight(5,round(length(files)/5),movnum);
    I = adapthisteq(vid{movnum}(1).cdata);
    if newcropyn==1
    if movnum==1
        [I, rect] = imcrop(I);
    else
        I = imcrop(I,rect);
    end
    else
        I = imcrop(I,rect);
    end
    %I = imcrop(I,[10*round(size(I,1)/100) round(10*size(I,1)/100) round(90*size(I,1)/100) round(70*size(I,2)/100)]);
    imagesc(I<mean(mean(I))-thresh*std(std(double(I))));colormap gray; axis off
end
%% analyze data (1)

thresh = 8;
Fs = 500;         % Hz
pulsestart = 100; % ms
duration = 100;   % ms
ws = 500;         % window size (pts)

sortorder = [11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29     1     2     3 4     5     6     7     8     9    10];

pulsestart = pulsestart/1e3;  % ms to s
duration = duration/1e3;

%cd '/Volumes/Promise Pegasus/Manual Backup/Lab/Videos/Zebrafish/High Speed/20150728/';
%files = dir('*.avi');

for j = 1:length(files)
    files(j).sortorder = sortorder(j);
end


for m = 1:length(files)
%    vid{m}=aviread(files(m).name);
tvec = (1:length(vid{m}))/Fs;
pulseind = findnearest(tvec,pulsestart);
durationind = findnearest(tvec,duration);

disp(['File #' num2str(m) '/' num2str(length(files))]);

%comments{m} = ['(' num2str(m) ') Filename:' files(m).name(1:end-4)];

 


for j = 1:length(vid{m})
    I = adapthisteq(vid{m}(j).cdata); % equalize
    I = imcrop(I,rect);
    %I = imcrop(I,[10*round(size(I,1)/100) round(10*size(I,1)/100) round(90*size(I,1)/100) round(70*size(I,2)/100)]); %crop
    vidmov{m}(j).cdata = I;
    q = I < mean(mean(I)) - thresh*std(std(double(I)));
    [a, b]=ind2sub(size(q),find(q==1));
    for k = 1:size(q,2)
        if isempty(find(b==k, 1))==0
            linevec(m,j,k) = mean(a(b==k));
            if k>1
                linevecangle(m,j,k-1) = atan2d(linevec(m,j,k) - linevec(m,j,k-1),1);
            end
        else
            linevec(m,j,k) = NaN;
            if k>1
                linevecangle(m,j,k-1) = NaN;
            end
        end
    end
    
    linevec(m,j,:) = smooth(linevec(m,j,:),ws);
    
    for k = 1:size(q,2)
        if isempty(find(b==k))==0            
            if k>1
                linevecangle(m,j,k-1) = atan2d(linevec(m,j,k) -linevec(m,j,k-1),1);
            end
        else            
            if k>1
                linevecangle(m,j,k-1) = NaN;
            end
        end
    end
    
    %linevecangle(m,j,:) = linevecangle(m,j,:)./mean(linevecangle(m,j,:));
    
    
    [U,S,V]=svd(squeeze(linevec(m,j,isnan(linevec(m,j,:))==0)));  % use SVD to find eigenvalues
    
    if size(V,1) >= 4 && size(V,2) >= 4
        a = V(1,4); b = V(2,4); % Choose eigenvector from V
        c = V(3,4); d = V(4,4); % with smallest eigenvalue
        xc(m,j) = -b/(2*a); yc(m,j) = -c/(2*a); % Find center and radius of the
        r(m,j) = sqrt(xc(m,j)^2+yc(m,j)^2-d/a); % circle, a*(x^2+y^2)+b*x+c*y+d=0
    end
    
    if mod(j,round(length(vid{m})/5))==0
        disp([num2str(j/length(vid{m})*100) '% complete']);
    end
end

disp('Finding eigenvalues')

% find covariance of the angles
    AngleCov = nancov(squeeze(linevecangle(m,:,:)));
    AngleCov(isnan(AngleCov)==1)=0;
    AngleEigenvalue_largest(m) = svds(AngleCov,1);
%S=svd(squeeze(linevecangle(m,j,isnan(linevecangle(m,j,:))==0))); 
%if isempty(S)==0
%lineveceig_smallest(m) = S(end);
%lineveceig_largest(m) = S(1);
%end
if isempty(squeeze(linevecangle(m,j,isnan(linevecangle(m,j,:))==0))) ==0
lineveceig_largest(m) = svds(squeeze(linevecangle(m,j,isnan(linevecangle(m,j,:))==0)));
end

disp('Calculating velocities');

tvecPRE = tvec(1:pulseind-1);
tvecPULSE = tvec(pulseind:pulseind+durationind);
tvecPOST = tvec(pulseind+durationind+1:end);

for k = 1:size(q,2)
    linevecvel(m,:,k) = squeeze(abs(gradient(tvec,linevec(m,:,k))));       % calculate velocity along the tail
    linevecvelPRE(m,:,k) = squeeze(abs(gradient(tvecPRE,linevec(m,1:pulseind-1,k))));
    linevecvelPULSE(m,:,k) = squeeze(abs(gradient(tvecPULSE,linevec(m,pulseind:pulseind+durationind,k))));
    linevecvelPOST(m,:,k) = squeeze(abs(gradient(tvecPOST,linevec(m,pulseind+durationind+1:end,k))));
end


linevecvel(m,isinf(linevecvel(m,:,:))) = NaN;
linevecvelPRE(m,isinf(linevecvelPRE(m,:,:))) = NaN;
linevecvelPULSE(m,isinf(linevecvelPULSE(m,:,:))) = NaN;
linevecvelPOST(m,isinf(linevecvelPOST(m,:,:))) = NaN;

VelMean(m) = (nanmean(nanmean(abs(linevecvel(m,:,:)))));
VelSTD(m) = (nanstd(nanstd(abs(linevecvel(m,:,:)))));
VelMax(m) = max(max(abs(linevecvel(m,:,:))));

VelMean_PRE(m) = (nanmean(nanmean(linevecvelPRE(m,:,:))));
VelSTD_PRE(m) = (nanstd(nanstd(linevecvelPRE(m,:,:))));

VelMean_PULSE(m) = (nanmean(nanmean(linevecvelPULSE(m,:,:))));
VelSTD_PULSE(m) = (nanstd(nanstd(linevecvelPULSE(m,:,:))));

VelMean_POST(m) = (nanmean(nanmean(linevecvelPOST(m,:,:))));
VelSTD_POST(m) = (nanstd(nanstd(linevecvelPOST(m,:,:))));

VelMeanLength(m,:) = squeeze(nanmean(abs(linevecvel(m,:,:)),2));
AngleLength(m,:) = squeeze(nanvar(abs(linevecangle(m,:,:)),[],2));
AngleVar(m) = var(AngleLength(m,:));
VelVarLength(m,:) =  squeeze(nanvar(abs(linevecvel(m,:,:)),[],2));

    disp('Complete.');
    
end
%% plot one of the vectors on an image
vid_ind=500;
movnum = 25;
I = adapthisteq(vid{movnum}(vid_ind).cdata);
I = imcrop(I,rect);
%I = imcrop(I,[10*round(size(I,1)/100) round(10*size(I,1)/100) round(90*size(I,1)/100) round(70*size(I,2)/100)]);
figure;
imshow(I);


hold on;
plot(squeeze(linevec(movnum,vid_ind,:)),'r','LineWidth',3);
    
    
%% create a video    
movnum=25;

path1=cd;

vid_out = [path1 '/Analysis/' comments{movnum}(15:end) '-analyzed.avi'];
writerObj = VideoWriter(vid_out);   % create writer object (Image Processing toolbox)
writerObj.Quality=100;
%writerObj.Width=sizeV(1)*2;writerObj.Height=sizeV(2)*2;
writerObj.FrameRate = 125;      % slowed to 0.25X

open(writerObj)                     % open the object

h2=figure;                          % set the figure for movie making
set(h2, 'Position', [100, 100, 1049, 895]);
sizeV = size(vid{movnum});    
    
for j = 1:length(vid{movnum})                   % loop through each frame
    I = adapthisteq(vid{movnum}(j).cdata);
    I = imcrop(I,rect);
    %I = imcrop(I,[10*round(size(I,1)/100) round(10*size(I,1)/100) round(90*size(I,1)/100) round(70*size(I,2)/100)]);
    imshow(I);colormap('gray');    % plot video frame
    hold on;
    title('0.25X Playback');        % CHECK this with writerObj.FrameRate!
    plot(squeeze(linevec(movnum,j,:)),'r','LineWidth',3);
    h(j) = getframe(h2);                % create object/snapshot of figure frame
    writeVideo(writerObj,h(j));         % write to video object
    if mod(j,round(length(vid{m})/5))==0
        disp([num2str(j/length(vid{m})*100) '% complete']);
    end
end

close(writerObj);                   % close video object
    
%% find eigenfish

for m = length(vid):-1:1
    disp(num2str(m));
    disp('calculating angles');
    for j = 1:length(vid{m})
        I = adapthisteq(vid{m}(j).cdata); % equalize
        I = imcrop(I,rect);
    	%I = imcrop(I,[10*round(size(I,1)/100) round(10*size(I,1)/100) round(90*size(I,1)/100) round(70*size(I,2)/100)]); %crop
        vidmov{m}(j).cdata = I;
        q = I < mean(mean(I)) - thresh*std(std(double(I)));
        [a, b]=ind2sub(size(q),find(q==1));
        
        linevec2 = linevec(:,:,1:5:end);
        for k = 1:size(linevec2,3)
            if isempty(find(b==k, 1))==0            
                if k>1
                    linevecangle(m,j,k-1) = atan2d(linevec2(m,j,k) -linevec2(m,j,k-1),1);
                end
            else            
                if k>1
                    linevecangle(m,j,k-1) = NaN;
                end
            end
        end
        %}
        
        linevecangle_norm(m,j,:) = linevecangle(m,j,:) - mean(linevecangle(m,j,:));
        
        clear thetacov
         %%for kk = 1:length(linevecangle(m,j,:))
            for ll = 1:length(linevecangle(m,j,:))
                thetacov(:,ll) = bsxfun(@times,bsxfun(@minus,linevecangle(m,j,:),mean(linevecangle(m,j,:))),bsxfun(@minus,linevecangle(m,j,ll),mean(linevecangle(m,j,:))));
                thetacov(:,ll) = thetacov(:,ll)./(length(linevecangle(m,j,:))-1);
               % thetacov(kk,ll)=(linevecangle(m,j,kk)-mean(linevecangle(m,j,:)))*(linevecangle(m,j,ll)-mean(linevecangle(m,j,:)))./(length(linevecangle(m,j,:))-1);
            if mod(ll,round(length(linevecangle(m,j,:))/2))==0
           %     disp(['Covariance for time point ' num2str(j) ': ' num2str(ll/length(linevecangle(m,j,:))*100) '% complete']);
            end
            end
            
        %[thetaU{m,j},thetaS{m,j},thetaV{m,j}] = svd(thetacov(isnan(thetacov)==0));
        if isempty(find(isnan(thetacov)==1))==1 
        [thetaU{m,j},thetaS{m,j},thetaV{m,j}] = eig(thetacov);
        
        
        % Extract eigenvectors for five largest eigenvalues
        eigvecA{m,j} = thetaV{m,j}(:,end);eigvalA(m,j) = thetaS{m,j}(end,end);
        eigvecB{m,j} = thetaV{m,j}(:,end-1);eigvalB(m,j) = thetaS{m,j}(end-1,end-1);
        eigvecC{m,j} = thetaV{m,j}(:,end-2);eigvalC(m,j) = thetaS{m,j}(end-2,end-2);
        eigvecD{m,j} = thetaV{m,j}(:,end-3);eigvalD(m,j) = thetaS{m,j}(end-3,end-3);
        eigvecE{m,j} = thetaV{m,j}(:,end-4);eigvalE(m,j) = thetaS{m,j}(end-4,end-4);
        
        curvature(m,j,:) = gradient(squeeze(linevecangle(m,j,:)));
        eigvecphaseAB(m,j,:) = atan2(-(eigvecB{m,j}-mean(eigvecB{m,j}))./std(eigvecB{m,j}),(eigvecA{m,j}-mean(eigvecA{m,j}))./std(eigvecA{m,j}));
        % calculate fractional variance for 1 to the maximum number of
        % eigenvalues, so fractionalVariance(3) is the fractional variance
        % for the three largest eigenvalues
        mn=1;
        for kk = length(thetaS{m,j}):-1:1
            fractionalVariance(m,j,mn) = sum(diag(thetaS{m,j}(1:kk,1:kk)))./sum(diag(thetaS{m,j}));
            mn=mn+1;
        end
        end
        if mod(j,round(length(vid{m})/5))==0
            disp([num2str(j/length(vid{m})*100) '% complete']);
        end
    end
end

