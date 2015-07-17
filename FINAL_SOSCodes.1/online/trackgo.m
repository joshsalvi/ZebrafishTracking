% -------------------------------------------------------------------------
% trackgo.m 
% 
% -------------------------------------------------------------------------
% Description:
%
% Main part of the track.m script. Once the background has been properly
% reconstructed, trackgo tracks the animal posture and position in a loop.
% -------------------------------------------------------------------------

function [] = trackgo()

global i threshold vidobj impro imraw finalbbox bboxref nextbboxref boxsize goodbackground

i=i+1;
% The running number and the time could be displayed 
%display(num2str(toc));
%display(num2str(i))

% Current frame snapshot 
myframe=getsnapshot(vidobj); 
% Image processing
removebgframe=(goodbackground)-(myframe);
itest=im2bw(removebgframe,threshold);
itestint=imcrop(itest,nextbboxref(i,:));
%
% Again, find the biggest object within the bbox.
% Not so crucial now but good to avoid bbox growing due to small dirts
[LO,num]=bwlabel(itestint,8);
Areas=regionprops(LO,'area');
area_val=[Areas.Area];
maxarea=max(area_val);
idxBig= find(maxarea == area_val);
it2=ismember(LO,idxBig);
%

% --
% B-boxes computation
% Crop the picture to focus on the biggest salient object
[r c]=find(it2);
maxc=max(c)+boxsize;
minc=min(c)-boxsize;
maxr=max(r)+boxsize;
minr=min(r)-boxsize;
bbox=[minc minr maxc-minc maxr-minr]; 
%
bboxref(i,:)=bbox; 
% Now box in the arena reference system (with respect to the entire frame)
absolutebbox = [bboxref(i,1)+nextbboxref(i,1) bboxref(i,2)+nextbboxref(i,2) bboxref(i,3) bboxref(i,4)];
% A reference to crop the next picture:
nextbbox=[absolutebbox(1)-boxsize absolutebbox(2)-boxsize absolutebbox(3)+2*boxsize absolutebbox(4)+2*boxsize];
if nextbbox(1)<0
   nextbbox(1)=0;
end
if nextbbox(2)<0
   nextbbox(2)=0;
end
%
% Finally, just use the current final good finalbbox
nextbboxref((i+1),:)=absolutebbox;
%nextbboxref((i+1),:)=nextbbox;


% ---
% Relevant data to save
%
% (A) The position of the bounding box in the absolute reference frame
finalbbox(i,:)=absolutebbox;
%
% (B) Output images:
% The raw cropped image
smallImage0=imcrop(myframe,absolutebbox);
% The processed and cropped image
smallImage=imcrop(it2,bbox);
% Efficient array format to save the images
imraw{i}=smallImage0;
impro{i}=smallImage;
% In case individual images are preferred:
% outname1 = sprintf('raw%06d.BMP',i);
% imwrite(smallImage0,outname1);
% outname2 = sprintf('pro%06d.BMP',i);
% imwrite(smallImage,outname2);

end

