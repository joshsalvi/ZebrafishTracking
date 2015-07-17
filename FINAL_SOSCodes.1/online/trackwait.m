% -------------------------------------------------------------------------
% trackwait.m 
% 
% -------------------------------------------------------------------------
% Description:
%
% The timer is waiting until the animal leaves the initial bounding box so
% that the whole background can be reconstructed and tracking can start.
% -------------------------------------------------------------------------

function [  ] = trackwait(  )

global x y l w iscrossing vidobj firstbg threshold

% The while loop is executed until iscrossing=2, which indicates that the 
% animal has left the initial bounding box and, accordingly the background
% reconstruction can be finalized

%--
while ~(iscrossing==2)

secondframe=getsnapshot(vidobj);
iscrossing=0;
removesecondframe= firstbg-secondframe; 
removebwsecondframe=im2bw(removesecondframe,threshold);

  % Rows
  for i=[x-2,x-1,x+l+1,x+l+2]
    for j=(y-2):(w+y+2)
        if removebwsecondframe(j,i)==1
           iscrossing=1;
           break
        end
    end
  end
  % Columns
  for j=[y-2,y-1,y+w+1,y+w+2]
    for i=(x-2):(l+x+2)
        if removebwsecondframe(j,i)==1
           iscrossing=1;
           break
        end
    end
  end
  
[LO,num]=bwlabel(removebwsecondframe,8);
Areas=regionprops(LO,'area');
area_val=[Areas.Area];
maxarea=max(area_val);
idxBig= find(maxarea == area_val);
it2=ismember(LO,idxBig);
[r c]=find(it2);
maxc=max(c);
minc=min(c);
maxr=max(r);
minr=min(r);

% Make sure there is nothing inside the little box and escape the while loop
if (iscrossing==0) &&  ((maxc>x+l) || (maxr>y+w) || (minc<x) || (minr<y))
     iscrossing=2;
end
%

end
%--

end
