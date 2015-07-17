% -------------------------------------------------------------------------
% trackbg.m 
% 
% -------------------------------------------------------------------------
% Description:
%
% Once the animal is outside the initial bounding box, we take a snapshot
% of the inner part so that a complete background reconstruction is achieved
% -------------------------------------------------------------------------

function [  ] = trackbg(  )

global vidobj x y l w goodbackground

secondbg= getsnapshot(vidobj);

% Replaceing pixels from "goodbackground" by those pixels in "secondbg" 
for ii=x:x+l
    for jj=y:y+w
        goodbackground(jj,ii)=secondbg(jj,ii);
    end
end

end
