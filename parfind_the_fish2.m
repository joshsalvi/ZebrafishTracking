function mid_pt_all = parfind_the_fish2(file)
    
    % turn  ON/OFF (1/0) visualization
    % can only do parallel computing with VIZ = 0
    VIZ = 0;
    
    % start parallel toolbox workers
    if ~VIZ
        matlabpool;
    end
    
    % movie file
    fishvideo = VideoReader(file);
    Nframes = fishvideo.NumberOfFrames;
    
    frame = read(fishvideo,1); 
    FrameSize = size(frame);
    
    imshow(frame); % show the first frame, manually click on the snout
    snout = ginput(1); % most anterior point on midline
    waist = ginput(1); % midline point at agar edge
    tip = ginput(1); % tail tip
    TL = sqrt((snout(1) - tip(1))^2 + (snout(2) - tip(2))^2); % total length snout to tip
    TL2 = sqrt((waist(1) - tip(1))^2 + (waist(2) - tip(2))^2); % total length waist to tip
    width = 0.045*10e2; % keep this value or tweak if needed, determines circle size
    R_step = TL/50; % step size between midline points
    num_caudal = round(TL2/R_step); % number of step points between waist and tip
    taper = 1*exp(-0.04*(1:num_caudal))+0.15;
    R_search = [width*taper, zeros(1,50)]; % determines how circle radii decrease
    mid_pt_all = zeros(2, num_caudal+1,500); % pre-allocate
    close;
    
    % beginning of automated tracking
    % loop through frames (using all cores on CPU)
    if VIZ
        
        for frame_num = 1:Nframes
            frame = read(fishvideo, frame_num);
            mid_pt = zeros(2, num_caudal+1); % pre-allocate
            frame_t = adaptivethreshold(frame, 100, 0.03, 0);
            
            imshow(frame_t);
            %pause;
            hold;
            
            mid_pt(:,1) = [snout(1); snout(2)]; % first center is first midline pt
            mid_pt(:,2) = [waist(1); waist(2)]; % second center is second midline pt

            done = 0; % not done
            i_mid = 2;
            while(~done)

                dir = mid_pt(:, i_mid) - mid_pt(:, i_mid-1); % points in the search direction
                theta_dir = atan2(dir(2),dir(1));
                [x_c,y_c] = make_half_circle(R_search(i_mid), [mid_pt(1,i_mid), mid_pt(2,i_mid)],theta_dir);
                c_p = [x_c; y_c; zeros(1, length(x_c))];
                
                % eliminate circle points out of the frame
                c_p(:, c_p(1, :) < 1) = []; 
                c_p(:, c_p(2, :) < 1) = [];
                c_p(:, c_p(1, :) > FrameSize(1)) = [];
                c_p(:, c_p(2, :) > FrameSize(2)) = [];

                plot(c_p(1, :), c_p(2, :),'b');
               
                % find points on the circle that are also on the fish
                for j = 1:length(c_p)
                    if frame_t(ceil(c_p(2,j)), ceil(c_p(1,j))) < 1
                        c_p(3, j) = 1;
                    end
                end

                c_p(:, c_p(3,:) == 0) = []; % get rid of circle points not on the fish

                plot(c_p(1, :), c_p(2, :), 'g');    
               
                if ~isempty(c_p) 
                    c_p = c_p(1:2, :);
                    c_p_theta = atan2(c_p(2, :) - mid_pt(2, i_mid), c_p(1, :) - mid_pt(1, i_mid));
                    mean_theta = circ_median(c_p_theta); % compute mean circle point, which is the new midline point
                    mid_pt(:, i_mid+1) = [real(R_step*exp(1i*mean_theta)); imag(R_step*exp(1i*mean_theta))] + mid_pt(:, i_mid); % new midline pt
                    i_mid = i_mid + 1;
                else  
                    done = 1; % done 
                end
            end         
        
            mid_pt = mid_pt(:, 1:num_caudal+1);

            plot(mid_pt(1, :), mid_pt(2, :),'or-');
            hold;
            pause(0.1);
            
            mid_pt_all(:, :, frame_num) = mid_pt;
        
        end
        
    else
        
        parfor frame_num = 1:Nframes
            frame = read(fishvideo,frame_num);
            mid_pt = zeros(2,num_caudal+1); % pre-allocate
            frame_t = adaptivethreshold(frame,100,0.02,0);

            mid_pt(:,1) = [snout(1);snout(2)]; % first center is first midline pt
            mid_pt(:,2) = [waist(1);waist(2)]; % second center is second midline pt

            done = 0; % not done
            i_mid = 2;
            while(~done)

                dir = mid_pt(:,i_mid) - mid_pt(:,i_mid-1); % points in the search direction
                theta_dir = atan2(dir(2),dir(1));
                [x_c,y_c] = make_half_circle(R_search(i_mid),[mid_pt(1,i_mid),mid_pt(2,i_mid)],theta_dir);
                c_p = [x_c;y_c;zeros(1,length(x_c))];

                % eliminate circle points out of the frame
                c_p(:,c_p(1,:) < 1) = []; 
                c_p(:,c_p(2,:) < 1) = [];
                c_p(:,c_p(1,:) > FrameSize(1)) = [];
                c_p(:,c_p(2,:) > FrameSize(2)) = [];

                % find points on the circle that are also on the fish
                for j = 1:length(c_p)
                    if frame_t(ceil(c_p(2,j)),ceil(c_p(1,j))) < 1
                        c_p(3,j) = 1;
                    end
                end

                c_p(:,c_p(3,:) == 0) = []; % get rid of circle points not on the fish

                if ~isempty(c_p) 
                    c_p = c_p(1:2,:);
                    c_p_theta = atan2(c_p(2,:) - mid_pt(2,i_mid),c_p(1,:) - mid_pt(1,i_mid));
                    mean_theta = circ_median(c_p_theta); % compute mean circle point, which is the new midline point
                    mid_pt(:,i_mid+1) = [real(R_step*exp(1i*mean_theta));imag(R_step*exp(1i*mean_theta))] + mid_pt(:,i_mid); % new midline pt
                    i_mid = i_mid + 1;
                else  
                    done = 1; % done 
                end
            end         
        
            mid_pt = mid_pt(:,1:num_caudal+1);  
            mid_pt_all(:,:,frame_num) = mid_pt;
        
        end
        
        matlabpool close;
    end
    
    % plot results overlayed on first frame
    figure;
    %imshow(fishvideo(1).cdata);
    hold;
    for k = 1:Nframes
        plot(mid_pt_all(1,:,k), mid_pt_all(2,:,k), 'r');
    end
    hold;
    
end

function [x,y] = make_half_circle(r,c,dir)
    p = r*exp(1i*(dir + pi*(-0.5:0.01:0.5)));
    x = real(p)+c(1);
    y = imag(p)+c(2);
end

function mean_theta = circ_median(theta)
    N = length(theta);
    mean_theta = atan2(sum(sin(theta))/N, sum(cos(theta))/N);
end

function bw=adaptivethreshold(IM,ws,C,tm)
%ADAPTIVETHRESHOLD An adaptive thresholding algorithm that seperates the
%foreground from the background with nonuniform illumination.
%  bw=adaptivethreshold(IM,ws,C) outputs a binary image bw with the local 
%   threshold mean-C or median-C to the image IM.
%  ws is the local window size.
%  tm is 0 or 1, a switch between mean and median. tm=0 mean(default); tm=1 median.
%
%  Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn)
%  at Tsinghua University, Beijing, China.
%
%  For more information, please see
%  http://homepages.inf.ed.ac.uk/rbf/HIPR2/adpthrsh.htm

    if (nargin<3)
        error('You must provide the image IM, the window size ws, and C.');
    elseif (nargin==3)
        tm=0;
    elseif (tm~=0 && tm~=1)
        error('tm must be 0 or 1.');
    end

    IM=mat2gray(IM);

    if tm==0
        mIM=imfilter(IM,fspecial('average',ws),'replicate');
    else
        mIM=medfilt2(IM,[ws ws]);
    end
    sIM=mIM-IM-C;
    bw=im2bw(sIM,0);
    bw=imcomplement(bw);
end
