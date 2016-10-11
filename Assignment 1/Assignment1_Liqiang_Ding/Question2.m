%% Initialize constants
Gaussian_size = 5;
SIGMA = 12;
Low_T = 0.08;
High_T = 0.16;

%% Load the Input Image
im = imread('step_snr1.gif');
im = double(im);
[row,col]=size(im);

%% Filter the input image with a Gaussian filter 
gaussian_filter = fspecial('gaussian', Gaussian_size, SIGMA);
im_gaussian = conv2(im, gaussian_filter, 'same');

%% Find directions and magnitude
[gaussian_x, gaussian_y] = gradient(gaussian_filter);
im_g_y = conv2(im_gaussian, gaussian_y, 'same');
im_g_x = conv2(im_gaussian, gaussian_x, 'same');
radian = atan2(im_g_y, im_g_x);

% Change from radian to degress
degree = radian*180/pi;

%% Assign the direction sector
degree_distance = zeros(row, col);
for i = 1  : row
    for j = 1 : col
        % if the degree is between 0 to 22.5 or between 157.5 to -157.5, 
        % assign the value of sector of 0
        if ((degree(i, j) > 0 ) && ...
            (degree(i, j) < 22.5) || ...
            (degree(i, j) > 157.5) && ...
            (degree(i, j) < -157.5))
        
            degree_distance(i, j) = 0;
        end
        % if the degree is between 22.5 to 67.5 or between -112.5 to -157.5, 
        % assign the value of sector of 45
        if ((degree(i, j) > 22.5) && ...
            (degree(i, j) < 67.5) || ...
            (degree(i, j) < -112.5) && ...
            (degree(i, j) > -157.5))
        
            degree_distance(i, j) = 45;
        end
        
        % if the degree is between 67.5 to 112.5 or between -67.5 to 112.5, 
        % assign the value of sector of 90
        if ((degree(i, j) > 67.5) && ...
            (degree(i, j) < 112.5) || ...
            (degree(i, j) < -67.5) && ...
            (degree(i, j) > 112.5))
        
            degree_distance(i, j) = 90;
        end
        
        % if the degree is between 112.5 to 157.5 or between -22.5 to -67.5, 
        % assign the value of sector of 135
        if ((degree(i, j) > 112.5) && ...
            (degree(i, j) <= 157.5) || ...
            (degree(i, j) < -22.5) && ...
            (degree(i, j) > -67.5))
         
            degree_distance(i, j) = 135;
        end
    end
end


%% Compute the magnitude of the edge using equation
im_g_mag = sqrt(im_g_x.^2 + im_g_y.^2);
figure(6);
imshow(im_g_mag, []);
title('Magnitude with Guassian Filter Size of 5 and Sigma of 12');

%% Non-maximal suppression
im_sup = zeros(row, col);
for i = 2  : (row-1)
    for j = 2 : (col-1)
        
        if (degree_distance(i, j) == 0)
        	if (im_g_mag(i, j) > im_g_mag(i, j - 1) && ...
                im_g_mag(i, j) > im_g_mag(i, j + 1))
                im_sup(i, j) = im_g_mag(i, j);
            else
                im_sup(i, j) = 0;
            end
        end

        if (degree_distance(i, j) == 45)
            if (im_g_mag(i, j) > im_g_mag(i + 1, j - 1) && ... 
                im_g_mag(i, j) > im_g_mag(i - 1, j + 1))
                im_sup(i, j) = im_g_mag(i, j);
            else
                im_sup(i, j) = 0;
            end
        end

        if (degree_distance(i, j) == 90)
            if (im_g_mag(i, j) > im_g_mag(i - 1, j) && ...
                im_g_mag(i, j) > im_g_mag(i + 1, j))
                im_sup(i, j) = im_g_mag(i, j);
            else
                im_sup(i, j) = 0;
            end
        end

        if (degree_distance(i, j) == 135)
            if (im_g_mag(i, j) > im_g_mag(i - 1, j - 1) && ...
                im_g_mag(i, j) > im_g_mag(i + 1, j + 1))
                im_sup(i, j) = im_g_mag(i, j);
            else
                im_sup(i, j) = 0;
            end
        end
    end
end

figure(7);
imshow(im_sup);
title('non-maximum suppressed with Guassian Filter Size of 5 and Sigma of 12');
%% Hysteresis Thresholding  
Th_L = Low_T * max(max(im_sup));
Th_H = High_T * max(max(im_sup));

T1=zeros(row,col);
T2=zeros(row,col);
for i = 1:row
    for j = 1:col     
        if (im_sup(i,j)>Th_H)
            T1(i,j)=1;
            T2(i,j)=1;
        elseif (im_sup(i,j)>Th_L)
            T1(i,j)=1;          
        else
        end
    end
end


indices=sub2ind(size(T1),find(T2));
im_thres=imfill(~T1,indices,4); 
im_thres=im_thres & T1;

% Print the T1 image
figure(8);
imshow(T1);
title('T1');

% Print the T2 image
figure(9);
imshow(T2);
title('T2');

% Print the result
figure(10);
imshow(im_thres);
title('Step Snr1 with Thres Lower of 0.05 and Thres High of 0.1');