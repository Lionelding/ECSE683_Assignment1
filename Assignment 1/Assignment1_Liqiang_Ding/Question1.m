%% Initialize Variables 
clear;
clc;
Sigma=1.5;
%Thres =0.75;
% For checkerboard.jpg
Thres =0.0005;
Thres_Discret=5;
Scaler=0.05;
W=5;

%% Filter image with Gaussian to reduce noise
I_RGB=imread('find_the_corners.jpg');
I_RGB=im2double(I_RGB);
I_GRAY_opt = rgb2gray(I_RGB);
I_GRAY_poly = rgb2gray(I_RGB);

figure(1);
imshow(I_GRAY_opt);
title('Input Image');

I_opt = imgaussfilt(I_GRAY_opt,Sigma);
I_poly = imgaussfilt(I_GRAY_poly,Sigma);
[row, col]=size(I_opt);

figure(2)
imshow(I_opt);
title('Blurred Input Image')


%% using a discret operator 
operator = [-1 0 1;-1 0 1;-1 0 1];
Ix_opt = conv2(I_opt, operator, 'same');
Iy_opt = conv2(I_opt, operator', 'same');

Ixy_opt=Ix_opt .* Iy_opt;

I_good_opt=zeros(row,col);
R_m_opt=zeros(row,col);

for x=1:(row-W),
   for y=1:(col-W),

       % Filter the gradient with Gaussian filter and Construct a window 
       Gaussian_Ix=imgaussfilt(Ix_opt,Sigma);
       Gaussian_Iy=imgaussfilt(Iy_opt,Sigma);
       Gaussian_Ixy=imgaussfilt(Ixy_opt,Sigma);  
       
       temp_Ix=Gaussian_Ix((x):(x+W),(y):(y+W));
       temp_Iy=Gaussian_Iy((x):(x+W),(y):(y+W));
       temp_Ixy=Gaussian_Ixy((x):(x+W),(y):(y+W));
       
       
       C_opt = [sum(sum(temp_Ix .^2)) sum(sum(temp_Ixy .^2)); sum(sum(temp_Ixy .^2)) sum(sum(temp_Iy .^2))];

       R_opt = (det(C_opt) - Scaler * (trace(C_opt) ^ 2));
       % Solve for ?s
       R_m_opt(x,y)=R_opt;


   end
end

% Label the Corner with and Remove the duplicated label using a window
for x=(1+5):(row-5),
   for y=(1+5):(col-5),
       if (R_m_opt(x,y)==max(max(R_m_opt(x-5:x+5,y-5:y+5))) && R_m_opt(x,y) > Thres_Discret )
            I_good_opt(x, y) = 125;
            I_GRAY_opt = insertMarker(I_GRAY_opt,[y x]);
       
       end
   end
end

figure (25)
imshow(I_good_opt);
title('Corner Detected with Discret Operators');

figure (26);
imshow(I_GRAY_opt);
title('Mapped to Orignial Image with Discret Operators');



%% using a least-squares approximation to estimate the values of a polynomial 
% h(x,y) = ax2 + bxy + cy2 + dx + ey + f at each coordinate. 

A=[4,4,4,-2,-2,1; ...
   1,2,4,-1,-2,1; ...
   0,0,4, 0,-2,1; ...
   1,-2,4,1,-2,1; ...
   4,-4,4,2,-2,1; ...
   4,2,1,-2,-1,1; ...
   1,1,1,-1,-1,1; ...
   0,0,1,0, -1,1; ...
   1,-1,1,1,-1,1; ...
   4,-2,1,2,-1,1; ...
   4,0,0, -2,0,1; ...
   1,0,0, -1,0,1; ...
   0,0,0,0,0,1;   ...
   1,0,0,1,0,1;   ...
   4,0,0,2,0,1;   ...
   4,-2,1,-2,1,1; ...
   1,-1,1,-1,1,1; ...
   0,0,1,0,1,1;   ...
   1,1,1,1,1,1;   ...
   4,2,1,2,1,1;   ...
   4,-4,4,-2,2,1; ...
   1,-2,4,-1,2,1; ...
   0,0,4,0,2,1;   ...
   1,2,4,1,2,1;   ...
   4,4,4,2,2,1    ...
    ];
R_m_poly=zeros(row,col);
I_good_poly=zeros(row,col);

for u=(1+2):(row-2)
    for v=(1+2):(col-2)
        
        % Compute the paramter of the polynomial function
        B_patch=I_GRAY_opt((u-2):(u+2),(v-2):(v+2));
        B = reshape(B_patch,[W^2,1]);
        
        parameter = (inv(transpose(A)*A))*transpose(A)*B;

        a=parameter(1);
        b=parameter(2);
        c=parameter(3);
        d=parameter(4);
        e=parameter(5);
        f=parameter(6);
        
    
        % Compute the derivative of polynomial function
        Ix_vector=[0,0,0,2*a,b,d]';
        Iy_vector=[0,0,0,b,2*c,e]';
 
        Ix_poly= A*Ix_vector;
        Iy_poly= A*Iy_vector;
        
        Ix_poly = reshape(Ix_poly,[W,W]);
        Iy_poly = reshape(Iy_poly,[W,W]);
        
        Ixy_poly=Ix_poly .* Iy_poly;
        
        % Filter the gradient with Gaussian filter and Construct a window 
        Gaussian_Ix_poly=imgaussfilt(Ix_poly,Sigma);
        Gaussian_Iy_poly=imgaussfilt(Iy_poly,Sigma);
        Gaussian_Ixy_poly=imgaussfilt(Ixy_poly,Sigma);        
        C_poly = [sum(sum(Gaussian_Ix_poly .^2)) sum(sum(Gaussian_Ixy_poly .^2)); sum(sum(Gaussian_Ixy_poly .^2)) sum(sum(Gaussian_Iy_poly .^2))];
        
        % Solve for lemda
        R_poly = (det(C_poly) - Scaler * (trace(C_poly) ^ 2));
        R_m_poly(u,v)=R_poly;

    end
end


% Label the Corner with and Remove the duplicated label using a window
for u=(1+5):(row-5),
   for v=(1+5):(col-5),
       if (R_m_poly(u,v)==max(max(R_m_poly(u-5:u+5,v-5:v+5))) && R_m_poly(u,v) > Thres )
            I_good_poly(u, v) = 125;
            I_GRAY_poly = insertMarker(I_GRAY_poly,[v u]);
       
       end
   end
end


figure (50)
imshow(I_good_poly);
title('Corner Detected with Least Square Method');

figure (51);
imshow(I_GRAY_poly);
title('Mapped to Orignial Image with Least Squared Method');


