%% Initializing Variables
SIGMA=0.5;
H_pos=0.05;
H_neg=-H_pos;
K_pos=0.05;
K_neg=-K_pos;

H_poly_pos=0.05;
H_poly_neg=-H_poly_pos;
K_poly_pos=0.05;
K_poly_neg=-K_poly_pos;

I_rgb=imread('problem3_angel-deg0.00.ppm');
I_rgb=im2double(I_rgb);
I_gray=rgb2gray(I_rgb);
I = imgaussfilt(I_gray,SIGMA);

z=load('problem3_angel_zdepth.mat');
z_depth = z.z_depth; % watch out the sign
[row, col]=size(z_depth);

%% Implementing operators for Ix, Iy, Ixx, Iyy, and Ixy using a 2D Gaussian as a basis. 
gaussian_filter = fspecial('gaussian', 5, SIGMA);
[Gx, Gy] = gradient(gaussian_filter);


[G_z_depth_x,G_z_depth_y] = gaussgradient(z_depth,SIGMA); %% not using z_depth2

[Gxx,Gxy]=gaussgradient(Gx,SIGMA);
[Gyx,Gyy]=gaussgradient(Gy,SIGMA);

[G_z_depth_xx, G_z_depth_xy] = gaussgradient(G_z_depth_x,SIGMA);
[G_z_depth_yx, G_z_depth_yy] = gaussgradient(G_z_depth_y,SIGMA);


%{
Ix = cat(3, Gx, G_z_depth_x);
Iy = cat(3, Gy, G_z_depth_y);
Ixx= cat(3, Gxx, G_z_depth_xx);
Iyy= cat(3, Gyy, G_z_depth_yy);
Ixy= cat(3, Gxy, G_z_depth_xy);
%}
%% Compute the H-K map for the angel data set

K_up=(G_z_depth_xx .*G_z_depth_yy)- (G_z_depth_xy .^2);
K_down=(1+ G_z_depth_x .^2 + G_z_depth_y .^2) .^2;

K=zeros(row,col);
for i=1:row
    for j=1:col
        K(i,j) = K_up(i,j)/K_down(i,j);
    end
end

H_up=((1+G_z_depth_x .^2).*G_z_depth_yy)-(2*G_z_depth_x .*G_z_depth_y.*G_z_depth_xy)+((1+G_z_depth_y .^2).*G_z_depth_xx);
H_down=((1+(G_z_depth_x .^2)+(G_z_depth_y .^2)) .^1.5);

H=zeros(row,col);
for i=1:row
    for j=1:col
        H(i,j) = H_up(i,j)/H_down(i,j);
    end
end

H=H*0.5;

figure(50);
imshow(K);
title('K map with operators')

figure(51);
imshow(H);
title('H map with operators')

S=zeros(col,row,3);
S=im2uint8(S);
%% Compute the shape label image, S, by assigning a shape label si
for i=1:row
    for j=1:col
        if ((abs(K(i,j))<K_pos)&&(abs(H(i,j))<H_pos))
            %plane 
            S(i,j)=0;
        elseif ((abs(K(i,j))<K_pos)&&(H(i,j)>H_pos))
            %concave cylindrical
            S(i,j,2)=50;           
        elseif ((abs(K(i,j))<K_pos)&&(H(i,j)<H_neg))
            %convex cylindrical
            S(i,j,2)=100; 
        elseif ((K(i,j)>K_pos)&&(H(i,j)>H_pos))
            %concave elliptic
            S(i,j,1)=150; 
        elseif ((K(i,j)>K_pos)&&(H(i,j)<H_neg))
            %convex elliptic
            S(i,j,3)=200; 
        elseif (K(i,j)<K_neg)
            %hyperbolic
            S(i,j,3)=250; 
        end
    end
end

figure(100)
imshow(S);
title('Shape Image using Operators')
%legend('plane','concave cylindrical','convex cylindrical','concave elliptic','convex elliptic','hyperbolic');


%% using a least-squares approximation to estimate the values of a polynomial 
% h(x,y) = ax2 + bxy + cy2 + dx + ey + f at each coordinate. 

% Sampling window size
W=5;
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

K_poly=zeros(row,col);
H_poly=zeros(row,col,1);

for u=1:(row-W)
    for v=1:(col-W)
        

    B_patch=z_depth(((u-1)+1):((u-1)+W),((v-1)+1):((v-1)+W));
    B = reshape(B_patch,[W^2,1]);
        
    
    parameter = (inv(transpose(A)*A))*transpose(A)*B;

    a=parameter(1);
    b=parameter(2);
    c=parameter(3);
    d=parameter(4);
    e=parameter(5);
    f=parameter(6);
    
    
        for x=((u-1)+1):((u-1)+W)
            for y=((v-1)+1):((v-1)+W)
                K_poly(x,y)=1000000*(-b^2+4*a*c)/(1+(d+2*a*x+b*y)^2+(e+b*x+2*c*y)^2)^2;
                H_poly(x,y)=1/2*(1+(d+2*a*x+b*y)^2+(e+b*x+2*c*y)^2)^1.5*(-2*b*(d+2*a*x+b*y)*(e+b*x+2*c*y)+2*c*(1+(d+2*a*x+b*y)^2)+2*a*(1+(e+b*x+2*c*y)^2));
            end
        end
    
    
    end
end



%K_new=im2uint8(K_poly);
figure(150);
imshow(K_poly);
title('K map with polynomial function');

figure(151);
imshow(H_poly);
title('H map with polynomial function');


S_poly=zeros(col,row,3);
S_poly=im2uint8(S_poly);
%% Compute the shape label image, S, by assigning a shape label si
for i=1:row
    for j=1:col
        if ((abs(K_poly(i,j))<K_poly_pos)&&(abs(H_poly(i,j))<H_poly_pos))
            %plane 
            S_poly(i,j)=0;
        elseif ((abs(K_poly(i,j))<K_poly_pos)&&(H_poly(i,j)>H_poly_pos))
            %concave cylindrical
            S_poly(i,j,2)=50;           
        elseif ((abs(K_poly(i,j))<K_poly_pos)&&(H_poly(i,j)<H_poly_neg))
            %convex cylindrical
            S_poly(i,j,2)=100; 
        elseif ((K_poly(i,j)>K_poly_pos)&&(H_poly(i,j)>H_poly_pos))
            %concave elliptic
            S_poly(i,j,1)=150; 
        elseif ((K_poly(i,j)>K_poly_pos)&&(H_poly(i,j)<H_poly_neg))
            %convex elliptic
            S_poly(i,j,3)=200; 
        elseif (K_poly(i,j)<K_poly_neg)
            %hyperbolic
            S_poly(i,j,3)=250; 
        end
    end
end


figure(200)
imshow(S_poly);
title('Shape Image using polynomial function');