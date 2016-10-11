I=imread('step_snr100.gif');
BW = edge(I,'Canny',[0.0503,0.1006],12);
imshow(BW);
title('Built-in Canny Detector on step snr100')



%{

for u=1:(row/W)
    for v=1:(col/W)
        

    B_patch=z_depth((W*(u-1)+1):(W*(u-1)+W),(W*(v-1)+1):(W*(v-1)+W));
    B = reshape(B_patch,[W^2,1]);
        
    
    parameter = (inv(transpose(A)*A))*transpose(A)*B;

    a=parameter(1);
    b=parameter(2);
    c=parameter(3);
    d=parameter(4);
    e=parameter(5);
    f=parameter(6);
    
    
        for x=(W*(u-1)+1):(W*(u-1)+W)
            for y=(W*(v-1)+1):(W*(v-1)+W)
                K_poly(x,y)=1000000*(-b^2+4*a*c)/(1+(d+2*a*x+b*y)^2+(e+b*x+2*c*y)^2)^2;
                H_poly(x,y)=1/2*(1+(d+2*a*x+b*y)^2+(e+b*x+2*c*y)^2)^1.5*(-2*b*(d+2*a*x+b*y)*(e+b*x+2*c*y)+2*c*(1+(d+2*a*x+b*y)^2)+2*a*(1+(e+b*x+2*c*y)^2));
            end
        end
    
    
    end
end


%}