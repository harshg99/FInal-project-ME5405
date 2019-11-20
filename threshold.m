% Thresholding an image
function t_image = threshold(imagex,varargin)


% --------------------------------
% Thresholds the image either based on input threshold value or 
% using Otsu's method
%--------------------
% Inputs:
% imagex    [uint8] image in uint8
% varargin  [int]   specifies threshold
%                   if absent thresholding done through otsu's method
%--------------------
% Outputs:
% t_image   [unit8] binary image


if(size(varargin)>=1)   
    t_image = zeros(size(imagex,1),size(imagex,2));
    t_image = (imagex > varargin{1});
    return
end

imghist=his(imagex,256,0);
maxerror=intmin;
thresh=[];
mu_t=sum(imghist.*[0:255]);
sig_t=sum(((0:255)-mu_t).*((0:255)-mu_t).*imghist);
for(k=1:254)
    w0=sum(imghist(1:k+1));
    w1=sum(imghist(k+2:256));
    mu_0=dot(imghist(1:k+1),(0:k))/w0;
    mu_1=dot(imghist(k+2:256),(k+1:255))/w1;
    sig_0=sum(((0:k)-mu_0).*((0:k)-mu_0).*(imghist(1:k+1)))/w0;
    sig_1=sum(((k+1:255)-mu_1).*((k+1:255)-mu_1).*(imghist(k+2:256)))/w1;
    sig_w=w0*sig_0+w1*sig_1;
    sig_b=w0*w1*((mu_0-mu_1)^2); 
    %eta=sig_b/sig_w;
    eta=sig_b;
    if(eta>maxerror)
        thresh=k;
        maxerror=eta;
    elseif(abs(eta-maxerror)<=0.01)
            thresh=[thresh,k];    
    end
end
    %thresh=sum(thresh)/length(thresh);
    thresh=max(thresh);
    t_image = zeros(size(imagex,1),size(imagex,2));
    t_image = (imagex > thresh);
    return
end
