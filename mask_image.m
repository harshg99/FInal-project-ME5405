function [img_fil]=mask_image(img,mask)

    
% Outputs a modfied image by convolving through a given a mask
% %--------------------
% Inputs:
% img       [uint8] image in uint8
% mask      [int]   specifies mask
%                   
%--------------------
% Outputs:
% img_fil   [double] processed rimage after convolving/correlating with mask
%
    
    
    %padding image with zerosto correct dimensions to get appropriate output
   
    size_i1=size(img);
    size_m=size(mask);
    %padding with image
    %img=[img(:,1:(size_m(2)-1)/2),img,img(:,size_i1(2)-(size_m(2)-1)/2:size_i1(2))];
    %img=[img(1:(size_m(1)-1)/2,:);img;img(size_i1(1)-(size_m(1)-1)/2:size_i1(1),:)];
    
    %padding with zeros
    img=[zeros(size(img,1),(size_m(2)-1)/2),img,zeros(size(img,1),(size_m(2)-1)/2)];
    img=[zeros((size_m(1)-1)/2,size(img,2));img;zeros((size_m(1)-1)/2,size(img,2))];
   
    %size_i1=size(img);
    img=im2double(img);
    img_fil=zeros(size_i1(1),size_i1(2));
    
    for(i=1:size(img_fil,1))
        for(j=1:size(img_fil,2))
            img_fil(i,j)=correl(img(i:(i+size_m(1)-1),j:(j+size_m(2)-1)),mask);
        end
    end

end