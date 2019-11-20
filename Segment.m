function [Image_Segments] = Segment(Labels,imagex,Min_size)
% Segment a grayscale image using the binary image
% Input: Labels - matrix of labels after component labelling
%        imagex - a grayscale image from which the segments are extracted
%        min_size - a minimum size for the segment (denoising)
% Output: Image_Segments - a cell of matrices, each matrix representing a segmented character



MAXI = max(Labels(:));
MINI = 1; % All labels start from 1
count = 1;
Labels = cast(Labels,'int8');
imagex = cast(imagex,'int8');
Full_Image_Segments = zeros(size(Labels));

for kk = MINI:MAXI
    I = (Labels == kk); Iint = int8(I);
    Full_Image_Segments(:,:,count) = Iint.*imagex; % Multiply the logical array with grayscale to extract the character
    count = count + 1; 
end


% Denoising to remove noise segments: commented out as component labelling already removes noisy segments


Image_Segments = cell(1,size(Full_Image_Segments,3));  

% Cropping: crop segments to fit the character and remove zero padding
for kk = 1:size(Full_Image_Segments,3)
    Image_Segments{kk} = Crops(Full_Image_Segments(:,:,kk));
end


end    
