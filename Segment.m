function [Image_Segments] = Segment(Labels,imagex,Min_size)
% Segment a grayscale image using the binary image
MAXI = max(Labels(:));
MINI = min(Labels(:));
count = 1;
Labels = cast(Labels,'int8');
imagex = cast(imagex,'int8');
Full_Image_Segments = zeros(size(Labels));
for kk = MINI:MAXI
    if (kk == 0) continue;
    else
        I = (Labels == kk); Iint = cast(I,'int8');
        New_Labels = Labels.*Iint;
        if (max(New_Labels(:)) == 0) continue; end 
        Full_Image_Segments(:,:,count) = New_Labels.*imagex;
        count = count + 1;      
    end
end


% Denoising to remove noise segments

count = 0;
for kk = 1:size(Full_Image_Segments,3)
    M = Full_Image_Segments(:,:,kk); M(M ~= 0) = 1;
    if (sum(M(:)) < Min_size)
        count = count + 1;
        check(count) = kk;
    end
end
if (count > 0)
    for kk = length(check):-1:1
        Full_Image_Segments(:,:,kk) = [];
    end
end

Image_Segments = cell(1,size(Full_Image_Segments,3));

% Cropping
for kk = 1:size(Full_Image_Segments,3)
    Image_Segments{kk} = Crops(Full_Image_Segments(:,:,kk));
end


end
        
    