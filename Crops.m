function [Cropped] = Crops(imagex)

% Returns a cropped version of imagex by removing all the background elements
% Imagex can be grayscale/binary 2D
% Cropped is of the same image type as Imagex (but lesser dimensions)


% To remove the 'background'
cond = (imagex ~= mode(imagex(:))); % Background defined as the mode of the pixel intensity distribution
iiarray = zeros(size(imagex,1));
jjarray = zeros(size(imagex,2));

%% Create the Window
% Scan by keeping row constant and move along columns
hcount = 0;
for ii = 1:size(cond,1)
    if (sum(cond(ii,:)) >= 1) 
        hcount = hcount + 1; 
        ii_array(hcount) = ii;
    end
end

% Scan by keeping column constant and move along rows
vcount = 0;
for ii = 1:size(cond,2)
    if (sum(cond(:,ii)) >= 1) 
        vcount = vcount + 1; 
        jj_array(vcount) = ii;
    end
end

Window = zeros(hcount,vcount);
Cropped = imagex(ii_array(1):ii_array(end),jj_array(1):jj_array(end));

end
