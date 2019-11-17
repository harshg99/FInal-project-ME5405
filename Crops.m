function [Cropped] = Crops(imagex)

% To remove the 'background'
cond = (imagex ~= mode(imagex(:)));
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
