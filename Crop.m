function [Cropped] = Crops(imagex,thold)
Cropped = imagex(imagex > thold);
cond = (imagex > thold);

hcount = 0;
for ii = 1:length(cond,1)
    if (sum(cond(ii,:)) >= 1)
        hcount = hcount + 1;
    end
end

vcount = 0;
for ii = 1:length(cond,2)
    if (sum(cond(:,ii)) >= 1)
        vcount = vcount + 1;
    end
end

Cropped = reshape(Cropped,hcount,vcount);
end