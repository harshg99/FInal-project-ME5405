function [Eroded] = Erosion(image,SE,icent,jcent)
% Erosion algorithm (mathematical morphology). 
% Input: Image (binary image), SE (structure element), icent and jcent (origin of the structure element)
% Output: Eroded = Image (-) SE where (-) is the erosion operator

L = size(image,1);
M = size(image,2);

image = cast(image,'int8');
Eroded = image;

% 'Quick and dirty' method to find the background/object in images (works
% for the samples provided: character and chip, but not a general method).
back_ground = mode(image(:));
object = 1 - back_ground;

%%
for ii = 1:L
for jj = 1:M
if (image(ii,jj) == 1)   
    % Run through the structure element
    for kk = 1:size(SE,1)
    for ll = 1:size(SE,2)
        if (SE(kk,ll) ~= 0)
            ipt = ii + (kk - icent);
            jpt = jj + (ll - jcent);
            
            % Erosion step
            if (ipt <= L && jpt <= M && ipt > 0 && jpt > 0)
                if (image(ipt,jpt) == back_ground)
                    Eroded(ii,jj) = back_ground;
                end
            else
                Eroded(ii,jj) = back_ground;
            end
        end
    end
    end
end
end
end

% Set tolerance to ensure binary image output (can remove if this is not
% critical)
tol = 0.5;
dil(abs(Eroded) < tol) = 0;
dil(abs(Eroded) >= tol) = 1;

end

