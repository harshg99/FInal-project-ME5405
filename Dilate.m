function [dil] = Dilate(image,SE,icent,jcent)
% Dilate (mathematical morphology) an image through the structure element with center as origin 
% Input: Image (binary image), SE (structure element), icent and jcent (origin of the structure element)
% Output: Eroded = Image (+) SE where (+) is the erosion operator

L = size(image,1);
M = size(image,2);

image = cast(image,'int8');
dil = image;
back_ground = mode(image(:));
object = 1 - back_ground;

%%
for ii = 1:size(image,1)
for jj = 1:size(image,2)
    if (image(ii,jj) == 1)
        for kk = 1:size(SE,1)
        for ll = 1:size(SE,2)
            if (SE(kk,ll) ~= 0)
                ipt = ii + kk - icent;
                jpt = jj + ll - jcent;

                % Dilation step
                if (ipt <= L && jpt <= M && ipt > 0 && jpt > 0)
                    dil(ipt,jpt) = object;
                end
            end
        end
        end
    end
end
end

% Set tolerance to ensure binary image output
tol = 0.5;
dil(abs(dil) < tol) = 0;
dil(abs(dil) >= tol) = 1;

end

