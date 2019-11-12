function [Skeleton] = Skeletonise(image,N,icorrect)
% Input: image (binary image), N is the number of iterations
% (corresponding to number of pixel depths to remove) and icorrect (1 or 0)
% determines whether image resolution needs to be improved prior to the
% morphological operation.
% Output: the full skeletonised image (binary)

% One application is able to reduce image to a pixel width from a max of 14 pixels to a width
% of max ~ 2 to 3 pixels (for chip image) and gives a one pixel thin image 
% for the characters image. 

% Once 'convergence' is reached, increasing the value of N does not affect the final output 
% since at that point there are no regions with pixel depths > N.

% Skeletonisation here is an application of the Lantuéjoul formula from, 
% (REF: Ch. Lantuéjoul, "Sur le modèle de
% Johnson-Mehl généralisé", Internal report of the Centre de Morph. Math., 
% Fontainebleau, France, 1977.)

% Formulation: if X is a binary image and (-) is an erosion and o is an opening operation,
% nB is a disk structural element with a distance (4/8 connectivity) of n
% units, and B is the unit disk structural element (i.e. nB where n = 1).
% Then S(n) = (X (-) nB) - (X (-) B)oB, and the output skeleton is the
% union of S(n) over all n (n = 1,...,N), i.e. S = U{S(n)} where U is the union symbol.

% Unit disk
B1 = [1 1 1;
    1 1 1;
    1 1 1];

%% Correction for low input resolution

% The correction is applied to the character image since the structure
% element in that case (3 by 3) is too big so image resolution is increased
% to 'fit' in the structure element
new_image = image;
if (icorrect == 1)
    Xarray = linspace(1,size(image,2),3*size(image,2)); 
    Yarray = linspace(1,size(image,1),3*size(image,1));
    [X,Y] = meshgrid(Xarray,Yarray);
    for ii = 1:size(X,1)
        for jj = 1:size(X,2)
            new_image(ii,jj) = image(round(Y(ii,jj)),round(X(ii,jj)));
        end
    end
end

%% 
% For now only works for neighbourhood of 8 (8-connectivity)
for n = 1:N
    B = ones(2*n+1,2*n+1);
    
    % Define the origin
    icent = (size(B,1)+1)/2;
    jcent = (size(B,2)+1)/2;
    
    % S_n = (X (-) nB) - (X (-) B)oB where (-) is an erosion and o is
    % opening. Here B is the structural element for radius = 1 and nB is
    % the structural element for a disk of radius = n. Note that an opening
    % AoB = (A (-) B) (+) B where (+) is a dilation and (-) an erosion.
    
    % Erosion
    Interior = Erosion(new_image,B,icent,jcent);
    
    % Opening
    Open1 = Erosion(Interior,B1,icent,jcent);
    Open2 = Dilate(Open1,B1,icent,jcent);
    
    S{1,n} = Interior - Open2;
    fprintf("Completed %i out of %i\n",n,N);
end

% The final skeletonised image is the union of all skeletons obtained from
% each disk

Skeleton = zeros(size(new_image)); 
Skeleton = cast(Skeleton,"int8");

for n = 1:N
    Sn = S{1,n};
    Skeleton = Skeleton + Sn;
    Skeleton(Skeleton > 1) = 1; %Ensures binary image
end

end



