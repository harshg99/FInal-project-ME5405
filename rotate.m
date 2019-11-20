function rotated_image = rotate(image,angles,interpol,constant)
% Rotates an image by angle degrees
% constant == 0: image that rotates outside is cropped out
% constant == 1: retains the image fully, changes matrix domain

image = cast(image,'uint8');
% Set image origin (where origin = image(originY,originX) in array notation
originX = 1; originY = 1;
maxX = size(image,2); maxY = size(image,1);

%% Evaluate the centroidal location
x = 0:1:size(image,2)-1; y = 0:1:size(image,1)-1;
%x = x - originX; y = y - originY;
image2 = double(image);
[X, Y] = meshgrid(x,y); total_mass = sum(image2(:));
X = double(X); Y = double(Y); 
dXmass = image2.*X; dXmass_sum = sum(dXmass(:)); 
dYmass = image2.*Y; dYmass_sum = sum(dYmass(:));

dXcent = round(dXmass_sum/total_mass);
dYcent = round(dYmass_sum/total_mass);

%% Rotation by angle degrees about centroid
L = max([size(image,1), size(image,2)]);

if (constant)
    x = -L:1:L; y = -L:1:L;
else
    x = 1:size(image,2); y = 1:size(image,1);
end

[X, Y] = meshgrid(x,y);

angles = angles*pi/180; % Convert to radians

% New coordinates at which the interpolation takes place
dXX = cos(angles)*X - sin(angles)*Y;
dYY = sin(angles)*X + cos(angles)*Y;

for ii = 1:size(dXX,1)
    for jj = 1:size(dXX,2)
        pointX = dXX(ii,jj); pointY = dYY(ii,jj);
        
        % For now define a 4 point neighbourhood. If bicubic - requires a
        % 16 point neighbourhood
        y1 = floor(pointY + dYcent); I1 = 0;
        y2 = ceil(pointY + dYcent); I2 = 0;
        x1 = floor(pointX + dXcent); I3 = 0;
        x2 = ceil(pointX + dXcent); I4 = 0;
        
        if (y1 > 0 && x1 > 0 && y1 <= maxY && x1 <= maxX)
            I1 = image(y1,x1);
        end
        if (y1 > 0 && x2 > 0 && y1 <= maxY && x2 <= maxX)
            I2 = image(y1,x2);
        end
        if (y2 > 0 && x1 > 0 && y2 <= maxY && x1 <= maxX)
            I3 = image(y2,x1);
        end
        if (y2 > 0 && x2 > 0 && y2 <= maxY && x2 <= maxX)
            I4 = image(y2,x2);
        end

        image_nhood = [I1, I2;
                     I3, I4];
 
        deta = pointX - floor(pointX);
        deps = pointY - floor(pointY); 
        
        switch interpol
            case "Nearest"
                if (deta < 0.5 && deps < 0.5)
                    out_image(ii,jj) = I1;
                elseif (deta < 0.5 && deps >= 0.5)
                    out_image(ii,jj) = I3;
                elseif (deta >= 0.5 && deps < 0.5)
                    out_image(ii,jj) = I2;
                else 
                    out_image(ii,jj) = I4;
                end
            case "Average"
                out_image(ii,jj) = mean(image_nhood(:));
            case "Bilinear"
                out_image(ii,jj) = image_nhood(1,1)*(1-deta)*(1-deps) +...
                    image_nhood(1,2)*(1-deps)*deta + image_nhood(2,1)*(1-deta)*deps+...
                    image_nhood(2,2)*deta*deps;
        end
    end
end

if (constant)
    rotated_image = Crops(out_image);
else
    rotated_image = out_image;
end

end







        

