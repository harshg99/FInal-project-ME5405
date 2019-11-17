function [Edges,Xderiv,Yderiv] = EdgeDet(imagex,thold,itype,icorrect)

%% Sobel derivative
% Function outputs the matrix of edges used for further analysis
% Cross-correlation with a specific gradient mask

switch itype
    case 1 % Sobel 3 x 3 derivative operator
        G_x = [-1 -2 -1;
            0 0 0;
            1 2 1];
        
        G_y = [-1 0 1;
            -2 0 2;
            -1 0 1];     
end

kern = size(G_x,1);

%% Correction for low input resolution

% The correction is applied to the character image since the structure
% element in that case (3 by 3) is too big so image resolution is increased
% to 'fit' in the structure element
new_image = imagex;
if (icorrect == 1)
    Xarray = linspace(1,size(imagex,2),3*size(imagex,2)); 
    Yarray = linspace(1,size(imagex,1),3*size(imagex,1));
    [X,Y] = meshgrid(Xarray,Yarray);
    for ii = 1:size(X,1)
        for jj = 1:size(X,2)
            new_image(ii,jj) = imagex(round(Y(ii,jj)),round(X(ii,jj)));
        end
    end
end

imagex = new_image;

%% Zero pad with a "bufferd" pixel thickness
bufferd = 2;
t_image = zeros(size(imagex,1)+2*bufferd,size(imagex,2)+2*bufferd);
t_image(bufferd+1:bufferd+size(imagex,1),bufferd+1:bufferd+size(imagex,2)) = imagex;

imagex = zeros(size(t_image,1),size(t_image,2));
imagex = t_image;

%%
% Populate the derivative matrices (\nabla f)
Xderiv = zeros(size(imagex));
Yderiv = zeros(size(imagex));

for ii = 2:size(imagex,1)-1
    for jj = 2:size(imagex,2)-1
        sumX = 0;
        sumY = 0;
        for kk = 1:kern^2
            i_star = floor((kk-1)/kern)+1; i_diff = i_star - (kern - 1);
            j_star = mod((kk-1),kern)+1; j_diff = j_star - (kern - 1);
            sumX = sumX + G_x(i_star,j_star)*imagex(ii-i_diff,jj-j_diff);
            sumY = sumY + G_y(i_star,j_star)*imagex(ii-i_diff,jj-j_diff);
        end
        Xderiv(ii,jj) = sumX;
        Yderiv(ii,jj) = sumY;
    end
end

Xderiv1 = abs(Xderiv);
Yderiv1 = abs(Yderiv);

% Edge detection using threshold
 
Edges = zeros(size(imagex));
%Edges(Xderiv1 > thold | Yderiv1 > thold) = 255;
Edges(sqrt(Xderiv1.^2 + Yderiv1.^2) > thold) = 255;

Xderiv(Xderiv1 < thold) = 0;
Yderiv(Yderiv1 < thold) = 0;

end
