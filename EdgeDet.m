
function [Edges,Xderiv,Yderiv] = EdgeDet(imagex,thold,itype)
% Function outputs the matrix of edges used for further analysis (Edges)

% INPUT
% Xderiv and Yderiv are gradients associated with the edges
% imagex: input image (grayscale: [uint8]);
% thold: threshold to determine edges
% itype: determines the type of derivative operator: Sobel (for now)

% OUTPUT
% Edges: Edges defined by magnitude of gradient greater than a certain thold (threshold) => grayscale
% Xderiv, Yderiv: derivatives w.r.t x and y matrices constructed separately.

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


Xderiv = mask_image(imagex,G_x);
Yderiv = mask_image(imagex,G_y);

Xderiv1 = abs(Xderiv);
Yderiv1 = abs(Yderiv);

% Edge detection using threshold
 
Edges = zeros(size(imagex));
%Edges(Xderiv1 > thold | Yderiv1 > thold) = 255;

Edges(sqrt(Xderiv1.^2 + Yderiv1.^2) > thold) = 255;

Xderiv(Xderiv1 < thold) = 0;
Yderiv(Yderiv1 < thold) = 0;

end
