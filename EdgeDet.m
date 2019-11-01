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

kern = size(G_x,1);

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
Edges(Xderiv1 > thold | Yderiv1 > thold) = 255;

Xderiv(Xderiv1 < thold) = 0;
Yderiv(Yderiv1 < thold) = 0;

end







