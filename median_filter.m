% Median filtering of an image
function m_image = median_filter(imagex,mask)

% INPUT
%
% imagex: grayscale image in 0 to 255
% mask            int (only 4 or 8)
%                 [int,int] specifies m by n matrix across which median is
%                 obtained
% OUTPUT
% m_image: median filtered image output (same type as imagex);



% INPUT
% nhood: [int8] Neighbourhood definition (4 or 8)
% imagex: grayscale image in 0 to 255

% OUTPUT
% m_image: median filtered image output (same type as imagex);

% nhood is the neighbourhood definition (4 or 8)
L = size(imagex,1);
M = size(imagex,2);

m_image = imagex;
if(length(mask)==1)
    nhood=mask;
    i_s=2;
    j_s=2;
    i_f=L-1;
    j_f=M-1;
elseif(length(mask)==2)
    nhood=0;
    %image padding boundaries
    m=mask(1);
    n=mask(2);
    bound_l=floor((n-1)/2);
    bound_r=ceil((n-1)/2);
    bound_u=floor((m-1)/2);
    bound_d=ceil((n-1)/2);
    
    %size_i1=size(imagex);
    %img=[zeros(:,1:bound_l),img,zeros(:,1:bound_r)];
    %img=[zeros(1:bound_u,:);img;zeros(1:bound_d,:)];
    i_s=1+bound_l;
    j_s=1+bound_u;
    i_f=L-bound_r;
    j_f=M-bound_d;
end
% Real-time updates
for ii = i_s:i_f
    for jj = j_s:j_f
        
        switch (nhood)
            case 4
                NH = [imagex(ii,jj-1),imagex(ii,jj+1);
                    imagex(ii-1,jj),imagex(ii+1,jj)];
            case 8
                NH = [imagex(ii-1,jj-1),imagex(ii-1,jj),imagex(ii-1,jj+1);
                    imagex(ii,jj-1),imagex(ii,jj),imagex(ii,jj+1);
                    imagex(ii+1,jj-1),imagex(ii,jj),imagex(ii,jj+1)];
            otherwise
                
        end
        if(length(mask)==2)
            NH = imagex((ii-bound_l):(ii+bound_r),(jj-bound_u):(jj+bound_d));
        end
            m_image(ii,jj) = median(NH(:));
    end
end
