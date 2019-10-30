% Median filtering of an image
function m_image = median_filter(imagex,nhood)

% nhood is the neighbourhood definition (4 or 8)
L = size(imagex,1);
M = size(imagex,2);

m_image = imagex;

% Real-time updates
for ii = 2:L-1
    for jj = 2:M-1
        
        switch nhood
            case 4
                NH = [imagex(ii,jj-1),imagex(ii,jj+1);
                    imagex(ii-1,jj),imagex(ii+1,jj)];
            case 8
                NH = [imagex(ii-1,jj-1),imagex(ii-1,jj),imagex(ii-1,jj+1);
                    imagex(ii,jj-1),imagex(ii,jj),imagex(ii,jj+1);
                    imagex(ii+1,jj-1),imagex(ii,jj),imagex(ii,jj+1)];
        end
        m_image(ii,jj) = median(NH(:));
    end
end



