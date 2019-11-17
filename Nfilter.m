% Nonlinear filter
function [imagex] = Nfilter(image,SE,icent,jcent,ctype)
imagex = image;
for ii = icent:size(image,1)-(icent-1)
    for jj = jcent:size(image,2)-(jcent-1)
        count = 1;
        nhood(count) = image(ii,jj);
        for kk =1:size(SE,1)
            for ll =1:size(SE,2)
                if (SE(kk,ll))
                    count = count + 1;
                    ipt = ii + kk - icent;
                    jpt = jj + ll - jcent;
                    nhood(count) = image(ipt,jpt);
                end
            end
        end
        
        switch ctype
            case "min"
                imagex(ii,jj) = min(nhood(:));
            case "max"
                imagex(ii,jj) = max(nhood(:));
            case "median"
                imagex(ii,jj) = median(nhood(:));
        end
        
    end
end

end
