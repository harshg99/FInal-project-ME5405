% Outputs the image histogram as a fraction (normalised pdf)

function h = his(image,levels)
h = zeros(1,levels);

for i = 1:size(image,1)
    for j = 1:size(image,2)
        h(image(i,j)+1) = h(image(i,j)+1) + 1;
    end
end
h = h/(size(image,1)*size(image,2));
end