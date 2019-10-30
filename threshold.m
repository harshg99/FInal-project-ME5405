% Thresholding an image
function t_image = threshold(imagex,thold)

t_image = zeros(size(imagex,1),size(imagex,2));

t_image = (imagex > thold);
end