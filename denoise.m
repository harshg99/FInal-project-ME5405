% De-noising an image
function n_image = denoise(imagex,ctype)

N = size(imagex,1);
M = size(imagex,2);

% Zero padding
pad_image = zeros(2*N,2*M);

for i = 1:size(imagex,1)
    for j = 1:size(imagex,2)
        pad_image(i,j) = imagex(i,j);
    end
end


% Convert to frequency domain
fftimage = fftshift(fft2(pad_image));

% Euclidean distance
cent = [N;M];
dist = zeros(2*N,2*M);
for i = 1:size(fftimage,1)
    for j = 1:size(fftimage,2)
        dist(i,j) = sqrt((i-cent(1))^2 + (j - cent(2))^2);
    end
end

% Create a (low-pass) filter over the dist matrix

filter = zeros(2*N,2*M);
k = 4.25;
var = (min([N M])/(k)).^2;
stdev = sqrt(var);

switch(ctype)
    case 'Ideal'
        filter = 1*(dist <= stdev);
    case 'Gaussian'
        filter = exp(-(dist.^2)/(2*var));
end

% Convolution
fftimage_filter = fftimage.*filter;

% Inverse FFT
n_image_pad = ifft2(ifftshift(fftimage_filter));

% Recover original (noise-free + blurred)
for i = 1:size(imagex,1)
    for j = 1:size(imagex,2)
        n_image(i,j) = n_image_pad(i,j);
    end
end

