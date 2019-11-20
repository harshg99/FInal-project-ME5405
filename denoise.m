% De-noising an image
function n_image = denoise(imagex,ctype,varargin)

% --------------------------------
% Denoises image using fourier transform
%--------------------
% Inputs:
% imagex    [uint8]    image in uint8
% ctype     [string]   specifies filter type (see filter type section)
% varargin  [      ]   depends on the filter type used(see filter type section)            
%--------------------
% Outputs:
% t_image   [unit8] binary image


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
        dist(i,j) = sqrt((i-cent(1))^2 + (N/M*(j - cent(2))^2));
    end
end

% Create a (low-pass) filter over the dist matrix

filter = zeros(2*N,2*M);

%var = (min([N M])/(4.25)).^2;
%stdev = sqrt(var);
%% Section: Filter Type
switch(length(varargin))
    case 0
    var=50;
    n=2;
    case 1
    var=varargin{1};
    n=2;
    case 2
    var =varargin{1};
    n=varargin{2}; 
end

switch(ctype)
    case 'Ideal'
        % varargin={'var'} var is variance
        filter = 1*(dist <= var);
    case 'Gaussian'
        % varargin={'var'} var is variance of gauss distribution
        filter = exp(-(dist.^2)/(2*var^2));
    case 'BLP'
        %varargin={'var'','n''}
        %var is cutoff frequencny
        %n is exponenent for Butterworth Low Pass
        filter = 1./(1+((dist/var).^(2*n)));
    case 'BRF'
    case 'BPF'
end
%%
%imshow(mat2gray(filter)*255,[0 255]);
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
n_image=uint8(abs(n_image));


