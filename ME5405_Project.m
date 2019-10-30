% ME5405: Computer Vision Project

clc;clear;close all;

%% Setting up the image matrix

% Read the images into a matrix
chip = imread('charact2.bmp');
fileID = fopen('charact1.txt','r');

% Characters matrix has 32 grayscale levels
charac = fscanf(fileID,'%s',[64,64]);
charac = charac.';

%charac_int = str2num(charac_txt)
check = charac > '9';
notcheck = ~check;
charac = charac - 65.*check - 26.*notcheck;

% Show both images
figure()
image(chip);
title(["Chip"]);

figure()
colormap(gray(32));
image(charac);
title(["Characters"]);

%% Image histogram

hist_chip_R = his(chip(:,:,1),256);
hist_chip_G = his(chip(:,:,2),256);
hist_chip_B = his(chip(:,:,3),256);
hist_charac = his(charac,32);
s1 = 0:1:length(hist_chip_R)-1;
s2 = 0:1:length(hist_charac)-1;

% Plot the frequency histograms
figure()
plot(s1,hist_chip_R);
title(["Histogram for chip: Red"]);
figure()
plot(s1,hist_chip_G);
title(["Histogram for chip: Blue"]);
figure()
plot(s1,hist_chip_B);
title(["Histogram for chip: Green"]);
figure()
plot(s2,hist_charac);
title(["Histogram for characters"]);

%% Denoising process

% Currently based on a Fourier denoising process (can also apply median
% filtering to avoid the coarsening effect).

chip(:,:,1) = abs(denoise(chip(:,:,1),'Gaussian'));
chip(:,:,2) = abs(denoise(chip(:,:,2),'Gaussian'));
chip(:,:,3) = abs(denoise(chip(:,:,3),'Gaussian'));

% Median filtering process
% chip(:,:,1) = median_filter(chip(:,:,1),8);
% chip(:,:,2) = median_filter(chip(:,:,2),8);
% chip(:,:,3) = median_filter(chip(:,:,3),8);


figure()
colormap(gray(255));
image(chip);
title(["Denoised chip"]);

%% Thresholding
t_charac = threshold(charac,16);
t_chip(:,:,1) = threshold(chip(:,:,1),116);
t_chip(:,:,2) = threshold(chip(:,:,2),116);
t_chip(:,:,3) = threshold(chip(:,:,3),116);


% Grayscale version of the rgb chip
R = 0.2989;
G = 0.5870;
B = 0.1140;
gray_chip(:,:) = R*chip(:,:,1) + G*chip(:,:,2) + B*chip(:,:,3);

t_chip_2 = threshold(gray_chip,120);

%% Segmentation process

% Background of the images
back_charac = mode(t_charac(:));
back_chip = mode(t_chip_2(:));

% Component labelling of thresholded images
[Labels_charac] = CompLabel(t_charac,8,back_charac);
[Labels_chip] = CompLabel(t_chip_2,8,back_chip);

% Segmentation
[Segment_charac] = Segment(Labels_charac,charac,10);
[Segment_chip] = Segment(Labels_chip,gray_chip,10);



%% Edge detection for grayscale images

t_chip_2_star = t_chip_2;
charac_star = charac;
num_trials = 1;
for ii = 1:num_trials
    [Edges_chip,Xderiv_chip,Yderiv_chip] = EdgeDet(t_chip_2_star,1,1); % Too much noise with graylevel
    [Edges_charac,Xderiv_charac,Yderiv_charac] = EdgeDet(charac_star,33,1);
    t_chip_2_star = Edges_chip;
    charac_star = Edges_charac;
end


%% Plots

figure()
colormap(gray(2));
image(t_charac);
title(["Thresholded characters"]);

figure()
colormap(gray(2));
image(t_chip);
title(["Thresholded chip: thresholding applied to each of R, G and B"]);

figure()
colormap(gray(2));
image(t_chip_2);
title(["Thresholded chip: thresholding applied to grayscale"]);

for kk = 1:length(Segment_charac)
    figure()
    colormap(gray(32));
    image(Segment_charac{1,kk});
    title(["Segment character"]);
end

for kk = 1:length(Segment_chip)
    figure()
    colormap(gray(255));
    image(Segment_chip{1,kk});
    title(["Segment chip"]);
end

figure()
colormap(gray(32))
image(Edges_charac);
title(["Edges of the character"]);

figure()
colormap(gray(255))
image(Edges_chip);
title(["Edges of the chip"]);









% Image segmentation
% Rotate characters of the image by -90 and 35
% Edge detection of characters
% Skeletonisation of characters
% Scale and display characters in a given sequence

