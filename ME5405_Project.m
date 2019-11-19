% ME5405: Computer Vision Project

clc;clear;close all;

%% Setting up the image matrix

% Read the images into a matrix
chip = imread('charact2.bmp');
fileID = fopen('charact1.txt','r');

% Characters matrix has 32 grayscale levels
charac = fscanf(fileID,'%s',[64,64]);
X = charac.';
charac = charac.';

%charac_int = str2num(charac_txt)
check = charac > '9';
notcheck = ~check;
charac = charac - 55.*check - 48.*notcheck;

% Show both images
figure()
imshow(chip);
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
% figure()
% plot(s1,hist_chip_R);
% title(["Histogram for chip: Red"]);
% figure()
% plot(s1,hist_chip_G);
% title(["Histogram for chip: Blue"]);
% figure()
% plot(s1,hist_chip_B);
% title(["Histogram for chip: Green"]);
% figure()
% plot(s2,hist_charac);
% title(["Histogram for characters"]);

% Grayscale version of the rgb chip
R = 0.2989;
G = 0.5870;
B = 0.1140;
gray_chip(:,:) = R*chip(:,:,1) + G*chip(:,:,2) + B*chip(:,:,3);

%% Denoising process

% Currently based on a Fourier denoising process (can also apply median
% filtering to avoid the coarsening effect).

chip(:,:,1) = abs(denoise(chip(:,:,1),'Gaussian'));
chip(:,:,2) = abs(denoise(chip(:,:,2),'Gaussian'));
chip(:,:,3) = abs(denoise(chip(:,:,3),'Gaussian'));
gray_chip = abs(denoise(gray_chip,'Gaussian'));


% % Median filtering process
% chip(:,:,1) = median_filter(chip(:,:,1),8);
% chip(:,:,2) = median_filter(chip(:,:,2),8);
% chip(:,:,3) = median_filter(chip(:,:,3),8);


figure()
colormap(gray(255));
image(chip);
title(["Denoised chip"]);

figure()
imshow(gray_chip,colormap(gray(255)));
title(["Denoised chip"]);

%% Contrast modification
% chip(:,:,1) = imadjust(chip(:,:,1));
% chip(:,:,2) = imadjust(chip(:,:,2));
% chip(:,:,3) = imadjust(chip(:,:,3));


%% Thresholding

t_charac = threshold(charac,16);
t_chip(:,:,1) = threshold(chip(:,:,1),115);
t_chip(:,:,2) = threshold(chip(:,:,2),115);
t_chip(:,:,3) = threshold(chip(:,:,3),115);

% Gray-scale thresholded chip used for further analysis
t_chip_2 = threshold(gray_chip,115);


%% Segmentation process

% Threshold for segmentation denoising
thres_charac = 10;
thres_chip = 300;

% Background of the images
back_charac = mode(t_charac(:));
back_chip = mode(t_chip_2(:));

% Component labelling of thresholded images

% Old version
% [Labels_charac] = CompLabel(t_charac,8,back_charac);
% [Labels_chip] = CompLabel(t_chip_2,8,back_chip);

[Labels_charac,~] = CompLabel(t_charac,8,0,thres_charac);
[Labels_chip,~] = CompLabel(t_chip_2,8,0,thres_chip);

% Segmentation
[Segment_charac] = Segment(Labels_charac,charac,thres_charac);
[Segment_chip] = Segment(Labels_chip,gray_chip,thres_chip);

%% Separating characters for the chip

buffer = 10; thold = 115;
for kk = 1:length(Segment_chip)
    col(kk) = size(Segment_chip{1,kk},2);
end

tot = length(Segment_chip);

for ii = 1:tot
    r = size(Segment_chip{1,ii},2);
    if (abs(r - max(col(:))) < abs(r - min(col(:)))) % Multiple characters in one segment (here assumes 2 characters)
        L = size(Segment_chip{1,ii},1);
        W = size(Segment_chip{1,ii},2);
        mid = floor(W/2);
        Sthold = threshold(Segment_chip{1,ii},thold);
        Slice = zeros(L,W);
        Slice(1:L,mid-buffer:mid+buffer) = Sthold(1:L,mid-buffer:mid+buffer); 
        count = 0;
        for ll = 1:50
            count = count + 1;
            Slice = Erosion(Slice,[1;1;1],2,1);
            if (~max(Slice(:,mid))) break; end
        end
        for ll = 1:count-1
             Slice = Dilate(Slice,[1;1;1],2,1);
        end
        Slice = Erosion(Slice,[1 1 1],1,2); %Trim edges
        Slice(1:L,1:mid-buffer) = 1; Slice(1:L,mid+buffer:W) = 1;
        Sthold = logical(Sthold).*logical(Slice);
%        Top = Sthold(1,:); Left = Sthold(:,1); Bottom = Sthold(end,:); Right = Sthold(:,end);
        Sthold(1,:) = 0; Sthold(end,:) = 0; Sthold(:,1) = 0; Sthold(:,end) = 0;
        [Labels_segment_chip] = CompLabel(Sthold,8,0,thres_chip);
        [Segment_segment_chip] = Segment(Labels_segment_chip,Segment_chip{1,ii},thres_chip);
        Segment_chip{1,ii} = Segment_segment_chip{1,1};
        Segment_chip{1,length(Segment_chip)+1} = Segment_segment_chip{1,2};
    end
end

%% Edge detection for thresholded images

t_chip_2_star = t_chip_2; %Thresholded chip
t_charac_star = t_charac; %Thresholded characters


% Using the thresholded version
for ii = 1:length(Segment_chip)
    tSegment_chip = threshold(Segment_chip{1,ii},115);
    [Edges_chip{1,ii},X_deriv_chip{1,ii},Y_deriv_chip{1,ii}] = EdgeDet(tSegment_chip,3,1,0); 
end

for jj = 1:length(Segment_charac)
    tSegment_charac = threshold(Segment_charac{1,jj},16);
    [Edges_charac{1,jj},X_deriv_charac{1,jj},Y_deriv_charac{1,jj}] = EdgeDet(tSegment_charac,3,1,1);
end



%% Rotation of segmented images
Rot_charac{1,1} = rotate(Segment_charac{1,1},90,"Bilinear",1);
Rot_charac = cell(1,length(Segment_charac));
Rot_chip = cell(1,length(Segment_chip));

% Rotation 1
angles = -90; % Rotation angle in degrees (+ve: anticlockwise, -ve: clockwise)

for kk = 1:length(Segment_charac)
    Rot_charac{1,kk} = rotate(Segment_charac{1,kk},angles,"Bilinear",1);
end

for kk = 1:length(Segment_chip)
    Rot_chip{1,kk} = rotate(Segment_chip{1,kk},angles,"Bilinear",1);
end

% Rotation 2
angles = 35; % Rotation angle in degrees (+ve: anticlockwise, -ve: clockwise)

for kk = 1:length(Segment_charac)
    Rot_charac{1,kk} = rotate(Rot_charac{1,kk},angles,"Bilinear",1);
end

for kk = 1:length(Segment_chip)
    Rot_chip{1,kk} = rotate(Rot_chip{1,kk},angles,"Bilinear",1);
end

%% Skeletonisation of images (~ one pixel thickness version of the images): iterative
iter = 4;

% Run 1
for kk = 1:length(Segment_charac)
    t_Segment_charac = threshold(Segment_charac{1,kk},16);
    Skele_charac{1,kk} = Skeletonise(t_Segment_charac,20,0,1);
end

for kk = 1:length(Segment_chip)
    t_Segment_chip = threshold(Segment_chip{1,kk},115);
    Skele_chip{1,kk} = Skeletonise(t_Segment_chip,20,0,1);
end

% Additional runs to get to one pixel
iter = iter - 1;
for ii = 1:iter
    for kk = 1:length(Segment_charac)
        Skele_charac{1,kk} = Skeletonise(Skele_charac{1,kk},20,0,0);
    end
    for kk = 1:length(Segment_chip)
        Skele_chip{1,kk} = Skeletonise(Skele_chip{1,kk},20,0,0);
    end
end

%% Plots

%% Thresholded images
figure()
colormap(gray(2));
image(t_charac);
title(["Thresholded characters"]);

% figure()
% colormap(gray(2));
% image(t_chip);
% title(["Thresholded chip: thresholding applied to each of R, G and B"]);

figure()
colormap(gray(2));
image(t_chip_2);
title(["Thresholded chip: thresholding applied to grayscale"]);

%% Segmented images

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

%% Edge detection

for ii = 1:length(Edges_chip)
    figure()
    colormap(gray(2))
    image(Edges_chip{1,ii});
    title(["Edge chip"]);
end

for ii = 1:length(Edges_charac)
    figure()
    colormap(gray(2))
    image(Edges_charac{1,ii});
    title(["Edge character"]);
end

%% Rotated segments
for kk = 1:length(Rot_charac)
    figure()
    colormap(gray(32));
    image(Rot_charac{1,kk});
    title(["Segment character: Rotated"]);
end

for kk = 1:length(Rot_chip)
    figure()
    colormap(gray(255));
    image(Rot_chip{1,kk});
    title(["Segment chip: Rotated"]);
end

%% Skeletonised images

for kk = 1:length(Skele_charac)
    figure()
    colormap(gray(2));
    image(Skele_charac{1,kk});
    title(["Segment character (Skeletonised)"]);
end

for kk = 1:length(Skele_chip)
    figure()
    colormap(gray(2));
    image(Skele_chip{1,kk});
    title(["Segment chip (Skeletonised)"]);
end

%% Scale and display characters in a given sequence
space = [30, 20];  % Spaces between each characters in x and y coordinates
new_size = [30, 24];  % Pixel sizes of each of the segmented character

% For image 1 character = 1A2B3C
size_char_segment = size(Segment_charac);
no_char = size_char_segment(1,2);
for i=1:no_char
    resized_img = scale(Segment_charac{1,i},new_size,"Bilinear",1)
    indiv_char(:,:,i) = resized_img
end
% Arrange charracter in the correct sequence/order
sqnced_char(:,:,1) = indiv_char(:,:,4);  %1
sqnced_char(:,:,2) = indiv_char(:,:,1);  %A
sqnced_char(:,:,3) = indiv_char(:,:,5);  %2
sqnced_char(:,:,4) = indiv_char(:,:,2);  %B
sqnced_char(:,:,5) = indiv_char(:,:,6);  %3
sqnced_char(:,:,6) = indiv_char(:,:,3);  %C

% Create zeros to place the sequence character: zero padding
Char_sequenced = zeros(2*space(1)+new_size(1), (no_char+1)*space(2)+no_char*new_size(2));
for i=1:no_char
    index = space(2)*i + new_size(2)*(i-1)
    Char_sequenced(space(1):space(1)+new_size(1)-1, index:index+new_size(2)-1) = sqnced_char(:,:,i)
end

figure;
colormap(gray(32));
image(Char_sequenced);
title("Scaled and sequenced character (32 gray scale)");

figure;
colormap(gray(2));
image(Char_sequenced);
title("Scaled and sequenced character (binary)");

% For image 2 = chip 7M2HD44780A00
size_chip_segment = size(Segment_chip);
no_chip = size_chip_segment(1,2);
new_size = [30,24]
for i=1:no_chip
    resized_img = scale(Segment_chip{1,i},new_size,"Bilinear",1);
    indiv_chip(:,:,i) = resized_img
end
% Arrange chip in the correct sequence/order
sqnced_chip(:,:,1) = indiv_chip(:,:,1);    %7
sqnced_chip(:,:,2) = indiv_chip(:,:,10);   %M
sqnced_chip(:,:,3) = indiv_chip(:,:,2);    %2
sqnced_chip(:,:,4) = indiv_chip(:,:,9);    %H
sqnced_chip(:,:,5) = indiv_chip(:,:,5);    %D
sqnced_chip(:,:,6) = indiv_chip(:,:,7);    %4
sqnced_chip(:,:,7) = indiv_chip(:,:,12);   %4
sqnced_chip(:,:,8) = indiv_chip(:,:,4);    %7
sqnced_chip(:,:,9) = indiv_chip(:,:,6);    %8
sqnced_chip(:,:,10) = indiv_chip(:,:,3);   %0
sqnced_chip(:,:,11) = indiv_chip(:,:,8);   %A
sqnced_chip(:,:,12) = indiv_chip(:,:,11);  %0
sqnced_chip(:,:,13) = indiv_chip(:,:,13);  %0

% Create zeros to place the sequence chip : zero padding
Chip_sequenced = zeros(2*space(1)+new_size(1), (no_chip+1)*space(2)+no_chip*new_size(2));
for i=1:no_chip
    index = space(2)*i + new_size(2)*(i-1)
    Chip_sequenced(space(1):space(1)+new_size(1)-1, index:index+new_size(2)-1) = sqnced_chip(:,:,i)
end

figure;
colormap(gray(255));
image(Chip_sequenced);
title("Scaled and sequenced character (255 gray scale)");

figure;
colormap(gray(2));
image(Chip_sequenced);
title("Scaled and sequenced character (binary)");