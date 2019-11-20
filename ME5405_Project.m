
% ME5405: Computer Vision Project

clc;clear;close all;
set(gcf,'color',[0.3 0.3 0.3]);
set(gca,'color','black');
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

hist_chip_R = his(chip(:,:,1),256,0);
hist_chip_G = his(chip(:,:,2),256,0);
hist_chip_B = his(chip(:,:,3),256,0);
hist_charac = his(charac,32,0);
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


% chip(:,:,1) = abs(denoise(chip(:,:,1),'Gaussian',100));
% chip(:,:,2) = abs(denoise(chip(:,:,2),'Gaussian',100));
% chip(:,:,3) = abs(denoise(chip(:,:,3),'Gaussian',100));

% BLP
chip(:,:,1) = abs(denoise(chip(:,:,1),'BLP',100,2));
chip(:,:,2) = abs(denoise(chip(:,:,2),'BLP',100,2));
chip(:,:,3) = abs(denoise(chip(:,:,3),'BLP',100,2));



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

t_charac = threshold(charac,0);
imshow(t_charac,[0 1]);
title(["Thresholded Image 1"]);

t_chip(:,:,1) = threshold(chip(:,:,1),115);
t_chip(:,:,2) = threshold(chip(:,:,2),115);
t_chip(:,:,3) = threshold(chip(:,:,3),115);

% Grayscale version of the rgb chip
R = 0.2990;
G = 0.5870;
B = 0.1140;
gray_chip(:,:) = R*chip(:,:,1) + G*chip(:,:,2) + B*chip(:,:,3);

t_chip_2 = threshold(gray_chip,118);
t_chip_3=t_chip_2;
for(j=1:2)
   t_chip_3=median_filter(t_chip_3,[5 5]);
end

figure();
imshow(t_chip_3,[0 1]);
title(["Thresholded Image 2"]);


%% Segmentation process

% Threshold for segmentation denoising
thres_charac = 10;
thres_chip = 300;

% Background of the images
back_charac = mode(t_charac(:));
back_chip = mode(t_chip_2(:));

% Component labelling of thresholded images

[Labels_charac,num_charac] = CompLabel(t_charac,8,back_charac);
figure();
imshow(Labels_charac+1,[0 num_charac+3]);
title(["Characters Image 1"]);

[Labels_chip,num_chip] = CompLabel(t_chip_3,8,back_chip,300);
figure();
imshow(Labels_chip+2,[0 num_chip+6]);
title(["Characters Image 2"]);



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

figure();
for ii = 1:length(Segment_charac)
    subplot(ceil((length(Segment_charac)/3)),3,ii);
    %colormap(gray(2))
    imshow(Segment_charac{1,ii}); 
end
a=axes;
a.Visible='off';
t=title(["Segmented Characters Image 1 "]);
t.Visible='on';

figure();
for ii = 1:length(Segment_chip)
    subplot(ceil((length(Segment_chip)/3)),3,ii);
    %colormap(gray(2))
    imshow(Segment_chip{1,ii}); 
end
a=axes;
a.Visible='off';
t=title(["Segmented Characters Image 2 "]);
t.Visible='on';

%% Edge detection for thresholded images

t_chip_2_star = Labels_chip>0; %Thresholded chip
t_charac_star = Labels_charac>0; %Thresholded characters

[Edges_charac_img,Xderiv_charac,Yderiv_charac]=EdgeDet(t_charac_star,3,1);
figure();
imshow(Edges_charac_img,[0 max(max(Edges_charac_img))]);
title(["Character Outlines Image 1"]);

[Edges_chip_img,Xderiv_chip,Yderiv_chip]=EdgeDet(t_chip_2_star,3,1);
figure();
imshow(Edges_chip_img,[0 max(max(Edges_chip_img))]);
title(["Character Outlines Image 2"]);

%% Using the thresholded version for each segmented character
for ii = 1:length(Segment_chip)
    tSegment_chip = threshold(Segment_chip{1,ii},115);
    [Edges_chip{1,ii},X_deriv_chip{1,ii},Y_deriv_chip{1,ii}] = EdgeDet(tSegment_chip,3,1); 
end

for jj = 1:length(Segment_charac)
    tSegment_charac = threshold(Segment_charac{1,jj},0);
    [Edges_charac{1,jj},X_deriv_charac{1,jj},Y_deriv_charac{1,jj}] = EdgeDet(tSegment_charac,3,1,3);
end

%TODO


h1=figure();

for ii = 1:length(Edges_charac)
    subplot(ceil(length(Edges_charac)/3),3,ii);
    %colormap(gray(2))
    imshow(Edges_charac{1,ii});
end
a=axes;
a.Visible='off';
t=title(["Edges Characters Image 1 "]);
t.Visible='on';

h2=figure();

for ii = 1:length(Edges_chip)
    subplot(ceil((length(Edges_chip)/3)),3,ii);
    %colormap(gray(2))
    imshow(Edges_chip{1,ii}); 
end
a=axes;
a.Visible='off';
t=title(["Edges Characters Image 2 "]);
t.Visible='on';


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

%TODO
h1=figure();

for ii = 1:length(Rot_charac)
    subplot(ceil(length(Rot_charac)/3),3,ii);
    %colormap(gray(2))
    imshow(Rot_charac{1,ii});
end
a=axes;
a.Visible='off';
t=title(["Rotated Characters Image 1 90 degree CW "]);
t.Visible='on';

h1=figure();

for ii = 1:length(Rot_chip)
    subplot(ceil(length(Rot_chip)/3),3,ii);
    %colormap(gray(2))
    imshow(Rot_chip{1,ii});
end
a=axes;
a.Visible='off';
t=title(["Rotated Characters Image 2 90 degree CW "]);
t.Visible='on';

% Rotation 2
angles = 35; % Rotation angle in degrees (+ve: anticlockwise, -ve: clockwise)

for kk = 1:length(Segment_charac)
    Rot_charac{1,kk} = rotate(Rot_charac{1,kk},angles,"Bilinear",1);
end

for kk = 1:length(Segment_chip)
    Rot_chip{1,kk} = rotate(Rot_chip{1,kk},angles,"Bilinear",1);
end

%TODO
h1=figure();

for ii = 1:length(Rot_charac)
    subplot(ceil(length(Rot_charac)/3),3,ii);
    %colormap(gray(2))
    imshow(Rot_charac{1,ii});
end
a=axes;
a.Visible='off';
t=title(["Rotated Characters Image 1 35 degree CCW from prior problem "]);
t.Visible='on';

h1=figure();

for ii = 1:length(Rot_chip)
    subplot(ceil(length(Rot_chip)/3),3,ii);
    %colormap(gray(2))
    imshow(Rot_chip{1,ii});
end
a=axes;
a.Visible='off';
t=title(["Rotated Characters Image 2 35 degree CCW from prior problem  "]);
t.Visible='on';

%% Skeletonisation of images (~ one pixel thickness version of the images): iterative
iter = 4;

% Run 1
for kk = 1:length(Segment_charac)
    t_Segment_charac = threshold(Segment_charac{1,kk},0);
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

h1=figure();

for ii = 1:length(Skele_charac)
    subplot(ceil(length(Skele_charac)/3),3,ii);
    %colormap(gray(2))
    imshow(Skele_charac{1,ii});
end
a=axes;
a.Visible='off';
t=title(["1 pixel thin image of characters in image 1 "]);
t.Visible='on';

h1=figure();

for ii = 1:length(Skele_chip)
    subplot(ceil(length(Skele_chip)/3),3,ii);
    %colormap(gray(2))
    imshow(Skele_chip{1,ii});
end
a=axes;
a.Visible='off';
t=title(["1 pixel thin image of characters in image 2  "]);
t.Visible='on';

%% Scale and display characters in a given sequence
space = [30, 20];  % Spaces between each characters in x and y coordinates
new_size = [30, 24];  % Pixel sizes of each of the segmented character

% For image 1 character = 1A2B3C
size_char_segment = size(Segment_charac);
no_char = size_char_segment(1,2);
for i=1:no_char
    resized_img = scale(Segment_charac{1,i},new_size,"Bilinear",1);
    indiv_char(:,:,i) = resized_img;
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
    index = space(2)*i + new_size(2)*(i-1);
    Char_sequenced(space(1):space(1)+new_size(1)-1, index:index+new_size(2)-1) = sqnced_char(:,:,i);
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
new_size = [30,24];
for i=1:no_chip
    resized_img = scale(Segment_chip{1,i},new_size,"Bilinear",1);
    indiv_chip(:,:,i) = resized_img;
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
    index = space(2)*i + new_size(2)*(i-1);
    Chip_sequenced(space(1):space(1)+new_size(1)-1, index:index+new_size(2)-1) = sqnced_chip(:,:,i);
end

figure;
colormap(gray(255));
image(Chip_sequenced);
title("Scaled and sequenced character (255 gray scale)");

figure;
colormap(gray(2));
image(Chip_sequenced);
title("Scaled and sequenced character (binary)");

