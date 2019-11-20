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
    Rot_charac{1,kk} = rotate(Rot_charac{1,kk},angles,"Nearest",1);
end

for kk = 1:length(Segment_chip)
    Rot_chip{1,kk} = rotate(Rot_chip{1,kk},angles,"Bilinear",1);
end

% Net rotation for verification
angles = -55; % Rotation angle in degrees (+ve: anticlockwise, -ve: clockwise)

for kk = 1:length(Segment_charac)
    Rot_charac2{1,kk} = rotate(Segment_charac{1,kk},angles,"Bilinear",1);
end

for kk = 1:length(Segment_chip)
    Rot_chip2{1,kk} = rotate(Segment_chip{1,kk},angles,"Bilinear",1);
end
