close all;clc;clear

% a)

% Read the input image
cells = imread('cells.jpg'); % Image is 'cells.jpg'

% ----hsv space k means segmentation:
hsv_cells = rgb2hsv(cells);
labeled_image = imsegkmeans(cells,5);
gray_cells = rgb2gray(cells);
[seg_image,~] = segmentLabels(gray_cells,labeled_image);
% imshow(seg_image)

%----isolating the nucleas:
mask = false(size(seg_image));
mask(labeled_image==1) = true;
imshow(mask)
mask = imopen(mask,strel('disk',3));
mask = imdilate(mask,strel('square',4));
% imshow(mask)


%----active contours to find the cells shape
cont = activecontour(cells,mask,400,'Chan-vese',ContractionBias=-0.85,SmoothFactor=2.8);
cont = imdilate(cont,strel("disk",5));

%----encircling the cells diameter:
circled_cells = DrawEdges(cells,cont,[0,0,255]);

%----calculating the mean diameter using major/minor axis and P/pi
stats = regionprops(cont,'Perimeter','MinorAxisLength','MajorAxisLength');
mean_diameter_1 = mean([stats.Perimeter])/pi;
mean_diameter_2 = mean([stats.MinorAxisLength]);
mean_diameter_3 = mean([stats.MajorAxisLength]);

mean_diameter = mean([mean_diameter_3,mean_diameter_2,mean_diameter_1]);

imshow(circled_cells)


% b)

%----isolating the ruler:
labeled_image = imsegkmeans(cells,4);
[seg_image,levels] = segmentLabels(gray_cells,labeled_image);

%----binarization and boundary extraction:
ruler_mask = false(size(labeled_image));
ruler_mask(seg_image == max(levels)) = 1;
ruler_mask = imopen(ruler_mask,strel("diamond",3));
ruler_stats = regionprops(ruler_mask,'BoundingBox');
box = floor(ruler_stats.BoundingBox);

%----converting the image into a signal:
y_0 = box(2);x_0 = box(1)+10 ;d_y = box(4);d_x = box(3)-15;
ruler = zeros(1,d_x);
for i = 1:d_x
    % subtracting the sum of each column from the height to get the height
    % of each tick
    col = ruler_mask(y_0:d_y+y_0-1,x_0-1+i);
    sum_col = sum(col);
    H_col = d_y - sum_col;
    ruler(i) =  H_col;
   
end
% figure;
% imshow(ruler_mask)
% figure;
% plot(ruler)
% 
% hold on
% 
%----calculating the distance between the ticks:
[pks,locs] = findpeaks(ruler,MinPeakHeight=3*max(ruler/4));
PPm =locs(2)-locs(1);
% scatter(locs,pks)
% hold off
mPP = 1/PPm;  % microns per pixel = 1/pixels per micron

mean_diameter_microns = mPP*mean_diameter;

fprintf('The mean diameter of the cells is %.2f [microns]. (%d pixels per micron)\n',mean_diameter_microns,PPm);


%% Functions


function [seg_image,levels] = segmentLabels(image, labeled_image)

    numValues = max(labeled_image(:));
    levels = round(linspace(0,255,numValues));
    labels = 1:numValues;
    seg_image = image;
    for i = 1:numValues
        seg_image(labeled_image==labels(i)) = levels(i);
    end
end



function overlay = DrawEdges(image, mask, edgeColor)

    % Validate edgeColor input
    if numel(edgeColor) ~= 3
        error('edgeColor must be a 1x3 vector representing RGB values.');
    end
    
    % Ensure edgeColor is in the range [0, 255]
    edgeColor = uint8(edgeColor);
    edgeColor = min(max(edgeColor, 0), 255);
    
    % Create the overlay image
    overlay = image;
    
    % Compute the edges of the mask
    edges = mask - imerode(mask, strel('disk', 3));
    
    % Apply the edge color to the overlay image
    overlay(:,:,1) = overlay(:,:,1) .* uint8(~edges) + edgeColor(1) * uint8(edges);
    overlay(:,:,2) = overlay(:,:,2) .* uint8(~edges) + edgeColor(2) * uint8(edges);
    overlay(:,:,3) = overlay(:,:,3) .* uint8(~edges) + edgeColor(3) * uint8(edges);

end
