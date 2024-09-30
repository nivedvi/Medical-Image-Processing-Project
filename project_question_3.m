%% Noah Mahagna 314994641 + Ruba Jabaren 208223552


% The explanation of question 3 steps (there is two part for the code first
% part is in the buttom of the code wheich it is to find how much pixels
% per micron and the second part is upove which is every thing asked in the
% question

% 1) reading the image.

% 2) image segmentation: by converted the image to grayscale then determined 
% threshold value to isolate dark nuclei.

% 3) label connected components for identifying and labeling distinct objects 
% in a binary image. 

% 4) Object Analysis: computes properties of image regions (area and centroid) 
% finding the two largest regions (assuming they are the nuclei we want).

% 5) Morphological Operations:by create a new binary image with only these
% two nuclei and createing a structuring element to dilates the nuclei to
% approximate the full cell area.
% 6) Boundary Detection: by using bwboundaries we can find the boundaries.

% 7) Visualization: we displays the original image and used plot to draws 
% the cell boundaries in blue inside loop.

% 8) Measurements: we created preallocate array to calculate diameters and 
% store them in it, by using poly2mask which converts boundary coordinates
% to a binary mask, and then we used pdist2 which computes pairwise distances 
% to find the maximum diameter.

% 9) Calculate the average diameter

clc; clear; format short; format compact

% 1)
img = imread('cells.jpg');
% 2)
gray_img = rgb2gray(img);
nuclei = gray_img < 40; 
% 3)
labeled = logical(nuclei);
% 4)
stats = regionprops(labeled, 'Area', 'Centroid');
[~, idx] = sort([stats.Area], 'descend');
nuclei_idx = idx(1:2);
% 5)
nuclei_img = ismember(labeled, nuclei_idx);
se = strel('disk', 15); 
cells = imdilate(nuclei_img, se);
% 6)
boundaries = bwboundaries(cells);
% 7)
imshow(img);
hold on;
for k = 1:length(boundaries)
    boundary = boundaries{k};
    plot(boundary(:,2), boundary(:,1), 'blue', 'LineWidth', 2);
end

% 8)
diameters_microns = zeros(1, length(boundaries)); 

for k = 1:length(boundaries)
    boundary = boundaries{k};
    mask = poly2mask(boundary(:,2), boundary(:,1), size(img,1), size(img,2));
    [y, x] = find(mask);
    pairwise_distances = pdist2([x, y], [x, y]);
    
    % maximum distance (diameter in pixels)
    diameter_pixels = max(pairwise_distances(:));
    
    % Convert to microns 
    % (based on the scale bar of the ruler  5 pixels = 1 micron)
    % we searched in google and found that 10 pixels = 1 micron, so we used
    % 10 10 pixels per 1 micron ( we hope it is OK)
    diameter_microns = diameter_pixels / 5;
    diameters_microns(k) = diameter_microns;
    fprintf('Cell %d diameter: %.2f microns\n', k, diameter_microns);
end

% 9)
average_diameter = mean(diameters_microns);
fprintf('Average cell diameter: %.2f microns\n', average_diameter);






%%

clc; clear; format short; format compact

% the ruler detection section 

% Read the image
img = imread('cells.jpg');

% Convert to grayscale if it's a color image
if size(img, 3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% Detect the yellow scale bar
% This assumes the scale bar is yellow. Adjust color thresholds as needed.
yellow_mask = img(:,:,1) > 200 & img(:,:,2) > 200 & img(:,:,3) < 100;

% Find the bounding box of the scale bar
stats = regionprops(yellow_mask, 'BoundingBox');
scale_bar = stats(1).BoundingBox;

% Get the length of the scale bar in pixels
scale_bar_length_px = scale_bar(3);

% Known length of the scale bar in microns
scale_bar_length_microns = 5; %  this based on the actual scale

% Calculate pixels per micron
pixels_per_micron =1 / (scale_bar_length_px / scale_bar_length_microns);

% we can use pixels_per_micron to measure other features in the image
% an examplefor the use of it: to measure the diameter of a cell:
% cell_diameter_px = [measure cell diameter in pixels]
% cell_diameter_microns = cell_diameter_px / pixels_per_micron;

% Display the result
disp(['1 micron is equal to ', num2str(pixels_per_micron), ...
    ' pixels in this image.']);






