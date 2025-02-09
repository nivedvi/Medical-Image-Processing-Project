clear; close all; clc

% a)

%----video read:
vid = VideoReader('ball.mp4');
ball_vid = read(vid);

%----frame enhancement:
frame1 = ball_vid(:,:,:,1);
frame1_en = HSVenhance(frame1);
diff_image = frame1_en(:,:,1)-frame1_en(:,:,3);

%----binarization:
bin_image = imopen(imbinarize(diff_image),strel('diamond',5));

stats = regionprops(bin_image,'Circularity','PixelIdxList');
[~,idx] = max([stats.Circularity]);
bin_mask = false(size(bin_image));
bin_mask(stats(idx).PixelIdxList) = true;
cont1 = bin_mask;

%----active contures object tracking:
[m,n,~,k] = size(ball_vid);
conbias = 0.12;smoothfact = 0.1;
bin_vid = zeros(m,n,k);
centroids = zeros(k,2);
ball_vid_detect = ball_vid;

for i = 1:k
    frame = ball_vid(:,:,:,i);
    [~,mask_frame] = contours(frame,cont1,conbias,smoothfact);
    bin_vid(:,:,i) = mask_frame;
    [area,cent] = extractCentroid(mask_frame);% centroids extraction for part b
    if i == 1
        box_area = area*4;
    end
    box = createBox(cent,[m,n],box_area);
    centroids(i,:) = cent;
    frame_detect = DrawEdges(frame,box,[0,255,0]);
    ball_vid_detect(:,:,:,i) = frame_detect;
    cont1 = imerode(mask_frame,strel('disk',3));
end

implay(ball_vid_detect);

% b) 

%----calculating the frequency of the bounces:
y_position = centroids(:,2);
time = vid.Duration;
[pks,locs] = findpeaks(y_position);
n_pks = length(pks);
frequency = n_pks*60/(time);
fprintf("The frequency of the ball is %.2f bounces per minute\n",frequency);
%% Functions

function binaryMask = createBox(centroid, imageSize, area)
    % Fixed line width
    lineWidth = 2;
    
    % Extract the centroid coordinates
    cx = centroid(1);
    cy = centroid(2);
    
    % Calculate the width and height of the bounding box assuming it's square
    sideLength = sqrt(area);
    halfSideLength = sideLength / 2;
    
    % Calculate the bounding box coordinates
    x = cx - halfSideLength;
    y = cy - halfSideLength;
    width = sideLength;
    height = sideLength;
    
    % Initialize the binary mask with the specified size
    m = imageSize(1);
    n = imageSize(2);
    binaryMask = false(m, n);
    
    % Calculate the indices for the bounding box edges
    xIndices = round(x):(round(x) + round(width) - 1);
    yIndices = round(y):(round(y) + round(height) - 1);
    
    % Ensure the indices are within the image bounds
    xIndices = max(1, min(n, xIndices));
    yIndices = max(1, min(m, yIndices));
    
    % Set the pixels for the edges of the bounding box with fixed line width
    for i = 0:lineWidth-1
        % Top edge
        binaryMask(max(1, yIndices(1)-i), xIndices) = true;
        binaryMask(min(m, yIndices(1)+i), xIndices) = true;
        
        % Bottom edge
        binaryMask(max(1, yIndices(end)-i), xIndices) = true;
        binaryMask(min(m, yIndices(end)+i), xIndices) = true;
        
        % Left edge
        binaryMask(yIndices, max(1, xIndices(1)-i)) = true;
        binaryMask(yIndices, min(n, xIndices(1)+i)) = true;
        
        % Right edge
        binaryMask(yIndices, max(1, xIndices(end)-i)) = true;
        binaryMask(yIndices, min(n, xIndices(end)+i)) = true;
    end

end

function [area,centroid] = extractCentroid(bin_img)
    stats = regionprops(bin_img,'Centroid','Area');
    [area,idx] = max([stats.Area]);
    centroid = stats(idx).Centroid;
end

function [cont,mask] = contours(image,first_mask,conbias,smoothfact)
  
    cont = activecontour(image,first_mask,300,"Chan-vese","ContractionBias",conbias,SmoothFactor=smoothfact);
    mask = imerode(cont,strel("disk",4));
    cont = cont-imerode(cont,strel("disk",3));
    cont = cat(3,cont,cont,cont);

end

function enhanced_image = HSVenhance(image)
    hsv_image = rgb2hsv(image);
    hsv_image(:,:,3) = histeq(hsv_image(:,:,3));
    enhanced_image = hsv2rgb(hsv_image);
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
