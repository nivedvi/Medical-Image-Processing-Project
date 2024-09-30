clc;clear;close all

% a) 
%----Loading the images - original and gray
hema_images = cell(1,8);
bw_images = cell(size(hema_images));
for i = 1:length(hema_images)
    hema_images{1,i} = imread(sprintf('hema%d.jpg',i));
    bw_images{1,i} = rgb2gray(hema_images{1,i});
end

%----Features and boundray extraction
hema_details = cell(4,8);
for i = 1:length(hema_images)
    image = double(bw_images{1,i});
    original = hema_images{1,i};
    [mask,details] = extractHemaKmeans(image,5); % k=5 kmeans segmentation
    hema_details{1,i} = details.Circularity;hema_details{2,i} = details.Solidity;hema_details{3,i} = details.Hu1;...
        hema_details{4,i} = details.Hu2; 
    overlay = DrawEdges(original,mask); % Creating red edges
    subplot(2,4,i); imshow(overlay)
    sub_title = sprintf('Hematoma Image %d',i);
    title(sub_title);
    
end

% b)
%----kmeans classification:
feature_matrix = cell2mat(hema_details);
reduced_features = pca(feature_matrix); % PCA dimensionallity reduction
[idx, C] =  kmeans(reduced_features,2,Replicates=10); % Kmeans classification
groups = unique(idx);
group_1=hema_images(idx==groups(1));
group_2=hema_images(idx==groups(2));
figure;
for i = 1:length(group_1)
    subplot(2,4,i);imshow(group_1{i});title(sprintf('Image %d - Subdural',i))
    subplot(2,4,i+4);imshow(group_2{i});title(sprintf('Image %d - Epidural',i+4))
end
sgtitle('Hematoma Classification')


%% Functions

function [hema,hema_details] = extractHemaKmeans(image,k)
    
% shape extraction and kmeans segmentation
    

    % initializing centroids to ensure repeatable segmentation
    [m,n] = size(image);
    image_reshape = image(:);
    minc = min(image_reshape);
    maxc = max(image_reshape);
    initialcenters = linspace(minc,maxc,k)';
    
    % segmentation
    [idx,~] = kmeans(image_reshape,k,'dist','sqeuclidean','start',initialcenters);
    clustered_image = reshape(idx,[m,n]);

    % extracting the desired segment
    mask = false(size(clustered_image));
    mask(clustered_image == 4) = true;
    mask = imclose(mask,strel("disk",2));
    mask = imopen(mask,strel("disk",2));
    stats = regionprops(mask,'Area','Circularity','Solidity','PixelIdxList');
    [~,idx] = max([stats.Area]); % extracting the hematoma
    pixels = stats(idx).PixelIdxList;
    largestObject = false(size(mask));
    largestObject(pixels) = true;
    largestObject = imfill(largestObject,'holes');
    hema = largestObject;
    hema_details = stats(idx);
    [hema_details.Hu1,hema_details.Hu2]= HuMoments(hema); % Hu moments
end

function [I1,I2] = HuMoments(binaryImage)
    F = logical(binaryImage);
    centroids = regionprops(F,'Centroid').Centroid;
    x_bar = centroids(1);y_bar = centroids(2);
    
    % Mu Calculation
    mu_00 = Mu_Moment(F,0,0,x_bar,y_bar);
    mu_11 =  Mu_Moment(F,1,1,x_bar,y_bar);
    mu_20 =  Mu_Moment(F,2,0,x_bar,y_bar);
    mu_02 =  Mu_Moment(F,0,2,x_bar,y_bar);
    I1 = (mu_02+mu_20)/(mu_00^2);
    I2 = ((mu_20-mu_02)^2 + 4*mu_11^2)/(mu_00)^4;
    
    
end

function Mu_pq = Mu_Moment(I,p,q,x_bar,y_bar)
      [m,n] = size(I);
      Mu_pq = 0;
      for i = 1:m
          for j = 1:n
              Mu_pq = Mu_pq + ((j-x_bar)^p)*((i-y_bar)^q)*I(i,j);
          end
      end
end

function overlay = DrawEdges(image,mask)

    overlay = image;
    edges = mask - imerode(mask,strel('disk',3));
    overlay(:,:,1) = overlay(:,:,1)+255*uint8(edges);
    overlay(:,:,2) = overlay(:,:,2)-255*uint8(edges);
    overlay(:,:,3) = overlay(:,:,3)-255*uint8(edges);


end
