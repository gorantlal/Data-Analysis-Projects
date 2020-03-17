clear all; close all; clc;

%CAM2_1
load('cam1_2.mat');numFrames1 = size(vidFrames1_2,4);
%creates a spatial filter, initializing

filter = zeros(480,640);
filter(213:405,296:417) = 1;
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));

for j = 1:numFrames1
   
    X = vidFrames1_2(:,:,:,j); 
    %imshow(X); drawnow
    
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid1 = rgb2gray(X); %turns to grayscale
    filt_vid1 = gray_vid1.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid1,0.98);
    %imshow(thresh)
    
    %finds all non-zero vectors
    indeces = find(thresh);
    
    %finds the matrix/vectors
    [Y, X] = ind2sub(size(thresh),indeces);
    
    %finds the centroid!
    mean_x1(j) = mean(X);
    mean_y1(j) = mean(Y);
    
end

%invert the mean_y1 because the axis of a picture is upside down
image_y = 480 - mean_y1;

%CAM2_2
load('cam2_2.mat');numFrames2 = size(vidFrames2_2,4);
%creates a spatial filter, initializing

filter = zeros(480,640);
filter(63:415,181:430) = 1;
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));

for j = 1:numFrames2
   
    X = vidFrames2_2(:,:,:,j); 
    %imshow(X); drawnow
    
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid2 = rgb2gray(X); %turns to grayscale
    filt_vid2 = gray_vid2.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid2,0.98);
    %imshow(thresh)
    
    %finds all non-zero vectors
    indeces = find(thresh);
    
    %finds the matrix/vectors
    [Y, X] = ind2sub(size(thresh),indeces);
    
    %finds the centroid!
    mean_x2(j) = mean(X);
    mean_y2(j) = mean(Y);
   
    
end


%invert the mean_y1 because the axis of a picture is upside down
image_y2 = 480 - mean_y2;



%CAM3_2
load('cam3_2.mat');numFrames3 = size(vidFrames3_2,4);
%creates a spatial filter, initializing

filter = zeros(480,640);
filter(184:338,264:510) = 1;
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));

for j = 1:numFrames3
   
    X = vidFrames3_2(:,:,:,j); 
    %imshow(X); drawnow
    
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid3 = rgb2gray(X); %turns to grayscale
    filt_vid3 = gray_vid3.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid3,0.93);
    %imshow(thresh)
    
    %finds all non-zero vectors
    indeces = find(thresh);
    
    %finds the matrix/vectors
    [Y, X] = ind2sub(size(thresh),indeces);
    
    %finds the centroid!
    mean_x3(j) = mean(X);
    mean_y3(j) = mean(Y);

end
image_y3 = 480-mean_y3;
%invert the mean_y1 because the axis of a picture is upside down
image_x3 = 640 - mean_x3;





plot(1:numFrames1,mean_x1); axis([0 360 0 400]);hold on;
plot(1:numFrames2,mean_x2); hold on; plot(1:numFrames3,400-mean_y3);
legend('1.2','2.2','3.2')
title('Case 2: Unsynced Original Data');
xlabel('Frame Number');
ylabel('XY Displacement (pixels)');




%%
plot(1:numFrames1,image_y); axis([0 360 0 400]); hold on;
plot(1:numFrames2,image_y2); hold on; plot(1:numFrames3,image_x3);
legend('1_1','2_1','3_1')

min1 = 13; min2 = 37; min3 = 16; width = 301; 
% 301
% 319
% 311
shift_mean1 = mean_y1(min1:min1+width);
shift_mean2 = mean_y2(min2:min2+width);
shift_mean3 = mean_x3(min3:min3+width);

set(0, 'DefaultLineLineWidth', 2);

figure();
subplot(2,2,1);
plot(1:302,480 - shift_mean1); hold on;
plot(1:302,480 - shift_mean2); hold on;
plot(1:302,640 - shift_mean3); hold on;
legend('1.2','2.2','3.2')
title('Case 2: Synced Original Data');
xlabel('Frame Number');
ylabel('Z Displacement (pixels)');

concat = [mean_x1(min1:min1+width);shift_mean1;mean_x2(min2:min2+width);...
    shift_mean2;mean_y3(min3:min3+width);shift_mean3];
%mean weights something less, normalize, average out teh center prosition


test1mat = [concat(1,:)-mean(concat(1,:));concat(2,:)-mean(concat(2,:));...
    concat(3,:)-mean(concat(3,:));concat(4,:)-mean(concat(4,:));...
    concat(5,:)-mean(concat(5,:));concat(6,:)-mean(concat(6,:))];


[u,s,v] = svd(test1mat/sqrt(301));

lambda = diag(s).^2; % produce diagonal variances
%Y= test1mat;
Y = u'*test1mat;
%Y = u'*test1mat;
sig= diag(s);
subplot(2,2,[2 4]);
plot (1:6 ,lambda/sum(lambda), 'mo ');
title (" Case 2: Energy Distribution");
xlabel ("Singular Values"); ylabel ("Energy");

subplot(2,2,3);
%plot(1:381,Y(1,:),1:381,Y(2,:),1:381,Y(3,:),1:381,Y(4,:),1:381,Y(5,:),1:381,Y(6,:));
plot(1:302,Y(1,:),1:302,Y(2,:),1:302,Y(3,:))
%,1:302,Y(4,:),1:302,Y(5,:),1:302,Y(6,:));
legend('PCA1','PCA2','PCA3');
title('Case 2: Principal Component Directions');
xlabel('Frame Number');
ylabel('Displacement (pixels)');