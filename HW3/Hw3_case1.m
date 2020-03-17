clear all; close all; clc; 

%CAM1_1
load('cam1_1.mat');numFrames1 = size(vidFrames1_1,4);
%creates a spatial filter, initializing the mean vectors
filter = zeros(480,640);
filter(200:430,300:400) = 1;
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));
 

set(0, 'DefaultLineLineWidth', 2);

for j = 1:numFrames1
    
    X = vidFrames1_1(:,:,:,j); 
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid1 = rgb2gray(X); %turns to grayscale
    filt_vid1 = gray_vid1.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid1,0.97);
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

% CAM2_1
load('cam2_1.mat'); numFrames2 = size(vidFrames2_1,4);
filter = zeros(480,640); 
filter(100:385,260:340) = 1;
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));
for j = 1:numFrames2
   
    X = vidFrames2_1(:,:,:,j); 
    filter_uint8 = uint8(filter);
    gray_vid2 = rgb2gray(X); 
    filt_vid2 = gray_vid2.*filter_uint8;
    
    %thresh = filt_vid2 > 250; %0.97
    %imshow(thresh); drawnow
    thresh = imbinarize(filt_vid2,0.97);
    indeces2 = find(thresh);
    [Y2, X2] = ind2sub(size(thresh),indeces2);
    mean_x2(j) = mean(X2);
    mean_y2(j) = mean(Y2);
    
end
image_y2 = 480-mean_y2;

% CAM3_1
load('cam3_1.mat'); numFrames3 = size(vidFrames3_1,4);
filter = zeros(480,640); 
filter(230:330,270:485) = 1; %manually found the filter
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));

for j = 1:numFrames3
   
    X = vidFrames3_1(:,:,:,j);
    %imshow(X); drawnow
    filter_uint8 = uint8(filter);
    gray_vid3 = rgb2gray(X); 
    filt_vid3 = gray_vid3.*filter_uint8; 
    thresh = imbinarize(filt_vid3,0.95);
    %thresh = filt_vid3 > 240; %0.95
    %imshow(thresh); drawnow
    indeces3 = find(thresh);
    [Y3, X3] = ind2sub(size(thresh),indeces3);
     mean_x3(j) = mean(X3);
     mean_y3(j) = mean(Y3);
end
%opposite axis
image_y3 = 480-mean_y3;


plot(1:numFrames1,mean_x1); axis([0 250 0 350]);hold on;
plot(1:numFrames2,mean_x2); hold on; plot(1:numFrames3,480-mean_y3);
legend('1.1','2.1','3.1')
title('Case 1: Unsynced Original Data');
xlabel('Frame Number');
ylabel('XY Displacement (pixels)');





%%
figure(1); %time!
image_x3 = 640-mean_x3;
plot(1:numFrames1,image_y); axis([0 numFrames1 0 350]); hold on;
plot(1:numFrames2,image_y2); hold on; plot(1:numFrames3,image_x3);
legend('1.1','2.1','3.1')

%crop the videos to match
min1 = 30; min2 = 39; min3 = 30; width = 196; 
shift_mean1 = mean_y1(min1:min1+width);
shift_mean2 = mean_y2(min2:min2+width);
shift_mean3 = mean_x3(min3:min3+width);

figure(2);
subplot(2,2,1)
plot(1:197,480 - shift_mean1); hold on;
plot(1:197,480 - shift_mean2); hold on;
plot(1:197,640 - shift_mean3); hold on;
legend('1.1','2.1','3.1')
title('Case 1: Synced Original Data');
xlabel('Frame Number');
ylabel('Z Displacement (pixels)');

%,shift_mean2,shift_mean3);
% figure(2); %positions
% subplot(3,1,1); plot(mean_x1, image_y,'.'); axis([0 640 0 480]);
% subplot(3,1,2); plot(mean_x2, image_y2,'.'); axis([0 640 0 480]);
% subplot(3,1,3); plot(mean_x3, image_y3,'.'); axis([0 640 0 480]);

concat = [mean_x1(min1:min1+width);shift_mean1;mean_x2(min2:min2+width);...
    shift_mean2;mean_x3(min3:min3+width);shift_mean3];
%mean weights something less, normalize, average out teh center prosition


test1mat = [concat(1,:)-mean(concat(1,:));concat(2,:)-mean(concat(2,:));...
    concat(3,:)-mean(concat(3,:));concat(4,:)-mean(concat(4,:));...
    concat(5,:)-mean(concat(5,:));concat(6,:)-mean(concat(6,:))];

[u,s,v] = svd(test1mat/sqrt(196));

lambda = diag(s).^2; % produce diagonal variances
%Y= test1mat;
Y = u'*test1mat;
%Y = u'*test1mat;
sig= diag(s);
subplot(2,2,[2 4])
plot (1:6 ,lambda/sum(lambda), 'mo ');
title (" Case 1: Energy Distribution");
xlabel ("Singular Values"); ylabel ("Energy");

subplot(2,2,3);
%plot(1:381,Y(1,:),1:381,Y(2,:),1:381,Y(3,:),1:381,Y(4,:),1:381,Y(5,:),1:381,Y(6,:));
plot(1:197,Y(1,:),1:197,Y(2,:),1:197,Y(3,:)); axis([0 200 -160 160]);
legend('PCA1','PCA2','PCA3');
title('Case 1: Principal Component Directions');
xlabel('Frame Number');
ylabel('Displacement (pixels)');