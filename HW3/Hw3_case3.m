clear all; close all; clc;

%CAM1_3
load('cam1_3.mat');numFrames1 = size(vidFrames1_3,4);
%creates a spatial filter, initializing

set(0, 'DefaultLineLineWidth', 2);

filter = zeros(480,640);
filter(214:427,268:400) = 1;
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));

for j = 1:numFrames1
   
    X = vidFrames1_3(:,:,:,j); 
    %imshow(X); drawnow
    
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid1 = rgb2gray(X); %turns to grayscale
    filt_vid1 = gray_vid1.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid1,0.9);
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


plot(mean_x1,image_y,'.'); axis([0 640 0 480])

%CAM1_3
load('cam2_3.mat');numFrames2 = size(vidFrames2_3,4);
%creates a spatial filter, initializing

filter = zeros(480,640);
filter(177:420,211:423) = 1;
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));

for j = 1:numFrames2
   
    X = vidFrames2_3(:,:,:,j); 
    %imshow(X); drawnow
    
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid2 = rgb2gray(X); %turns to grayscale
    filt_vid2 = gray_vid2.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid2,0.9);
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

%CAM1_3
load('cam3_3.mat');numFrames3 = size(vidFrames3_3,4);
%creates a spatial filter, initializing

filter = zeros(480,640);
filter(154:337,245:495) = 1;
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));

for j = 1:numFrames3
   
    X = vidFrames3_3(:,:,:,j); 
    %imshow(X); drawnow
    
    filter_uint8 = uint8(filter); %converts filter to an uint8 type
    gray_vid3 = rgb2gray(X); %turns to grayscale
    filt_vid3 = gray_vid3.*filter_uint8; %applies the spatial filter
    
    %thresh = filt_vid1 > 250; %0.97
    %could also binarize it and drop the tolerance
    thresh = imbinarize(filt_vid3,0.9);
    %imshow(thresh)
    
    %finds all non-zero vectors
    indeces = find(thresh);
    
    %finds the matrix/vectors
    [Y, X] = ind2sub(size(thresh),indeces);
    
    %finds the centroid!
    mean_x3(j) = mean(X);
    mean_y3(j) = mean(Y);
    
end
   
%invert the mean_y1 because the axis of a picture is upside down
image_y3 = 480 - mean_y3;
image_x3 = 640 - mean_x3;



plot(1:numFrames1,mean_x1); axis([0 260 0 480]);hold on;
plot(1:numFrames2,mean_x2); hold on; plot(1:numFrames3,480-mean_y3);
legend('1.3','2.3','3.3')
title('Case 3: Unsynced Original Data');
xlabel('Frame Number');
ylabel('XY Displacement (pixels)');



%%
figure(); 
plot(1:numFrames1,image_y); axis([0 360 0 400]); hold on;
plot(1:numFrames2,image_y2); hold on; plot(1:numFrames3,image_x3);
title('Z');
legend('1_3','2_3','3_3')

figure();
plot(1:numFrames1,mean_x1); axis([0 360 0 400]); hold on;
plot(1:numFrames2,mean_x2); hold on; plot(1:numFrames3,image_y3);
title('XY');
legend('1_3','2_3','3_3')


min1 = 38; min2 = 26; min3 = 31; width = 200; 
% 200
% 255
% 206
shift_mean1 = mean_y1(min1:min1+width);
shift_mean2 = mean_y2(min2:min2+width);
shift_mean3 = mean_x3(min3:min3+width);

figure(); subplot(2,2,1);
plot(1:201,480 - shift_mean1); hold on;
plot(1:201,480 - shift_mean2); hold on;
plot(1:201,640 - shift_mean3); hold on;
axis([0 200 80 350]);
legend('1.3','2.3','3.3')
title('Case 3: Synced Original Data');
xlabel('Frame Number');
ylabel('Z Displacement (pixels)');

concat = [mean_x1(min1:min1+width);shift_mean1;mean_x2(min2:min2+width);...
    shift_mean2;mean_y3(min3:min3+width);shift_mean3];
%mean weights something less, normalize, average out teh center prosition


test1mat = [concat(1,:)-mean(concat(1,:));concat(2,:)-mean(concat(2,:));...
    concat(3,:)-mean(concat(3,:));concat(4,:)-mean(concat(4,:));...
    concat(5,:)-mean(concat(5,:));concat(6,:)-mean(concat(6,:))];

[u,s,v] = svd(test1mat/sqrt(200));

lambda = diag(s).^2; % produce diagonal variances
%Y= test1mat;
Y = u'*test1mat;
%Y = u'*test1mat;
sig= diag(s);
subplot(2,2,[2 4]);
plot (1:6 ,lambda/sum(lambda), 'mo ');
title (" Case 3: Energy Distribution");
xlabel ("Singular Values"); ylabel ("Energy");

subplot(2,2,3);
%plot(1:381,Y(1,:),1:381,Y(2,:),1:381,Y(3,:),1:381,Y(4,:),1:381,Y(5,:),1:381,Y(6,:));
plot(1:201,Y(1,:),1:201,Y(2,:),1:201,Y(3,:))
%,1:201,Y(4,:),1:201,Y(5,:),1:201,Y(6,:));
axis([0 200 -100 100]);
legend('PCA1','PCA2','PCA3');
title('Case 3: Principal Component Directions');
xlabel('Frame Number');
ylabel('Displacement (pixels)');

