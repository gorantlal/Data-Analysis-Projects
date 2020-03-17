clear; close all; clc;

load Testdata

L=15; % spatial domain
n=64; % Fourier modes

x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x; %spatial domain discretization
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k); %freq components of fft
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); %makes a coordinate system

ave_sig = zeros(n,n,n); %creates a 3D array of 0s

%Averaging the spectrum
for j=1:20 
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    ut = fftn(Un);
    ave_sig = ave_sig + ut;
end
ave_sig = abs((fftshift(ave_sig)))./20;

%finding the max value, finding corresponding indices
[max_amp,index] = max(ave_sig(:));
[center_x,center_y,center_z] = ind2sub([n,n,n],index);

%utilizing fake indices to find the correct indicis
center_kx = Kx(center_x,center_y,center_z);
center_ky = Ky(center_x,center_y,center_z);
center_kz = Kz(center_x,center_y,center_z);

% Guassian Filter
filter = exp(-0.2*((Kx - center_kx).^2 + (Ky - center_ky).^2 + (Kz - center_kz).^2));

pos_x = zeros(20,1);
pos_y = zeros(20,1);
pos_z = zeros(20,1);

%applying the filter on each time point
for j=1:20
    
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    ut = fftshift(fftn(Un));
    filtered_sig = filter.*ut;
    og = (ifftn(filtered_sig));
    
    [max_sig,Index] = max(Un(:)); %finding max position
  
    
    [center_X,center_Y,center_Z] = ind2sub(size(Un),Index);
    pos_x(j) = X(center_X,center_Y,center_Z);
    pos_y(j) = Y(center_X,center_Y,center_Z);
    pos_z(j) = Z(center_X,center_Y,center_Z);
    
    
%     [max_sig,Index] = max(og(:)); %finding max position
%     [max_sig_1,Index_1] = max(Un(:))
%     
%     [center_X,center_Y,center_Z] = ind2sub(size(og),Index);
%     pos_x(j) = X(center_X,center_Y,center_Z);
%     pos_y(j) = Y(center_X,center_Y,center_Z);
%     pos_z(j) = Z(center_X,center_Y,center_Z);
end

%mapping the trajectory
plot3(pos_x,pos_y,pos_z,'-.','Linewidth',3); grid on;
title('Position of Marble Over Time Course', 'FontSize',12);
xlabel('X direction');ylabel('Y direction');zlabel('Z direction');
labels = {'t=1','t=2','t=3','t=4','t=5','t=6','t=7','t=8','t=9'...
    't=10','t=11','t=12','t=13','t=14','t=15','t=16','t=17','t=18'...
    't=19','t=20'};
text(pos_x,pos_y,pos_z,labels,'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','FontSize',12)

%where accoustic wave should be focused
lastx = pos_x(end);
lasty = pos_y(end);
lastz = pos_z(end);
