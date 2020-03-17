clear all; close all; clc;
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
%plot((1:length(y))/Fs,y);
%xlabel('Time [sec]'); ylabel('Amplitude');
%title('Mary had a little lamb (piano)');
p8 = audioplayer(y,Fs); playblocking(p8);
%piano_freq = fftshift(fft(y));
%setting up the domain
L=length(y)/Fs; n=length(y);
t2=linspace(0,L,n+1); t=t2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1];
ks=fftshift(k);
y=y';
%creating the slide
tslide=0:.1:length(y)/Fs;
Sgt_spec = zeros(length(tslide),length(y));
a = 100;
for j=1:length(tslide)
%guassian window
g=exp(-a*(t-tslide(j)).^2);
Sg=g.*y;
Sgt=fft(Sg);
%filtering around fundamental frequency
%[center_strength,index] = max(abs(Sgt));
%center_freq = k(index);
%filter = exp(-0.001*((k - center_freq).^2));
%Sgt=filter.*Sgt;
Sgt_spec(j,:) = fftshift(abs(Sgt));
end
figure(1);
pcolor(tslide,ks/(2*pi),Sgt_spec.')
set(gca,'Ylim',[150 500],'Fontsize',10)
title('Piano with Overtones','FontSize',14)
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
shading interp
colormap(hot)
colorbar;
clear all; close all; clc;
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
%plot((1:length(y))/Fs,y);
%xlabel('Time [sec]'); ylabel('Amplitude');
%title('Mary had a little lamb (recorder)');
p8 = audioplayer(y,Fs); playblocking(p8);
%rec_freq = fftshift(fft(y));
%setting up the domain
L=length(y)/Fs; n=length(y);
t2=linspace(0,L,n+1); t=t2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1];
ks=fftshift(k);
y=y';
%creating the slide
tslide=0:.1:length(y)/Fs;
Sgt_spec = zeros(length(tslide),length(y));
a = 100;
for j=1:length(tslide)
%guassian window
g=exp(-a*(t-tslide(j)).^2);
Sg=g.*y;
Sgt=fft(Sg);
%filtering around fundamental frequency
[center_strength,index] = max(abs(Sgt));
center_freq = k(index);
filter = exp(-0.001*((k - center_freq).^2));
Sgt=filter.*Sgt;
Sgt_spec(j,:) = fftshift(abs(Sgt));
end
figure(1);
pcolor(tslide,ks/(2*pi),Sgt_spec.')
set(gca,'Ylim',[650 1200],'Fontsize',10)
title('Recorder without Overtones','FontSize',14)
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
shading interp
colormap(hot)
colorbar;