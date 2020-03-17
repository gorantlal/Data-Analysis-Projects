clear all; close all; clc; 

load handel
v = y'; figure(1)
plot((1:length(v))/Fs,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest, v(n)');

p8 = audioplayer(v,Fs);
playblocking(p8);
signal_trans = fftshift(fft(v));

tslide=0:0.1:length(v)/Fs; %sampling
Sgt_spec = zeros(length(tslide),length(v));
%indow = 0.1; %window width
new_array = zeros(length(v),1);
L=length(v)/Fs;
n=length(v);
t2=linspace(0,L,n+1);
t=t2(1:n); 
k=(2*pi/L)*[0:(n+1)/2-1 -n/2:-1]; 
ks=fftshift(k);
height = 1;
count =1;
%Shannon Filter

for a = [0.01 0.1 1 10] %Varying window length

    for j=1:length(tslide)

        new_array = zeros(length(v),1); 
        down = tslide(j)-(a/2);
        up = tslide(j)+(a/2);

        if down < 0 %creates a pulse
           B = t <= up;
           new_array(B') = height.*v(B');

        else
           B = t>=down & t<=up;
           new_array(B') = height.*v(B');

        end
        Sgt=fft(new_array);
        Sgt_spec(j,:) = fftshift(abs(Sgt));
    end
    
    %Making a plot
    subplot(2,2,count); hold on;
    sgtitle('Spectrogram of Varying Shannon Window Lengths',...
        'FontSize',16,'FontWeight','bold');
    pcolor(tslide,ks/(2*pi),Sgt_spec.'), 
    shading interp
    title(a, 'FontSize',14);
    ylabel('Frequency (Hz)', 'FontSize',12);
    xlabel('Time (sec)','FontSize',12);
    set(gca,'Ylim',[-2000 2000],'Xlim',[0 8.9],'Fontsize',10)
    colormap(hot)
    colorbar;
    count = count + 1;
    
end
