%Loading the music
load handel
v = y';
figure(1); plot((1:length(v))/Fs,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest, v(n)');
p8 = audioplayer(v,Fs); playblocking(p8);
signal_trans = fftshift(fft(v));

%Setting up the domains
L=length(v)/Fs; n=length(v);
t2=linspace(0,L,n+1); t=t2(1:n);
k=(2*pi/L)*[0:(n+1)/2-1 -n/2:-1];
ks=fftshift(k);

figure(2) %The frequency plot
plot(ks,abs(signal_trans)/max(abs(signal_trans)));

%a = 100; defining the window width
count = 1;
tslide=0:.1:L; %list of tau values to sample
Sgt_spec = zeros(length(tslide),length(v));

for a = [1 10 100 1000] %Varying window length
    for j=1:length(tslide) %creating the spectrogram
    g=exp(-a*(t-tslide(j)).^2);
    Sg=g.*v;
    Sgt=fft(Sg);
    Sgt_spec(j,:) = fftshift(abs(Sgt));
    end

 %Making a plot
subplot(2,2,count); hold on;
sgtitle('Spectrogram of Varying Gaussian Window Lengths',...
'FontSize',16,'FontWeight','bold');
pcolor(tslide,ks/(2*pi),Sgt_spec.'),
shading interp
 title(a, 'FontSize',14);
 ylabel('Frequency (Hz)', 'FontSize',12);
xlabel('Time (sec)','FontSize',12);
 set(gca,'Ylim',[-2500 2500],'Xlim',[0 8.9],'Fontsize',10)
colormap(hot)
 colorbar;
 count = count + 1;

end
