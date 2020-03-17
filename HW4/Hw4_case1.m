% AMATH 482 Homework 4
% Lahari Gorantla
clear all; close all; clc; 

% Test Case 1
% Choose 3 different bands of different genres

%% Uploading the music

num_song = 30;
[flume] = getSong('flume', num_song);
[hall] = getSong('hall', num_song);
[nat] = getSong('nat', num_song);

%% Applying the Spectrogram with a Gabor function

ds = 4; %downsampling
L = 5; % 5 seconds for each audio clip
fs = 44100; fs_d = 44100/ds; %frequency electronic music was sampled at
n = fs_d * L; % number of Fourier modes
t2=linspace(0,L,n+1); t=t2(1:n); % time domain discretization
k=(1/L)*[0:(n-1)/2 (-n+1)/2:-1]; ks=fftshift(k); % freq domain

%obtaining spectro data through a function
[Sgt_flume] = mySpectro(flume, 50, 0.1, L, n, num_song, t);
[Sgt_hall] = mySpectro(hall, 50, 0.1, L, n, num_song, t);
[Sgt_nat] = mySpectro(nat, 50, 0.1, L, n, num_song, t);

%% Taking the SVD

[U,S,V] = svd([Sgt_flume, Sgt_hall, Sgt_nat],'econ'); 
%% Plotting the Singular Values

figure(1);
semilogy(diag(S),'ko','Linewidth',2);
title('SVD Decomposition');
xlabel('Singular Values');
ylabel('Energy');

%% LDA!

feature = 20;
[U,S,V,threshold_1, threshold_2,w,low_p,med_p,high_p] = dc_trainer(Sgt_flume,Sgt_hall,Sgt_nat,feature);

%% Plotting

figure(1);
semilogy(diag(S),'ko','Linewidth',2)
title('SVD Decomposition: Flume, Hall & Oats, Nat King Cole');
xlabel('Singular Values');
ylabel('Energy');

%plot dog/cat projections onto w
figure(2); subplot(2,1,1);
plot(low_p,zeros(length(low_p)),'ob','Linewidth',2); hold on
plot(med_p,1.*ones(length(med_p)),'or','Linewidth',2); hold on
plot(high_p, 2.*ones(length(high_p)),'og','Linewidth',2); hold on
title('Projections onto w basis');
ylim([0 2.5])

%making a histogram of overlap with thresholds
subplot(2,1,2);
histogram(low_p,length(low_p)) 
hold on 
histogram(med_p,length(med_p))
histogram(high_p,length(high_p)) 
plot([threshold_1 threshold_1],[0 18],'r','LineWidth',2)
plot([threshold_2 threshold_2],[0 18],'b','LineWidth',2)
title('Statistics on LDA Basis Projection: Flume, Hall, Nat')
set(gca,'Fontsize',12)
legend('Flume', 'Hall & Oats','Nat King Cole','Threshold 1','Threshold 2')

%% Classifying Code

[testset] = getSong('test1', 30);
[test_spect] = mySpectro(testset, 50, 0.1, L, n, num_song, t);

TestNum = 30;
TestMat = U'*test_spect;  % PCA projection
pval = w'*TestMat;  % LDA projection

hiddenlabels = zeros(30,1);
hiddenlabels(1:10) = 0;
hiddenlabels(11:20) = 1;
hiddenlabels(21:30) = 2;

label = zeros(30,1);
for i = 1:30
    
    if pval(1,i) > threshold_1 && pval(1,i) < threshold_2
        label(i) = 1;
        
    elseif pval(1,i) > threshold_1 && pval(1,i) > threshold_2
        label(i) = 2;
        
    else
        label(i) = 0;
        
    end
end

% add the success rate
count = 0;

for i = 1:30
    if hiddenlabels(i) == label(i)
        count = count + 1;
    end
end

sucRate = count/TestNum;
disp(sucRate);

%%
set(0, 'DefaultLineLineWidth', 2);

subplot(1,2,1)
plot(low_p,zeros(length(low_p)),'ob','Linewidth',2); hold on
plot(med_p,1.*ones(length(med_p)),'or','Linewidth',2); hold on
plot(high_p, 2.*ones(length(high_p)),'og','Linewidth',2); hold on
plot(pval(1:10), zeros(length(1:10)),'ok','Linewidth',2); hold on;
plot(pval(11:20), 1.*ones(length(1:10)),'ok','Linewidth',2); hold on;
plot(pval(21:30),2.*ones(length(1:10)),'ok','Linewidth',2); hold on;
plot([threshold_1 threshold_1],[0 18],'r','LineWidth',2); hold on;
plot([threshold_2 threshold_2],[0 18],'b','LineWidth',2); hold on;
%legend('BB','BTS','GG')
title('Projections Onto W Basis');
ylim([0 3.5])

subplot(1,2,2);
histogram(low_p,length(low_p)) 
hold on 
histogram(med_p,length(med_p))
histogram(high_p,length(high_p)) 
plot([threshold_1 threshold_1],[0 18],'r','LineWidth',2)
plot([threshold_2 threshold_2],[0 18],'b','LineWidth',2)
title('Statistics on LDA Basis Projection: Flume, Hall, Nat')
set(gca,'Fontsize',12)
legend('Flume', 'Hall & Oats','Nat King Cole','Threshold 1','Threshold 2')



%% Uploading music function
% goal: build a function
% name of the song and the number associated with it
function [music_array] = getSong(song_name,num_song)
    for i = 1:num_song
        file_name = sprintf('%s_%d.wav', song_name, i);
        [y, fs] = audioread(file_name);
        ds_song = downsample(y,4);
        music_array(i,:) = ds_song(:,1);
    end
end

%% My Spectrogram function

function [sg] = mySpectro(v, width, tau, L, n, num_song, t) 
    a = width; tslide=0:tau:L; 
    Spec_g = zeros(length(tslide),floor(n)); % store filtered frequency data
    sg = []; % store spectrogram of all 30 songs a row vectors
    for song = 1:num_song 
        for j=1:length(tslide)
            % Gabor filter function / window
            g = exp(-a*(t-tslide(j)).^2);
            filtdat = g.*v(song,:); % apply filter to signal (mutiplication in time domain) 
            sgt = fft(filtdat); 
            Spec_g(j,:) = fftshift(abs(sgt)); % We don't want to scale it
        end

        row_vec = [];
        for i = 1:length(tslide)
            row_vec = [row_vec Spec_g(i, :)];
        end
        sg =[sg; row_vec];
    end
    sg = sg';
end

%% LDA

function [U,S,V,threshold_1, threshold_2,w,low_p,med_p,high_p] = dc_trainer(set1,set2,set3,feature);

    % number of columns
    n1 = size(set1,2); n2 = size(set2,2); n3 = size(set3,2);
    
    [U,S,V] = svd([set1,set2,set3],'econ');
    
    proj = S*V'; % projection onto principal components
    U = U(:,1:feature);
    proj1 = proj(1:feature,1:n1);
    proj2 = proj(1:feature,n1+1:n1+n2);
    proj3 = proj(1:feature,n1+n2+1:n1+n2+n3);
    
    mean1 = mean(proj1,2);
    mean2 = mean(proj2,2);
    mean3 = mean(proj3,2);
    
    Sw = 0; % within class variances
    
    for k=1:n1
        Sw = Sw + (proj1(:,k)-mean1)*(proj1(:,k)-mean1)';
    end
    for k=1:n2
        Sw = Sw + (proj2(:,k)-mean2)*(proj2(:,k)-mean2)';
    end
    for k=1:n3
        Sw = Sw + (proj3(:,k)-mean3)*(proj3(:,k)-mean3)';
    end
    
    % between class variance
    bigMean = (mean1+mean2+mean3)/3;
    Sb1 = (mean1-bigMean)*(mean1-bigMean)';
    Sb2 = (mean2-bigMean)*(mean2-bigMean)';
    Sb3 = (mean3-bigMean)*(mean3-bigMean)';
    Sb = (Sb1 + Sb2 + Sb3)/3;
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    v_1 = w'*proj1; 
    v_2 = w'*proj2;
    v_3 = w'*proj3;
    
    result = [v_1;v_2;v_3];
    
    sort_mean_1 = mean(v_1);
    sort_mean_2 = mean(v_2);
    sort_mean_3 = mean(v_3);    
    [~,I] = sort([sort_mean_1, sort_mean_2, sort_mean_3]);
    low_vals = I(1);
    med_vals = I(2);
    high_vals = I(3);
          
    low_p = result(low_vals,:); 
    low_p = sort(low_p);
    med_p = result(med_vals,:); 
    med_p = sort(med_p);
    high_p = result(high_vals,:); 
    high_p = sort(high_p);
      
    t1 = length(low_p);
    t2 = 1;
    while low_p(t1)>med_p(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    threshold_1 = (low_p(t1)+med_p(t2))/2;
    
    t2 = length(med_p);
    t3 = 1;
    while med_p(t2)>high_p(t3)
        t2 = t2-1;
        t3 = t3+1;
    end
    threshold_2 = (med_p(t2)+high_p(t3))/2;   
end
