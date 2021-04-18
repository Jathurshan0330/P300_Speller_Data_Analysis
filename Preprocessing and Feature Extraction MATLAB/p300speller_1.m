%% Read Data
path = 'C:\Users\Jathu\Desktop\P300 Speller\p300\p300\';
data_num = 1;
name = strcat('S',int2str(data_num));
file = strcat(name,'.mat');
filename = strcat(path, file);

P300_speller = load(filename);
fs = P300_speller.fs;
data = P300_speller.y;
trig = P300_speller.trig;
num_channels = 8;
chan = 8;
%% Plot Frequency Response
ca_filtered = data;  %84

ca_fft=fft(ca_filtered(:,chan));
[n,p]= size(ca_filtered(:,:));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of  Raw EEG')
xlabel('f (Hz)')
ylabel('|P1(f)|')

[Pxx_relax,F_relax] = periodogram(ca_filtered(:,chan),hamming(length(ca_filtered(:,chan))),length(ca_filtered(:,chan)),fs);
Pxx_relax = 10*log10(Pxx_relax);
figure;
plot(F_relax,Pxx_relax)
title("Power Spectra of Raw EEG Signal");

%% Common average Filtering
data = com_avg(data,num_channels); 

[Pxx_relax,F_relax] = periodogram(data(:,chan),hamming(length(data(:,chan))),length(data(:,chan)),fs);
Pxx_relax = 10*log10(Pxx_relax);
figure;
plot(F_relax,Pxx_relax)
title("Power Spectra of Common Average Filtered Signal");

ca_fft=fft(data(:,chan));
[n,p]= size(data(:,:));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of  Common Average Filtered EEG')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% Spectral Whitening

data_whitened = data;

for i =1:num_channels
        
    data_whitened(:,i)=whitening_Jathu(data_whitened,10,i) ;
    
end
for i =1:num_channels
     data_whitened(:,i)= high_pass(data_whitened(:,i),0.5,fs) ;
     data_whitened(:,i)= low_pass(data_whitened(:,i),35,fs);              %200 for standford data
     %data_filtered(:,i)=filter(bandpass,data_filtered(:,i));
end
[Pxx_relax,F_relax] = periodogram(data_whitened(:,chan),hamming(length(data_whitened(:,chan))),length(data_whitened(:,chan)),fs);
Pxx_relax = 10*log10(Pxx_relax);
figure;
plot(F_relax,Pxx_relax)
title("Power Spectra of Spectral Whitened Signal");

ca_fft=fft(data_whitened(:,chan));
[n,p]= size(data_whitened(:,:));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of  Spectral Whitened Signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Bandpass filtering

for i =1:num_channels
     data(:,i)= high_pass(data(:,i),0.5,fs) ;
     data(:,i)= low_pass(data(:,i),35,fs);              %200 for standford data
     %data_filtered(:,i)=filter(bandpass,data_filtered(:,i));
end

[Pxx_relax,F_relax] = periodogram(data(:,chan),hamming(length(data(:,chan))),length(data(:,chan)),fs);
Pxx_relax = 10*log10(Pxx_relax);
figure;
plot(F_relax,Pxx_relax)
title("Power Spectra of Bandpass Filtered (0.5 - 30 Hz)Signal");

ca_fft=fft(data(:,chan));
[n,p]= size(data(:,:));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of  Bandpass Filtered (0.5 - 30 Hz)Signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Plotting trigger data

figure;
plot(trig);
title('Trigger Signal')


%% Plot pos targets
[n,m] = size(trig);
count_pos = 0;
count_neg =0;
pos_trails_onset = [];
neg_trails_onset = [];
for j = 1:n
    if trig(j,1) == 1
        count_pos = count_pos + 1;
        pos_trails_onset(count_pos,1) = j;
    end
    
    if trig(j,1) == -1
        count_neg = count_neg + 1;
        neg_trails_onset(count_neg,1) = j;
    end
    
end


%% checking EEG plots
window = data_whitened(pos_trails_onset(4,1)-24:pos_trails_onset(4,1)+25*7-1,:);
figure;
for i=1:num_channels
    plot(window(:,i)+(i*0.18))
    hold on
end
legend('1','2','3','4','5','6','7','8');
title(" EEG for Targert +1 in range of pre -100ms to post +700ms Spectral whithened")


% window = data(neg_trails_onset(4,1)-24:neg_trails_onset(4,1)+25*7-1,:);
% figure;
% for i=1:num_channels
%     plot(window(:,i)+(i*50))
%     hold on
% end
% legend('1','2','3','4','5','6','7','8');
% title(" EEG for Targert -1 in range of pre -100ms to post +700ms")
