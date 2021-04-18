%% Load data
ECoG_Handpose=load("ECoG_Handpose.mat");
data=ECoG_Handpose.y;
data=transpose(data);
Fs=1200;
%% %% Apply a high pass filter for all channels
data_filtered=data;
Fs=1200;
for i =2:61
    data_filtered(:,i)=high_pass(data_filtered(:,i),0.5,Fs) ;
    data_filtered(:,i)=low_pass(data_filtered(:,i),200,Fs);
end

% data_filtered(:,2:61)=high_pass(data_filtered(:,2:61),0.5,Fs) ;
% data_filtered(:,2:61)=low_pass(data_filtered(:,2:61),200,Fs);

%% %% Common average the filter 
ca_filtered=com_avg(data_filtered,60);

% for i =2:61
%     ca_filtered(:,i)=high_pass(ca_filtered(:,i),0.5,Fs) ;
%     ca_filtered(:,i)=low_pass(ca_filtered(:,i),200,Fs);
% end

ca_fft=fft(ca_filtered(:,2));
[n,p]= size(ca_filtered(:,:));
P2 = abs(ca_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of  Common Average Filtered ECoG')
xlabel('f (Hz)')
ylabel('|P1(f)|')


data_fft=fft(data_filtered(:,2));
[n,p]= size(data_filtered(:,:));
P2 = abs(data_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of ECoG')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%% Notch filter
notch_filtered=ca_filtered;

for i = 1:6
    Fc=[50*i-0.2,50*i+0.2];
    Fc
    [b,a] = butter(3,Fc/(Fs/2),'stop');

    %freqz(b,a);
    for i =2:61
        notch_filtered(:,i)=filtfilt(b,a,notch_filtered(:,i)) ;

    end
end

notch_fft=fft(notch_filtered(:,61));
[n,p]= size(notch_filtered(:,:));
P2 = abs(notch_fft/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of Notch filtered ECoG')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% 
% 
