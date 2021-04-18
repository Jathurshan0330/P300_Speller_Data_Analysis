path = 'C:\Users\Jathu\Desktop\P300 Speller\p300\p300\';
do_whitening = 1;
window_size = 7;   %7 for 800ms     5 for 600ms


for data_num = 1:5
    data_num =4
    name = strcat('S',int2str(data_num));
    file = strcat(name,'.mat');
    filename = strcat(path, file)
    P300_speller = load(filename);
    fs = P300_speller.fs;
    data = P300_speller.y;
    trig = P300_speller.trig;
    num_channels = 8;
    chan = 8;
    
    %Common Average filtering
    data = com_avg(data,num_channels);
    
    %Spectral Whitening if necessary
    if do_whitening == 1
        for i =1:num_channels
            data(:,i)=whitening_Jathu(data,10,i) ;
        end
        x = "Data is whithened"
    end
    
    %Band Pass Filtering
    for i =1:num_channels
         data(:,i)= high_pass(data(:,i),0.5,fs) ;
         data(:,i)= low_pass(data(:,i),32,fs);              %200 for standford data
         %data_filtered(:,i)=filter(bandpass,data_filtered(:,i));
    end
    
    
    [Pxx_relax,F_relax] = periodogram(data(:,chan),hamming(length(data(:,chan))),length(data(:,chan)),fs);
    Pxx_relax = 10*log10(Pxx_relax);
    figure;
    plot(F_relax,Pxx_relax)
    title("Power Spectra of Preprocessed Signal");

    ca_fft=fft(data(:,chan));
    [n,p]= size(data(:,:));
    P2 = abs(ca_fft/n);
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(n/2))/n;
    figure;
    plot(f,P1)
    title('Single-Sided Amplitude Spectrum of Preprocessed Signal')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')

   
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
    
    [n,m] = size(pos_trails_onset);
    pos_data = [];
    
    for sam = 1:n
        window = data(pos_trails_onset(sam,1)-25:pos_trails_onset(sam,1)+(25*window_size)-1,:);
        for win = 1:window_size-2
            start = (win-1)*25 +1;
            [Pxx,F] = periodogram(window(start:start+99,:),hamming(length(window(start:start+99,:))),length(window(start:start+99,:)),fs);
            Pxx = 10*log10(Pxx);
            pos_data(win,1:8,sam) = mean(Pxx(1:2,:),1);
            pos_data(win,9:16,sam) = mean(Pxx(3:4,:),1);
            pos_data(win,17:24,sam) = mean(Pxx(4:6,:),1);
            pos_data(win,25:32,sam) = mean(Pxx(7:13,:),1);
        end
       % break
    end
    
    [n,m] = size(neg_trails_onset);
    neg_data = [];
    
    for sam = 1:n
        window = data(neg_trails_onset(sam,1)-25:neg_trails_onset(sam,1)+(25*window_size)-1,:);
        for win = 1:window_size-2
            start = (win-1)*25 +1;
            [Pxx,F] = periodogram(window(start:start+99,:),hamming(length(window(start:start+99,:))),length(window(start:start+99,:)),fs);
            Pxx = 10*log10(Pxx);
            neg_data(win,1:8,sam) = mean(Pxx(1:2,:),1);
            neg_data(win,9:16,sam) = mean(Pxx(3:4,:),1);
            neg_data(win,17:24,sam) = mean(Pxx(4:6,:),1);
            neg_data(win,25:32,sam) = mean(Pxx(7:13,:),1);
        end
    end
    
    
    pos_m = mean(pos_data(1:2,:,:),1);
    pos_data_mean=mean(pos_data,3);
    pos_data_norm = pos_data - pos_m;
    pos_data_norm_m  = mean(pos_data_norm,3);

    neg_m = mean(neg_data(1:2,:,:),1);
    neg_data_mean=mean(neg_data,3);
    neg_data_norm = neg_data - neg_m;
    neg_data_norm_m  = mean(neg_data_norm,3);

    figure;
    f=1:num_channels*4;
    t=-0.1:0.1:0.8;
    imagesc(t,f,pos_data_mean);   %temporal_rock(:,:,1)
    grid on
    colormap('jet')
    title("Feature Map for Target before Normalization")
    xlabel('Time(s)','FontSize',12,'FontWeight','bold')
    ylabel('Channels','FontSize',12,'FontWeight','bold')

    figure;
    f=1:num_channels*4;
    t=-0.1:0.1:0.8;
    imagesc(t,f,neg_data_mean);   %temporal_rock(:,:,1)
    grid on
    colormap('jet')
    title("Feature Map for Non-Target before Normalization")
    xlabel('Time(s)','FontSize',12,'FontWeight','bold')
    ylabel('Channels','FontSize',12,'FontWeight','bold')

    figure;
    f=1:num_channels*4;
    t=-0.1:0.1:0.8;
    imagesc(t,f,pos_data_norm_m);   %temporal_rock(:,:,1)
    grid on
    colormap('jet')
    title("Feature Map for Target after Normalization")
    xlabel('Time(s)','FontSize',12,'FontWeight','bold')
    ylabel('Channels','FontSize',12,'FontWeight','bold')

    figure;
    f=1:num_channels*4;
    t=-0.1:0.1:0.8;
    imagesc(t,f,neg_data_norm_m);   %temporal_rock(:,:,1)
    grid on
    colormap('jet')
    title("Feature Map for Non-Target after Normalization")
    xlabel('Time(s)','FontSize',12,'FontWeight','bold')
    ylabel('Channels','FontSize',12,'FontWeight','bold')

    
    
    
    
    
    name_1 = strcat("pos_data_S",int2str(data_num));
    name_1 = strcat(name_1,"_800ms");
    name_1 = strcat(name_1,"_Spectral_whitened");
    save(name_1, "pos_data_norm");
    
    name_1 = strcat("neg_data_S",int2str(data_num));
    name_1 = strcat(name_1,"_800ms");
    name_1 = strcat(name_1,"_Spectral_whitened");
    save(name_1, "neg_data_norm");
    
    break
end



%% 
start = min(pos_trails_onset(1,1),neg_trails_onset(1,1))
window = data(start-25:start+(25*window_size)-1,:);
relax = [];
for win = 1:window_size-2
    start = (win-1)*25 +1;
    [Pxx,F] = periodogram(window(start:start+99,:),hamming(length(window(start:start+99,:))),length(window(start:start+99,:)),fs);
    Pxx = 10*log10(Pxx);
    relax(win,1:8) = mean(Pxx(1:2,:),1);
    relax(win,9:16) = mean(Pxx(3:4,:),1);
    relax(win,17:24) = mean(Pxx(4:6,:),1);
    relax(win,25:32) = mean(Pxx(7:13,:),1);
end
relax_m = mean(relax,1);





%% 
pos_m = mean(pos_data(1:2,:,:),1);
pos_data_mean=mean(pos_data,3);

pos_data_norm  = mean(pos_data - pos_m,3);

neg_m = mean(neg_data(1:2,:,:),1);
neg_data_mean=mean(neg_data,3);
neg_data_norm  = mean(neg_data - neg_m,3);

figure;
f=1:num_channels*4;
t=-0.1:0.1:0.8;
imagesc(t,f,pos_data_mean);   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("Feature Map for Target before Normalization")
xlabel('Time(s)','FontSize',12,'FontWeight','bold')
ylabel('Channels','FontSize',12,'FontWeight','bold')

figure;
f=1:num_channels*4;
t=-0.1:0.1:0.8;
imagesc(t,f,neg_data_mean);   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("Feature Map for Non-Target before Normalization")
xlabel('Time(s)','FontSize',12,'FontWeight','bold')
ylabel('Channels','FontSize',12,'FontWeight','bold')

figure;
f=1:num_channels*4;
t=-0.1:0.1:0.8;
imagesc(t,f,pos_data_norm);   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("Feature Map for Target after Normalization")
xlabel('Time(s)','FontSize',12,'FontWeight','bold')
ylabel('Channels','FontSize',12,'FontWeight','bold')

figure;
f=1:num_channels*4;
t=-0.1:0.1:0.8;
imagesc(t,f,neg_data_norm);   %temporal_rock(:,:,1)
grid on
colormap('jet')
title("Feature Map for Non-Target after Normalization")
xlabel('Time(s)','FontSize',12,'FontWeight','bold')
ylabel('Channels','FontSize',12,'FontWeight','bold')
