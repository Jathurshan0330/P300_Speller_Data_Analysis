path = 'C:\Users\Jathu\Desktop\P300 Speller\p300\p300\';
do_whitening = 1;
window_size = 7;   %7 for 800ms     5 for 600ms


for data_num = 1:5
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
    pos_stat = [];
    
    for sam = 1:n
        window = data(pos_trails_onset(sam,1)-25:pos_trails_onset(sam,1)+(25*window_size)-1,:);
        [Pxx,F] = periodogram(window,hamming(length(window)),length(window),fs);
        Pxx = 10*log10(Pxx);
        pos_stat(1,:,sam) = mean(Pxx(1:4,:),1);
        pos_stat(2,:,sam) = mean(Pxx(4:7,:),1);
        pos_stat(3,:,sam) = mean(Pxx(8:12,:),1);
        pos_stat(4,:,sam) = mean(Pxx(13:25,:),1);
        
        pos_stat(5,:,sam) = mean(window,1);
        pos_stat(6,:,sam) = max(window);
        pos_stat(7,:,sam) = min(window);
        pos_stat(8,:,sam) = rms(window);
    end
    
    [n,m] = size(neg_trails_onset);
    neg_stat = [];
    
    for sam = 1:n
        window = data(neg_trails_onset(sam,1)-25:neg_trails_onset(sam,1)+(25*window_size)-1,:);
        [Pxx,F] = periodogram(window,hamming(length(window)),length(window),fs);
        Pxx = 10*log10(Pxx);
        neg_stat(1,:,sam) = mean(Pxx(1:4,:),1);
        neg_stat(2,:,sam) = mean(Pxx(4:7,:),1);
        neg_stat(3,:,sam) = mean(Pxx(8:12,:),1);
        neg_stat(4,:,sam) = mean(Pxx(13:25,:),1);
        
        neg_stat(5,:,sam) = mean(window,1);
        neg_stat(6,:,sam) = max(window);
        neg_stat(7,:,sam) = min(window);
        neg_stat(8,:,sam) = rms(window);
    end
    
    name_1 = strcat("pos_data_stat_S",int2str(data_num));
    name_1 = strcat(name_1,"_800ms");
    name_1 = strcat(name_1,"_Spectral_whitened");
    save(name_1, "pos_stat");
    
    name_1 = strcat("neg_data_stat_S",int2str(data_num));
    name_1 = strcat(name_1,"_800ms");
    name_1 = strcat(name_1,"_Spectral_whitened");
    save(name_1, "neg_stat");
    
    %break
end