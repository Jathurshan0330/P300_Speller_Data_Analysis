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
    
    d_1 = [];  %0.5-4Hz
    d_2 = [];  %4 - 7.5 Hz
    d_3 = [];  %7.5 - 13.5Hz
    d_4 = [];  %13.5 - 30.5 Hz
    %Band Pass Filtering
    for i =1:num_channels
         %d_1(:,i)= high_pass(data(:,i),0.5,fs) ;
         %d_1(:,i)= low_pass(data(:,i),4,fs);              %200 for standford data
         d_1(:,i)= low_pass(data(:,i),4,fs); 
         
%          d_2(:,i)= high_pass(data(:,i),4,fs) ;
%          d_2(:,i)= low_pass(data(:,i),7.5,fs);
         d_2(:,i)= band_pass(data(:,i),4.2,7.5,fs); 
%          d_3(:,i)= high_pass(data(:,i),7.5,fs) ;
%          d_3(:,i)= low_pass(data(:,i),13.5,fs);
         d_3(:,i)= band_pass(data(:,i),7.5,13.5,fs); 
%          d_4(:,i)= high_pass(data(:,i),13.5,fs) ;
%          d_4(:,i)= low_pass(data(:,i),30.5,fs);
         d_4(:,i)= band_pass(data(:,i),13.5,30.5,fs);
         %data_filtered(:,i)=filter(bandpass,data_filtered(:,i));
    end
    
    chan = 1;
    [Pxx_relax,F_relax] = periodogram(d_1(:,chan),hamming(length(d_1(:,chan))),length(d_1(:,chan)),fs);
    Pxx_relax = 10*log10(Pxx_relax);
    figure;
    plot(F_relax,Pxx_relax)
    title("Power Spectra of Preprocessed Signal d1");

    [Pxx_relax,F_relax] = periodogram(d_2(:,chan),hamming(length(d_2(:,chan))),length(d_2(:,chan)),fs);
    Pxx_relax = 10*log10(Pxx_relax);
    figure;
    plot(F_relax,Pxx_relax)
    title("Power Spectra of Preprocessed Signal d2");
    
    [Pxx_relax,F_relax] = periodogram(d_3(:,chan),hamming(length(d_3(:,chan))),length(d_3(:,chan)),fs);
    Pxx_relax = 10*log10(Pxx_relax);
    figure;
    plot(F_relax,Pxx_relax)
    title("Power Spectra of Preprocessed Signal d3");
    
    [Pxx_relax,F_relax] = periodogram(d_4(:,chan),hamming(length(d_4(:,chan))),length(d_4(:,chan)),fs);
    Pxx_relax = 10*log10(Pxx_relax);
    figure;
    plot(F_relax,Pxx_relax)
    title("Power Spectra of Preprocessed Signal d4");

   
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
        pos_data(:,1:8,sam) = d_1(pos_trails_onset(sam,1)-25:pos_trails_onset(sam,1)+(25*window_size)-1,:);
        pos_data(:,9:16,sam) = d_2(pos_trails_onset(sam,1)-25:pos_trails_onset(sam,1)+(25*window_size)-1,:);
        pos_data(:,17:24,sam) = d_3(pos_trails_onset(sam,1)-25:pos_trails_onset(sam,1)+(25*window_size)-1,:);
        pos_data(:,25:32,sam) = d_4(pos_trails_onset(sam,1)-25:pos_trails_onset(sam,1)+(25*window_size)-1,:);
    end
    
    [n,m] = size(neg_trails_onset);
    neg_data = [];
    
    for sam = 1:n
        neg_data(:,1:8,sam) = d_1(neg_trails_onset(sam,1)-25:neg_trails_onset(sam,1)+(25*window_size)-1,:);
        neg_data(:,9:16,sam) = d_2(neg_trails_onset(sam,1)-25:neg_trails_onset(sam,1)+(25*window_size)-1,:);
        neg_data(:,17:24,sam) = d_3(neg_trails_onset(sam,1)-25:neg_trails_onset(sam,1)+(25*window_size)-1,:);
        neg_data(:,25:32,sam) = d_4(neg_trails_onset(sam,1)-25:neg_trails_onset(sam,1)+(25*window_size)-1,:);
    end
    
    name_1 = strcat("pos_data_S",int2str(data_num));
    name_1 = strcat(name_1,"_800ms");
    name_1 = strcat(name_1,"_Spectral_whitened");
    save(name_1, "pos_data");
    
    name_1 = strcat("neg_data_S",int2str(data_num));
    name_1 = strcat(name_1,"_800ms");
    name_1 = strcat(name_1,"_Spectral_whitened");
    save(name_1, "neg_data");
    
    break
end