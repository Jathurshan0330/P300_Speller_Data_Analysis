function x=band_pass(y,Fp1,Fp2,Fs) 
    Fn = Fs/2;                                              % Nyquist Frequency
    Wp = [Fp1/Fn Fp2/Fn];                                       % Theta Passband
    Ws = [Wp(1)*0.9,Wp(2)*1.1];   %0.4                              % Bandstop Frequencies
    Rp = 10;                                                % Passband Ripple
    Rs = 40;                                                % Stopband Ripple
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                          % Determine Optimal Order   %%buttord
    [b,a] = cheby2(n,Rs,Ws);                                % Transfer Function Coefficients %%butter(n,Wp,"high")   
    [sos,g] = tf2sos(b,a);                      % Second-Order-Section For Stability
    x=filtfilt(sos,g,y);
end 