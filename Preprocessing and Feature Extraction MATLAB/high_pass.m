function x=high_pass(y,Fp,Fs) 
    Fn = Fs/2;                                              % Nyquist Frequency
    Wp = Fp/Fn;                                       % Theta Passband
    Ws = Wp *0.4;   %0.4                              % Bandstop Frequencies
    Rp = 10;                                                % Passband Ripple
    Rs = 40;                                                % Stopband Ripple
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                          % Determine Optimal Order   %%buttord
    [b,a] = cheby2(n,Rs,Ws,"high");                                % Transfer Function Coefficients %%butter(n,Wp,"high")   
    [sos,g] = tf2sos(b,a);                      % Second-Order-Section For Stability
    x=filtfilt(sos,g,y);
end 