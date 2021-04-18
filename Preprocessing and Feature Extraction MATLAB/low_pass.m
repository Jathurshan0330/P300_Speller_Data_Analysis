function x=low_pass(y,Fp,Fs) 
    Fn = Fs/2;                                              % Nyquist Frequency
    Wp = Fp/Fn;                                       % Theta Passband
    Ws = Wp *1.1;                                 % Bandstop Frequencies
    Rp = 10;                                                % Passband Ripple
    Rs = 40;                                                % Stopband Ripple
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Determine Optimal Order
    [b,a] = cheby2(n,Rs,Ws);                                % Transfer Function Coefficients
    [sos,g] = tf2sos(b,a);                      % Second-Order-Section For Stability
    x=filtfilt(sos,g,y);
end 