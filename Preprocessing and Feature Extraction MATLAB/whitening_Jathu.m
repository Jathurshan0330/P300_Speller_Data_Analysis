function y_out=whitening_Jathu(data,order,channel) 
y = data(:,channel);
A = lpc(y,order);
Filt = [0 -A(2:end)];
est_y = filter(Filt,1,y);
y_out = y - est_y;
end
    