function y_out=whitening(data,order,channel) 
y_n = data(channel,:);
A = aryule(y_n,order);
y_out = zeros(1,(size(y_n,2)+order));
for i = 0:order
    y_p = zeros(1,(size(y_n,2)+i));
    y_p(i+1:end) = y_n;
    y_p = [y_p zeros(1,order-i)];
    y_p = A(1,i+1)*y_p;
    y_out = y_out + y_p;
end
y_out = y_out(1,1:(size(y_n,2)));
%y_out = y_n-y_out;
end
