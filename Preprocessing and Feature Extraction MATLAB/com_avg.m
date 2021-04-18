function y=com_avg(x,no_of_channel)
    y=x;
    x1=sum(x,2)/no_of_channel;
    %[n,m]=size(x);
    size(x1);
    y=x-x1;
    
    
end
