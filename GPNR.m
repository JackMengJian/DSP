function [X_Rx,carrier_phase]=GPNR(X_Rx,X_Tx,avg_window)
N=length(X_Rx);
for kk=1:N-avg_window
    sum1=0;
    for j=kk:kk+avg_window
        sum1=sum1+X_Rx(j)*conj(X_Tx(j));
        carrier_phase(kk)=angle(sum1);
    end
end

for kk=N-2*avg_window+1:N-avg_window
    sum1=0;
    for j=kk:kk+avg_window
           sum1=sum1+X_Rx(j)*conj(X_Tx(j));
        carrier_phase(kk+avg_window)=angle(sum1);
    end
end

carrier_phase=unwrap(carrier_phase);
X_Rx=X_Rx.*exp(-1i*carrier_phase);