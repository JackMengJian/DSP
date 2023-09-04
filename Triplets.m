clear all;
close all;
for o = -2 : 10
    file1 = ['E:\Transformer\30G\streamX_', num2str(o), 'dBm_100km_12_Rx.mat'];
    file2 = ['E:\Transformer\30G\streamY_', num2str(o), 'dBm_100km_12_Rx.mat'];
    load(file1);
    load(file2);
    cc = 30000;
    for k = 1001:(1000+cc)
        s=1;
        L=25;
        for i = -L:L
            for j = -L:L
                if abs(i)*abs(j) <= 15
                    aa=X_Rx(k+i)*X_Rx(k+j)*conj(X_Rx(k+i+j))+Y_Rx(k+j)*conj(Y_Rx(k+i+j))*X_Rx(k+i);
                    Rtr(k-1000,s+1)=imag(aa);
                    Rtr(k-1000,s)=real(aa);
                    s=s+2;
                end 
            end
        end
    end
    triplet = Rtr;
    file3 = ['E:\Transformer\30G\Triplets\2515\', num2str(o), 'dBm_100km_12_triplet.mat'];
    save(file3,"triplet");
end