function  [X_Rx,carrier_phase]=BPS(X_Rx,M,B,BPS_window,switch_case,phase_origin)
    N       =length(X_Rx);
    %M :调制阶数
    %B：测试相位个数
    %BPS_window：BPS窗口长度
    %switch_case:1 高斯距离 2 曼哈顿距离  3 切比雪夫距离
    %phase_orgin  0
    %carrier_phase:跟踪相位
    phB     =(-B/2:B/2-1)/B*pi/2;%(b-B/2)/B
    carrier_phase=zeros(1,N);
    switch switch_case
        case 1
            X_Rx_tmp=kron(X_Rx,ones(B,1)).';
            phB_tmp=kron(phB,ones(N,1));
            difB=X_Rx_tmp.*exp(1j*phB_tmp);
            temp = qamdemod(difB,M,'gray');
            Tx_de=qammod(temp,M,'gray');
            dis=(abs(difB-Tx_de)).^2;
        case 2
            X_Rx_tmp=kron(X_Rx,ones(B,1)).';
            phB_tmp=kron(phB,ones(N,1));
            difB=X_Rx_tmp.*exp(1j*phB_tmp);
            temp = qamdemod(difB,M,'gray');
            Tx_de=qammod(temp,M,'gray');
            dis_i=abs(real(difB)-real(Tx_de));
            dis_q=abs(imag(difB)-imag(Tx_de));
            dis=dis_i+dis_q;
        case 3
            X_Rx_tmp=kron(X_Rx,ones(B,1)).';
            phB_tmp=kron(phB,ones(N,1));
            difB=X_Rx_tmp.*exp(1j*phB_tmp);
            temp = qamdemod(difB,M,'gray');
            Tx_de=qammod(temp,M,'gray');
            dis_i=abs(real(difB)-real(Tx_de));
            dis_q=abs(imag(difB)-imag(Tx_de));
            dis=max(dis_i,dis_q);
        case 4
            X_Rx_tmp=kron(X_Rx,ones(B,1)).';
            phB_tmp=kron(phB,ones(N,1));
            difB=X_Rx_tmp.*exp(1j*phB_tmp);
            dis_i=abs(abs(abs(abs(abs(real(difB))-8)-4)-2)-1);
            dis_q=abs(abs(abs(abs(abs(imag(difB))-8)-4)-2)-1);
            dis=dis_i+dis_q;
    end
            
 
    
    for kk =1:BPS_window
    dis_all  = sum(dis(kk:kk+BPS_window,:));
    [~,bb]=min(dis_all);
    carrier_phase(kk)=phB(bb);
    end
    
    for kk = BPS_window+1:N-BPS_window
    dis_all  = sum(dis(kk-BPS_window:kk+BPS_window,:));
    [~,bb]=min(dis_all);
    carrier_phase(kk)=phB(bb);
    end
    
    for kk =N-2*BPS_window+1:N-BPS_window
    dis_all  = sum(dis(kk:kk+BPS_window,:));
    [~,bb]=min(dis_all);
    carrier_phase(kk+BPS_window)=phB(bb);
    end
    carrier_phase=-carrier_phase;
    
    %%%%%%%相位解绕%%%%%%
    %%%%%%%%%%%%%%%%%%%%
   switch M
       case 16
             carrier_phase=phase_unwrap_2(carrier_phase,phase_origin);
       case 64 
            carrier_phase=phase_unwrap_2(carrier_phase,phase_origin);
       case 32
             carrier_phase=phase_unwrap_3(carrier_phase,phase_origin);
             carrier_phase=phase_unwrap_3(carrier_phase,phase_origin);
   end
%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%  %%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%
X_Rx=X_Rx.*exp(-1i*carrier_phase);
end

