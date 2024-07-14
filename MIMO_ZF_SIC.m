clc;
clear all;
close all;
tic
%% initialistion
N=10^5;
SNR_dB=0:25;
SNR=10.^(0.1*SNR_dB);
L=length(SNR_dB);
no_Rx=2;
no_Tx=2;
mod=input("Enter the modulation scheme (1-->BPSK) (2-->QPSK)");
no_err=zeros(1,L);
no_err_zf=zeros(1,L);

%% BPSK modulation

for k=1:L
    if(mod==1)
        s=bpsk_mod(N);
    else
        s=qpsk_mod(N);
    end

    s1=kron(ones(no_Rx,1),s);
    s2=reshape(s,1,no_Tx,N/no_Tx);
    h=sqrt(1/2)*(randn(no_Rx,no_Tx,N/no_Tx)+1i*randn(no_Rx,no_Tx,N/no_Tx));
    n=sqrt(1/2)*(randn(no_Rx,N/no_Rx)+1i*randn(no_Rx,N/no_Rx));
    d=zeros(no_Rx,N/no_Tx);
    for kk=1:N/no_Tx
        d(:,kk)=h(:,:,kk)*flip(s2(:,:,kk).'); 
    end
    d=d+n*sqrt(1/SNR(k));
    %at receiver
    h2=zeros(no_Rx,no_Tx,N/no_Tx); % to store inv(H^H*H)
    h3=zeros(no_Rx,no_Tx,N/no_Tx);% % to store (H^H)
    %below for loop is to get inv(H^H*H)
    for j=1:N/no_Tx
        temp=h(:,:,j)'*h(:,:,j);
        h3(:,:,j)=h(:,:,j)';
        h2(:,:,j)=inv(temp);%to get inv(H^H*H)
    end
    y3=zeros(no_Rx,N/no_Tx);
   det=[];
   det_zf=[];
    for pp=1:N/no_Tx
        yt=h2(:,:,pp)*h3(:,:,pp)*d(:,pp);
        y3(:,pp)=h2(:,:,pp)*h3(:,:,pp)*d(:,pp); %rx symbols as [y1;y2];
        if(mod==1)
            if(real(yt(2))>0)
                det2=1;
            else
                det2=-1;
            end
            if(real(yt(1))>0)
                det1=1;
            else
                det1=-1;
            end
        else
            det2=qpsk_decision(yt(2)); %detemining x2;
            det1=qpsk_decision(yt(1));
        end
            det=[det det2]; %estimated y2
            det_zf=[det_zf det2 det1];
        
     end

    
    
    %now need to extract h12 and h22 from h
    he=h(:,2,:);
    % now we need to multipy h12*estimated y2  and h22*estimated y2 
    dety2=kron(det,ones(no_Rx,1));
    dety2=reshape(dety2,no_Rx,1,N/no_Tx);
    det_inter=he.*dety2; %interference due to x2; ie [h12*x2; h22*x2]
    det_inter=squeeze(det_inter);
    y_x1=d-det_inter;%subracting the interfernce of x2 from recevied symbol;
     % using MRC to decode x1; get the h11 and h21 
    h_mrc=h(:,1,:);
    h_mrc=squeeze(h_mrc);
    yhat=zeros(1,N/no_Tx);
    det_x1=[];
    yhat = sum(conj(h_mrc).*y_x1,1)./sum(h_mrc.*conj(h_mrc),1);
    for a=1:N/no_Tx 
        if(mod==1)
            if(real(yhat(a))>0)
                det11=1;
            else
                det11=-1;
            end
        else 
            det11=qpsk_decision(yhat(a)); %yhat is not scalar
        end
        det_x1=[det_x1 det11];
    end
    %combining x1 and x2
    y=[];
    for t=1:N/no_Tx
        y=[y det(t) det_x1(t)];
    end
 
    no_err(k) = size(find([s- y]),2);
    no_err_zf(k)=size(find([s- det_zf]),2);
     
     
end
simBer_sic = no_err/N; % simulated ber
simber_zf=no_err_zf/N;
SNRLin = 10.^(SNR_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./SNRLin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./SNRLin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p));

figure
semilogy(SNR_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(SNR_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
semilogy(SNR_dB,simBer_sic,'mo-','LineWidth',2);
semilogy(SNR_dB,simber_zf,'p-','LineWidth',2);
axis([0 25 10^-5 0.5])
grid on
if(mod==1)
    legend('theory (nTx=2,nRx=2, ZF)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, ZF-SIC)','sim (nTx=2, nRx=2, only ZF)');
    xlabel('SNR in dB');
    ylabel('Bit Error Rate');
    title('BER for BPSK modulation with 2x2 MIMO and ZF-SIC equalizer (Rayleigh channel)');
else
    legend('theory (nTx=2,nRx=2, ZF)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, ZF-SIC)','sim (nTx=2, nRx=2, only ZF)');
    xlabel('SNR in dB');
    ylabel('symbol error rate(SER)');
    title('SER for QPSK modulation with 2x2 MIMO and ZF-SIC equalizer (Rayleigh channel)');

end

toc

function [sp]=bpsk_mod(N)
    sp=zeros(1,N);
    for i=1:N
        if(rand<0.5)
            sp(i)=1;
        else
            sp(i)=-1;
        end
    end
end
function [sq]=qpsk_mod(N)
    sq=zeros(1,N);
    for i=1:N
        temp=rand;
        if(temp<0.25)
            sq(i)=1+1i;
        elseif(temp<0.5)
            sq(i)=-1+1i;
        elseif(temp<0.75)
            sq(i)=-1-1i;
        else
            sq(i)=1-1i;
        end
    end
end
function det=qpsk_decision(yt)
        m1=(abs(yt-(1+1i)))^2;
        m2=(abs(yt-(-1+1i)))^2;
        m3=(abs(yt-(-1-1i)))^2;
        m4=(abs(yt-(1-1i)))^2;
        m=min([m1 m2 m3 m4]);
        if(m==m1)
            det=1+1i;
        elseif(m==m2)
            det=-1+1i;
        elseif(m==m3)
            det=-1-1i;
        else
            det=1-1i;
        end
end
 

