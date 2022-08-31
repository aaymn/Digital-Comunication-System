%%%% 128 Qam modulation
clear all;
close all;

M=128; %% Symbols number
N=log2(M); % Number of Bits Per Symbol
n= 896*10; %no of generated bits
symbols_no=n/7;
Bit_stream= randi([0,1],1,n); 

Symbolm=reshape(Bit_stream,length(Bit_stream)/N,N);

Eb_No_db=[-10 -5 0 5 10 15 20];
Eb_No=[];
snr_db=[];
snr1=[];
Es=[];
Eavg=[];
Eb=[];
a=[];

SER=[];

% SNR=1000; %SNR=30dB
% SNR=10; %SNR=10dB%
% SNR=3.16227766; %SNR=5dB
% SNR=1; %SNR=0dB
% SNR=0.5; %SNR=-3dB

SNR=[1000 10 3.16227766 1 0.5];
ar=[];

p=0; %p=0 plots the constellation diagrams using the given SNR, while p=1 plots BER vs Eb/No

if p==0
    ar=SNR;
else
    ar=Eb_No_db;
end

No=2e-9;

for l=1:length(ar)
    yy=[];
    
    Eb_No(l)=power(10,(Eb_No_db(l)/10));
    snr1(l)=Eb_No(l)*N/2;  %N/2 = spectral efficiency


if p==0
    %Etot(l)=4*((4*symsum((power((((2*n1)-1)*a),2)),n1,1,8))+(8*symsum((power((((2*n1)-1)*a),2)),n1,1,4)));
    Eavgs(l)=(SNR(l)*No*2);
    a(l)=sqrt(Eavgs(l)/106);
    Eb(l)=Eavgs(l)/N;
else
    Eavgs(l)=(snr1(l)*No*2);
    a(l)=sqrt(Eavgs(l)/106);
    Eb(l)=Eavgs(l)/N;

end


for i=1:1:symbols_no
    x=Symbolm(i,[4:7]);
    y=Symbolm(i,[1:3]);
    
    if x==[0 0 0 0]
        s(i)=15*a(l);
    elseif x==[0 0 0 1]
        s(i)=13*a(l);
    elseif x==[0 0 1 1]
        s(i)=11*a(l);
    elseif x==[0 0 1 0]
        s(i)= 9*a(l);
    elseif x==[0 1 1 0]
        s(i)= 7*a(l);
    elseif x==[0 1 1 1]
        s(i)= 5*a(l);
    elseif x==[0 1 0 1]
        s(i)= 3*a(l);
    elseif x==[0 1 0 0]
        s(i)= 1*a(l);
    elseif x==[1 1 0 0]
        s(i)= -1*a(l);
     elseif x==[1 1 0 1]
        s(i)= -3*a(l);
     elseif x==[1 1 1 1]
        s(i)= -5*a(l);
     elseif x==[1 1 1 0]
        s(i)= -7*a(l);
     elseif x==[1 0 1 0]
        s(i)= -9*a(l);
     elseif x==[1 0 1 1]
        s(i)= -11*a(l);
     elseif x==[1 0 0 1]
        s(i)= -13*a(l);
     elseif x==[1 0 0 0]
        s(i)= -15*a(l);
    end
    
    if y==[0 0 0]
        q(i)=1j*7*a(l);
    elseif y==[0 0 1]
        q(i)=1j*5*a(l);
    elseif y==[0 1 1]
        q(i)=1j*3*a(l);
    elseif y==[0 1 0]
        q(i)=1j*1*a(l);
    elseif y==[1 1 0]
        q(i)=1j*-1*a(l);
    elseif y==[1 1 1]
        q(i)=1j*-3*a(l);
    elseif y==[1 0 1]
        q(i)=1j*-5*a(l);
    elseif y==[1 0 0]
        q(i)=1j*-7*a(l);
    
    end
    
sig=s(i)+q(i);
yy=[yy sig];
end

 
noise=[];

whitenoise = (randn(1,symbols_no)+ sqrt(-1)*randn(1,symbols_no))*sqrt(No/2); 

sig_with_noise=yy+whitenoise;
     

if p==0
    scatterplot(yy);
    xlabel('In-phase ({\Phi}1)');
    ylabel('Quadrature ({\Phi}2)');
    title ({'Constellation diagram of the signals';'at the input of the AWGN channel'});

    scatterplot(sig_with_noise);
    xlabel('In-phase ({\Phi}1)');
    ylabel('Quadrature ({\Phi}2)');
    title ({'Constellation diagram of the signals';' at the output of the AWGN channel'});
     
end     

sn=[];
qn=[];

for i=1:length(sig_with_noise)
    sn(i)=real(sig_with_noise(i));
    qn(i)=imag(sig_with_noise(i));
end


%%DEMAPER

s2=[];
q2=[];
ret_sig=[];


for i=1:length(sig_with_noise)
    min_len=sqrt((power((sn(i)-s(1)),2))+(power((qn(i)-imag(q(1))),2)));
    for j=1:length(yy)
       len=sqrt((power((sn(i)-s(j)),2))+(power((qn(i)-imag(q(j))),2)));
        if len<=min_len
            min_len=len;
            s2(i)=s(j);
            q2(i)=q(j);
            ret_sig(i)=s2(i)+q2(i);


        end
    end
end     



%SER
error=0;
for i=1:length(ret_sig)
    if ret_sig(i)~=yy(i)
        error=error+1;
    end
end

SER(l)=error/symbols_no;

x2=[];
y2=[];
symbolm2=[];

% SER THEORETICAL
Pes1(l)=(2*(1-(1/16))*qfunc((2*a(l))/(sqrt(2*No))))+(2*(1-(1/8))*qfunc((2*a(l))/(sqrt(2*No))));
Pes2(l)=(2*(1-(1/16))*qfunc((2*a(l))/(sqrt(2*No))))+(2*(1-(1/8))*qfunc((2*a(l))/(sqrt(2*No))))-(2*((2*(1-(1/16))*qfunc((2*a(l))/(sqrt(2*No)))))*(2*(1-(1/8))*qfunc((2*a(l))/(sqrt(2*No)))));


% BER THEORETICAL
Peb1(l)=Pes1(l)/N;
Peb2(l)=Pes2(l)/N;

for i=1:1:symbols_no
    
     if s2(i)==15*a(l)
        x2=[0 0 0 0];
    elseif s2(i)==13*a(l)
        x2=[0 0 0 1];
    elseif s2(i)==11*a(l)
        x2=[0 0 1 1];
    elseif s2(i)== 9*a(l)
        x2=[0 0 1 0];
    elseif s2(i)==7*a(l)
        x2=[0 1 1 0];
    elseif s2(i)== 5*a(l)
       x2=[0 1 1 1] ;
    elseif s2(i)== 3*a(l)
        x2=[0 1 0 1];
    elseif s2(i)== 1*a(l)
        x2=[0 1 0 0];
    elseif s2(i)== -1*a(l)
        x2=[1 1 0 0];
     elseif s2(i)== -3*a(l)
        x2=[1 1 0 1];
     elseif s2(i)== -5*a(l)
        x2=[1 1 1 1];
     elseif s2(i)== -7*a(l)
        x2=[1 1 1 0];
     elseif s2(i)== -9*a(l)
        x2=[1 0 1 0];
     elseif s2(i)== -11*a(l)
        x2=[1 0 1 1];
     elseif s2(i)== -13*a(l)
        x2=[1 0 0 1];
     elseif s2(i)== -15*a(l)
       x2=[1 0 0 0] ;
    end
    
    if q2(i)==1j*7*a(l)
        y2=[0 0 0];
    elseif q2(i)==1j*5*a(l)
        y2=[0 0 1];
    elseif q2(i)==1j*3*a(l)
       y2=[0 1 1];
    elseif q2(i)==1j*1*a(l)
        y2=[0 1 0];
    elseif q2(i)==1j*-1*a(l)
        y2=[1 1 0];
    elseif q2(i)==1j*-3*a(l)
        y2=[1 1 1];
    elseif q2(i)==1j*-5*a(l)
        y2=[1 0 1];
    elseif q2(i)==1j*-7*a(l)
        y2=[1 0 0];
    
    end    
    
symbolm2(i,[4:7])=x2;
symbolm2(i,[1:3])=y2;

end
   

%BER

error2=0;
for j=1:symbols_no
    for i=1:7
        if symbolm2(j,i)~= Symbolm(j,i)
        error2 = error2 + 1;
        end
    end


end

BER(l) = error2/n;


end

if p==1
    
 Eb_No=-10:5:20;
     
 figure()
 semilogy(Eb_No,BER);
 xlim([-10 20]);
 xlabel('Eb/No(dB)');
 ylabel('BER');
 title ('BER vs. Eb/No');
 

 figure()
 semilogy(Eb_No,Peb1);
 xlim([-10 20]);
 xlabel('Eb/No(dB)');
 ylabel('BER');
 title ('BER vs. Eb/No (Theoretical)');
 

 figure()
 semilogy(Eb_No,Peb2);
 xlim([-10 20]);
 xlabel('Eb/No(dB)');
 ylabel('BER');
 title ('BER vs. Eb/No (Theoretical)');
 
  figure()
 semilogy(Eb_No,Pes1);
 xlim([-10 20]);
 xlabel('Eb/No(dB)');
 ylabel('SER');
 title ('SER vs. Eb/No (Theoretical)');
 

 figure()
 semilogy(Eb_No,Pes2);
 xlim([-10 20]);
 xlabel('Eb/No(dB)');
 ylabel('SER');
 title ('SER vs. Eb/No (Theoretical)');
 
end



