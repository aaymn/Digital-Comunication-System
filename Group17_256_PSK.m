%256-PSK
clear all;
close all;

%decleration
M=256;
N=log2(M);
x=8000;
data = rand(1,x);
mat=[];
SNR_dB=[30 10 5 0 -3];
SNR=10.^(SNR_dB./10);
Eb_No_dB=[-10 -5 0 5 10 15 20];
Eb_No=10.^(Eb_No_dB./10);
whitenoise = randn(1,x/8)+1i*randn(1,x/8);
BER=[]; rng=(0); norm=[]; index=[]; Rx_dec=[];

%for loop to take the random stream and arrange it in an array
for f = 1:(x/N)      %each symbol is 8 bits
    i=1;             % i should be initialized each round 
    for t=(8*f)-7:(8*f)
        mat(f,i)=data(t); 
        i=i+1;
    end           
end

%stage 1 (change gray code to binary code)
g=mat>0.5; %gray array 
b(:,1) = g(:,1); %copying the 1st  gray bit as it is, to the resultant vector
for j = 2:size(g,2) %i is the column , size(g,2) means the size of the second dimension , size of the column
    b(:,j) = xor( b(:,j-1), g(:,j) ); %XORing the following bits as per the procedure 
end 

%stage 2 (get (real+i imaginary) from the 8 bit of each row)
 for r = 1:size(b,1) %for loop to pass on each row
     for c = 1:size(b,2) %for loop to pass on each column
         z(r)=( (power(2,0))*b(r,8) )+( (power(2,1))*b(r,7) )+( (power(2,2))*b(r,6) )+( (power(2,3))*b(r,5) )+( (power(2,4))*b(r,4) )+((power(2,5))*b(r,3))+((power(2,6))*b(r,2))+( (power(2,7))*b(r,1) );
         y(r)=cos( ((2*pi)/(256))*z(r) ) + 1i*(sin( ((2*pi)/(256))*z(r )) ); %the output
      end
 end
 scatterplot(y);
 title('Constellation Diagram Tx');
 xlabel('In-phase ({\Phi}1)');
 ylabel('Quadrature ({\Phi}2)');
 
 %stage 3 adding noise AWGN 
for i=1:length(SNR)
     n=y+(whitenoise/sqrt(SNR(i)./4));
         scatterplot(n);
         title('Constellation Diagram After Channel');
         xlabel('In-phase ({\Phi}1)');
         ylabel('Quadrature ({\Phi}2)');
end 

%stage 4 demapper
for h=1:length(Eb_No)
     n=y+(whitenoise/sqrt(Eb_No(h)*4));   %take Eb/No in account
 for j=1:length(n)                    %loop for noisy symbols
   for i=1:length(y)                    %loop for Tx symbols
       dis=y(i)- n(j);                
       norm(i)=abs(dis);
   end
   g=min(norm);                   %min distance calculation btn Rx symbol and Tx symbol " region determination"
   for k=1:length(n)
   if (norm(k)==min(norm))
   index(j)=k;                     %index for tx symbol achieve min distance
   end
end
 end
 q=ceil(length(index)/256);
 for i=1:length(index)
     Rx_dec(i)=(angle(y(index(i))))*(256/360)*(180/pi);      %Rx symbol region decleration 
     for p= 1:q
      if(Rx_dec(i)>(q*256))
                 Rx_dec(i)= Rx_dec(i)-(q*256);
      end
      q=q-1;
     end
     if (Rx_dec(i)<0)
        Rx_dec(i)= Rx_dec(i)+256;
     end
 end
        Rx_dec=round(Rx_dec); 
        row_num=de2bi(Rx_dec,8,'left-msb');   % row_num is the Rx symbols
        count=0;
 for i=1:size(row_num,1)
     for j=1:size(row_num,2)                     % simulation BER
         BER_Rx(i,j)=abs(row_num(i,j)-b(i,j));       
         if (BER_Rx(i,j)>0)
             count=count+1;
         end 
     end
 end
 Err_Bits(h)=count/(size(BER_Rx,1)*size(BER_Rx,2)) ;
 BER(h)=0.25*qfunc(sqrt(16*(Eb_No(h)))*sin(pi/256));  % theortical BER
end

figure(7)
semilogy(Eb_No_dB,Err_Bits);
hold on;
semilogy(Eb_No_dB,BER);
ylim([1e-2 1]);
title('Bit Error Rate for 256-PSK Modulation');
legend('Simulated','Theortical');
xlabel('Eb/No (in dB)'); ylabel('BER');
grid on;
hold off;