close all;
M = 16;                 % 16-QAM
n = log2(M);            %bit to QAM
%N = 96e4;                %data length
NoS = 10;
BER = zeros(1,NoS);
step = 5;
snr = 0:step:step*NoS;      %SNR = step*NoS  ---- for testing

%%
%Input Data   ---------------------------------------------------
                       %bitstream
%D_in = randi([0 1],1,N);
img = imread('Lenna.png');
gray = rgb2gray(img);
adj_img = imadjust(gray, [0.3,0.7],[]);
bw_img = imbinarize(adj_img);
ima_Tx = double(bw_img);
size_img = size(ima_Tx,1);
D_in = reshape(ima_Tx,1,[]);
length_data = length(D_in);

a = ceil(length_data/(48*n)) * (48*n) - length_data;
D_in = [D_in zeros(1,a)];
N = length(D_in);

%%
%Modulation------------------------------------------------------
D_in_Matrix = reshape(D_in,[], n);          
D_in_Symbol = bi2de(D_in_Matrix);

D_in_mod = qammod(D_in_Symbol,M);              %16-QAM?
N_mod = N/n;

cols = ceil(N_mod/48);                      %Number of subcarriers
%-----------------------------------------------------------------

%%
%Serial to Parallel---------------------------------------------
D_in_P = reshape(D_in_mod,48,[]);
%-----------------------------------------------------------------

%%
%Insert Pilot and Zeros symbols------------------------------------
D_in_P_inserted = zeros(64, cols);
D_in_P_inserted([8,22,46,60],:)       = 3+3i;

D_in_P_inserted(2:7,:)       = D_in_P(1:6,:);
D_in_P_inserted(9:21,:)      = D_in_P(7:19,:);
D_in_P_inserted(23:27,:)     = D_in_P(20:24,:);
D_in_P_inserted(39:45,:)     = D_in_P(25:31,:);
D_in_P_inserted(47:59,:)     = D_in_P(32:44,:);
D_in_P_inserted(61:64,:)     = D_in_P(45:48,:);
%-----------------------------------------------------------------

%%
%IFFT-------------------------------------------------------------
D_in_P_ifft = ifft(D_in_P_inserted, 64); 
%-----------------------------------------------------------------

%%
%Guard Interval Insertion-----------------------------------------
D_Tx_P = zeros(80, cols);
D_Tx_P(17:80,:) = D_in_P_ifft(:,:); 
D_Tx_P(1:16,:) = D_in_P_ifft(49:64,:);
%-----------------------------------------------------------------

%%
% Parallel to Serial----------------------------------------------
D_Tx = reshape(D_Tx_P,[],1);

%-----------------------------------------------------------------

%%
%                           RECEIVER
%                           --------
for index = 1:NoS
    
    %%
    %Add White Gaussian Noise - AWGN ---------------------------------
     D_Rx = awgn(D_Tx,snr(index));
     %D_Rx = D_Tx;

    %-----------------------------------------------------------------
    
    %%
    %Serial to Parallel ----------------------------------------------
    D_Rx_P = reshape(D_Rx,80,[]);

    %-----------------------------------------------------------------
    
    %%
    %Remove Guard Interval--------------------------------------------
    D_Rx_NoG = D_Rx_P(17:80,:); 
    %-----------------------------------------------------------------

    %%
    %FFT-------------------------------------------------------------
    D_Rx_fft =fft(D_Rx_NoG, 64); 
    %-----------------------------------------------------------------

    %%
    %Pilot Extraction-------------------------------------------------
    D_Rx_Pilot = zeros(4,cols);
    D_Rx_Pilot(1,:) = D_Rx_fft(8,:);
    D_Rx_Pilot(2,:) = D_Rx_fft(22,:);
    D_Rx_Pilot(3,:) = D_Rx_fft(46,:);
    D_Rx_Pilot(4,:) = D_Rx_fft(60,:);
    %-----------------------------------------------------------------

    %%
    %Remove Pilot + Zeros Signal--------------------------------------
    D_Rx_Remove = zeros(48,cols);

    D_Rx_Remove(1:6,:)      = D_Rx_fft(2:7,:);
    D_Rx_Remove(7:19,:)     = D_Rx_fft(9:21,:);
    D_Rx_Remove(20:24,:)    = D_Rx_fft(23:27,:);
    D_Rx_Remove(25:31,:)    = D_Rx_fft(39:45,:);
    D_Rx_Remove(32:44,:)    = D_Rx_fft(47:59,:);
    D_Rx_Remove(45:48,:)    = D_Rx_fft(61:64,:);
    %-----------------------------------------------------------------

    %%
    %Channel Estimation-----------------------------------------------
    D_Rx_Pilot_Abs = abs(D_Rx_Pilot);
    D_Rx_Pilot_est = D_Rx_Pilot_Abs./abs(3+3j);
    channel_est = mean(mean(D_Rx_Pilot_est),2);
    %-----------------------------------------------------------------

    %%
    %Data Calculation-------------------------------------------------
    D_Rx_cal = D_Rx_Remove./channel_est;
    %-----------------------------------------------------------------

    %%
    %Parallel to Serial---------------------------------------------             
    D_Rx_S = reshape(D_Rx_cal,[],1);
    %-----------------------------------------------------------------

    %%
    %Demodulation-----------------------------------------------------
    D_Rx_demod = qamdemod(D_Rx_S,M);              %16-QAM?

    D_Rx_Matrix = de2bi(D_Rx_demod);
    D_out = reshape(D_Rx_Matrix,1, []);
    %-----------------------------------------------------------------

    %%
    %BER Calculate----------------------------------------------------
    N_error = 0;
    
    for i = 1:N
        if D_out(i) ~= D_in(i)
            N_error = N_error + 1;
        end
    end
    BER(index) = N_error/N;
    %-----------------------------------------------------------------
    
end

%  semilogy(snr,BER,'o--');
%  grid on;
%  xlabel('SNR(dB)');
%  ylabel('BER');

D_out(:,(N-a+1):N) = [] ;
ima_Rx = reshape(D_out,size_img,[]);

subplot(121);
imshow(ima_Tx);
title('Tx image');

subplot(122);
imshow(ima_Rx);
title('Rx image'); 