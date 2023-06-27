%% Main development script for the IEEE 802.11ad PHY simulator
% used in the PRIN 2017 Liquid Edge project

% The main reference is the standard:
%
% "IEEE Standard for Information technology--Telecommunications and information
% exchange between systems--Local and metropolitan area networks--Specific
% requirements-Part 11: Wireless LAN Medium Access Control (MAC) and Physical
% Layer (PHY) Specifications Amendment 3: Enhancements for Very High Throughput
% in the 60 GHz Band," in IEEE Std 802.11ad-2012 (Amendment to IEEE Std 802.11-2012,
% as amended by IEEE Std 802.11ae-2012 and IEEE Std 802.11aa-2012), pp.1-628,
% 28 Dec. 2012, doi: 10.1109/IEEESTD.2012.6392842.


%% clear workspace
clear
clc

%% %%%%%%%%%%%%
% parameters %%
%%%%%%%%%%%%%%%

if exist("parms.mat", 'file')

  % defined outside of this file
  load parms.mat
  delete("parms.mat")

else
  % defined inside of this file

  % standard variable parameters
  MCS = 12; % SC MCS index (1-5=BPSK, 6-9=QPSK, 10-12=16QAM)
  L_PS = 1500; % N_PSDU, Length of data octets (1–262143)

  % simulation variable parameters
  ftypetx = 'RRC'; % TX filtering type: RRC or GAU
  %ftypetx = 'GAU'; % TX filtering type: RRC or GAU
  ftyperx = 'RRC'; % RX filtering type: RRC or GAU
  %ftyperx = 'GAU'; % RX filtering type: RRC or GAU
  M_o = 1; % channel oversampling factor (i.e., "analog" signal rate)
  numpacks = 5000; % number of packets in each cycle

  %%% Select channel model below %%%
  %h_ch = 1; chname = 'id'; % ideal channel
  %h_ch = loadchan_but('../data/measured_channels/CIR2', M_o*ieee80211ad.Fc); chname = 'BUT2med'; % measured channel Brno
  h_ch = variablechannel(RandStream('mt19937ar', 'Seed', 823), M_o); chname = 'rand823'; % random sparse channel

  %SNRdBlist = -6:1:8; % average SNR for MCS1-5 - coarse
  %SNRdBlist = -4:0.5:2; % average SNR for MCS1-5 - fine
  %SNRdBlist = -1:0.25:1; % average SNR for MCS1-5 - ultrafine
  %SNRdBlist = -2:1:10; % average SNR for MCS6-9 - coarse
  %SNRdBlist = -2:0.5:4; % average SNR for MCS6-9 - fine
  %SNRdBlist = 6:1:20; % average SNR for MCS10-12 - coarse
  SNRdBlist = 6:0.5:20; % average SNR for MCS10-12 - fine
  %SNRdBlist = 7:0.25:9; % average SNR for MCS10-12 - ultrafine
  %SNRdBlist = -2:1:24;
  usepreamble = true; % generate (and transmit) the preamble part
  useheader = false; % generate (and transmit) the header part
  usedata = true; % generate (and transmit) the data part
  chestimate = 6; % estimate the channel:
                  % 1 = perfect channel, 4 = LXC, 5 = PI,
                  % 6 = SXC, 7 = SMS
  equalize = 2; % equalize the signal:
                % 0 = recovers delay/amplitude/phase of the largest tap,
                % 1 = ZFEQ, 2 = MMSEQ
  ibicancel = 1; % IBI cancellation (0 = disabled, 1 = active)
  chshorten = 0; % channel shortening (0 = disabled, 1 = active)
  chshalpha = 0; % unused parameter
  llrmethod = 1; % LLR update method: 0 = conventional, 1 =  proposed
  savepath = 'results/'; % path where figs and pdfs are saved

end

% debugging parameters
debugplot = true; % draw plots for debugging
debugplot2 = true; % draw plots for debugging
debugpow = false; % draw more plots of PSDs and print powers along several points of the RX chain
debugber = false; % print BERs more frequently

% utility parameters
N_fft = 512;
threshold = 1e9;

%% %%%%%%%%%%%%%
% error check %%
%%%%%%%%%%%%%%%%

% adjustments
if length(SNRdBlist) > 1
  debugplot = false;
  debugpow = false;
else
  numpacks = 1;
end

%% cycle on the SNR
bpskberdet = [];
qpskberdet = [];
qam16berdet = [];
SNRdBdetlist = nan(size(SNRdBlist));
SNRdBdet2list = nan(size(SNRdBlist));
coded_ber_list = nan(size(SNRdBlist));
coded_per_list = nan(size(SNRdBlist));
uncoded_ber_list = nan(size(SNRdBlist));
avg_niters_list = nan(size(SNRdBlist));
num_coded_bit_errors_list = nan(size(SNRdBlist));
for SNRdB = SNRdBlist

  %% cycle on the packets
  num_uncoded_bit_errors = 0;
  num_uncoded_bits = 0;
  uncoded_ber = nan;
  num_coded_bit_errors = 0;
  num_coded_bits = 0;
  coded_ber = nan;
  num_coded_packet_errors = 0;
  num_coded_packets = 0;
  coded_per = nan;
  avg_niters = 0;
  tic;
  txrxnode = ieee80211ad(MCS, ftypetx, ftyperx, M_o, usepreamble, useheader, ...
    usedata, debugplot, debugpow, debugber, chestimate, h_ch, equalize, ...
    threshold, ibicancel, chshorten, chshalpha, llrmethod, SNRdB);
  [sigma2_w, p_S, p_I, p_N, p_X, p_E, p_S1, p_S2, p_I1, p_I2, p_N1, p_N2] = txrxnode.calcpowers(db2pow(SNRdB));
  disp(['S=' num2str(mean(p_S)) ', S1=' num2str(mean(p_S1)) ', S2=' num2str(mean(p_S2)) ', I=' num2str(mean(p_I)) ', N=' num2str(mean(p_N)) ', S+I=' num2str(mean(p_S + p_I)) ', N+I=' num2str(mean(p_N + p_I)) ...
    ', S+I+N=' num2str(mean(p_S + p_I + p_N)) ', S/(N+I)=' num2str(pow2db(mean(p_S)/mean(p_I + p_N))) ...
    'dB, (S+I)/N=' num2str(pow2db(mean(p_S + p_I)/mean(p_N))) 'dB, S/N=' num2str(pow2db(mean(p_S)/mean(p_N))) 'dB'])
  SNRdBdetlist(SNRdBlist == SNRdB) = pow2db(mean(p_S)/(mean(p_N)/2 + mean(p_I)/2));
  SNRdBdet2list(SNRdBlist == SNRdB) = pow2db(mean(p_S)/(mean(p_N)/2 + mean(p_I)/2));
  for np = 1:numpacks

    % random information bits
    b = randi([0 1], 8*L_PS, 1, 'int8');

    % generate a frame
    [x_up, Lambda, r_b_th, r_b_act] = txrxnode.transmit(b);

    % channel effect
    e_up = conv(x_up, h_ch);

    % add WG noise
    %sigma2_w = 1;
    w = sqrt(sigma2_w/2)*(randn(size(e_up)) + 1i*randn(size(e_up)));
    y_up = 1*e_up + 1*w;

    if debugpow
      figure(301)
      clf
      plot(0:N_fft - 1, pow2db(pwelch(x_up, N_fft, 1)*2*pi), 'DisplayName', 'x')
      hold on
      plot(0:N_fft - 1, pow2db(pwelch(e_up, N_fft, 1)*2*pi), 'DisplayName', 'e')
      plot(0:N_fft - 1, pow2db(pwelch(y_up, N_fft, 1)*2*pi), 'DisplayName', 'y')
      hold off
      legend('show')
      set(gca, 'XLim', [0 (N_fft-1)])
      grid
      fprintf('x(≈%.4f)=%.4f\n', p_X, var(x_up, 1))
      fprintf('e(≈%.4f)=%.4f\n', p_E, var(e_up, 1))
      disp(['w(≈' num2str(sigma2_w) ')=' num2str(var(w, 1))])
      disp(['y(≈' num2str(p_E + sigma2_w) ')=' num2str(var(y_up, 1))])
    end

    % receive the frame
    [bhat, Lambdahat, avgni] = txrxnode.receive(y_up, L_PS, sigma2_w, p_S, p_I, p_N, SNRdB, p_S1, p_S2, p_I1, p_I2, p_N1, p_N2);
    avg_niters = avg_niters + avgni;

    if debugplot
      figure(4)
      clf
      [X_up, f] = pwelch(x_up, blackman(1024), 128, 1024, M_o*txrxnode.Fc, 'centered', 'psd');
      plot(f/1e9, pow2db(X_up), 'DisplayName', 'TX spectrum')
      hold on
      toplevdB = pow2db(max(X_up));
      plot([-8 -3.06 -2.7 -1.2 -0.94 +0.94 +1.2 +2.7 +3.06 +8], toplevdB + [-30 -30 -22 -17 0 0 -17 -22 -30 -30], 'm', 'DisplayName', 'TX mask')
      hold off
      grid
      xlim(M_o*txrxnode.Fc*([-0.5 0.5])/1e9)
      xlabel('Frequency (GHz)')
      ylabel('PSD (dB)')
      title(ftypetx)
      legend('show', 'Location', 'SouthEast')
    end

    %% error rates

    % uncoded
    num_uncoded_bit_errors = num_uncoded_bit_errors + sum(Lambda(:) ~= int8(Lambdahat(:) < 0));
    num_uncoded_bits = num_uncoded_bits + numel(Lambda);
    uncoded_ber = num_uncoded_bit_errors / num_uncoded_bits;

    % coded
    num_coded_bit_errors = num_coded_bit_errors + sum(b ~= bhat);
    num_coded_bits = num_coded_bits + numel(b);
    coded_ber = num_coded_bit_errors / num_coded_bits;
    num_coded_packet_errors = num_coded_packet_errors + any(b ~= bhat);
    num_coded_packets = num_coded_packets + 1;
    coded_per = num_coded_packet_errors / num_coded_packets;

    if debugber || np == numpacks || np == 1 || ~mod(np, 40)
      % print
      disp([int2str(toc) 's, NP=' int2str(np) ', SNR=' num2str(SNRdB) 'dB, BERc=' num2str(coded_ber, '%.1e') ...
        ', BR=' num2str(r_b_act/1e6, '%.0f') '(' num2str(r_b_th/1e6, '%.0f') ')Mbps, ERRc=' int2str(num_coded_bit_errors)])
    end

    % check when to stop
    if ~mod(np, 40) && ((np >= 400 && num_coded_bit_errors >= 20000) || (toc >= 2400))
      break
    end

  end

  %% collect stats for plotting
  coded_ber_list(SNRdBlist == SNRdB) = coded_ber;
  coded_per_list(SNRdBlist == SNRdB) = coded_per;
  uncoded_ber_list(SNRdBlist == SNRdB) = uncoded_ber;
  avg_niters_list(SNRdBlist == SNRdB) = avg_niters / np;
  num_coded_bit_errors_list(SNRdBlist == SNRdB) = num_coded_bit_errors;

  resnum = llrmethod*10000 + chestimate*1000 + equalize*500 + ibicancel*100 + chshorten * 50 + MCS;

  if debugplot2
    figure(resnum)
    clf
    semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/ieee80211ad.L_CW), coded_ber_list, 'o-', ...
      'DisplayName', ['BER cod. (sim.) MCS' int2str(MCS)])
    hold on
    semilogy(SNRdBlist(1:length(uncoded_ber_list)), uncoded_ber_list, '*', ...
      'DisplayName', ['BER unc. (sim.) MCS' int2str(MCS)])
    gamma = db2pow(SNRdBlist(1):0.1:SNRdBlist(end));
    gammadet = db2pow(SNRdBdetlist);
    gammadet2 = db2pow(SNRdBdet2list);

    if txrxnode.nu == 1
      bpskber = qfunc(sqrt(2*gamma));
      semilogy(pow2db(gamma), bpskber, 'r-', 'DisplayName', ['BER unc. (th.) AWGN MCS' int2str(MCS)])
      bpskberdet = qfunc(sqrt(gammadet));
      semilogy(SNRdBlist, bpskberdet, 'r--', 'DisplayName', ['BER unc. (th.) MCS' int2str(MCS)])
    elseif txrxnode.nu == 2
      qpskber = qfunc(sqrt(gamma));
      semilogy(pow2db(gamma), qpskber, 'r-', 'DisplayName', ['BER unc. (th.) AWGN MCS' int2str(MCS)])
      qpskberdet = 0.5*qfunc(sqrt(gammadet/2)) + 0.5*qfunc(sqrt(gammadet2/2));
      semilogy(SNRdBlist, qpskberdet, 'r--', 'DisplayName', ['BER unc. (th.) MCS' int2str(MCS)])
    elseif txrxnode.nu == 4
      qam16ber = 0.75*qfunc(sqrt(0.2*gamma));
      semilogy(pow2db(gamma), qam16ber, 'r-', 'DisplayName', ['BER unc. (th.) AWGN MCS' int2str(MCS)])
      qam16berdet = 0.375*qfunc(sqrt(0.2*gammadet/2)) + 0.375*qfunc(sqrt(0.2*gammadet2/2));
      semilogy(SNRdBlist, qam16berdet, 'r--', 'DisplayName', ['BER unc. (th.) MCS' int2str(MCS)])
    end

    semilogy(SNRdBlist(1:length(coded_per_list)), coded_per_list, '+-', ...
      'DisplayName', ['FER (th.) MCS' int2str(MCS)])
    semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/ieee80211ad.L_CW), min(1, coded_ber_list*L_PS*8), ...
      '^-', 'DisplayName', ['FER (sim.) MCS' int2str(MCS)])
    hold off
    grid
    legend('show', 'Location', 'Southwest', 'FontSize', 8)
    xlabel('SNR (dB)')
    ylabel('Error rate')
    ylim([1e-6 1])
    title(['\rm' ftypetx '-' ftyperx ', MCS' int2str(MCS) ', chest=' int2str(chestimate) ...
      ', eq=' int2str(equalize) ', chsh=' int2str(chshorten) ', ibic=' int2str(ibicancel) ', ' chname])
    savefig(gcf, [savepath int2str(resnum) '.fig'], 'compact')
    exportgraphics(gcf, [savepath int2str(resnum) '.pdf'])
    exportgraphics(gcf, [savepath int2str(resnum) '.png'])
    drawnow
  end

  save([savepath int2str(resnum) '.mat'], ...
    'SNRdBdetlist', 'SNRdBdet2list', ...
    'uncoded_ber_list', 'coded_ber_list', 'coded_per_list', 'avg_niters_list', ...
    'bpskberdet', 'qpskberdet', 'qam16berdet', 'r_b_act', 'r_b_th', ...
    'txrxnode', 'MCS', 'L_PS', 'ftypetx', 'ftyperx', 'M_o', 'numpacks', 'h_ch', ...
    'SNRdBlist', 'usepreamble', 'useheader', 'usedata', 'chestimate', 'equalize', 'ibicancel', ...
    'savepath', 'debugplot', 'debugplot2', 'debugpow', 'debugber', 'chshorten', ...
    'chname', 'num_coded_bit_errors_list', 'llrmethod')

end


% load measured channel from MAT file
function hp = loadchan_but(matname, fs)
S = load([matname '.mat']);
fs_orig = 10e9;
h_orig = squeeze(S.CIR_mat(1,1,:));
t_orig = (0:length(h_orig) - 1)/fs_orig;
Tmax = length(h_orig) / fs_orig;
K = ceil(fs_orig/fs);
fs_i = fs*K;
h_i = zeros(K*ceil(Tmax*fs_i/K), 1);
h_i(round(t_orig*fs_i + 1)) = h_orig;
hp = sum(reshape(h_i, K, []), 1).';
hp = hp/sqrt(hp'*hp);
while abs(hp(1)) < 1e-2
  hp(1) = [];
end
hp(129:end) = [];
hp = hp/sqrt(hp'*hp);
end

% generates a random channel, possibly different each time
function [h, n_i] = variablechannel(strm, M)
SC = 20;
hpos = upsample(randerr(1, 128, SC, strm), M);
n_i = find(hpos) - 1;
h = sqrt(1/2)*(randn(strm, length(hpos), 1) + 1i*randn(strm, length(hpos), 1)).*(hpos.*exp(-0.01*(0:length(hpos) - 1)))';
end