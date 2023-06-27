% SIMS for the paper figures
clear
clc
cd ..

L_PS = 1500; % N_PSDU, Length of data octets (1–262143)
ftypetx = 'RRC'; % TX filtering type: RRC or GAU
ftyperx = 'RRC'; % RX filtering type: RRC or GAU
M_o = 1; % channel oversampling factor (i.e., "analog" signal rate)
numpacks = 500; % number of packets in each cycle

%%% Select channel model below %%%
%h_ch = 1; chname = 'id'; % ideal channel
%h_ch = loadchan_but('../data/measured_channels/CIR2', M_o*ieee80211ad.Fc); chname = 'BUT2med'; % measured channel Brno
h_ch = variablechannel(RandStream('mt19937ar', 'Seed', 823), M_o); chname = 'rand823'; % random sparse channel

usepreamble = true; % generate (and transmit) the preamble part
useheader = false; % generate (and transmit) the header part
usedata = true; % generate (and transmit) the data part
chestimate = 7; % estimate the channel:
                % 1 = perfect channel, 4 = LXC, 5 = PI,
                % 6 = SXC, 7 = SMS
equalize = 0; % equalize the signal:
              % 0 = recovers delay/amplitude/phase of the largest tap,
              % 1 = ZFEQ, 2 = MMSEQ
ibicancel = 0; % IBI cancellation (0 = disabled, 1 = active)
chshorten = 0; % channel shortening (0 = disabled, 1 = active)
chshalpha = 0; % unused parameter
savepath = 'results/but2_'; % path where figs and pdfs are saved
llrmethod = 1; % LLR update method: 0 = conventional, 1 =  proposed

save common_parms.mat

%%

for MCS = [1 5 9 12]

  if MCS>=1&&MCS<=5
    sl = -6:1:19;  % BPSK
  elseif MCS>=6&&MCS<=9
    sl=-2:1:24; % QPSK
  elseif MCS>=10&&MCS<=12
    sl=4:1:24; % 16QAM
  end

  load common_parms.mat
  equalize = 2;
  ibicancel = 0;
  chshorten = 0;
  SNRdBlist = sl;
  save parms.mat
  main

  load common_parms.mat
  equalize = 2;
  ibicancel = 1;
  chshorten = 0;
  SNRdBlist = sl;
  save parms.mat
  main

  load common_parms.mat
  equalize = 2;
  ibicancel = 0;
  chshorten = 1;
  SNRdBlist = sl;
  save parms.mat
  main
end

%% final draw
clear
clc
figure(807)
clf
hv = "off";
spath = 'results/but2_'; % path where figs and pdfs are saved
load([spath '18001.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'bs-.', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
hold on
load([spath '18101.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'b*--', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18051.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'bv:', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18005.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'rs-.', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18105.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'r*--', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18055.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'rv:', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18009.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'ms-.', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18109.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'm*--', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18059.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'mv:', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18012.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'gs-.', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18112.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'g*--', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
load([spath '18062.mat'])
semilogy(SNRdBlist(1:length(coded_ber_list)) - 0*10*log10(txrxnode.L_MW/txrxnode.L_CW), coded_ber_list, 'gv:', ...
  'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
  'HandleVisibility', hv)
if hv == "on"
  legend('show', 'Location', 'eastoutside', 'FontSize', 8, 'Numcolumns', 1)
  title(['\rm' ftypetx '-' ftyperx ', ' chname ', chest = ' int2str(chestimate)])
else
  plot(NaN, NaN, 'b', 'DisplayName', "MCS1")
  plot(NaN, NaN, 'r', 'DisplayName', "MCS5")
  plot(NaN, NaN, 'm', 'DisplayName', "MCS9")
  plot(NaN, NaN, 'g', 'DisplayName', "MCS12")
  plot(NaN, NaN, 'ks-.', 'DisplayName', "MMSEQ")
  plot(NaN, NaN, 'k*--', 'DisplayName', "MMSEQ+IBIC")
  plot(NaN, NaN, 'kv:', 'DisplayName', "MMSEQ+CS")
  legend('show', 'Location', 'SouthWest', 'FontSize', 6, 'Numcolumns', 3)
  title('')
end
hold off
grid
xlabel('SNR (dB)')
ylabel('BER')
ylim([1e-6 1])
xlim([-6 24])
basename = [spath 'chest' int2str(chestimate)];
exportgraphics(gcf, [basename '.pdf'])
exportgraphics(gcf, [basename '.png'])
savefig(gcf, [basename '.fig'], 'compact')

%% cleanup
delete("common_parms.mat")
cd batch

%% load measured channel from MAT file
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

% generates a random channel different each time
function [h, n_i] = variablechannel(strm, M)
SC = 20;
hpos = upsample(randerr(1, 128, SC, strm), M);
n_i = find(hpos) - 1;
h = sqrt(1/2)*(randn(strm, length(hpos), 1) + 1i*randn(strm, length(hpos), 1)).*(hpos.*exp(-0.01*(0:length(hpos) - 1)))';
end
