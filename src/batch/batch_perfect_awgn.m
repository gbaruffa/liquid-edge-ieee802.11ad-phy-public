% SIMS for the paper figures
clear
clc
cd ..

% common parameters
L_PS = 1500; % N_PSDU, Length of data octets (1–262143)
ftypetx = 'RRC'; % TX filtering type: RRC or GAU
ftyperx = 'RRC'; % RX filtering type: RRC or GAU
M_o = 1; % channel oversampling factor (i.e., "analog" signal rate)
numpacks = 500; % number of packets in each cycle

%%% Select channel model below %%%
h_ch = 1; chname = 'id'; % ideal channel
%h_ch = loadchan_but('../data/measured_channels/CIR2', M_o*ieee80211ad.Fc); chname = 'BUT2med'; % measured channel Brno
%h_ch = variablechannel(RandStream('mt19937ar', 'Seed', 823), M_o); chname = 'rand823'; % random sparse channel

usepreamble = false; % generate (and transmit) the preamble part
useheader = false; % generate (and transmit) the header part
usedata = true; % generate (and transmit) the data part
chestimate = 1; % estimate the channel:
                % 1 = perfect channel, 4 = LXC, 5 = PI,
                % 6 = SXC, 7 = SMS
equalize = 0; % equalize the signal:
              % 0 = recovers delay/amplitude/phase of the largest tap,
              % 1 = ZFEQ, 2 = MMSEQ
ibicancel = 0; % IBI cancellation (0 = disabled, 1 = active)
chshorten = 0; % channel shortening (0 = disabled, 1 = active)
chshalpha = 0; % unused parameter
llrmethod = 0; % LLR update method: 0 = conventional, 1 =  proposed
savepath = 'results/awgn2_'; % path where figs and pdfs are saved

save common_parms.mat

%%
MCSlist = 1:12;
for MCS = MCSlist

  if MCS>=1&&MCS<=5
    sl = -6:0.5:6;  % BPSK
  elseif MCS>=6&&MCS<=9
    sl=-2:0.5:9; % QPSK
  elseif MCS>=10&&MCS<=12
    sl=4:0.5:15; % 16QAM
  end

  load common_parms.mat
  SNRdBlist = sl;
  save parms.mat
  main

end

%% final draw
clear
clc
figure(901)
clf
hv = "off";
spath = 'results/awgn2_'; % path where figs and pdfs are saved
MCSlist = 1:12;
stylist = {'bo-' 'bs-' 'b^-' 'bv-' 'bh-' 'rs-' 'r^-' 'rv-' 'rh-' 'gs-' 'g^-' 'gv-'};
for MCS = 1:12
  load([spath int2str(1000+MCS) '.mat'])
  semilogy(SNRdBlist(1:length(coded_ber_list)), coded_ber_list, stylist{MCS}, ...
    'DisplayName', ['M' int2str(MCS) 'e' int2str(equalize) 'i' int2str(ibicancel) 's' int2str(chshorten)], ...
    'HandleVisibility', hv)
  hold on
end
if hv == "on"
  legend('show', 'Location', 'northeast', 'FontSize', 6, 'Numcolumns', 1)
  title(['\rm' ftypetx '-' ftyperx ', ' chname ', chest = ' int2str(chestimate)])
else
  plot(NaN, NaN, 'b', 'DisplayName', "BPSK")
  plot(NaN, NaN, 'r', 'DisplayName', "QPSK")
  plot(NaN, NaN, 'g', 'DisplayName', "16QAM")
  plot(NaN, NaN, 'bo', 'DisplayName', "1/2 \rho=2")
  plot(NaN, NaN, 'ks', 'DisplayName', "1/2")
  plot(NaN, NaN, 'k^', 'DisplayName', "5/8")
  plot(NaN, NaN, 'kv', 'DisplayName', "3/4")
  plot(NaN, NaN, 'kh', 'DisplayName', "13/16")
  legend('show', 'Location', 'northeast', 'FontSize', 6, 'Numcolumns', 3)
  title('')
end
hold off
grid
xlabel('SNR (dB)')
ylabel('BER')
ylim([1e-6 1])
xlim([-6 13])
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
