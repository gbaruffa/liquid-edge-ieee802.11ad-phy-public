%% test several channel estimation methods
clear
clc

%% parameters
L_b = 128; 
L_s_spi = 128; 
L_s_sms = 64; 
N_S = 17;
N_b = 9;
L_fft = 512;
N_pack = 200;
L_h_est = 128;
%snrdb_list = 0:2:40; 
snrdb_list = 0:0.25:40;
%k_u_list = [0.5 1 2 4 8 16 32];
k_u_list = unique([10.^(-3:0.05:2.5) 1 4 32]);

%%% Select channel model below %%%
% h=[1; zeros(L_b - 1,1)];% ideal channel
%h = loadchan_but('../../data/measured_channels/CIR2', 1.7600e+09); chname = 'BUT2med'; % measured channel Brno
[h, n_i] = variablechannel(RandStream('mt19937ar', 'Seed', 823), 1); chname = 'rand823';% random sparse channel

%% initializations

% Ga128
g_a = [
  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1  1  1  1  ...
  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 ...
  -1 -1  1  1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 ...
  1 -1 -1  1
  ].';

% Gb128
g_b = [
  -1 -1  1  1  1  1  1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 ...
  -1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 ...
  -1 -1  1  1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  ...
  1 -1 -1  1
  ].';

s = [
  repmat(g_a, N_S - 1, 1)
  -g_a
  ];

c = [
  -g_b
  -g_a
  +g_b
  -g_a
  -g_b
  +g_a
  -g_b
  -g_a
  -g_b
  ];

S(:, 1) = [ 0 -1  0 -1  0 +1  0 -1  0].';
S(:, 2) = [-1  0 +1  0 -1  0 -1  0 -1].';

x = [
  s
  c
  ];

Chat = toeplitz(c, [c(1) s(end:-1:end - (L_h_est - 1) + 1).']);
pinvChat = pinv(Chat);
C1hat = toeplitz(x, [x(1) zeros(1, L_h_est - 1)]);
pinvC1hat = pinv(C1hat);
Ga = toeplitz([g_a; zeros((N_b - 1)*L_b, 1)], [g_a(1) zeros(1, N_b*L_b - 1)]);
Gb = toeplitz([g_b; zeros((N_b - 1)*L_b, 1)], [g_b(1) zeros(1, N_b*L_b - 1)]);
J_a = kron([S(1:8, 1).' 0], eye(L_b));
J_b = kron([S(1:8, 2).' 0], eye(L_b));

epsilon2_pi_list = zeros(size(snrdb_list));
epsilon2_stf_list = zeros(size(snrdb_list));
epsilon2_sxc_list = zeros(size(snrdb_list));
epsilon2_spi_list = zeros(size(snrdb_list));
epsilon2_sms_listlist = zeros(length(k_u_list), length(snrdb_list));
epsilon2_sms_list = zeros(size(snrdb_list));
epsilon2_lxc_list = zeros(size(snrdb_list));

tic
for snrdb = snrdb_list

  disp([num2str(toc, "%.1f") 's, snrdb=' num2str(snrdb)])
  hhat_story = zeros(L_b, 0);

  for np = 1:N_pack

    %% auxiliary channel parameters
    hnorm0 = sum(abs(h) > eps);
    hf = fft(h, L_fft);
    L_h = length(h);

    %% signal goes through the channel
    e = conv(x, h);

    %% add noise
    sigma2_w = var(e, 1)/10^(snrdb/10);
    w = sqrt(sigma2_w/2)*(randn(size(e)) + 1i*randn(size(e)));
    y = e + w;

    % receive and select samples
    p = y(N_S*L_b + (1:N_b*L_b));
    p1 = y(1:(N_S + N_b)*L_b);

    %% PI
    hhat_pi = pinvChat*p;
    hfhat_pi = fft(hhat_pi, L_fft);
    epsilon2_pi = channelNMSE(h, hhat_pi);
    epsilon2_pi_list(snrdb_list == snrdb) = epsilon2_pi_list(snrdb_list == snrdb) + epsilon2_pi;

    %% STF
    hhat_stf = pinvC1hat*p1;
    hfhat_stf = fft(hhat_stf, L_fft);
    epsilon2_stf = channelNMSE(h, hhat_stf);
    epsilon2_stf_list(snrdb_list == snrdb) = epsilon2_stf_list( snrdb_list == snrdb) + epsilon2_stf;

    %% LXC
    hhat_lxc = (1/(L_b*N_b))*Chat'*p;
    hfhat_lxc = fft(hhat_lxc, L_fft);
    epsilon2_lxc = channelNMSE(h, hhat_lxc);
    epsilon2_lxc_list(snrdb_list == snrdb) = epsilon2_lxc_list(snrdb_list == snrdb) + epsilon2_lxc;

    %% SXC
    hhat_sxc = (1/(8*L_b))*(J_a*Ga' + J_b*Gb')*p;
    hfhat_sxc = fft(hhat_sxc, L_fft);
    epsilon2_sxc = channelNMSE(h, hhat_sxc);
    epsilon2_sxc_list(snrdb_list == snrdb) = epsilon2_sxc_list(snrdb_list == snrdb) + epsilon2_sxc;

    %% SPI
    hhat_in = hhat_sxc;
    
    % Please select only one of the three options below by uncommenting
    % the three corresponding lines below

    % 1: non genie-aided (no a-priori knowledge of the channel)
    [~, highestIdxs] = maxk(abs(hhat_in), L_s_spi); 

    % 2: partial genie-aided (knowledge of number of largest taps but not their position)
    %L_s_spi = length(n_i); [~, highestIdxs] = maxk(abs(hhat_in), L_s_spi); 

    % 3: full genie-aided (full knowledge of number of largest taps and their position)
    %L_s_spi = length(n_i); highestIdxs = n_i.' + 1;  

    E = zeros(L_s_spi, L_h_est);
    E(sub2ind([L_s_spi L_h_est], 1:length(highestIdxs), highestIdxs.')) = 1; 

    hhat_spi = E.'*pinv(Chat*E.')*p;
    hf_hat_spi = fft(hhat_spi, L_fft);
    epsilon2_spi = channelNMSE(h, hhat_spi);
    epsilon2_spi_list(snrdb_list == snrdb) = epsilon2_spi_list(snrdb_list == snrdb) + epsilon2_spi;

    %% SMS
    for k_u = k_u_list
      hhat_in = hhat_sxc;


      % Please select only one of the three options below by uncommenting
      % the three corresponding lines below

      % 1: non genie-aided (no a-priori knowledge of the channel)
      [~, highestIdxs] = maxk(abs(hhat_in), L_s_sms);

      % 2: partial genie-aided (knowledge of number of largest taps but not their position)
      %L_s_sms = length(n_i); [~, highestIdxs] = maxk(abs(hhat_in), L_s_sms);

      % 3: full genie-aided (full knowledge of number of largest taps and their position)
      %L_s_sms = length(n_i); highestIdxs = n_i.' + 1;

      E = zeros(L_s_sms, L_h_est); 
      E(sub2ind([L_s_sms L_h_est], 1:length(highestIdxs), highestIdxs.')) = 1; 

      Uhat_in = diag(conj(hhat_in).*hhat_in);
      Gamma = E*Uhat_in*Chat';
      hhat_sms = (1/(k_u*sigma2_w))*E.'*(eye(L_s_sms) - Gamma*Chat*E.'/(Gamma*Chat*E.' + k_u*sigma2_w*eye(L_s_sms)))*Gamma*p;
      hfhat_sms = fft(hhat_sms, L_fft);
      epsilon2_sms = channelNMSE(h, hhat_sms);
      epsilon2_sms_listlist(k_u_list == k_u, snrdb_list == snrdb) = epsilon2_sms_listlist(k_u_list == k_u, snrdb_list == snrdb) + epsilon2_sms;
    end

    % exhaustive search of the best ku
    [~, kuminIdx] = min(epsilon2_sms_listlist(:, snrdb_list == snrdb));
    epsilon2_sms_list(snrdb_list == snrdb) = epsilon2_sms_listlist(kuminIdx, snrdb_list == snrdb);

    if np == N_pack
      %%
      figure(11)
      clf
      subplot(211)
      plot(0:L_fft - 1, abs(hf), 'DisplayName', 'Original')
      hold on
      plot(0:L_fft-1, abs(hfhat_pi), 'DisplayName', 'PI')
      plot(0:L_fft-1, abs(hfhat_stf), 'DisplayName', 'STF')
      plot(0:L_fft-1, abs(hfhat_lxc), 'DisplayName', 'LXC')
      plot(0:L_fft-1, abs(hfhat_sxc), 'DisplayName', 'SXC')
      plot(0:L_fft-1, abs(hf_hat_spi), 'DisplayName', 'SPI')
      plot(0:L_fft-1, abs(hfhat_sms), 'DisplayName', 'SMS')
      hold off
      xlim([0 (L_fft - 1)])
      xlabel('Frequency bin \itk\rm')
      ylabel('Magnitude')
      legend('show', 'NumColumns', 2)
      grid
      f1t1h = title(['\rmChannel frequency response, SNR=' num2str(snrdb) 'dB, k_u=' num2str(k_u)]);
      subplot(212)
      stem(0:L_h - 1, abs(h), 'DisplayName', 'Original')
      hold on
      stem(0:L_h_est - 1, abs(hhat_pi), 'DisplayName', 'PI', 'Marker', 'v')
      stem(0:L_h_est - 1, abs(hhat_stf), 'DisplayName', 'STF', 'Marker', 'x')
      stem(0:L_b - 1, abs(hhat_sxc), 'DisplayName', 'SXC', 'Marker', '<')
      stem(0:L_h_est - 1, abs(hhat_lxc), 'DisplayName', 'LXC', 'Marker', 's')
      stem(0:L_h_est - 1, abs(hhat_spi), 'DisplayName', 'SPI', 'Marker', '>')
      stem(0:L_h_est - 1, abs(hhat_sms), 'DisplayName', 'SMS', 'Marker', '^')
      hold off
      xlim([0 max([L_h L_b L_h_est])])
      xlabel('Time lag \itn\rm')
      ylabel('Magnitude')
      legend('show', 'NumColumns', 2)
      grid
      f1t2h = title(['\rmCIR, \rm||\bfh\rm||_0=' int2str(hnorm0)]);
    end
  end
end

% average the values
epsilon2_pi_list = epsilon2_pi_list/N_pack;
epsilon2_stf_list = epsilon2_stf_list/N_pack;
epsilon2_sxc_list = epsilon2_sxc_list/N_pack;
epsilon2_spi_list = epsilon2_spi_list/N_pack;
epsilon2_sms_listlist = epsilon2_sms_listlist/N_pack;
epsilon2_sms_list = epsilon2_sms_list/N_pack;
epsilon2_lxc_list = epsilon2_lxc_list/N_pack;
[~, kuminIdx_list] = min(epsilon2_sms_listlist);

%%
figure(11)
set([f1t1h f1t2h], 'Visible', 'off')
exportgraphics(gcf, "../results/chest1.pdf")
exportgraphics(gcf, "../results/chest1.png")
savefig(gcf, '../results/chest1.fig', 'compact')
set([f1t1h f1t2h], 'Visible', 'on')

%%
figure(21)
clf
semilogy(snrdb_list, epsilon2_pi_list, 'DisplayName', 'PI')
hold on
semilogy(snrdb_list, epsilon2_stf_list, 'DisplayName', 'STF')
semilogy(snrdb_list, epsilon2_sxc_list, 'DisplayName', 'SXC')
semilogy(snrdb_list, epsilon2_lxc_list, 'DisplayName', 'LXC')
semilogy(snrdb_list, epsilon2_spi_list, 'DisplayName', ['SPI \itL_s\rm=' int2str(L_s_spi)])
semilogy(snrdb_list, epsilon2_sms_list, 'DisplayName', ['SMS \itL_s\rm=' int2str(L_s_sms)], 'LineStyle', '--')
hold off
grid
xlabel("SNR (dB)")
ylabel('NMSE \it\epsilon^2')
legend('show', 'NumColumns', 2)
exportgraphics(gcf, "../results/chest2.pdf")
exportgraphics(gcf, "../results/chest2.png")
savefig(gcf, '../results/chest2.fig', 'compact')

%%
if length(k_u_list) > 1
  figure(31)
  clf
  surf(snrdb_list, k_u_list, epsilon2_sms_listlist)
  set(gca, 'ZScale', 'log', 'YScale', 'log')
  xlabel("SNR (dB)")
  ylabel('k_u')
  zlabel('NMSE \epsilon^2')
  view([60 30])
  hold on
  plot3(snrdb_list, k_u_list(kuminIdx_list), epsilon2_sms_listlist(sub2ind(size(epsilon2_sms_listlist), kuminIdx_list, 1:length(snrdb_list))), 'r-', ...
    'LineWidth', 2)
end

%%
if length(k_u_list) > 1
  [xData, yData] = prepareCurveData( snrdb_list, k_u_list(kuminIdx_list) );
  ft = fittype('poly2');
  [fitresult, gof] = fit( xData, yData, ft );
  k_u_opt = fitresult.p1*(snrdb_list).^2 + fitresult.p2*(snrdb_list) + fitresult.p3; % for channel (58)
  figure(4)
  clf
  plot(snrdb_list, k_u_list(kuminIdx_list), 'DisplayName', chname)
  hold on
  plot(snrdb_list, k_u_opt, 'r--', 'DisplayName', num2str([fitresult.p1 fitresult.p2 fitresult.p3]))
  hold off
  grid on
  xlabel('SNR (dB)')
  ylabel('\itk_u^*')
  legend('show')
  exportgraphics(gcf, "../results/chest3.pdf")
  exportgraphics(gcf, "../results/chest3.png")
  savefig(gcf, '../results/chest3.fig', 'compact')
end

%% generates a random channel (possibly) different each time
function [h, n_i] = variablechannel(strm, M)
SC = 20;
hpos = upsample(randerr(1, 128, SC, strm), M);
n_i = find(hpos) - 1;
h = sqrt(1/2)*(randn(strm, length(hpos), 1) + 1i*randn(strm, length(hpos), 1)).*(hpos.*exp(-0.01*(0:length(hpos) - 1)))';
end

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

%% NMSE of channel
function e = channelNMSE(h, h_est)
lmax = max(length(h), length(h_est));
e = sum(abs([h; zeros(lmax - length(h), 1)] - [h_est; zeros(lmax - length(h_est), 1)]).^2) / sum(abs(h).^2);
end
