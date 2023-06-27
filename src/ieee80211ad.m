classdef ieee80211ad < handle
  % A class representing an IEEE 802.11ad modem.

  properties (Access = private)
    pnIC
    kappa
    pnseqrho2
    N_CBPB
    rho
    r_C
    g_a % Ga_128
    g_b % Gb_128
    gamma_a % Ga_64
    L_G
    Tc
    MCS = -1
    ldpcencoderD
    r_STF
    r_CEF
    ldpcencoderH
    v
    L_v
    L_B
    h_ch
    s_preamble
    ldpcdecoder
    Ga_128_L
    Gb_128_L
    Chat
    h
    L_h
    hf
    q
    h_sh
    vf
    hf_sh
    H_tx
  end

  properties (Access = public)
    ftypetx
    ftyperx
    M_o
    h_tx
    h_rx
    L_MW
    nu
    L_D
    usepreamble
    useheader
    usedata
    chestimate
    equalize
    threshold
    ibicancel
    chshorten
    chshalpha
    debugplot
    debugber
    debugpow
    llrmethod
  end

  properties (Constant)

    % LDPC parity matrices
    L_CW = 672 % codeword size tab. 21-19 p. 473
    Z = 42 % submatrix size
    Fs = 2640e6 % Hz, sample rate
    Fc = 2640e6*2/3 % sp/s, sample rate of SC
    % Rate-1/2 LDPC code matrix H = 336 rows x 672 columns, Z = 42, tab. 21.6 p. 452
    P_12 = [
      40 -1 38 -1 13 -1  5 -1 18 -1 -1 -1 -1 -1 -1 -1
      34 -1 35 -1 27 -1 -1 30  2  1 -1 -1 -1 -1 -1 -1
      -1 36 -1 31 -1  7 -1 34 -1 10 41 -1 -1 -1 -1 -1
      -1 27 -1 18 -1 12 20 -1 -1 -1 15  6 -1 -1 -1 -1
      35 -1 41 -1 40 -1 39 -1 28 -1 -1  3 28 -1 -1 -1
      29 -1  0 -1 -1 22 -1  4 -1 28 -1 27 -1 23 -1 -1
      -1 31 -1 23 -1 21 -1 20 -1 -1 12 -1 -1  0 13 -1
      -1 22 -1 34 31 -1 14 -1  4 -1 -1 -1 13 -1 22 24
      ]
    P_58 = [
      20 36 34 31 20  7 41 34 -1 10 41 -1 -1 -1 -1 -1
      30 27 -1 18 -1 12 20 14  2 25 15  6 -1 -1 -1 -1
      35 -1 41 -1 40 -1 39 -1 28 -1 -1  3 28 -1 -1 -1
      29 -1  0 -1 -1 22 -1  4 -1 28 -1 27 24 23 -1 -1
      -1 31 -1 23 -1 21 -1 20 -1  9 12 -1 -1  0 13 -1
      -1 22 -1 34 31 -1 14 -1  4 -1 -1 -1 -1 -1 22 24
      ]
    P_34 = [
      35 19 41 22 40 41 39  6 28 18 17  3 28 -1 -1 -1
      29 30  0  8 33 22 17  4 27 28 20 27 24 23 -1 -1
      37 31 18 23 11 21  6 20 32  9 12 29 -1  0 13 -1
      25 22  4 34 31  3 14 15  4 -1 14 18 13 13 22 24
      ]
    P_1316 = [
      29 30  0  8 33 22 17  4 27 28 20 27 24 23 -1 -1
      37 31 18 23 11 21  6 20 32  9 12 29 10  0 13 -1
      25 22  4 34 31  3 14 15  4  2 14 18 13 13 22 24
      ]
    N_fft = 512
    L_b = 128 % length of a single preamble Golay sequence
    xheaderlen = 1024
  end

  methods (Access = public)

    % Constructs an IEEE 802.11ad modem
    function obj = ieee80211ad(MCS, ftypetx, ftyperx, M_o, usepreamble, useheader, usedata, debugplot, debugpow, debugber, ...
        chestimate, h_ch, equalize, threshold, ibicancel, chshorten, chalpha, llrmethod, snrdb)

      obj.Tc = 1/obj.Fc; % s, SC chip time
      % Golay sequences 21.11 (p. 490)
      obj.g_a = obj.A_k(7, 128 - (1:128)', [1 8 2 4 16 32 64], [-1 -1 -1 -1 +1 -1 -1]); % Ga_128 is A_7(128-n), tab. 21-24
      obj.g_b = obj.B_k(7, 128 - (1:128)', [1 8 2 4 16 32 64], [-1 -1 -1 -1 +1 -1 -1]); % Gb_128 is B_7(128-n), tab. 21-25
      obj.gamma_a = obj.A_k(6, 64 - (1:64)', [2 1 4 8 16 32], [+1 +1 -1 -1 +1 -1]); % Ga_64 is A_6(64-n), tab. 21-26
      obj.pnIC = [1 0 0 1 0 0 0]; % PN initial conditions x1 ... x7
      obj.kappa = comm.PNSequence('Polynomial', [1 0 0 0 1 0 0 1], ...
        'VariableSizeOutput', true, 'MaximumOutputSize', [50000 1], ...
        'Mask', 7, 'InitialConditions', obj.pnIC);
      obj.kappa.release();
      obj.kappa.InitialConditions = [1 1 1 1 1 1 1];
      obj.kappa.reset();
      obj.pnseqrho2 = int8(obj.kappa(20000));
      obj.kappa.release();
      obj.kappa.InitialConditions = obj.pnIC;
      obj.kappa.reset();

      [obj.nu, obj.N_CBPB, obj.rho, obj.r_C] = obj.MCS_order(MCS);% N_CBPS, N_CBPB, rho, R
      obj.L_D = obj.N_CBPB/obj.nu;
      obj.L_G = length(obj.gamma_a);
      obj.L_B = obj.L_D + obj.L_G;

      % prepare LDPC encoder
      if obj.r_C == 1/2
        obj.ldpcencoderD = ldpcEncoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_12));
      elseif obj.r_C == 5/8
        obj.ldpcencoderD = ldpcEncoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_58));
      elseif obj.r_C == 3/4
        obj.ldpcencoderD = ldpcEncoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_34));
      elseif obj.r_C == 13/16
        obj.ldpcencoderD = ldpcEncoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_1316));
      else
        error('TODO')
      end

      obj.ldpcencoderH = ldpcEncoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_34));

      % prepare LDPC decoder
      if obj.r_C == 1/2
        obj.ldpcdecoder = ldpcDecoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_12));
      elseif obj.r_C == 5/8
        obj.ldpcdecoder = ldpcDecoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_58));
      elseif obj.r_C == 3/4
        obj.ldpcdecoder = ldpcDecoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_34));
      elseif obj.r_C == 13/16
        obj.ldpcdecoder = ldpcDecoderConfig(ldpcQuasiCyclicMatrix(obj.Z, obj.P_1316));
      else
        error('TODO')
      end

      if M_o == 1
        % ideal filters without oversampling
        obj.h_tx = 1;
        obj.h_rx = 1;
        span = 0;
      else
        % depends on the type of TX filter
        if strcmp(ftypetx, 'RRC')
          % root raised cosine
          beta = 0.25; % Filter rolloff
          span = 16; % filter span in symbols (at least 16)
          obj.h_tx = rcosdesign(beta, span, M_o, 'sqrt').';
        elseif strcmp(ftypetx, 'GAU')
          % Gaussian
          bt = 0.3; % 3-dB bandwidth-symbol time product
          span = 3;
          obj.h_tx = gaussdesign(bt, span, M_o).';
        else
          error('Invalid TX filter type')
        end
        % depends on the type of RX filter
        if strcmp(ftyperx, 'RRC')
          % root raised cosine
          beta = 0.25; % Filter rolloff
          span = 16; % filter span in symbols (at least 16)
          obj.h_rx = rcosdesign(beta, span, M_o, 'sqrt').';
        elseif strcmp(ftyperx, 'GAU')      % shorten
          % Gaussian
          bt = 0.3; % 3-dB bandwidth-symbol time product
          span = 3;
          obj.h_rx = gaussdesign(bt, span, M_o).';
        else
          error('Invalid RX filter type')
        end
      end

      % frequency response of the filters
      [obj.H_tx, ~] = freqz(obj.h_tx, 1, 512, 'whole', M_o*obj.Fc);

      % calculate global filter
      h_trx = conv(obj.h_tx, obj.h_rx);
      [H_trx, f] = freqz(h_trx, 1, 512, 'whole', M_o*obj.Fc);

      if debugplot
        figure(1)
        clf
        subplot 211
        stem(0:length(obj.h_tx)-1, obj.h_tx, 'DisplayName', 'tx')
        hold on
        stem(0:length(h_trx)-1, h_trx, 'DisplayName', 'tx*rx')
        hold off
        grid on
        xlim([0 length(h_trx)])
        legend('show')
        xlabel('n')
        ylabel('Amplitude')
        title(ftypetx)
        subplot 212
        plot(f/1e9, mag2db(abs(obj.H_tx)), 'DisplayName', 'tx')
        hold on
        plot(f/1e9, mag2db(abs(H_trx)), 'DisplayName', 'tx*rx')
        legend('show')
        grid on
        xlim([0, M_o*obj.Fc/2/1e9])
        xlabel('Frequency (GHz)')
        ylabel('Magnitude (dB)')
        title(ftypetx)
        figure(61)
        clf
        plot(f/1e9, (abs(obj.H_tx)), 'DisplayName', 'tx')
        hold on
        plot(f/1e9, (abs(H_trx)), 'DisplayName', 'tx*rx')
        legend('show')
        grid on
      end

      % Common preamble 21.3.6 (p. 448)

      % Short Training field (STF) 21.3.6.2: composed of 16 repetitions of sequences Ga_128(n)
      % of length 128 defined in 21.11, followed by a single repetition of –Ga_128(n)
      obj.r_STF = [repmat(obj.g_a, 16, 1); -obj.g_a] .* exp(1i * pi * (0:2175)' / 2); % 2176 chips

      % Channel Estimation field (CEF) 21.3.6.3
      Gu_512 = [-obj.g_b; -obj.g_a; +obj.g_b; -obj.g_a];
      Gv_512 = [-obj.g_b; +obj.g_a; -obj.g_b; -obj.g_a];
      obj.r_CEF = [Gu_512; Gv_512; -obj.g_b] .* exp(1i * pi * (0:1151)' / 2); % SC 1152 chips
      obj.s_preamble = [obj.r_STF; obj.r_CEF];

      if chestimate == 2
        % long DFT
        LDFT = repmat(fft(r_CEF(1:1024)), L_o, 1);
      elseif chestimate == 3
        % short DFT
        SDFT = repmat(fft(mean(reshape(r_CEF(1:1024), 512, 2), 2)), L_o, 1);
      elseif chestimate == 4
        % long correlation
        r_CE_L = upsample(r_CEF, L_o);
      elseif chestimate == 5
        % pseudoinverse
        r_CE_L_mat = toeplitz(upsample(r_CEF, L_o), [r_CEF(1) zeros(1, L_o*128 - 1)]);
        r_CEF_L_mat_pinv=(r_CE_L_mat'*r_CE_L_mat)\r_CE_L_mat';
      elseif chestimate == 6
        % short correlation
        obj.Ga_128_L = upsample(obj.g_a .* exp(1i * pi * (0:127)' / 2), 1);
        obj.Gb_128_L = upsample(obj.g_b .* exp(1i * pi * (0:127)' / 2), 1);
      elseif chestimate == 7
        % SMS
        obj.Ga_128_L = upsample(obj.g_a .* exp(1i * pi * (0:127)' / 2), 1);
        obj.Gb_128_L = upsample(obj.g_b .* exp(1i * pi * (0:127)' / 2), 1);
        N_S = 17;
        s = [
          repmat(obj.Ga_128_L, N_S - 1, 1)
          -obj.Ga_128_L
          ];
        c = [
          -obj.Gb_128_L
          -obj.Ga_128_L
          +obj.Gb_128_L
          -obj.Ga_128_L
          -obj.Gb_128_L
          +obj.Ga_128_L
          -obj.Gb_128_L
          -obj.Ga_128_L
          -obj.Gb_128_L
          ];
        obj.Chat = toeplitz(c, [c(1) s(end:-1:end - (128 - 1) + 1).']);
      end


      obj.h = upfirdn(conv(obj.h_tx, h_ch), obj.h_rx, 1, M_o);
      obj.L_h = length(obj.h);
      obj.hf = fft([obj.h; 0], obj.N_fft);

      % error check
      if chestimate >= 2 && ~usepreamble
        error('Preamble must be active to estimate the channel')
      end
      if equalize && ~chestimate
        error('Choose an estimation method to perform equalization')
      end

      obj.MCS = MCS;
      obj.ftypetx = ftypetx;
      obj.ftyperx = ftyperx;
      obj.M_o = M_o;
      obj.usepreamble = usepreamble;
      obj.useheader = useheader;
      obj.usedata = usedata;
      obj.debugplot = debugplot;
      obj.debugber = debugber;
      obj.debugpow = debugpow;
      obj.chestimate = chestimate;
      obj.h_ch = h_ch;
      obj.equalize = equalize;
      obj.threshold = threshold;
      obj.ibicancel = ibicancel;
      obj.chshorten = chshorten;
      obj.chshalpha = chalpha;
      obj.llrmethod = llrmethod;

      % shortening
      obj.L_v = 260;
      if obj.chshorten
        % calculate shortener
        obj.v = obj.shortening_the_channel(obj.h, 56, obj.L_v, sqrt(sum(abs(obj.h_ch).^2))*db2pow(-snrdb));
      else
        obj.v = 1;
      end
      obj.L_v = length(obj.v);

      % shorten
      obj.vf = fft([obj.v; 0], obj.L_B);
      obj.h_sh = conv(obj.h, obj.v);
      obj.hf_sh = fft([obj.h_sh; 0], obj.L_B);

      % calc noise power
      s2w = obj.calcpowers(db2pow(snrdb));

      % equalization
      if obj.equalize == 1
        % ZF-FDE
        obj.q = 1 ./ obj.hf_sh;
        qbadidx = abs(obj.q) > obj.threshold;
        obj.q(qbadidx) = obj.threshold*angle(obj.q(qbadidx));
      elseif obj.equalize == 2
        % MMSE-FDE
        mmse_scale = mean(abs(obj.hf_sh).^2 ./ (abs(obj.hf_sh).^2 + s2w*abs(obj.vf).^2));
        obj.q = (conj(obj.hf_sh) ./ (abs(obj.hf_sh).^2 + s2w*abs(obj.vf).^2)) / mmse_scale;
      else
        % recover delay, amplitude and phase
        [amp, dly] = max(abs(obj.h_sh));
        dly = dly - 1;
        phs = angle(mean(obj.h_sh));
        obj.q = exp(2i*pi*(0:obj.L_B - 1).'*dly/obj.L_B)*exp(-1i*phs)/amp;
      end

    end

    % Generate transmitted samples
    function [x_up, Lambda, r_b_th, r_b_act] = transmit(obj, b)

      % data length
      L_PS = length(b)/8;

      % Data field payload of the PSDU 21.6.3.2 p. 473
      N_W = ceil((8*L_PS)/((obj.L_CW/obj.rho) * obj.r_C)); % TODO: check if BRP packet
      L_DP = N_W*obj.L_CW * obj.r_C/obj.rho - 8*L_PS; % N_DATA_PAD
      N_B = ceil(N_W*obj.L_CW / obj.N_CBPB); % N_BLKS
      L_BP = N_B*obj.N_CBPB - N_W*obj.L_CW; % N_BLK_PAD
      obj.L_MW = obj.L_CW * obj.r_C / obj.rho; % L_CWD
      r_b_th = (obj.L_D/obj.L_B)*obj.nu*obj.Fc*obj.r_C/obj.rho;
      r = exp(1i*(pi/2)*(0:N_B*obj.L_B + obj.L_G - 1)');

      % data padding and scrambling
      obj.kappa.release();
      obj.kappa.InitialConditions = obj.pnIC;
      obj.kappa.reset();
      B_MW = [reshape(bitxor([b; zeros(L_DP, 1, 'int8')], int8(obj.kappa(8*L_PS + L_DP))), obj.L_MW, N_W)
        zeros((obj.rho - 1)*obj.L_MW, N_W)];

      % channel encoding 21.6.3.2.3.3 p. 473
      B_CW = zeros(obj.L_CW, N_W, 'int8');

      % encode
      for n = 1:N_W
        B_CW(:, n) = ldpcEncode(B_MW(:, n), obj.ldpcencoderD);
        if obj.rho == 2
          % repetition code
          B_CW(obj.L_MW + (1:obj.L_MW), n) = bitxor(B_CW(1:obj.L_MW, n), obj.pnseqrho2(1:obj.L_MW)); % repeat
        end
      end

      % padding
      b_CW = [B_CW(:); int8(obj.kappa(L_BP))];

      % mapping
      Lambda = reshape(b_CW, obj.nu, N_B*obj.L_D);
      if obj.nu == 1
        % BPSK
        mu = (2*double(Lambda(1, :)) - 1).';
      elseif obj.nu == 2
        % QPSK
        mu = (exp(-1i*pi/4)/sqrt(2)) * ((2*double(Lambda(1, :)) - 1) + 1i*(2*double(Lambda(2, :)) - 1)).';
      elseif obj.nu == 4
        % 16QAM
        mu = (1/sqrt(10)) * ( (4*double(Lambda(1, :)) - 2) - (2*double(Lambda(1, :)) - 1) .* (2*double(Lambda(2, :)) - 1) ).' ...
           + (1i/sqrt(10)) * ( (4*double(Lambda(3, :)) - 2) - (2*double(Lambda(3, :)) - 1) .* (2*double(Lambda(4, :)) - 1) ).';
      else
        error('TODO')
      end

      % guard intervals
      M = reshape(mu, obj.L_D, N_B);
      D = [M; repmat(obj.gamma_a, 1, N_B)];

      % modulate
      x_data = r(1:obj.L_B*N_B + obj.L_G).*[obj.gamma_a; D(:)];

      if obj.debugplot
        figure(3)
        clf
        plot(x_data(1:20:end), 'o', 'DisplayName', 'tx')
        xlabel('I')
        ylabel('Q')
        grid
        axis equal
        legend('show')
      end

      % Generation of samples
      x_preamble = obj.generate_preamble();
      x_header = obj.generate_header();

      % put together all the pieces
      x = [x_preamble; x_header; x_data];

      % calculate actual bitrate
      r_b_act = 8*L_PS/(obj.Tc*length(x));

      % interpolate and shape the signal
      x_up = upfirdn(x, obj.h_tx, obj.M_o, 1);

    end

    % Receive samples
    function [bhat, Lambdahat, avg_niters] = receive(obj, y_up, L_PS, sigma2_w, p_S, p_I, p_N, snrdb, p_S1, p_S2, p_I1, p_I2, p_N1, p_N2)

      % filter the signal and decimate
      y = upfirdn(y_up, obj.h_rx, 1, obj.M_o);

      % Data field payload of the PSDU 21.6.3.2 p. 473
      N_W = ceil((8*L_PS)/((obj.L_CW/obj.rho) * obj.r_C)); 
      N_B = ceil(N_W*obj.L_CW / obj.N_CBPB); % N_BLKS
      obj.L_MW = obj.L_CW * obj.r_C / obj.rho; % L_CWD
      r = exp(1i*(pi/2)*(0:N_B*obj.L_B + obj.L_G - 1)');

      % Preamble part
      if obj.usepreamble

        % extract preamble parts
        rhat_CEF = y(length(obj.r_STF) + (1:length(obj.r_CEF)));

        % types of estimation
        if obj.chestimate == 2

          % long DFT estimation
          LDFT_r = fft(rhat_CEF(1:L_o*1024));
          hfhat = LDFT_r ./ LDFT;
          hhat = ifft(hfhat);

        elseif obj.chestimate == 3

          % short DFT estimation
          SDFT_r = fft(mean(reshape(rhat_CEF(1:L_o*1024), L_o*512, 2), 2));
          hfhat = SDFT_r ./ SDFT;
          hhat = ifft(hfhat);

        elseif obj.chestimate == 4

          % long correlation estimation
          [corr, corr_lags] = xcorr(rhat_CEF, r_CE_L);
          hhat = L_o*corr(corr_lags < L_o*512 & corr_lags >= 0)/length(r_CE_L);

        elseif obj.chestimate == 5

          % pseudoinverse estimation
          hhat = r_CEF_L_mat_pinv * rhat_CEF;

        elseif obj.chestimate == 6

          % short correlation estimation
          corr = xcorr(rhat_CEF, obj.Ga_128_L);
          midx = (length(corr) - 1)/2 + 1;
          alpha = reshape(corr(midx:midx + 128*8 - 1), 128, 8);
          alpha = (-alpha(:, 2) - alpha(:, 4) + alpha(:, 6) - alpha(:, 8))/(4*128);
          corr = xcorr(rhat_CEF, obj.Gb_128_L);
          beta = reshape(corr(midx:midx + 128*8 - 1), 128, 8);
          beta = (-beta(:, 1) + beta(:, 3) - beta(:, 5) - beta(:, 7))/(4*128);
          hhat = (alpha(1:128) + beta(1:128))/2;

        elseif obj.chestimate == 7

          % short correlation estimation as starting point
          corr = xcorr(rhat_CEF, obj.Ga_128_L);
          midx = (length(corr) - 1)/2 + 1;
          alpha = reshape(corr(midx:midx + 128*8 - 1), 128, 8);
          alpha = (-alpha(:, 2) - alpha(:, 4) + alpha(:, 6) - alpha(:, 8))/(4*128);
          corr = xcorr(rhat_CEF, obj.Gb_128_L);
          beta = reshape(corr(midx:midx + 128*8 - 1), 128, 8);
          beta = (-beta(:, 1) + beta(:, 3) - beta(:, 5) - beta(:, 7))/(4*128);
          hhat = (alpha(1:128) + beta(1:128))/2;
          hhat1 = hhat;

          % refinement of hhat
          hhat = [hhat; min(abs(hhat(hhat ~= 0)))*randn(max(0, 128 - length(hhat)), 1)];
          L_s_sms = 128;
          [~, highestIdxs] = maxk(abs(hhat), L_s_sms); % non genie-aided
          E = zeros(L_s_sms, 128);
          E(sub2ind([L_s_sms 128], 1:length(highestIdxs), highestIdxs.')) = 1; % use sparsity
          Uhat_in = diag(conj(hhat).*hhat);
          Gamma = E*Uhat_in*obj.Chat';
          % set best k_u for sparse channel
          if snrdb < 0
            k_u = 5; 
          elseif snrdb > 40
            k_u = 40; 
          else
            k_u = (7/240) *(snrdb).^2 + (-7/24)*(snrdb) + 5; 
          end
          % % set best k_u for Brno channel
          % if snrdb < 0
          %   k_u = 7/2; 
          % elseif snrdb > 36
          %   k_u = 0; 
          % else
          %   k_u = (19/9360) *(snrdb).^2 + (-797/4680)*(snrdb) + 7/2; 
          % end
          hhat = (1/(k_u*sigma2_w))*E.'*(eye(L_s_sms) - Gamma*obj.Chat*E.'/(Gamma*obj.Chat*E.' + k_u*sigma2_w*eye(L_s_sms)))*Gamma*rhat_CEF;

        end

      end

      % select ideal channel estimate
      if obj.chestimate == 1
        hhat = obj.h;
      end
      hfhat = fft([hhat; 0], obj.N_fft);

      % shortening
      if obj.chshorten
        % calculate shortener
        if isempty(obj.v)
          obj.v = obj.shortening_the_channel(hhat, 56, 188, sqrt(sum(abs(obj.h_ch).^2))*db2pow(-snrdb));
        end
      else
        obj.v = 1;
      end

      % shorten
      obj.vf = fft([obj.v; 0], obj.L_B);
      hhat_sh = conv(hhat, obj.v);
      hfhat_sh = fft([hhat_sh; 0], obj.L_B);
      y_sh = conv(y, obj.v);

      if obj.debugplot
        figure(91)
        clf
        subplot 211
        plot(0:length(hfhat) - 1, abs(hfhat), 'DisplayName', ['est' int2str(obj.chestimate)])
        hold on
        plot(0:length(obj.hf) - 1, abs(obj.hf), '--', 'DisplayName', 'th')
        if obj.chshorten == 1
          plot(0:length(hfhat_sh) - 1, abs(hfhat_sh), '-.', 'DisplayName', 'estsh')
        end
        xlim([0 length(obj.hf) - 1])
        grid on
        legend('show')
        ylabel('Magnitude')
        xlabel('Frequency bin \itk')
        title('CTF')
        subplot 212
        stem(0:length(hhat) - 1, abs(hhat), 'x', 'DisplayName', ['est' int2str(obj.chestimate)])
        hold on
        stem(0:length(obj.h) - 1, abs(obj.h), 'DisplayName', 'th')
        if obj.chshorten == 1
          stem(0:length(hhat_sh) - 1, abs(hhat_sh), '^', 'DisplayName', 'estsh')
        end
        grid on
        xlim tight
        legend('show')
        ylabel('Magnitude')
        xlabel('Time lag \itn')
        title('CIR')
      end

      % Data part
      if obj.usedata

        % select data only without first GI
        Dhat = reshape(y_sh(obj.usepreamble*length(obj.s_preamble) + obj.useheader*obj.xheaderlen + obj.L_G + (1:N_B*obj.L_B)), ...
          obj.L_B, N_B).';

        % initialize
        Dhat_eq = zeros(N_B, obj.L_B);
        Mhat = zeros(N_B, obj.L_D);
        Lhat_h = length(hhat);
        if obj.ibicancel > 0 && obj.usepreamble
          iota = conv([obj.g_a(64 + 1:end).' -obj.g_b.' -obj.g_a.' -obj.g_b.' zeros(1, obj.L_G)] .* r(1:obj.L_B).', hhat);
          iota = [iota(obj.L_B + (1:Lhat_h - 1)) zeros(1, obj.L_B - (Lhat_h - 1))];
        else
          iota = zeros(1, obj.L_B);
        end

        % do all the blocks
        for ii = 1:N_B

          % successive iterations
          for ibicycle = 1:obj.ibicancel + 1

            % IBI cancellation, just once
            if ibicycle == 1
              dhat_noibi = Dhat(ii, :) - iota;
            end

            % equalization
            if obj.equalize == 1

              % ZF-FDE
              obj.q = 1 ./ hfhat_sh;
              qbadidx = abs(obj.q)>obj.threshold;
              obj.q(qbadidx) = obj.threshold*angle(obj.q(qbadidx));
              Dhat_eq(ii, :) = ifft(obj.q.'.*fft(dhat_noibi));

            elseif obj.equalize == 2

              % diagonal MMSE-FDE
              mmse_scale = mean(abs(hfhat_sh).^2 ./ (abs(hfhat_sh).^2 + sigma2_w*abs(obj.vf).^2)); 
              obj.q = (conj(hfhat_sh) ./ (abs(hfhat_sh).^2 + sigma2_w*abs(obj.vf).^2)) / mmse_scale;
              Dhat_eq(ii, :) = ifft(obj.q.'.*fft(dhat_noibi));

            else

              % recover delay, amplitude and phase
              [amp, dly] = max(abs(hhat_sh));
              dly = dly - 1;
              phs = angle(mean(hhat_sh));
              obj.q = exp(2i*pi*(0:obj.L_B - 1).'*dly/obj.L_B)*exp(-1i*phs)/amp;
              Dhat_eq(ii, :) = ifft(obj.q.'.*fft(dhat_noibi));

            end

            % de-pi/2, guard removal
            Mhat(ii, :) = r(1:obj.L_D)'.*Dhat_eq(ii, 1:obj.L_D);

            if obj.ibicancel > 0

              % regenerate
              if obj.nu == 1
                % BPSK
                tmp = sign(real(Mhat(ii, :)));
              elseif obj.nu == 2
                % QPSK
                tmp = Mhat(ii, :)*exp(1i*pi/4)*sqrt(2);
                tmp = (sign(real(tmp)) + 1i*sign(imag(tmp)))*exp(-1i*pi/4)/sqrt(2);
              else
                % 16QAM
                tmp = Mhat(ii, :)*sqrt(10);
                tmp = (...
                  (4*double(real(tmp) >= 0) - 2) - (2*double(real(tmp) >= 0) - 1).*(2*double(abs(real(tmp)) < 2) - 1) ...
                  + 1i*((4*double(imag(tmp) >= 0) - 2) - (2*double(imag(tmp) >= 0) - 1).*(2*double(abs(imag(tmp)) < 2) - 1)) ...
                  )/sqrt(10);
              end

              % remove old ISI
              if ibicycle > 1
                dhat_noibi = dhat_noibi - iota;
              end

              % calculate new ISI
              iota = conv([tmp zeros(1, obj.L_G)].*r(1:obj.L_B).', hhat);
              iota = [iota(obj.L_B + (1:min(obj.L_b, obj.L_h) - 1)) zeros(1, obj.L_B - (min(obj.L_b, obj.L_h) - 1))];

              % add new ISI
              dhat_noibi = dhat_noibi + iota;

            end

          end

        end

        if obj.debugpow
          figure(301)
          hold on
          plot(0:obj.N_fft - 1, pow2db(pwelch(Mhat(:), obj.N_fft, 1)*2*pi), 'DisplayName', 'M\^')
          hold off
          legend('show')
          disp(['Mhat(≈S+I+N)=' num2str(var(Mhat(:), 1))])
        end

        if obj.debugplot
          figure(3)
          hold on
          plot(Dhat(1:20:end), 'k.', 'DisplayName', 'rx')
          plot(Mhat(1:20:end), 'g*', 'DisplayName', 'eq')
          hold off
          axis square
          legend off
          title([])
          set(gca, 'FontSize', 13)
          xlim([-2 2])
          ylim([-2 2])
        end

        % channel decoding
        Lambdahat = zeros(obj.nu, N_B*obj.L_D);

        Mhat1 = Mhat.';
        if obj.llrmethod == 1
          sigma2 = (p_N + p_I)./(2*p_S); % PROPOSED - variance of the real or imaginary part
        else
          sigma2 = (p_N)./(2*(p_S + p_I)); % CONVENTIONAL - alternative imprecise method for LLR
        end
        if obj.rho == 1

          % rho=1, no repetition
          if obj.nu == 1
            % BPSK
            if obj.llrmethod == 1
              tmp = -2*real(Mhat1./sqrt(p_S))./sigma2; % OPTIMUM
            else
              tmp = -2*real(Mhat1./sqrt(p_S+p_I))./sigma2; % INFIMUM
            end
            Lambdahat(1, :) = tmp(:);
          elseif obj.nu == 2
            % QPSK
            if obj.llrmethod == 1
              tmp = -2*real(exp(1i*pi/4)*Mhat1./sqrt(p_S))./sigma2; % OPTIMUM
              Lambdahat(1, :) = tmp(:);
              tmp = -2*imag(exp(1i*pi/4)*Mhat1./sqrt(p_S))./sigma2; % OPTIMUM
              Lambdahat(2, :) = tmp(:);
            else
              tmp = -2*real(exp(1i*pi/4)*Mhat1./sqrt(p_S+p_I))./sigma2; % INFIMUM
              Lambdahat(1, :) = tmp(:);
              tmp = -2*imag(exp(1i*pi/4)*Mhat1./sqrt(p_S+p_I))./sigma2; % INFIMUM
              Lambdahat(2, :) = tmp(:);
            end
          elseif obj.nu == 4
            % 16QAM
            if obj.llrmethod == 1
              tmp = sqrt(10./p_S).*Mhat1; % OPTIMUM
            else
              tmp = sqrt(10./(p_S+p_I)).*Mhat1; % INFIMUM
            end
            tmp2 = repmat(sigma2, 1, N_B);
            Lambdahat(1, :) = log(exp(-(real(tmp(:)) + 3).^2./(20*tmp2(:))) + exp(-(real(tmp(:)) + 1).^2./(20*tmp2(:)))) ...
              - log(exp(-(real(tmp(:)) - 1).^2./(20*tmp2(:))) + exp(-(real(tmp(:)) - 3).^2./(20*tmp2(:))));
            Lambdahat(2, :) = log(exp(-(real(tmp(:)) + 3).^2./(20*tmp2(:))) + exp(-(real(tmp(:)) - 3).^2./(20*tmp2(:)))) ...
              - log(exp(-(real(tmp(:)) - 1).^2./(20*tmp2(:))) + exp(-(real(tmp(:)) + 1).^2./(20*tmp2(:))));
            Lambdahat(3, :) = log(exp(-(imag(tmp(:)) + 3).^2./(20*tmp2(:))) + exp(-(imag(tmp(:)) + 1).^2./(20*tmp2(:)))) ...
              - log(exp(-(imag(tmp(:)) - 1).^2./(20*tmp2(:))) + exp(-(imag(tmp(:)) - 3).^2./(20*tmp2(:))));
            Lambdahat(4, :) = log(exp(-(imag(tmp(:)) + 3).^2./(20*tmp2(:))) + exp(-(imag(tmp(:)) - 3).^2./(20*tmp2(:)))) ...
              - log(exp(-(imag(tmp(:)) - 1).^2./(20*tmp2(:))) + exp(-(imag(tmp(:)) + 1).^2./(20*tmp2(:))));
          end

          bhat_CW = Lambdahat(:);
          Bhat_CW = reshape(bhat_CW(1:obj.L_CW*N_W), obj.L_CW, N_W);

        elseif obj.rho == 2

          % rho == 2, repetition (only for BPSK 1/2 -> MCS 1)
          tmp_muhat = Mhat1(:);
          tmp_muhat_b = reshape(tmp_muhat(1:obj.L_CW*N_W), obj.L_CW, N_W);
          tmp_sigma2 = repmat(sigma2, N_B, 1);
          tmp_sigma2_b = reshape(tmp_sigma2(1:obj.L_CW*N_W), obj.L_CW, N_W);
          tmp_Ps = repmat(p_S, N_B, 1);
          tmp_Ps3b = reshape(tmp_Ps(1:obj.L_CW*N_W), obj.L_CW, N_W);

        
          if obj.llrmethod == 1
            % PROPOSED
            sigma2_1 = (p_N1 + p_I1)./(2*p_S1); 
            sigma2_2 = (p_N2 + p_I2)./(2*p_S2); 
          else
            % CONVENTIONAL
            sigma2_1 = (p_N1)./(2*(p_S1 + p_I1)); 
            sigma2_2 = (p_N2)./(2*(p_S2 + p_I2)); 
          end

          Bhat_CW(1:obj.L_MW, :) = -2*real(((tmp_muhat_b(1:obj.L_MW, :) + (-1).^double(obj.pnseqrho2(1:obj.L_MW)).*tmp_muhat_b(obj.L_MW + (1:obj.L_MW), :))/2)./sqrt(repmat([p_S1 p_S2], 1, N_W/2)))./repmat([sigma2_1 sigma2_2], 1, N_W/2);
          Bhat_CW(obj.L_MW + (1:obj.L_MW), :) = +Inf;
          Bhat_CW((obj.rho*obj.L_MW + 1):obj.L_CW, :) = -2*real(tmp_muhat_b((obj.rho*obj.L_MW + 1):obj.L_CW, :))./(sqrt(tmp_Ps3b((obj.rho*obj.L_MW + 1):obj.L_CW, :)).*tmp_sigma2_b((obj.rho*obj.L_MW + 1):obj.L_CW, :));
          Lambdahat(:) = -sign(tmp_muhat);

        end

        % At high SNR sometimes NaN are generated in the LLR calculation
        Bhat_CW(isnan(Bhat_CW)) = 0;

        % actual decode
        tot_niters = 0;
        Bhat_MW = zeros(obj.rho*obj.L_MW, N_W, 'int8');
        for ii = 1:N_W
          [Bhat_MW(:, ii), niters] = ldpcDecode(Bhat_CW(:, ii), obj.ldpcdecoder, 50);
          tot_niters = tot_niters + niters;
        end
        avg_niters = tot_niters/N_W;

        % remove zeros
        Bhat_MW(obj.L_MW + (1:(obj.rho - 1)*obj.L_MW), :) = [];

        % original bits are retrieved
        obj.kappa.release();
        obj.kappa.InitialConditions = obj.pnIC;
        obj.kappa.reset();
        bhat = bitxor(Bhat_MW((1:8*L_PS).'), int8(obj.kappa(8*L_PS)));

      end

    end

    % Calculate all the required powers
    function [s2w, p_S, p_I, p_N, p_X, p_E, p_S1, p_S2, p_I1, p_I2, p_N1, p_N2] = calcpowers(obj, snr_bsh)

      % filtered signal power at TX
      p_X = mean(abs(obj.H_tx).^2*1);

      % signal gone through the channel
      H_ch = fft(obj.h_ch,length(obj.H_tx));
      p_E = mean(abs(obj.H_tx .* H_ch).^2*1);

      % transmitted signal power
      s2d = 1;

      % signal power at the shortener input
      Sshin = s2d * (obj.h' * obj.h);

      % The noise power at the shortener input is
      %      Nshin = s2w
      % Thus the calculated SNR before the shortener is
      %      snr_bsh = Sshin / Nshin
      % Substituting, we solve for the needed noise power
      %      snr_bsh = Sshin / s2w
      s2w = Sshin / snr_bsh;

      if nargout > 1
        F=fft(eye(obj.L_B))/sqrt(obj.L_B);
        R=diag(exp(1i*(pi/2)*(0:obj.L_B-1).'));
        P=[eye(obj.L_D) zeros(obj.L_D,obj.L_G)];
        % for the signal
        H = toeplitz([obj.h;zeros(obj.L_G+3*obj.L_B-1,1)],[obj.h(1) zeros(1,obj.L_G+3*obj.L_B-1)]);
        V = toeplitz([obj.v;zeros(obj.L_h+obj.L_G+3*obj.L_B-1-1,1)],[obj.v(1) zeros(1,obj.L_h+obj.L_G+3*obj.L_B-1-1)]);
        K_i = [zeros(obj.L_B,obj.L_G) zeros(obj.L_B) eye(obj.L_B) zeros(obj.L_B,obj.L_h+(3-2)*obj.L_B-1+obj.L_v-1)];
        N_i = P*R'*F'*diag(obj.q)*F*K_i*V;
        B_i = N_i*H;
        B_iim1 = B_i(:, obj.L_G + 0*obj.L_B + (1:obj.L_D));
        B_ii = B_i(:, obj.L_G + 1*obj.L_B + (1:obj.L_D));
        N_ii = N_i(:, obj.L_G + 1*obj.L_B + (1:obj.L_D));
        TTX = [
          eye(obj.L_CW/4) zeros(obj.L_CW/4)
          diag((-1).^double(obj.pnseqrho2(1:obj.L_CW/4))) zeros(obj.L_CW/4)
          zeros(obj.L_CW/2)
          zeros(obj.L_CW/4) eye(obj.L_CW/4)
          zeros(obj.L_CW/4) diag((-1).^double(obj.pnseqrho2(1:obj.L_CW/4)))
          zeros(obj.L_CW/2)
          ];
        TRX = 0.5*[
          eye(obj.L_CW/4) diag((-1).^double(obj.pnseqrho2(1:obj.L_CW/4))) zeros(obj.L_CW/4,obj.L_CW*3/2)
          zeros(obj.L_CW/4,obj.L_CW) eye(obj.L_CW/4) diag((-1).^double(obj.pnseqrho2(1:obj.L_CW/4))) zeros(obj.L_CW/4,obj.L_CW/2)
          ];
        Bbig_ii = TRX*[
          B_ii zeros(obj.L_D, obj.L_D*2)
          zeros(obj.L_D) B_ii zeros(obj.L_D)
          zeros(obj.L_D, obj.L_D*2) B_ii
          ]*TTX;
        Dbig_ii = Bbig_ii.*eye(size(Bbig_ii));
        pbig_S = s2d*diag(Dbig_ii*Dbig_ii');
        p_S1 = pbig_S(1:obj.L_CW/4);
        p_S2 = pbig_S(obj.L_CW/4+1:end);
        ISIbig = s2d*diag((Bbig_ii - Dbig_ii)*(Bbig_ii - Dbig_ii)');
        ISI1 = ISIbig(1:obj.L_CW/4);
        ISI2 = ISIbig(obj.L_CW/4+1:end);
        Bbig_iim1 = TRX*[
          B_iim1 zeros(obj.L_D, obj.L_D*2)
          zeros(obj.L_D) B_iim1 zeros(obj.L_D)
          zeros(obj.L_D, obj.L_D*2) B_iim1
          ];
        IBIbig = s2d*diag(Bbig_iim1*Bbig_iim1');
        IBI1 = IBIbig(1:obj.L_CW/4);
        IBI2 = IBIbig(obj.L_CW/4+1:end);
        p_I1 = ISI1 + (1 - obj.ibicancel)*IBI1;
        p_I2 = ISI2 + (1 - obj.ibicancel)*IBI2;
        Nbig_ii = TRX*[
          N_ii zeros(obj.L_D, obj.L_D*2)
          zeros(obj.L_D) N_ii zeros(obj.L_D)
          zeros(obj.L_D, obj.L_D*2) N_ii
          ];
        pbig_N = s2w*diag(Nbig_ii*Nbig_ii');
        p_N1 = pbig_N(1:obj.L_CW/4);
        p_N2 = pbig_N(obj.L_CW/4+1:end);
        D_ii = B_ii.*eye(size(B_ii));
        % results
        p_S = s2d*diag(D_ii*D_ii');
        ISI = s2d*diag((B_ii - D_ii)*(B_ii - D_ii)');
        IBI = s2d*diag(B_iim1*B_iim1');
        p_I = ISI + (1 - obj.ibicancel)*IBI;
        p_N = s2w*diag(N_ii*N_ii'); % TODO include also previous block noise

      end

    end

  end
  
  methods (Access = private)

    % Generate the preamble
    function x_preamble = generate_preamble(obj)

      if obj.usepreamble
        % assemble the preamble for SC
        x_preamble = obj.s_preamble;
      else
        x_preamble = zeros(0, 1);
      end
    end

    % Generate the header
    function x_header = generate_header(obj)

      if obj.useheader

        % 64 bits
        header = zeros(64, 1, 'int8');

        % Scrambler Initialization: bits X1–X7 of the initial scrambler state
        header(1:7) = obj.pnIC;

        % MCS: index into the Modulation and Coding Scheme table
        header(8:12) = de2bi(obj.MCS, 5, 'right-msb');

        % Length: number of data octets in the PSDU. Range 1–262143
        header(13:30) = de2bi(L_PS, 18, 'right-msb');

        % Additional PPDU: contains a copy of the parameter ADD-PPDU from the TXVECTOR.
        % A value of 1 indicates that this PPDU is immediately followed by another
        % PPDU with no IFS or preamble on the subsequent PPDU. A value of 0 indicates
        % that no additional PPDU follows this PPDU
        header(31) = 0;

        % Packet Type: See definition of Packet Type field in Table 21-11.
        header(32) = 0;

        % Training Length: Corresponds to the TXVECTOR parameter TRN-LEN.
        % If the Beam Tracking Request field is 0, the Training Length field indicates
        % the length of the training field. The use of this field is defined in 21.10.2.2.3.
        % A value of 0 indicates that no training field is present in this PPDU.
        % If the Beam Tracking Request field is 1 and the Packet Type field is 1,
        % the Training Length field indicates the length of the training field.
        % If the Packet Type field is 0, the Training Length field indicates the
        % length of the training field requested for receive training.
        header(33:37) = de2bi(0, 5, 'right-msb');

        % Aggregation: Set to 1 to indicate that the PPDU in the data portion of
        % the packet contains an A-MPDU; otherwise, set to 0.
        header(38) = 0;

        % Beam Tracking Request: Corresponds to the TXVECTOR parameter BEAM_TRACKING_REQUEST.
        % Set to 1 to indicate the need for beam tracking (9.35.7); otherwise, set to 0.
        % The Beam Tracking Request field is reserved when the Training Length field is 0.
        header(39) = 1;

        % Last RSSI: Contains a copy of the parameter LAST_RSSI from the TXVECTOR.
        % When set to 0, this field is reserved and ignored by the receiver.
        % The value is an unsigned integer:
        % Values of 2 to 14 represent power levels (–71+value×2) dBm.
        % A value of 15 represents a power greater than or equal to -42 dBm.
        % A value of 1 represents a power less than or equal to –68 dBm.
        % Value of 0 indicates that the previous packet was not received a SIFS
        % period before the current transmission.
        header(40:43) = de2bi(0, 4, 'right-msb');

        % Turnaround: As defined in Table 21-1.
        header(44) = 0;

        % Reserved: set to 0, ignored by the receiver
        header(45:48) = 0;

        % HCS: Header check sequence
        header(1:64) = crcGen(logical(header(1:48)));

        % Encode and modulate header 21.6.3.1.4.p. 472
        ob = header;
        LH = length(ob);
        d1s = zeros(size(ob), 'int8');

        % scrambling
        obj.kappa.reset();
        d1s(1:7) = ob(1:7);
        d1s(8:LH) = bitxor(ob(8:end), int8(obj.kappa(LH - 7)));

        % channel encoding
        pc = ldpcEncode([d1s; zeros(504 - LH, 1, 'int8')], obj.ldpcencoderH);
        cs1 = [pc(1:LH); pc(505:664); pc(673:end)];

        % scrambling, again
        obj.kappa.release();
        obj.kappa.InitialConditions = [1 1 1 1 1 1 1];
        obj.kappa.reset();
        cs2 = bitxor([pc(1:LH); pc(505:656); pc(665:end)], int8(obj.kappa(224)));

        % modulate
        N_CBPB_H = 448;
        k = (0:N_CBPB_H + obj.L_G - 1)';
        x_header = [
          [obj.gamma_a; (+2*double([cs1; cs2]) - 1)] .* exp(1i * pi * k / 2);
          [obj.gamma_a; (-2*double([cs1; cs2]) + 1)] .* exp(1i * pi * k / 2)
          ];
      else
        x_header = zeros(0, 1);
      end
    end

    % For Golay sequences
    function akn = A_k(obj, k, n, D_k, W_k)
      if k == 0
        akn = double(n == 0);
      else
        akn = W_k(k) * obj.A_k(k - 1, n, D_k, W_k) + obj.B_k(k - 1, n - D_k(k), D_k, W_k);
      end
    end

    % For Golay sequences
    function bkn = B_k(obj, k, n, D_k, W_k)
      if k == 0
        bkn = double(n == 0);
      else
        bkn = W_k(k) * obj.A_k(k - 1, n, D_k, W_k) - obj.B_k(k - 1, n - D_k(k), D_k, W_k);
      end
    end

    % MCS code rate and modulation, tab 21-18 p. 472, tab 21-20 p. 476
    function [ncbps, ncbpb, rho, r_C, dr] = MCS_order(~, MCS)

      dr = nan;
      switch MCS
        case 1
          ncbps = 1; % BPSK
          ncbpb = 448; % coded bits per block
          rho = 2; % repetition
          r_C = 1/2; % code rate R
        case 2
          ncbps = 1; % BPSK
          ncbpb = 448; % coded bits per block
          rho = 1; % no repetition
          r_C = 1/2; % code rate R
        case 3
          ncbps = 1; % BPSK
          ncbpb = 448; % coded bits per block
          rho = 1; % no repetition
          r_C = 5/8; % code rate R
        case 4
          ncbps = 1; % BPSK
          ncbpb = 448; % coded bits per block
          rho = 1; % no repetition
          r_C = 3/4; % code rate R
        case 5
          ncbps = 1; % BPSK
          ncbpb = 448; % coded bits per block
          rho = 1; % no repetition
          r_C = 13/16; % code rate R
        case 6
          ncbps = 2; % QPSK
          ncbpb = 896; % coded bits per block
          rho = 1; % no repetition
          r_C = 1/2; % code rate R
        case 7
          ncbps = 2; % QPSK
          ncbpb = 896; % coded bits per block
          rho = 1; % no repetition
          r_C = 5/8; % code rate R
        case 8
          ncbps = 2; % QPSK
          ncbpb = 896; % coded bits per block
          rho = 1; % no repetition
          r_C = 3/4; % code rate R
        case 9
          ncbps = 2; % QPSK
          ncbpb = 896; % coded bits per block
          rho = 1; % no repetition
          r_C = 13/16; % code rate R
        case 10
          ncbps = 4; % 16-QAM
          ncbpb = 1792; % coded bits per block
          rho = 1; % no repetition
          r_C = 1/2; % code rate R
        case 11
          ncbps = 4; % 16-QAM
          ncbpb = 1792; % coded bits per block
          rho = 1; % no repetition
          r_C = 5/8; % code rate R
        case 12
          ncbps = 4; % 16-QAM
          ncbpb = 1792; % coded bits per block
          rho = 1; % no repetition
          r_C = 3/4; % code rate R
        otherwise
          error('Invalid MCS index')
      end
    end

    % Channel shortening
    function v = shortening_the_channel(obj, h, nu, L_v, s2w)
      if length(h) <= nu
        v = 1;
        return
      end
      lambdamin = +Inf;
      for Delta = 0:10
        H = toeplitz([h; zeros(L_v - 1, 1)], [h(1) zeros(1, L_v - 1)]);
        H_win = H(Delta + (0:nu - 1) + 1, :);
        H_wall = [H((0:Delta - 1) + 1, :); H(Delta + nu + 1:end, :)];
        A = H_wall'*H_wall;
        B = H_win'*H_win;
        C = s2w*eye(size(A));
        [V, D] = eig(A, B);
        d = diag(D);
        d(isinf(d)) = +Inf;
        [~, lambdaminIdx] = min(real(d));
        vnow = V(:, lambdaminIdx);
        if abs(d(lambdaminIdx)) < abs(lambdamin)
          lambdamin = d(lambdaminIdx);
          v = vnow;
        end
      end
    end

  end


end
