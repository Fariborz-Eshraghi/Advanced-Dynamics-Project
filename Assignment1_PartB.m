% Part B - 1st Assignment - Dynamics
clc;
clear all;
close all;

%% 1) Identify the Modal Parameters (Natural Frequencies, Damping Ratios, and Mode Shapes) -> first two axial modes
% Load the data

df = 0.333; % Resolution (Hz)
%cd 'D:\Origins\Polimi\ADMS (Advance Dynamics of Mechanical Systems)\Assignments\First\Prat B';
d = load ("Data.mat");
freq = d.freq;
frf = d.frf;
cohe = d.cohe;

% Ploting the experimental data
H   = cell(1,12);          % pre-allocate cell array
mag = cell(1,12);
Co = cell(1,12);

figure ('Name','Original Data'); % Opening the figure to draw
for i = 1 : 12
    H {i} = frf (:,i);           % Complex FRF  (m/s^2/N)
    mag{i} = abs(H{i});
    Co {i} = cohe (:,i);
    
    % Ploting the FRF
    subplot (3,1,1), semilogy(freq, abs(H{i}));
    hold on;
    ylabel('|H|  (m s^{-2}/N)'), grid on

    % Ploting the Phase
    subplot (3,1,2), plot(freq, rad2deg(angle(H{i})));
    hold on;
    ylabel('phase  (rad)'), xlabel('Frequency (Hz)')

    % Ploting the Coherrence
    subplot(3,1,3), plot(freq, Co {i});
    hold on;
    xlabel('Frequency (Hz)'), ylabel('Coherence'), ylim([0 1]), grid on
end


%% Calculating the estimated part
% Add points depending on preferences
points = frf;
F_resp = d.freq;
n_modes = 2;
n_samples = length(frf(1,:));
F_resp_size = size (points,1);
FRFs = zeros (n_samples, F_resp_size);
f_ranges = {};

FRFs = points.'; 

% Initial guesses (may change implementation)
% These are independent of n_samples
fnGuess   = [670, 1610];                % Hz
w_guess   = 2*pi*fnGuess;        % rad/s

zeta_guess = 0.006;

% These are 1 per sample
A_guess = 0.08;
Rl_guess = 0;
Rh_guess = 0;

% Estimate Parameters
p1 = repmat([A_guess, Rl_guess, Rh_guess], 1, n_samples);

f1 = [1, 1300]; % Hz
f2 = [850, 1800]; % Hz

for i=1:n_modes
    p0(i, :) = [w_guess(i), zeta_guess, p1];
    nPar (i) = numel(p0(i, :));

end

lb   = -Inf(1,nPar(1));      ub =  Inf(1,nPar(1));   % full length
lb(1) =   0;              ub(1) =  Inf;        % w_i  ≥ 0
lb(2) =   0;              ub(2) =  0.30;       % 0  ≤ ζ ≤ 0.30 (example)


% lb(4:3:end) = -1e-0; ub(4:3:end) =  1e-0;     % Rl constrained
% lb(5:3:end) = 0;  ub(5:3:end) = 0;  % Lock Rh = 0

opts = optimoptions('lsqnonlin','Display','iter');
for i= 1:n_modes
    p_est(i, :) = lsqnonlin(@(p)errfunc(FRFs,f1(i) ,f2(i) ,F_resp,p), ...
                  p0 (i, :), lb, ub, opts);
end

%% From here
Gnum_jki = {};

for i =1: n_modes
    % Searching for the closest frequency points to the margines
    [~, idx_f1(i)] = min(abs(F_resp - f1(i)));
    [~, idx_f2(i)] = min(abs(F_resp - f2(i)));
    
    Gnum = zeros(1, length(F_resp));
    Gnum_jk = zeros(n_samples, length(F_resp));

    % Estimated w_nat and zeta
    w_i_est    = p_est(i, 1);
    zeta_i_est = p_est(i, 2);
    
    for f = 1:length(F_resp)
        freq = F_resp(f);
        omega = 2*pi*freq;
        for idx_r = 1:n_samples
            r = 3 + 3*(idx_r-1);
            % Sample dependent parameters
            A_i_est   = p_est(i, r);
            Rl_est    = p_est(i, r+1);
            Rh_est    = p_est(i, r+2);
    
            Gnum_jk(idx_r, f) = freq_resp_numerical(omega, w_i_est, zeta_i_est, A_i_est, ...
            Rl_est, Rh_est);
        end
    end
    Gnum_jki{i} = Gnum_jk;
end

%% 2) Check the quality of the identification comparing the identified FRFs and the experimental ones.
Rls = p_est(4:3:end);
Rhs = p_est(5:3:end);

for i = 1:n_samples
    Magg = abs(H{i});
    % Plot magnitude for point i
    figure (i);
    subplot(2,1,1);
    semilogy(F_resp, Magg, 'b', 'LineWidth', 1.5); % Experimental magnitude
    grid on;
    hold on;

    % Plot phase
    subplot(2,1,2);
    plot(F_resp, rad2deg(angle(H{i})), 'LineWidth', 1.5);  % Experimental Phase
    grid on;
    hold on;%  Experimental
    
    % Semilogy(F_resp, abs(Gnum),'or');
    for g = 1:n_modes
        w = p_est(g, 1);
        f_search = w/(2*pi);

        [~, idx_w_nat] = min(abs(F_resp - f_search));
        marg = 500;
        i1 = idx_w_nat - marg;
        i2 = idx_w_nat + marg;
  

        % Select cell corresponding to correct mode estimation
        Gnums = Gnum_jki{g};
        Gnumm = Gnums(i,:);
        
        subplot(2,1,1);
        semilogy(F_resp(i1:i2), abs(Gnumm(i1:i2)),'or'); % Estimated magnitude
        xlabel('Frequency (Hz)');
        %xlim([3 6]);
        %ylim([1e-4 1e-1]);
        ylabel('|G| (m/N)');
        
        subplot(2,1,2);
        % plot(F_resp, rad2deg(angle(Gnum)),'or');
        plot(F_resp(i1:i2), rad2deg(angle(Gnumm(i1:i2))),'or'); % Estimated Phase
        % ylim([-270 270])
        xlabel('Frequency (Hz)');
        % xlim([3 6]);
        ylabel('Phase (rad)');
            
    end
end


%% Part 3: Plot a diagram showing the identified mode shapes with the indication of the corresponding natural frequencies and damping ratios.
% Measuring grid with regular angular spacing of 15°   -> define an angular spatial domain 
% Polar symmetry of the system                         -> Polarplot Matlab function

for g = 1: n_modes
    for i = 1: n_samples; theta1(i) = 15 * (2*pi/360) *(i-1); theta2(i) = 15 * (2*pi/360)*(i + n_samples); end   % Defining the 12th theta
    for i = 1: n_samples; rho1 (i) = p_est (g, 3*i); end           % Having the A values!
    M = max (abs(rho1));
    rho1 = rho1./M;
    rho2 = flip(rho1); % Normalizing by the maximum one

    rho1 = rho1 + 10*ones(size(theta1));
    rho2 = rho2 + 10*ones(size(theta2));
    figure;
    polaraxes;
    polarplot (theta1, rho1,'or-', 'LineWidth',1.5);     % Ploting the diagram
    hold on;
    polarplot(theta2, rho2, 'ob-', 'LineWidth',1.5)
    
    theta = [theta1 theta2];
    
    hold on; % Adding the 
    % Add a dashed line to show the reference zero point
    polarplot(theta, 10*ones(size(theta)), 'k--', 'LineWidth', 1);
    ax = gca;
    ax.ThetaZeroLocation = 'bottom';
    ax.ThetaDir = 'counterclockwise';
    ax.ThetaTick = 0:15:360;
    %ax.RLim = [0, 5];
end


















%% Annex: Functions

function Gnum = freq_resp_numerical(omega, w_i, zeta_i, A_i, Rl, Rh)    
    i = sqrt(-1);
    resonant = A_i/(-omega^2 + i*2*zeta_i*w_i*omega + w_i^2);
    low_freq = Rl/omega^2;
    high_freq = Rh;

    Gnum = resonant + low_freq + high_freq;
end

function epsilon = errfunc(FRFs, f1, f2, F_resp, p)
    n = size(FRFs, 1);  % Number of sampled FRFs (rows)

    % Frequency range indices
    [~, idx_f1] = min(abs(F_resp - f1));
    [~, idx_f2] = min(abs(F_resp - f2));

    F = F_resp(idx_f1:idx_f2);  % Frequency vector in selected range
    m = length(F);              % Number of frequency points in range

    % Estimated global parameters
    w_i    = p(1);
    zeta_i = p(2);

    % Preallocate residual vector
    epsilon = zeros(2 * n * m, 1);
    idx = 1;

    for idx_r = 1:n
        r = 3 + 3*(idx_r - 1);

        % Extract sample-specific parameters
        A_i = p(r);
        Rl  = p(r + 1);
        Rh  = p(r + 2);

        % Extract and slice experimental FRF
        Gexp = FRFs(idx_r, idx_f1:idx_f2);
        scale = max(abs(Gexp));  % or norm(Gexp)

        for kk = 1:m
            f = F(kk);
            omega = 2 * pi * f;

            Gnum = freq_resp_numerical(omega, w_i, zeta_i, A_i, Rl, Rh);

            epsilon(idx)   = real(Gexp(kk) - Gnum)/scale;
            epsilon(idx+1) = imag(Gexp(kk) - Gnum)/scale;

            idx = idx + 2;
        end
    end
end








%% Block 2 Locate the first two axial modes (initial estimates)
%for i = 1:12                     % 12 response positions
%    Hcol  = frf(:,i);            % complex FRF column (vector)
 %   mag   = abs(Hcol);           % magnitude of that FRF
%
 %   % --- do your peak search on this column ------------------
  %  [pks,locs] = findpeaks(mag,'MinPeakProminence',1e-3);
   % % keep the two highest peaks below 1000 Hz
    %below1k    = freq(locs) <= 1000;
%    [~,order]  = sort(pks(below1k),'descend');
 %   top2Idx    = locs(below1k);
  %  top2Idx    = top2Idx(order(1:2));       % indices of the two peaks
   % fnGuess    = freq(top2Idx);             % f_n1(0), f_n2(0)
   % ResGuess   = Hcol(top2Idx);             % complex residues
    % ----------------------------------------------------------



    % … store fnGuess, ResGuess if you need them later
% end



