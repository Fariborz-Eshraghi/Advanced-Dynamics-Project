% Part B - 1st Assignment - Dynamics
clc;
clear all;

%% 1) Identify the Modal Parameters (Natural Frequencies, Damping Ratios, and Mode Shapes) -> first two axial modes
% Load the data

df = 0.333; % Resolution (Hz)
cd 'D:\Origins\Polimi\ADMS (Advance Dynamics of Mechanical Systems)\Assignments\First\Prat B';
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


Omega1 = 660;
Omega2 = 1600;


% Add points depending on preferences
points = frf;
F_resp = d.freq;
n_samples = size (points,2);
F_resp_size = size (points,1);
FRFs = zeros (n_samples, F_resp_size);
f_ranges = {};

FRFs = points.'; 

% Initial guesses (may change implementation)
% These are independent of n_samples
fnGuess   = Omega1;                % Hz
w_guess   = 2*pi*fnGuess;        % rad/s

zeta_guess = 0.01;

% These are 1 per sample
A_guess = 0.08;
Rl_guess = 0;
Rh_guess = 0;

p1 = repmat([A_guess, Rl_guess, Rh_guess], 1, n_samples);
p0 = [w_guess, zeta_guess, p1];

% Estimate Parameters
f1 = 1; %Hz
f2 = 1000; %Hz
% p_est = lsqnonlin(@(p) errfunc(FRFs, f1, f2, F_resp, p), p0, [], [], []);

nPar = numel(p0);
lb   = -Inf(1,nPar);      ub =  Inf(1,nPar);   % full length
lb(1) =   0;              ub(1) =  Inf;        % w_i  ≥ 0
lb(2) =   0;              ub(2) =  0.30;       % 0  ≤ ζ ≤ 0.30 (example)

opts = optimoptions('lsqnonlin','Display','iter');
p_est = lsqnonlin(@(p)errfunc(FRFs,f1,f2,F_resp,p), ...
                  p0, lb, ub, opts);


% Searching for the closest frequency points to the margines
[~, idx_f1] = min(abs(F_resp - f1));
[~, idx_f2] = min(abs(F_resp - f2));

Gnum = zeros(1, length(F_resp));
Gnum_jk = zeros(n_samples, length(F_resp));

% Estimated w_nat and zeta
w_i_est    = p_est(1);
zeta_i_est = p_est(2);

for f = 1:length(F_resp)
    freq = F_resp(f);
    omega = 2*pi*freq;
    for idx_r = 1:n_samples
        r = 3 + 3*(idx_r-1);
        % Sample dependent parameters
        A_i_est    = p_est(r);
        Rl_est     = p_est(r+1);
        Rh_est     = p_est(r+2);

        Gnum_jk(idx_r, f) = Gnum_jk(idx_r, f) + freq_resp_numerical (omega, w_i_est, zeta_i_est, A_i_est, ...
        Rl_est, Rh_est);
    end
 
end


%% 2) Check the quality of the identification comparing the identified FRFs and the experimental ones.

Gnumm = cell(1,12); 
for i = 1:12

    Magg = abs(H{i});            %  Experimental
    Gnumm = Gnum_jk(i,:);        %  Estimated Data
    
    % Plot magnitude for point i
    figure (i);
    subplot(2,1,1);
    semilogy(F_resp, Magg, 'b', 'LineWidth', 1.5); % Experimental magnitude
    grid on;
    hold on;

    % Semilogy(F_resp, abs(Gnum),'or');
    semilogy(F_resp(idx_f1:idx_f2), abs(Gnumm(idx_f1:idx_f2)),'or'); % Estimated magnitude
    xlabel('Frequency (Hz)');
    %xlim([3 6]);
    %ylim([1e-4 1e-1]);
    ylabel('|G| (m/N)');
    title('Magnitude Estimation');

    % Plot phase
    subplot(2,1,2);
    plot(F_resp, rad2deg(angle(H{i})), 'LineWidth', 1.5);  % Experimental Phase
    grid on;
    hold on;
    
    % plot(F_resp, rad2deg(angle(Gnum)),'or');
    plot(F_resp(idx_f1:idx_f2), rad2deg(angle(Gnumm(idx_f1:idx_f2))),'or'); % Estimated Phase
    % ylim([-270 270])
    xlabel('Frequency (Hz)');
    % xlim([3 6]);
    ylabel('Phase (rad)');
    title('Phase Estimation');
end


%% Part 3: Plot a diagram showing the identified mode shapes with the indication of the corresponding natural frequencies and damping ratios.
% Measuring grid with regular angular spacing of 15°   -> define an angular spatial domain 
% Polar symmetry of the system                         -> Polarplot Matlab function

for i = 1: 12; theta (i) = 15 * (2*pi/360) *(i-1); end   % Defining the 12th theta
for i = 13: 24; theta (i) = 15 * (2*pi/360) *(i); end   % Defining the mirrored ones
for i = 1: 12; rho1 (i) = p_est (3*i); end           % Having the A values!
M = max (abs(rho1));
rho2 = rho1./M; % Normilizing by the maximum one
rho = [rho2, flip(rho2)]; % Mirroring the numbers
polarplot (theta,rho);


%% Annex: Functions

function Gnum = freq_resp_numerical(omega, w_i, zeta_i, A_i, Rl, Rh)    
    i = sqrt(-1);
    resonant = A_i/(-omega^2 + i*2*zeta_i*w_i*omega + w_i^2);
    low_freq = Rl/omega^2;
    high_freq = Rh;

    Gnum = resonant + low_freq + high_freq;
end

function epsilon = errfunc(FRFs, f1, f2, F_resp, p)
    n = length(FRFs(:,1)');  % Number of sampled FRFs

    [~, idx_f1] = min(abs(F_resp - f1));
    [~, idx_f2] = min(abs(F_resp - f2));

    F = F_resp(idx_f1:idx_f2);
    F = F.';  
    m = length(F);  % Number of frequency points
    
    % Estimated w_nat and zeta
    w_i    = p(1);
    zeta_i = p(2);

    epsilon = zeros(2 * n * m, 1);
    idx = 1;

    for idx_r = 1:n
        r = 3 + 3*(idx_r-1);
        Gexp = FRFs(idx_r,:);
        
        % Sample dependent parameters
        A_i    = p(r);
        Rl     = p(r+1);
        Rh     = p(r+2);
        for fval = F 
            [~, pos] = min(abs(F_resp - fval));
            omega = fval*2*pi;
            Gnum = freq_resp_numerical(omega, w_i, zeta_i, A_i, Rl, Rh);
            
            epsilon(idx)   = real(Gexp(pos) - Gnum);   % Real part
            epsilon(idx+1) = imag(Gexp(pos) - Gnum);   % Imaginary part

            % Update index for the next residual
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



