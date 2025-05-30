clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL properties of the beam
L  = 1200e-3; %m - length
h  = 8e-3; %m thickness
b  = 40e-3; %m width
rho = 2700; %kg/m3 density
E  = 68e9;  %Pa Youngs modulus

m = b*h*rho; % kg/m
J = b*h^3/12;

%%%%
%% External force (Subject to change)
fo  = 500; %N - Amplitude
f_f = 5; %Hz - Frequency 
%%%%

%%%
%% Set frequency range and frequency resolution
df   = 1.0000e-04; %Hz
fmin = 0; %Hz
fmax = 200; %Hz
freq_size = (fmax - fmin)/df;
%%%

%% Set space domain and resolution
n_points = 1000;
x = linspace(0, L, n_points);
dx = L/n_points;


F=0:df:fmax;
w=F.*2.*pi;

c = (m/(E*J))^0.25; 

H=@(w)     [    1                                        0                        1                                     0;
                0                                        1                        0                                     1;
                -cos(w^(1/2)*c*L)                 -sin(w^(1/2)*c*L)               cosh(w^(1/2)*c*L)      sinh(w^(1/2)*c*L);
                sin(w^(1/2)*c*L)                  -cos(w^(1/2)*c*L)               sinh(w^(1/2)*c*L)      cosh(w^(1/2)*c*L)];


for i=1:length(w)
    dets(i)=det(H(w(i)));
end


%% Find Natural Frequencies
[pks, locs] = findpeaks(-abs(dets), 'MinPeakProminence',0.1);
f_nat = locs.*df;
w_nat = locs.*df.*2*pi;

figure(1), box on
semilogy(F,abs(dets),'-b')
hold on, grid on, xlabel('f [Hz]')

plot(F(locs),abs(dets(locs)),'or')

%% Find Analytical Mode Shapes
for i=1:length(w_nat)
   H_i = H(w_nat(i));
   H_ihat = H_i(2:end, 2:end );
   N_ihat = H_i(2:end, 1);
   
   X_ihat = -inv(H_ihat)*N_ihat;
   X_i = [1, X_ihat'];
   xj_n{i} = X_i;
   xj{i} = X_i;
end

xj = reshape(cell2mat(xj),4,[]);
xj_n = reshape(cell2mat(xj_n),4,[]);

n_modes = length(w_nat);
phi_matrix = zeros(length(x), n_modes);
phi_matrix_normalized = zeros(length(x), n_modes);
gamma = sqrt(w_nat)*c;

figure(2)

for j=1:n_modes
   A = xj_n(1,j);
   B = xj_n(2,j);
   C = xj_n(3,j);
   D = xj_n(4,j);
   
   phi_matrix(:,j) = A*cos(gamma(j)*x) + B*sin(gamma(j)*x) + C*cosh(gamma(j)*x) + D*sinh(gamma(j)*x);
   phi_matrix_normalized(:,j) = phi_matrix(:,j)/max(abs(phi_matrix(:,j)));

   plot(x, phi_matrix_normalized(:,j),'LineWidth',1.5);
   xlabel('Position (m)')
   ylabel('Displacement (m)')
   title('Normalized Mode Shapes')
   hold on
   grid on
end


%% Analytical Frequency Response

% Set Parameters
zeta_i = 0.01;
f1 = fmin + df;
f2 = fmax;
F_resp_size = 10000;

% indices over the span of vector x ([0, 1000])
pos_i = 166;
pos_k = 1000;

% Compute Frequency Response
[G, F_resp]  = freq_resp(phi_matrix_normalized, w_nat, pos_i, pos_k, m, zeta_i, x, ...
                         f1, f2, F_resp_size);

% Compute magnitude and phase
FRF_magnitude = abs(G);
FRF_phase = rad2deg(angle(-G)); % convert to degrees

% Plot magnitude
figure;
subplot(2,1,1);
semilogy(F_resp, FRF_magnitude, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('|G| (m/N)');


% Plot phase
subplot(2,1,2);
plot(F_resp, FRF_phase, 'r', 'LineWidth', 1.5);
grid on;
ylim([-270 270])
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');

%% Modal Identification

% Add points depending on preferences, need to change how to implement it
% for given data set in part B
points = [[166, 1000]; [400, 1000]; [600, 1000]; [950, 1000]];
n_samples = length(points);
FRFs = zeros(n_samples, F_resp_size);
f_ranges = {};

% "Experimental" frequency responses
for k=1:n_samples
    [FRFs(k,:), f_ranges{k}] = freq_resp(phi_matrix_normalized, w_nat, points(k,1), points(k,2), m, zeta_i, x, f1, f2, F_resp_size);
end

% Initial guesses (may change implementation)
% These are independent of n_samples
'MODAL IDENTIFICATION: Initial Guesses'
mode_number = input('Which mode are you trying to identify?\n n = ');
w_guess = input('Initial guess for natural frequency (Hz):\n w_guess = ')*2*pi;
zeta_guess = input('Initial guess for damping ratio (0-1 values):\n zeta_guess = ');

% These are 1 per sample
A_guess = input('Initial Guess for Resonance Gain:\n A_guess = ');
Rl_guess = 0;
Rh_guess = 0;

p1 = repmat([A_guess, Rl_guess, Rh_guess], 1, n_samples);
p0 = [w_guess, zeta_guess, p1];

% Estimate Parameters
f1 = input('Lower Frequency (Hz) for estimation:\n f1 = '); %Hz
f2 = input('Upper Frequency (Hz) for estimation:\n f2 = '); %Hz
p_est = lsqnonlin(@(p) errfunc(FRFs, f1, f2, F_resp, p), p0, [], [], []);

% Estimated w_nat and zeta
w_i_est    = p_est(1);
zeta_i_est = p_est(2);

fprintf("Estimated Natural frequency Value for Mode %d: %.2f rad/s (%.2f Hz)\n", mode_number, w_i_est, w_i_est/(2*pi));
fprintf("Estimated damping ratio for Mode %d: zeta = %.2f\n", mode_number, zeta_i_est);


Gnum = zeros(1, length(F_resp));
Gnum_jk = zeros(n_samples, length(F_resp));

for f = 1:length(F_resp)
    freq = F_resp(f);
    omega = 2*pi*freq;
    for idx_r = 1:n_samples
        r = 3 + 3*(idx_r-1);
        % Sample dependent parameters
        A_i_est    = p_est(r);
        Rl_est     = p_est(r+1);
        Rh_est     = p_est(r+2);

        Gnum_jk(idx_r, f) = Gnum_jk(idx_r, f) + freq_resp_numerical(omega, w_i_est, zeta_i_est, A_i_est, ...
        Rl_est, Rh_est);
    end
 
end

Gnum_1 = Gnum_jk(1,:);
%Gnum_1 = Gnum_1(idx_f1:idx_f2);

w = w_i_est;
f_search = w/(2*pi);

[~, idx_w_nat] = min(abs(F_resp - f_search));
marg = 500;
if idx_w_nat - marg < 1
    idx_f1 = 1;
else
    idx_f1 = idx_w_nat - marg;
end

if idx_w_nat + marg > length(F_resp)
    idx_f2 = F_resp(end);
else
    idx_f2 = idx_w_nat + marg;
end

% Plot magnitude for point 1
figure;
semilogy(F_resp, FRF_magnitude, 'b', 'LineWidth', 1.5);
grid on;
hold on;
%semilogy(F_resp, abs(Gnum),'or');
semilogy(F_resp(idx_f1:idx_f2), abs(Gnum_1(idx_f1:idx_f2)),'or');
xlabel('Frequency (Hz)');
%xlim([3 6]);
%ylim([1e-4 1e-1]);
ylabel('|G| (m/N)');
title('Magnitude Estimation');


% Plot phase
figure
plot(F_resp, FRF_phase, 'LineWidth', 1.5);
grid on;
hold on;
%plot(F_resp, rad2deg(angle(Gnum)),'or');
plot(F_resp(idx_f1:idx_f2), rad2deg(angle(-Gnum_1(idx_f1:idx_f2))),'or');
ylim([-270 270])
xlabel('Frequency (Hz)');
%xlim([3 6]);
ylabel('Phase (rad)');
title('Phase Estimation');

%% Mode Shape Comparison

% Here we use our estimation around the 1st natural frequency, so we have
% to compare with the 1st mode shape obtained in parts 1 and 2

A1_est = p_est(3:3:end);
phi_exp = real(A1_est);

% Find maximum experimental value and corresponding index
[max_val, max_pos] = max(A1_est);
max_idx = points(max_pos, 1);

% Normalize Experimental Mode Shape
phi_exp_norm = phi_exp/max_val;

x_positions = points(:,1).*dx;

figure;
plot(x, (-phi_matrix(:,mode_number)/phi_matrix(max_idx, mode_number)),'LineWidth',1.5);
hold on;
grid on;
plot(x_positions, -phi_exp_norm, 'or', 'DisplayName', 'Experimental');
xlabel('Position (m)');
ylabel('Displacement');
legend('Model', 'Identified');
title('Mode Shape Comparison');

%% Annex: Functions
function [Gexp, F_resp] = freq_resp(phi_matrix, w_nat, pos_i, pos_k, m, ...
                            zeta_i, x, f1, f2, F_resp_size)
    Gexp = [];
    j = sqrt(-1);
    F_resp = linspace(f1, f2, F_resp_size);
    n_modes = length(w_nat);
    
    for f = F_resp
        sum = 0;
        for i=1:n_modes
            m_i = trapz(x, m.*phi_matrix(:,i).^2);
            omega = f*2*pi;
            w_i =  w_nat(i);
            sum = sum + (phi_matrix(pos_i,i)*phi_matrix(pos_k,i)/m_i)/( ...
                -omega^2 + j*2*zeta_i*w_i*omega + w_i^2);
        end
        
        Gexp(end+1) = sum;
    end
    
end

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



