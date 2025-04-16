clc
clear all
close all

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
fmax = 200; %Hz
%%%
%% Set space domain and resolution
n_points = 1000;
x = linspace(0, L, n_points);

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
[pks, locs] = findpeaks(-abs(dets), 'MinPeakProminence',1);
f_nat = locs.*df;
w_nat = locs.*df.*2*pi;

figure(1), box on
semilogy(F,abs(dets),'-b')
hold on, grid on, xlabel('f [Hz]')

plot(F(locs),abs(dets(locs)),'or')

phi = cell(1, length(w_nat));
phi_n = cell(1, length(w_nat));


%% Find Mode Shapes
for i=1:length(w_nat)
   H_i = H(w_nat(i));
   H_ihat = H_i(2:end, 2:end );
   N_ihat = H_i(2:end, 1);
   
   X_ihat = -inv(H_ihat)*N_ihat;
   X_i = [1, X_ihat'];
   xj_n{i} = X_i/max(abs(X_i));
   xj{i} = X_i;
end

xj = reshape(cell2mat(xj),4,[]);
xj_n = reshape(cell2mat(xj_n),4,[]);

n_modes = length(w_nat);
phi_matrix = zeros(length(x), n_modes);
gamma = sqrt(w_nat)*c;

figure(2) 

for j=1:n_modes
   A = xj_n(1,j);
   B = xj_n(2,j);
   C = xj_n(3,j);
   D = xj_n(4,j);
   
   phi_matrix(:,j) = A*cos(gamma(j)*x) + B*sin(gamma(j)*x)+C*cosh(gamma(j)*x)+ D*sinh(gamma(j)*x);
   plot(x, phi_matrix(:,j),'LineWidth',1.5);
   xlabel('Position (m)')
   ylabel('Displacement (m)')
   title('Normalized Mode Shapes')
   hold on
end


%% Frequency Response
zeta_i = 0.01;
F_resp = linspace(0, fmax, 10000);

% indices over the span of vector x ([0, 1000])
pos_i = 166;
pos_k = 1000;
x_i = L*pos_i/n_points;
x_k = L*pos_k/n_points;

G = [];
vals=zeros(4,2);
i = sqrt(-1);


for f=1:length(F_resp)
    sum = 0;
    for j=1:n_modes
        m_i = trapz(length(phi_matrix), m.*phi_matrix(:,j).^2);
        sum = sum + phi_matrix(pos_i, j)*phi_matrix(pos_k, j)/m_i/(-(F_resp(f)*2*pi)^2+2*i*zeta_i*w_nat(j)*(F_resp(f)*2*pi)+w_nat(j)^2);
    end
    
    G(f) = sum;
end

% Compute magnitude and phase
FRF_magnitude = abs(G);
FRF_phase = angle(G); % convert to degrees

% Plot magnitude
figure;
subplot(2,1,1);
semilogy(F_resp, FRF_magnitude, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('|FRF|');
title('Frequency Response Function - Magnitude');

% Plot phase
subplot(2,1,2);
plot(F_resp, FRF_phase, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (°)');
title('Frequency Response Function - Phase');