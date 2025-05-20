clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw structure
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble mass and stiffness matrices
[M,K] = assem(incid,l,m,EA,EJ,gamma,idb);

%% Free Coordinate Submatrices
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);

MCF = M(ndof+1:end, 1:ndof);
KCF = K(ndof+1:end, 1:ndof);

%% Point 2) Eigenmodes/shapes
[modes, omega2] = eig(MFF\KFF);
omega = sqrt(diag(omega2));
% Sort frequencies in ascending order put
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;
% Sort mode shapes in ascending order 
modes = modes(:,i_omega);

%% Damping Matrix
alfa = 0.1;
beta = 2e-4;

C = alfa*M + beta*K;
CFF = C(1:ndof,1:ndof);
CCF = C(ndof+1:end, 1:ndof);

%% Force Vector
F0 = zeros(ndof, 1);
F0(idb(41,2)) = 1; % We have a vertical force (idx 2) on node 21 (row 21)

% Frequency range
om = (0:0.01:20)*2*pi;

for j = 1:length(om)
    i = sqrt(-1);
    A = -om(j)^2*MFF + i*om(j)*CFF + KFF;

    % Displacement FRF
    X(:,j) = A\F0;

    % Acceleration FRF
    Xpp(:,j) = -om(j)^2*X(:,j);

    % Reaction Forces FRF
    rr(:,j) = (-om(j)^2*MCF + i*om(j)*CCF + KCF)*X(:,j);
end

%% Point 3) Plotting FRF
% Transfer Function Vertical Displacement in A
idx_Ay = idb(41,2);
G_Ay = X(idx_Ay,:);

% Transfer Function Vertical Displacement in B
idx_By = idb(13,2);
G_By = X(idx_By,:);

% Plot Displacement FRFs
figure
subplot(2,1,1)
semilogy(om, abs(G_Ay), 'LineWidth', 1.5);
ylabel('abs(FRF)')
title('Displacement and Acceleration FRFs');
hold on
subplot(2,1,2)
plot(om, angle(G_Ay), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')
hold on

subplot(2,1,1)
semilogy(om, abs(G_By), 'LineWidth', 1.5);
ylabel('abs(FRF)')
title('Vertical Displacement FRFs');
legend(['Point A';'Point B'])
hold on
subplot(2,1,2)
plot(om, angle(G_By), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')
hold on


%% Mode Shape Plots
scale_factor = -2;
figure
for i = 1:2
    mode = modes(:,i);
    subplot(2,1,i)
    diseg2(mode,scale_factor,incid,l,gamma,posiz,idb,xy)
    title(sprintf('Mode Shape %d, freq %.2f rad/s', [i, omega(i)]));
end

%% Point 4 Modal Superposition
% Modal Matrices
ii = 1:2;
Phi = modes(:,ii);
Mmod = Phi'*MFF*Phi;
Kmod = Phi'*KFF*Phi;
Cmod = Phi'*CFF*Phi;
Fmod = Phi'*F0;

for ii = 1:length(om)
    i = sqrt(-1);
    xx_mod(:,ii) = (-om(ii)^2*Mmod + i*om(ii)*Cmod + Kmod)\Fmod;
end

xx_m = Phi*xx_mod;

% Transfer Function of Modal Superimposition
idx_By = idb(13,2);
G_By_superimposed = xx_m(idx_By,:);

figure
subplot(2,1,1)
semilogy(om, abs(G_By), 'LineWidth', 1.5);
ylabel('abs(FRF)')
title('FRF: Vertical Displacement in B vs Vertical force in A');
hold on
subplot(2,1,2)
plot(om, angle(G_By), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')
hold on

subplot(2,1,1)
semilogy(om, abs(G_By_superimposed), 'LineWidth', 1.5);
ylabel('abs(FRF)')
legend({'FEM'; 'Modal'})
hold on
subplot(2,1,2)
plot(om, angle(G_By_superimposed), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')
