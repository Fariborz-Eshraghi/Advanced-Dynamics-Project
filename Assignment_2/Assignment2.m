clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------Point 1) FE Model Definition---------------------------
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure;

% Draw structure
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);

% Assemble mass and stiffness matrices
[M,K] = assem(incid,l,m,EA,EJ,gamma,idb);

% Free Coordinate Submatrices
MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);

MCF = M(ndof+1:end, 1:ndof);
KCF = K(ndof+1:end, 1:ndof);

%% --------------------- Point 2) Eigenmodes/shapes------------------------
[modes, omega2] = eig(MFF\KFF);
omega = sqrt(diag(omega2));

% Sort frequencies in ascending order
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;

% Sort mode shapes in ascending order 
modes = modes(:,i_omega);

%----------------- Mode Shape Plots up to 3rd mode-------------------------
scale_factor = -2;
figure
for i = 1:3
    mode = modes(:,i);
    subplot(3,1,i)
    diseg2(mode,scale_factor,incid,l,gamma,posiz,idb,xy)
    title(sprintf('Mode Shape %d, freq %.2f rad/s', [i, omega(i)/(2*pi)]));
end

% Damping Matrix
alfa = 0.1;
beta = 2e-4;

C = alfa*M + beta*K;
CFF = C(1:ndof,1:ndof);
CCF = C(ndof+1:end, 1:ndof);

%% ----------- Point 3) FRF of Force applied on point A -------------------
F0 = zeros(ndof, 1);
F0(idb(41,2)) = 1; % We have a vertical force (idx 2) on node 21 (row 21)

% ---------------------- Generate FRF arrays-------------------------------
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

%----------------------- Obtain FRFs -------------------------------------
% Transfer Function Vertical Displacement in A
idx_Ay = idb(41,2);
G_Ay = X(idx_Ay,:);

% Transfer Function Vertical Displacement in B
idx_By = idb(13,2);
G_By = X(idx_By,:);

% --------------------- Plot Displacement FRFs ---------------------------
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
title('Point 3) Vertical Displacement FRFs');
legend(['Point A';'Point B'])
hold on
subplot(2,1,2)
plot(om, angle(G_By), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')
hold on

%% ----------------Point 4) Modal Superposition-----------------------------
% Modal Matrices
ii = 1:3;
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
G_By_superimposed = xx_m(idx_By,:);

%-------------------Modal Superimposition Plots----------------------------
figure
subplot(2,1,1)
semilogy(om, abs(G_By), 'LineWidth', 1.5);
ylabel('abs(FRF)')
title('Point 4) FRF: Vertical Displacement in B vs Vertical force in A');
hold on
subplot(2,1,2)
plot(om, angle(G_By), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')
hold on

subplot(2,1,1)
semilogy(om, abs(G_By_superimposed), 'LineWidth', 1);
ylabel('abs(FRF)')
legend({'FEM'; 'Modal'})
hold on
subplot(2,1,2)
plot(om, angle(G_By_superimposed), 'LineWidth', 1);
xlabel('Frequency (rad/s)')
ylabel('Phase (rad/s)')

%% -------------Point 5) Static Response of Weight and Deflection-------------------------
acc_0 = zeros(ndof, 1);
acc_0(idb(:,2)) = -9.81;
% acc_0(2:3:end) = -9.81;

F_grav = M*acc_0;

xF = KFF\F_grav(1:ndof);
figure
diseg2(xF,80,incid,l,gamma,posiz,idb,xy)
title('Point 5) Static Deflection');

[max_deflection, idx] = max(abs(xF))

%% -----------Point 6) Moving load time history--------------------------------------------

F = 50000; %[N]
v_M = 2; %[m/s]
dt = 0.01; %[s]
dx = v_M*dt;
Fn_global = [];

T_rot = @(c,s) [
    c  s  0   0   0   0;
   -s  c  0   0   0   0;
    0  0  1   0   0   0;
    0  0  0   c   s   0;
    0  0  0  -s   c   0;
    0  0  0   0   0   1];

for idx = 13:40
    if idx == 13
        a_vals = 0:dx:l(idx);  % include the initial point only once
    else
        a_vals = dx:dx:l(idx); % skip a = 0 to avoid double-counting at shared node
    end

    c = cos(gamma(idx));
    s = sin(gamma(idx));
    T = T_rot(c,s);

    for i = 1:length(a_vals)
        a = a_vals(i);
        b = l(idx) - a;

        % Global force (vertical)
        F_global = [0; -F];
        R = [c s; -s c];
        F_local = R * F_global;
        Pu = F_local(1);
        Pv = F_local(2);

        % Equivalent nodal local forces
        f1 = Pv * b^2 * (3*a + b) / l(idx)^3;
        m1 = Pv * a * b^2 / l(idx)^2;
        f2 = Pv * a^2 * (3*b + a) / l(idx)^3;
        m2 = -Pv * a^2 * b / l(idx)^2;

        fa1 = Pu * b / l(idx);
        fa2 = Pu * a / l(idx);

        Fn_local = [fa1; f1; m1; fa2; f2; m2];

        % Transform to global
        Fn_local_rot = T' * Fn_local;

        % Assemble into global vector
        Temp = zeros(size(M,1),1);
        Temp(incid(idx,:)) = Fn_local_rot;

        % Store
        Fn_global = [Fn_global Temp];
    end
end

Fn_global(:,end+2/dt) = 0;
Fn_global(ndof+1:end,:) = [];

t_vect = (1:size(Fn_global,2))*dt;

figure
plot(t_vect, Fn_global')
ylabel('Nodal forces')
xlabel('Time [s]')

Qn_global = Phi'*Fn_global;
figure
plot(t_vect, Qn_global)
grid on
ylabel('Q')
xlabel('Time [s]')
legend({'Mode 1', 'Mode 2', 'Mode 3'})

n_modes = 3;
X0 = zeros(1,2*n_modes);
[t, X] = ode45(@(t,X) odefn(t,X, t_vect,Qn_global,Mmod,Kmod,Cmod), t_vect, X0);

% Step 1: Get modal displacements only
q_modal = X(:, 1:3);  % 1589 × 3

% Step 2: Reconstruct physical DOFs
TH = q_modal * Phi';  % 1589 × 221
TH = TH';             % 221 × 1589

% Step 3: Get displacement at DOF 41, Y direction
X_displ = TH(idb(41,2), :);  % row vector

% Step 4: Plot
figure(1000)
plot(t_vect, X_displ), grid on
xlabel('Time [s]')
ylabel('y_A [m]')
title('Point 6) Vertical Displacement History Point A, Load: 50kN, Velocity: 2m/s')

%% -----------Point 7) Proposal for reducing the static deformation--------------------------------------------
total_mass1 = l*m';


[file_i2,xy2,nnod2,sizee2,idb2,ndof2,incid2,l2,gamma2,m2,EA2,EJ2,posiz2,nbeam2,pr2] = loadstructure;

% Draw structure
dis_stru(posiz2,l2,gamma2,xy2,pr2,idb2,ndof2);

% Assemble mass and stiffness matrices
[M2,K2] = assem(incid2,l2,m2,EA2,EJ2,gamma2,idb2);

% Free Coordinate Submatrices
MFF2 = M2(1:ndof2,1:ndof2);
KFF2 = K2(1:ndof2,1:ndof2);

MCF2 = M2(ndof2+1:end, 1:ndof2);
KCF2 = K2(ndof2+1:end, 1:ndof2);

%-------------------------------------------------

total_mass2 = l2*m2';

acc_02 = zeros(ndof2, 1);
acc_02(idb2(:,2)) = -9.81;
% acc_0(2:3:end) = -9.81;

F_grav2 = M2*acc_02;

xF2 = KFF2\F_grav2(1:ndof2);
figure
total_deformation2 = diseg3(xF2,80,incid2,l2,gamma2,posiz2,idb2,xy2);
title('Point 7) Static Deflection proposal');

%--------------------
perc_mass = total_mass2/total_mass1
prtg_of_improvement = (total_deformation1-total_deformation2)/total_deformation1


%% ODE FUNCTIONS FOR POINT 6)

function dXdt = odefn(t, X, t_vect, Qn_global, Mmod, Kmod, Cmod)
    % Interpolate the modal force Qn_global at time t
    Q = interp1(t_vect, Qn_global', t, 'linear', 0)';  % Linear interp, 0 if out of bounds

    n = length(X)/2;
    q = X(1:n);       % Modal displacements
    q_dot = X(n+1:end); % Modal velocities

    % Compute derivatives
    dq = q_dot;
    dq_dot = Mmod \ (Q - Cmod*q_dot - Kmod*q);

    % Return concatenated derivative vector
    dXdt = [dq; dq_dot];
end


%% DISSEG MODIFICATION FOR POINT 7) %%

function sum_def = diseg3(mode,scale_factor,incidenze,l,gamma,posiz,idb,xy)
    sum_def = 0;
    % Checking consistency input data
    [n_el,nc]=size(incidenze);
    
    if length(posiz) ~= n_el
        sprintf('Error: number of nodes in posit matrix differs from number of elements')
        return
    end
    
    n_gdl=length(mode);
    
    
    hold on
    % looping on the finite elements
    for k=1:n_el
    % building the nodal displacements vector of each element in the global
    % reference frame
        for iri=1:6
            if incidenze(k,iri) <= n_gdl
            xkG(iri,1)=mode(incidenze(k,iri));
            else
            xkG(iri,1)=0.;
            end
        end
    % Applying scale factor
        xkG=scale_factor*xkG;
    % Global to Local reference frame rotation
        lambda = [ cos(gamma(k)) sin(gamma(k)) 0. 
                  -sin(gamma(k)) cos(gamma(k)) 0.
	                    0.      0.     1. ] ;
        Lambda = [ lambda     zeros(3)
                  zeros(3)   lambda      ] ;
        xkL=Lambda*xkG;
    
    % Computing the axial (u) and transversal (w) displacements by means of
    % shape functions
        csi=l(k)*[0:0.05:1];
        fu=zeros(6,length(csi));
        fu(1,:)=1-csi/l(k);
        fu(4,:)=csi/l(k);
        u=(fu'*xkL)';
        fw=zeros(6,length(csi));
        fw(2,:)=2*(csi/l(k)).^3-3*(csi/l(k)).^2+1;
        fw(3,:)=l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k));
        fw(5,:)=-2*(csi/l(k)).^3+3*(csi/l(k)).^2;
        fw(6,:)=l(k)*((csi/l(k)).^3-(csi/l(k)).^2);
        w=(fw'*xkL)';
        
        sum_def = sum_def + norm(u+w*1j);
    
    
    % Local to global transformation of the element's deformation
        xyG=lambda(1:2,1:2)'*[u+csi;w];
        undef=lambda(1:2,1:2)'*[csi;zeros(1,length(csi))];
     % Plotting deformed and undeformed element's shape
        plot(undef(1,:)+posiz(k,1),undef(2,:)+posiz(k,2),'b--')
        plot(xyG(1,:)+posiz(k,1),xyG(2,:)+posiz(k,2),'b')
    end
    
    % Looping the nodes
    n_nodi=size(idb,1);
    xkG=[];
    for k=1:n_nodi
        for ixy=1:2
            if idb(k,ixy) <= n_gdl
            xkG(k,ixy)=mode(idb(k,ixy));
            else
            xkG(k,ixy)=0.;
            end
        end
    end
    xkG=scale_factor*xkG;
    xyG=xkG+xy;
    plot(xy(:,1),xy(:,2),'b.')
    plot(xyG(:,1),xyG(:,2),'bo')
    
    grid on
    box on
    axis equal

end