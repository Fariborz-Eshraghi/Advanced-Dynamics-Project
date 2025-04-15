clc
clear all
close all

%PHYSICAL properties of the beam
L  = 1200e-3; %m - length
h  = 8e-3; %m thickness
b  = 40e-3; %m width
rho = 2700; %kg/m3 density
E  = 68e9;  %Pa Youngs modulus

m = b* h* rho; % kg/m
J = b* h^3/12;

%%%%
% External force (Subject to change)
fo  = 500; %N - Amplitude
f_f = 5; %Hz - Frequency 
%%%%

%%%
% Set frequency range and frequency resolution
df   = 1.0000e-04; %Hz
fmax = 200; %Hz
%%%
% Set space domain and resolution
dx = 0.01;%m
x  = 0:dx:L;

F=0:df:fmax;
w=F.* 2.* pi;

c = (m/(E* J))^0.25;
G = w.^(1/2)* c; % gamma 

H=@(w)     [    1                                        0                        1                                     0;
                0                                        1                        0                                     1;
                -cos(w^(1/2)* c* L)                 -sin(w^(1/2)* c* L)               cosh(w^(1/2)* c* L)      sinh(w^(1/2)* c* L);
                sin(w^(1/2)* c* L)                  -cos(w^(1/2)* c* L)               sinh(w^(1/2)* c* L)      cosh(w^(1/2)* c* L)];


for i=1:length(w)
    dets(i)=det(H(w(i)));
end

figure(1), box on
semilogy(F,abs(dets),'-b')
hold on, grid on, xlabel('f [Hz]')

i_nat=[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
    end
end

plot(F(i_nat),abs(dets(i_nat)),'or') 