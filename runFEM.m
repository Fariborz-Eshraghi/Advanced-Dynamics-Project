clear all;
close all;
clc;
[file_i, xy, nnod, sizee, idb, ndof, incid, l, gamma, m, EA, EJ, posiz, nbeam, pr] = loadstructure;
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);
[M,K] = assem(incid,l,m,EA,EJ,gamma,idb);
freeDofs = sort(idb(idb <= ndof));
%assert(numel(freeDofs)==ndof);
MFF= M(freeDofs, freeDofs);
KFF = K(freeDofs, freeDofs);
%MFF = M(1:ndof,1:ndof);
%KFF = K(1:ndof,1:ndof);
[modes,omega2] = eig(MFF\KFF);
omega = sqrt(diag(omega2));
[omega, i_omega] = sort(omega);
freq0 = omega/(2*pi) ;
modes = modes(:, i_omega);
disp(freq0(1:4));
for i= 2:5
    figure;
    diseg2( -modes(:,i),1,incid,l,gamma,posiz,idb,xy);
end