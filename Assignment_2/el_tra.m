function [mG,kG] = el_tra (l,m,EA,EJ,gamma)

% mass matrix in the local reference frame
mL = m*l*[1./3.  0.          0.          1./6.  0.          0.
          0.     13./35.     11.*l/210.  0.     9./70.      -13*l/420.
          0.     11.*l/210.  l^2/105.    0.     13*l/420.   -l^2/140.
          1./6.  0.          0.          1./3.  0.          0.
          0.     9./70.      13*l/420.   0.     13./35.     -11.*l/210.
          0.     -13*l/420.  -l^2/140.   0.     -11.*l/210. l^2/105.   ] ;

% stiffness matrix in the local reference frame
% contribution due to the axial deformation
kL_ax = EA/l* [ 1 0 0 -1 0 0
                0 0 0  0 0 0 
                0 0 0  0 0 0 
                -1 0 0  1 0 0 
                0 0 0  0 0 0 
                0 0 0  0 0 0 ] ; 

% contribution due to bending deformation
kL_fl = EJ * [ 0.    0.       0.      0.    0.       0.     
               0.  +12./l^3  6./l^2   0.  -12./l^3  6./l^2
               0.   6./l^2  +4./l     0.   -6./l^2  +2./l
               0.    0.       0.      0.    0.       0. 
               0.  -12./l^3  -6./l^2  0.  +12./l^3  -6./l^2
               0.   6./l^2  +2./l    0.   -6./l^2  +4./l    ] ;
           
kL = kL_ax+kL_fl ;

% matrix transformation from local to global reference frame
% rotation matrix 3x3 (single node DoF)
lambda = [ cos(gamma) sin(gamma) 0. 
          -sin(gamma) cos(gamma) 0.
              0.         0.      1.] ;

% rotation matrix 6x6 (element DoF)
Lambda = [ lambda     zeros(3)
           zeros(3)   lambda      ] ;

mG = Lambda' * mL * Lambda ;
kG = Lambda' * kL * Lambda ;