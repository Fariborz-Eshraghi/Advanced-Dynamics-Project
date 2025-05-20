function diseg2(mode,scale_factor,incidenze,l,gamma,posiz,idb,xy)

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
