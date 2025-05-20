function [M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb)

% Checking consistency input data
[n_el,nc]=size(incidenze);
if nc ~= 6
    disp('Error: number of columns of incidence matrix different from 6')
    return
end
if length(l) ~= n_el
    sprintf('Error: number of elements in l differenc from n')
    return
end
if length(m) ~= n_el    
    sprintf('Error: number of elements in m differenc from number of elements')
    return
end
if length(EA) ~= n_el
    sprintf('Error: number of elements in EA differenc number of elements')
    return
end
if length(EJ) ~= n_el
    sprintf('Error: number of elements in EJ differenc number of elements')
    return
end
if length(gamma) ~= n_el
    sprintf('Error: number of elements in alpha differenc number of elements')
    return
end

if min(min(incidenze)) ~= 1    
    sprintf('Error: dof sequence does not start from 1')
    return
end

% Limiting the total dof of the structure
n_dof=max(max(idb));
% if n_dof > 100
%     sprintf('Errore: dof > 100, not allowed!')
%     return
% end

% Assembling matrices M and K
M=zeros(n_dof,n_dof);
K=zeros(n_dof,n_dof);
for k=1:n_el
    [mG,kG] = el_tra(l(k),m(k),EA(k),EJ(k),gamma(k));
    for iri=1:6
        for ico=1:6
            i1=incidenze(k,iri);
            i2=incidenze(k,ico);
            M(i1,i2)=M(i1,i2)+mG(iri,ico);
            K(i1,i2)=K(i1,i2)+kG(iri,ico);
        end
    end
end