function [file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure

% Loads *.inp file
disp(' ');
file_i=input(' Name of the input file *.inp (without extension) for the structure to be analyzed = ','s') ;
disp(' ')
% Check input file
if exist([file_i '.inp'])~=2
    disp(['File ',file_i,'.inp does not exist' ])
    file_i=[];
    return
end
nprova=file_i;

% File opening
eval(['fid_i=fopen(''',file_i,'.inp'',''r'');']);

% Reading card NODES
findcard(fid_i,'*NODES')
% line=scom(fid_i)
% finewhile=findstr(line,'*ENDNODES')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDNODES'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i %f %f')';
        iconta=iconta+1;
        if iconta ~=tmp(1)
        disp('Errore: nodi non numerati in ordine progressivo')
        break  
        end
        ivinc(iconta,:)=tmp(2:4);
        xy(iconta,:)=tmp(5:6);
    end
end
nnod=iconta;
disp(['Number structure nodes ',int2str(nnod)])
sizee=sqrt((max(xy(:,1))-min(xy(:,1)))^2+(max(xy(:,2))-min(xy(:,2)))^2);
% End of reading card NODES

% Building IDB matrix
idof=0;
for i=1:nnod
    for j=1:3
        if ivinc(i,j) == 0
            idof=idof+1;
            idb(i,j)=idof;
        end
    end
end
ndof=idof;
disp(['Number of structure d.o.f. ',int2str(ndof)])
for i=1:nnod
    for j=1:3
        if ivinc(i,j) == 1
            idof=idof+1;
            idb(i,j)=idof;
        end
    end
end


% reading card BEAMS
findcard(fid_i,'*PROPERTIES')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDPROPERTIES'));
    if finewhile == 1
        tmp=sscanf(line,'%i %f %f %f')';
        iconta=iconta+1;
        properties(iconta,:) = tmp(2:4);
    end
end
disp(['Number of properties ',int2str(iconta)])

% reading card BEAMS
findcard(fid_i,'*BEAMS')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDBEAMS'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i')';
        iconta=iconta+1;
        incid(iconta,:)=tmp(2:3);
        pr(1,iconta) = tmp(4);
        m(1,iconta)  = properties(tmp(4),1);
        EA(1,iconta) = properties(tmp(4),2);
        EJ(1,iconta) = properties(tmp(4),3);
        l(1,iconta)=sqrt((xy(incid(iconta,2),1)-xy(incid(iconta,1),1))^2+(xy(incid(iconta,2),2)-xy(incid(iconta,1),2))^2);
        gamma(1,iconta)=atan2(xy(incid(iconta,2),2)-xy(incid(iconta,1),2),xy(incid(iconta,2),1)-xy(incid(iconta,1),1));
        incidenze(iconta,:)=[idb(incid(iconta,1),:) idb(incid(iconta,2),:)];
        posiz(iconta,:)=xy(incid(iconta,1),:);
    end
end
nbeam=iconta;
disp(['Number of beam FE ',int2str(nbeam)])





fclose(fid_i);