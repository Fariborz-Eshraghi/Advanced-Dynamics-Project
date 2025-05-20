function line=scom(fid)
line=fgetl(fid);
tmp=sscanf(line,'%s',1);
while(tmp=='!')
    if feof(fid)
      error('Impossible to find uncommented string')
   end
 line=fgetl(fid);
 tmp=sscanf(line,'%s',1);
%  disp('in scom')
%  pause
end 

