function tmp=findcard(fid,card)
% Searchimg, within the file specified by fid, and place the file reading
% pointer on the next row. Reteurns 1 if the string is found, 0 otherwise
maxiter=1e5;

frewind(fid);
for i=1:maxiter
   if feof(fid)
      error(['The following string is not found: ' card])
   end
   riga=fgets(fid);
   if ~isempty(findstr(riga,card))
      return
   end
end

error(['The following string is not found: ' card])

