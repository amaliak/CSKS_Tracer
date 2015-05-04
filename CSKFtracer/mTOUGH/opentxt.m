function opentxt(filename)

   fprintf('You have requested file: %s\n', filename);

   wh = which(filename);
   if exist(filename, 'file') == 2
     fprintf('Opening in MATLAB Editor: %s\n', filename);
     edit(filename);
   elseif ~isempty(wh)
     fprintf('Opening in MATLAB Editor: %s\n', wh);
     edit(wh);
   else
     warning('MATLAB:fileNotFound', ...
             'File was not found: %s', filename);
   end

end 