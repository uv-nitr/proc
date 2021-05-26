function [msg,cal] = read_opus_cal(file)

% Data in CalFile: [ DateTime ENO3_(PathLength) ESW_(CalTemp) Reference ]
%
%        file: ''
%        info: [0x0 struct]
%         wvl: [1x0 double]
%        eno3: [1x0 double]
%         esw: [1x0 double] %referring to bromide seasalt, ebr
%         ref: [1x0 double] %referring to reference intensities in DIwater
%    esw_temp: []
%



msg = '';

z = zeros(1,0);
d = struct;

cal = struct( 'file' , { '' } , ...
              'info' , { d([]) } , ...
              'wvl'  , { z } , ...
              'eno3' , { z } , ...
              'esw'  , { z } , ...
              'ref'  , { z } , ...
              'esw_temp' , { [] } );

if nargin < 1
   return
end

if ~chkstr(file,1)

    msg = sprintf('Input must be a String for CalFile.');
    
elseif ~( exist(file,'file') == 2 )

    msg = sprintf('CalFile "%s" doesn''t exist.',file);

end

if ~isempty(msg)
    return
end


f = which(file);
if ~isempty(f)
    file = f;
end

[p,n,e] = fileparts(file);


if strcmp(lower(e),'.csv')

   [m,c] = rd_csv(file); 

   if ~isempty(m)
       msg = sprintf('Error call RD_CSV to read CSV-CalFile "%s".\n%s',file,m);
       return
   end

   if isempty(c)
       msg = sprintf('Empty Data in CSV-CalFile "%s".\n%s',file);
       return
   end
      
   h = fieldnames(c);

   c = struct2cell(c);

   c = cat(2,c{:});


else

   try
       [c,s,t] = xlsread(file);
   catch
       msg = sprintf('Error read XLS-CalFile "%s".\n%s',file,lasterr);
       return

   end


   if ~all( size(c,2) == [ size(s,2) size(t,2) ] )
       msg = sprintf('Invalid Data in XLS-CalFile "%s".\n%s',file, ...
               'Number of Columns of Data and Header Mismatch');
       return
   end

   h = t(1,:);

end

if ~( ( size(c,1) > 1 ) & ( size(c,2) >= 4 ) )
    msg = sprintf('Invalid data in CalFile "%s", 4 columns required.',file);
    return
end

ct = h{3};
ii = findstr(ct,'_');
if ~isempty(ii)
    ct = ct((max(ii)+1):end);  % Remaining String after "_"
end

[m,cal_temp] = str2vec(ct);    % Retrieve Number from String

if ~( isempty(m) & ( prod(size(cal_temp)) == 1 ) )
    msg = sprintf('Can''t retrieve SW-CalTemp from 3. Header "%s" in CalFile "%s".', ...
                   h{3},file);
    return
end

c = permute(c,[2 1]);

cal = struct( 'file' , { file } , ...
              'info' , { dir(file) } , ...
              'wvl'  , { c(1,:) } , ...
              'eno3' , { c(2,:) } , ...
              'esw'  , { c(3,:) } , ...
              'ref'  , { c(4,:) } , ...
              'esw_temp' , { cal_temp } );
