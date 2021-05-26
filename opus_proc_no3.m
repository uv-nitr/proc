function [msg,dat] = opus_proc_no3(mode,dfile,cfile,varargin)

% PROC_NO3 Processing of OPUS Nitrate Measurements
%
% [ Msg, DataStruct ] = opus_proc_no3( Mode, SRC , CAL , STP , Int_NO3, Int_BSL )
%
% Mode:  string: OPUS
%
% SRC: Raw DataSource OPUS Light and Dark Intensity Measurements
%
% CAL: Calibration file: wavelength, eno3, ebr, cal T, reference intensity
%
% STP timeseries: [ MTime Sal Temp [Press] ]
% STP single sample: Sal or [ Sal Temp ] or [ Sal Temp Press ]
%     defaults: Sal = 35, Temp = InstrumentTemp or 20, Press = 0,
%
% STP can be a Structure with fields 'stp' and 'coeff'
%     STP.stp = Sal or [ Sal Temp ] or [ Sal Temp Press ] or [ MTime Sal Temp [Press] ]
%     STP.coeff = CoefficientStructure for fit of Seawaterabsorption versus Temperature in CALC_NO3
%     STP.coeff.wl  Reference Wavelength
%     STP.coeff.cf|cp1|ce1   see CALC_NO3
%
% Int_NO3: Wavelength Interval for NO3-fit
% Int_BSL: Wavelength Interval for Baseline(CDOM)-fit
%
% Processing Steps: 
%
%   (1) Check for CalStructure or CalFile, 
%          read CalFile
%
%   (2) Check for DataStructure (Raw or Pre) or DataFileList 
%         
%
%   (3) Basic check of CTD-Data "STP" 
%         valid: [Sal] or [Sal Temp] or [ Sal Temp Press] or [ MTime Sal Temp Press ]
%
%   [4] Sort Merged Data by Time, DarkRecords before LightRecords (at same Timestamp)
%
% Save RawMatFile here
%
%   [5] Find Raw-Dark values corresponding to Raw-Light samples, extract Light samples only
%       Correct Light samples with Dark values, 'raw_intens' --> 'intens'
%
%   (6) Match CTD-records "STP" with Light samples, clip data to time range of CTD-data
%
%
%   [7] Bromide Correction, Baseline(CDOM)-Fit and no3-Fit with
%   proc_calc_no3

Nin = nargin;
vin = varargin;

dat = struct;
dat = dat([]);

if Nin < 2
   msg = 'Not enough Input Arguments';
   return
elseif ~chkstr(mode,1)
   msg = 'Mode must be a non-empty strimng: OPUS.';
   return
end
   
mode = lower(mode);


% Functions to read data by File-extension
% OPUS: DAT or CSV

mergerep = '@EXT@';  
                                  
mergefcn = sprintf('rd_%s_%s',mode,mergerep);


calfcn = sprintf('rd_%s_cal',mode);


opus_wvl = linspace(190,394,256);

% Seconds to Match RawDark with RawLight

acc = [1 3 10];


sal  = 35; % default Salinity
temp = 20; % default Temperature
press = 0;

wrg = [ 180 400 ];  % Default Interval (Range) for Wavelengths

%-----------------------------------------------------------------------

msg = cell(0,1);
dat = struct;
dat = dat([]);

%-----------------------------------------------------------------------
% Basic Check of Inputs

if Nin < 3, cfile = ''; end

Nin = Nin - 3;

% Check for StringSelection-Input for MergeData
str = {};
if Nin > 0
   if iscell(vin{1})
      str = vin{1};
      if iscellstr(str), str = {str}; end
      vin = vin(2:Nin);
      Nin = Nin - 1;
   end
end

% Check for RawMatfile-Input to save RawData
rawfile = -1;
if Nin > 0
   if chkstr(vin{1},0)
      rawfile = vin{1};
      vin = vin(2:Nin);
      Nin = Nin - 1;
      if ~isempty(rawfile)
          [p,n,e] = fileparts(rawfile);
          if ~strcmp(lower(e),'.mat')
              rawfile = fullfile(p,[n e '.mat']);
          end
      end
   end
end

% Check for PreMatfile-Input to save PreData
prefile = -1;
if Nin > 2
   if chkstr(vin{3},0)
      prefile = vin{3};
      vin = vin([1 2 (4:Nin)]);
      Nin = Nin - 1;
      if ~isempty(prefile)
          [p,n,e] = fileparts(prefile);
          if ~strcmp(lower(e),'.mat')
              prefile = fullfile(p,[n e '.mat']);
          end
      end
   end
end

if Nin < 1; stp = []; else, stp = vin{1}; end
if Nin < 2; int = []; else, int = vin{2}; end
if Nin < 3; int_bsl = []; else, int_bsl = vin{3}; end
if Nin < 4; int_no3 = []; else, int_no3 = vin{4}; end

auto_int = isempty(int);

%-----------------------------------------------------------------------
   
if ~( isempty(int_no3) | ( isnumeric(int_no3) & isequal(size(int_no3,2),[2]) ) )
    msg = cat(1,msg,{'Interval for no3-Fit must be a 2-element numeric.'});
end

if ~( isempty(int_bsl) | ( isnumeric(int_bsl) & isequal(size(int_bsl),[1 2]) ) )
    msg = cat(1,msg,{'Interval for Baseline-Fit must be a 2-element numeric.'});
end


%***********************************************************************
% (1) Check CalFile or CalData
%***********************************************************************

m = '';

[mm,def] = feval(calfcn);

defld = fieldnames(def);

if isempty(cfile)

   cal = def;

elseif isstruct(cfile) & ( prod(size(cfile)) == 1 )

   if isequal(fieldnames(cfile),defld)
      cal   = cfile;
      cfile = 'InputStructure';
   else
      m = sprintf('CalStructure if Input must be created by "%s".',upper(calfcn));
   end

elseif ~chkstr(cfile,1)

    m = 'Input for CalFile must be a string or single element structure.';

else

   [m,cal] = feval(calfcn,cfile);
   
   if ~isempty(m)
       m = sprintf('Error call %s to read CalFile "%s".\n%s', ...
                       upper(calfcn),cfile,m);
   end

end  % isempty(cfile)


cal_ok = isempty(m);
if cal_ok
   cal_ok = ( isstruct(cal) & ( prod(size(cal)) == 1 ) );
   if cal_ok
      cal_ok = ~isempty(cal.wvl);
   end
end

if ~isempty(m)
    msg = cat(1,msg,{m});
end

%***********************************************************************
% (2) Check first Data for MAT-File or Structure
%***********************************************************************

is_mat = chkstr(dfile,1);
if is_mat
   [p,n,e] = fileparts(dfile);
   if isempty(e)
      f = fullfile(p,[n '.mat']);
      if exist(f,'file') == 2
         file = f;
         e = '.mat';
      end
   end   
   is_mat = strcmp(lower(e),'.mat');
end

chk_struct = ( is_mat | ( isstruct(dfile) & ( prod(size(dfile)) == 1 ) ) );

dtext = dfile;
if ~chkstr(dtext,1)
    dtext = 'FileList';
end

is_raw = 1;
is_pre = 0;

% Data from Input-Structure or MAT-File must have following Fields

raw_fld = { 'mtime' 'wvl' 'raw_intens' 'light_flag' 'drk_avg' 'cal' };
pre_fld = { 'mtime' 'wvl' 'intens' 'stp' 'cal' };

m = '';

if chk_struct

   if is_mat

      try
        dat = load('-mat',dfile);
      catch
          m = sprintf('Error load MAT-DataFile "%s".\n%s',dtext,lasterr);
      end

   else

      dat   = dfile;
      dtext = 'Structure';

   end

   if isempty(m)

      f = fieldnames(dat);

      for ff = raw_fld
          is_raw = any(strcmp(f,ff{1}));
          if ~is_raw
              break
          end
      end

      if ~is_raw
          for ff = pre_fld
              is_pre = any(strcmp(f,ff{1}));
              if ~is_pre
                  break
              end
          end
      end

      if ~( is_raw | is_pre )
          rf = sprintf('%s , ',raw_fld{:});
          pf = sprintf('%s , ',pre_fld{:});

          m = sprintf('Invalid Data in "%s", missing fieldnames:\n RawData: %s\n or\nPreData: %s', ...
                         dtext,rf(1:(end-3)),pf(1:(end-3)));
      end

   end

   if isempty(m) & isfield(dat,'cal')

      if ~isequal(fieldnames(cal),fieldnames(dat.cal))
          m = sprintf('CalStructure in Data of "%s" must be created by %s.', ...
                       dtext,upper(calfcn));
      elseif ~cal_ok
          cal = dat.cal;
          cal_ok = ~isempty(cal.wvl);
      end

   end

else

   [ok,c] = chkcstr(dfile,0);
   if ok
      [p,n,e] = fileparts(lower(c{1}));
      if ~isempty(e) & ( size(e,2) > 1 )
          mergefcn = strrep( mergefcn , mergerep , e(2:end) );
          if ~exist(mergefcn,'file') == 2;
              m = sprintf('Can''t find MergeFcn "%s".',mergefcn);
          end
      else
          m = sprintf('Can''t identify MergeFcn by Extension of "%s".',c{1})
      end
   else
      m = 'Input for FileList must be a String or Cellarray of Strings.';
   end

   if isempty(m)

      mergefcn = {mergefcn};


      [m,dat] = mergedata(mergefcn,dfile,str);
      if ~isempty(m)
          m = sprintf('Error call MERGEDATA(%s) for List "%s".\n%s', ...
                      upper(mergefcn{1}),dtext,m);
      end

   end

   if isempty(m) & ~cal_ok
      cal_ok = isfield(dat,'cal');
      if cal_ok
         cal = dat.cal;
         cal_ok = ( isstruct(cal) & ( prod(size(cal)) == 1 ) );
         if cal_ok
            cal_ok = isequal(fieldnames(dat.cal),defld);
            if cal_ok
               dcal = dat.cal;
               cal_ok = ~isempty(dat.cal.wvl);
            else
      m = sprintf('Invalid CalStructure returned by MERGEDATA(%s) for List "%s".', ...
                              upper(mergefcn{1}),dtext);
            end
         end
      end
   end

end

if isempty(m)
   if size(dat.wvl,1) > 1
      if ~( sum(sum(abs(diff(dat.wvl,1,1)),1),2) == 0 )
           m = sprintf('WaveLengths in Data of "%s" are different.\n%s',dtext);
      else
           dat.wvl = dat.wvl(1,:);
      end
   end
end

if ~isempty(m)
    msg = cat(1,msg,{m});
end

if ~cal_ok & isempty(msg)

     cal = def;

      nn = NaN * ones(1,size(dat.wvl,2));
 
      cal.wvl  = dat.wvl;
      cal.eno3 = nn;
      cal.esw  = nn;
      cal.ref  = nn;
      cal.esw_temp = NaN;

end

%***********************************************************************
% (3) Basic check of CTD-Data
%***********************************************************************

m = '';

coeff = [];    % Coefficients for SW-Absorption determined in calc_no3

if isstruct(stp) 
   if ~( ( prod(size(stp)) == 1 ) & ...
         ( any(strcmp(fieldnames(stp),'coeff')) | any(strcmp(fieldnames(stp),'stp')) ) )
       m = 'StructureInput for STP must be a single element Structure with fields "coeff" and "stp".';
   else
      if isfield(stp,'coeff')
         cf = getfield(stp,'coeff');
         if isstruct(cf) & ( prod(size(cf)) == 1 )
            coeff = cf;
         elseif ~isempty(cf)
            m = 'Value for field "coeff" in STP-Structure must be a single element Structure, see CALC_NO3';
         end
      end
      if isfield(stp,'stp')
         stp = getfield(stp,'stp');
      else
         stp = [];
      end
   end
end

if isempty(m) & ( prod(size(stp)) > 1 )

   ncol = [ 2  3 ] + ( size(stp,1) > 1 );
   if ~( isnumeric(stp) & any( size(stp,2) == ncol ) & ...
         ( prod(size(stp)) == size(stp,1)*size(stp,2) ) )
       m = sprintf('Invalid STP-Data, Numeric with %.0f or %.0f Columns required.',ncol);
   end

   if isempty(m) & ( size(stp,2) == ncol(1) )
       % No pressure-Column, use default Pressure
       stp = cat( 2 , stp , press+0*stp(:,1) );
   end

end

if isempty(m) & ~isempty(stp)
   dat(1).stp = stp;
end

if isstruct(dat) & ~isempty(dat)
   dat(1).coeff = coeff;
end

if ~isempty(m)
    msg = cat(1,msg,{m});
end

%-----------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if nargout == 0
       error(msg)
    end
    return
end

%***********************************************************************
% [4] Sort Data by Time
%***********************************************************************

scl = 24 * 3600 * 100;  % Scale Time to Full Hundreds of Seconds

m0 = round( min(dat.mtime) * scl );
mt = round( ( dat.mtime - m0/scl ) * scl );

if ~all( diff(mt) >= 0 )

    [h,si] = sort( mt + dat.light_flag/10 , 1 , 'ascend');  % Sort Dark before Light
  
    nn = size(mt,1);
    nw = size(dat.wvl,2);

    ff = fieldnames(dat);
    for f = ff(:)'
        v = getfield(dat,f{1});
        if ( isnumeric(v) | islogical(v) | iscell(v) ) 
           if ( size(v,1) == nn ) & any( size(v,2) == [ 1  nw ] )
              dat = setfield(dat,f{1},v(si,:));
           end
        end
    end

    dat.mtime = ( mt(si) + m0 ) / scl ;

end
    
%***********************************************************************
% Save Raw-Data
%***********************************************************************

if ~isequal(rawfile,-1)

    dat.cal = cal;

    if isempty(rawfile) & ( isempty(prefile) | isequal(prefile,-1) )
       fprintf(1,'Return RAW-Data only. Further Inputs following STP are ignored.\n')
       return
    end

    fprintf(1,'Save RawData to "%s" ... ',rawfile);

    try
       save(rawfile,'-mat','-struct','dat');
       fprintf(1,'done\n');
    catch
       msg = sprintf('Can''t save RawData to "%s".\n%s',rawfile,lasterr);
       fprintf(1,'error\n');
    end

    if ~isempty(msg)
        if ( nargout == 0 )
           error(msg)
        end
        return
    end

end

%***********************************************************************
% (5) Find Wavelengths of Data and Calibration within Interval
%***********************************************************************

if isempty(int)
   int = wrg;
   int(1) = max(int(1),min(dat.wvl));
   int(1) = max(int(1),min(cal.wvl));
   int(2) = min(int(2),max(dat.wvl));
   int(2) = min(int(2),max(cal.wvl));
end

ii_dat = ( ( int(1) <= dat.wvl(1,:) ) & ( dat.wvl(1,:) <= int(2) ) );

if ~any(ii_dat)
    msg = cat(1,msg,{'Can''t find WavelLengths for Data within Interval'});
elseif all(ii_dat)
    ii_dat = [];
else
    ii_dat = find(ii_dat);
end

ii_cal = ( ( int(1) <= cal.wvl ) & ( cal.wvl <= int(2) ) );

if ~any(ii_cal)
    msg = cat(1,msg,{'Can''t find WavelLengths for Calibration within Interval'});
elseif all(ii_cal)
    ii_cal = [];
else
    ii_cal = find(ii_cal);
end


if isempty(msg)

   dw = dat.wvl; if ~isempty(ii_dat), dw = dw(ii_dat); end
   cw = cal.wvl; if ~isempty(ii_cal), cw = cw(ii_cal); end

   ok = ( prod(size(dw)) == prod(size(cw)) );

   if ok
      ok = all( abs( dw - cw ) <= 1e-3 );
   end

   if ~ok
       msg = cat(1,msg,{'WavelLengths of Data and Calibration don''t match within Interval.'});
   end

end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if nargout == 0
       error(msg)
    end
    return
end


%***********************************************************************
if is_raw
%***********************************************************************
% (6) Get RawDarkValues corresponding to RawLightMeasurements
% (7) Reduce Raw-Data to LightMeasurements and Wavelength Interval
%***********************************************************************

   is_lgt = ( dat.light_flag == 1 );

   if ~any(is_lgt)
       msg = 'No Light Records in Data.';
   elseif all(is_lgt)
       msg = 'No Dark Records in Data.';
   end

   if ~isempty(msg)
       return
   end

   mt = dat.mtime;
   sc = ( mt - mt(1) ) * 24 * 3600;           % Seconds since first record

   ii_lgt = find( is_lgt); nl = sum( is_lgt); % Index and Number of LightRecords
   ii_drk = find(~is_lgt); nd = sum(~is_lgt); % Index and Number of DarkRecords

   % TimeVector for Light- and Dark-Records
   lgt_sc = sc(ii_lgt); lgt_sc = lgt_sc(:); lgt_sc = permute(lgt_sc,[2 1]);
   drk_sc = sc(ii_drk); drk_sc = drk_sc(:);

   % Matrice of TimeDifferences, ascending over Columns
   dt = drk_sc(:,ones(1,nl)) - lgt_sc(ones(nd,1),:); 

   % Row-Index at which Difference is ZERO+acc or close to ZERO
   i1 = sum( dt <= 0+acc(1) , 1 );

   i1 = i1 + ( i1 == 0 );
   
   ii = i1 + size(dt,1) * ( 0 : (size(dt,2)-1) );
  

   jj = ( abs(dt(ii)) > acc(2) );

   if any(jj)
      jj = find(jj);
      [dm,i1(jj)] = min(abs(dt(:,jj)),[],1);
      if any( dm > acc(3) )
         warning('Large time differences for Dark vs. Light')
      end
   end

   raw_drk_avg = dat.drk_avg(ii_drk(i1));
   raw_drk_mt  = mt(ii_drk(i1));

   %-----------------------------------------------------------------------
   % Reduce Structure to LightRecords at selected Wavelengths

   nn = size(mt,1);
   nw = size(dat.wvl,2);

   ff = fieldnames(dat);
    
   for f = ff(:)'
       v = getfield(dat,f{1});
       if ( isnumeric(v) | islogical(v) | iscell(v) ) 
          sz = size(v);
          if ( sz(1) == nn ) & any( sz(2) == [ 1  nw ] )
             v = v(ii_lgt,:);
          end
          if ~isempty(ii_dat) & ( sz(2) == nw ) & any( sz(1) == [ 1  nn ] )
             v = v(:,ii_dat);
          end
          dat = setfield(dat,f{1},v);
       end
   end

   %-----------------------------------------------------------------------
   % Clean up, Add Dark-Corrected Intensity, Calibration and CTD-Data

   dat = rmfield(dat,'drk_avg');

   nr = size(dat.wvl,2);  % Reduced Numbre of WaveLengths

   dat.raw_drk_avg = raw_drk_avg;
   dat.raw_drk_mt  = raw_drk_mt;

   dat.intens = dat.raw_intens - raw_drk_avg(:,ones(1,nr));

%-----------------------------------------------------------------------
elseif ~isempty(ii_dat)
%-----------------------------------------------------------------------
% (7) Reduce PreData to Wavelength-Interval
%-----------------------------------------------------------------------

   nn = size(mt,1);
   nw = size(dat.wvl,2);

   ff = fieldnames(dat);
    
   for f = ff(:)'
       v = getfield(dat,f{1});
       if ( isnumeric(v) | islogical(v) | iscell(v) ) 
          sz = size(v);
          if ( sz(2) == nw ) & any( sz(1) == [ 1  nn ] )
             dat = setfield(dat,f{1},v(:,ii_dat));
          end
       end
   end


end % is_raw

%***********************************************************************
if ~isempty(ii_cal)
%***********************************************************************
% (8) Reduce CalData to Wavelength-Interval
%***********************************************************************

   nw = size(cal.wvl,2);

   ff = fieldnames(cal);
    
   for f = ff(:)'
       v = getfield(cal,f{1});
       if ( isnumeric(v) | islogical(v) | iscell(v) ) 
          sz = size(v);
          if ( sz(2) == nw )
             cal = setfield(cal,f{1},v(:,ii_cal));
          end
       end
   end


end 


dat.cal = cal;


%***********************************************************************
% (9) Match CTD-Data with Records
%***********************************************************************

mt = dat.mtime;
nn = size(mt,1);
nw = size(dat.wvl,2);

if isempty(stp)
   if isfield(dat,'stp')
      stp = dat.stp;
   end
end

if prod(size(stp)) == 1
   sal = stp;
   stp = [];
end

if isempty(stp)

   switch mode
       case 'opus'
             fld = 'Temperature';
   end

   if isfield(dat,fld)
      tp = getfield(dat,fld);
   else
      warning(sprintf('Can''t find Temperature Field "%s" in %s-Data.',fld,upper(mode)));
      tp = temp;
   end

   stp = cat( 2 , sal+0*tp , tp , press+0*tp );

   if size(stp,1) > 1
      stp = cat( 2 , mt , stp );
   end
   
end


if size(stp,1) == 1

   stp = stp( ones(nn,1) , : );

elseif ~isequal(size(stp),[nn 3]) 

   ct  = stp(:,1);     % CTD-Time
   stp = stp(:,2:end); % [ S T P ]

   % Check if CTD-Data have same TimeVector
   ok = ( size(ct,1) == nn );
   if ok
      ok = all( abs( ct - mt ) < 1e-6 );
   end

   if ~ok
 
       isn = any(isnan([ct stp]),2);
       if any(isn)
          warning('STP Data contain NaN-records.')
          ii  = find(~isn);
          ct  = ct(ii);
          stp = stp(ii);
       end

       dt = round( (ct-ct(1)) * 24 * 3600 * 100 );  

       dt = diff(dt,1,1);

       if ~( all( dt > 0 ) | all( dt < 0 ) )
           warning('STP Data should be monotonic in time.')
            [ct,ix] = uniqued(ct,1+i);
                stp = stp(ix,:);
       elseif ct(1) > ct(2)
           ct  =  ct(end:-1:1);
           stp = stp(end:-1:1,:);
       end

     
       % Expand Start- and End-Time of CTD-Data by 1/1000 sec
       ct([1 end]) = ct([1 end]) + 1e-3/24/3600 * [ -1 ; 1 ];

       ii = ( ( ct(1) <= mt ) & ( mt <= ct(end) ) );

       if ~any(ii)

          msg = 'No Spectra Data found within Range of STP-Data.';
          return

       elseif ~all(ii)

          warning('Clip Spectra Data to range of CTD-Data')

          ii = find(ii);

          ff = fieldnames(dat);

          for f = ff(:)'
              v = getfield(dat,f{1});
              if ( isnumeric(v) | islogical(v) | iscell(v) ) 
                 if ( size(v,1) == nn ) & any( size(v,2) == [ 1  nw ] )
                    dat = setfield(dat,f{1},v(ii,:));
                 end
              end
          end

          mt = dat.mtime;

       end

       stp = interp1(ct,stp,mt);
   
   end
end

dat.stp = stp;

%***********************************************************************
% Save Pre-Data
%***********************************************************************

if ~isequal(prefile,-1)

    if isempty(prefile)
       fprintf(1,'Return Pre-Processed Data only.\n')
       return
    end

    fprintf(1,'Save PreData to "%s" ... ',prefile);

    try
       save(prefile,'-mat','-struct','dat');
       fprintf(1,'done\n');
    catch
       msg = sprintf('Can''t save PreData to "%s".\n%s',prefile,lasterr);
       fprintf(1,'error\n');
    end

    if ~isempty(msg)
        if nargout == 0
           error(msg)
        end
        return
    end

end

if isempty(cfile) & isnan(cal.esw_temp)
   return
end

%***********************************************************************
% (10) Calculate Absorbance, Baseline- and NO3-Fit
%***********************************************************************

coeff = dat.coeff;

if isstruct(coeff) & ( prod(size(coeff)) == 1 )
    eswt = coeff;
    eswt.esw_temp = cal.esw_temp;
else
    eswt = cal.esw_temp;
end
    
dat.fit = opus_calc_no3( dat.stp, dat.wvl, dat.intens , ...
               cal.ref, cal.esw_temp, cal.esw, cal.eno3 , ...
                    int_no3 , int_bsl );


%***********************************************************************
