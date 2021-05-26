function [out,stp] = opus_calc_no3(stp,wvl,intens,ref,esw_temp,esw,eno3,int_no3,int_bsl,scl,M)

% OPUS_CALC_NO3  Calculate in-situ Nitrate Concentration from spectrophotometric measurements
%
% STP       [ N by 3 ] CTD-Data: [ Sal Temp Press ]  
% WVL       [ 1 by M ] WaveLengths                   
% Intens    [ N by M ] In-situ Intensity Measurements, Dark-corrected
%
% Ref       [ 1 by M ] Reference Intensities measured in deionized water
%                        Dark-corrected, also called Waterbase spectra
%
% ebr_temp  [ 1 by 1 ] CalTemp [degC] for ebr
% ebr       [ 1 by M ] Sea Water Molar Extinction Coefficients from Calibration
% Eno3      [ 1 by M ] Nitrate Molar Extinction Coefficients from Calibration
% Eno3      [ K by M ] Multiple Molar Extinction Coefficients from Calibration
%
%
% Int_NO3   [ MinWvl MaxWvl ] WaveLength-Intervall to determine no3 
%                              by linear Regression, default: [ 217  240 ]
%Int_NO3 : interval (wavelength range) for the NO3 fit
% Int_BSL  [ MinWvl MaxWvl ] WaveLength-Intervall to determine CDOM Offset Correction
%                              by linear Regression, default: empty, See note below !
%
%int_bsl : interval (wavelength range) for the CDOM-related baseline fit
% linear regression between absorbance and wavelength [240 260]
% Output Structure:
% 
%     ebr_in_situ: [ N by M ]  in-situ Salt Water Extinction, Temperature and Pressure corrected
%     abs_in_situ: [ N by M ]  in-situ Absorbance (Beer Lambert Law)
%         abs_cor: [ N by M ]  Sea Salt Corrected Absorbance: abs_in_situ - ebr_in_situ * Salitinty
%         int_no3: [ 1 by 2 ]  Interval for Wavelengths for NO3-fit
%         ind_no3: [ 1 by O ]  Index Vector for Wavelengths in Intervall of Int_NO3
%         wvl_no3: [ 1 by O ]  Wavelengths in Intervall of Int_NO3
%         int_bsl: [ 1 by 2 ]  Interval for Wavelengths for separate BaselineFit
%         ind_bsl: [ 1 by P ]  Index Vector for Wavelengths in Intervall of Int_BSL
%         wvl_bsl: [ 1 by P ]  Wavelengths in Intervall of Int_BSL
%             bsl: [ N by 2 ]  [ Offset Slope ] Coefficients for Baseline-Fit (CDOM Offset Correction)
%             no3: [ N by K ]  no3 value, linear Regression Coefficient for Eno3
%     abs_bsl_fit: [ N by M ]  Baseline-Fit: bsl(1) + bsl(2) * WVL
%     abs_bsl_cor: [ N by M ]  Baseline Corrected abs_in_situ: abs_in_situ - abs_bsl_fit
%     abs_no3_fit: [ N by M by K ]  NOX-Fit: no3 * Eno3
%    abs_diff_fit: [ N by M ]  Difference of abs_cor - abs_bsl_fit - sum(abs_no3_fit,3)
%    abs_diff_rms: [ N by 1 ]  Root Mean Square of Difference
%
% Note: An empty Value for INT_BSL performs the Baseline and NO3 Regression in a single step,
%        i.e. BaselineOffset, BaselineSlope and NO3-Slope in a single Matrice Fit,
%       applied to the Sea Salt Corrected Absorbance (abs_cor) within the Intervall of Int_NO3
%
%       A nonempty Value for INT_BSL applies the Baseline Regression (Offset and Slope) first to 
%        the Corrected Absorbance (abs_cor) within the Intervall of Int_BSL, and in a second step 
%        the no3-fit (Slope) to the Corrected Absorbance reduced by the Baseline Fit 
%       (abs_cor - abs_bsl_fit) within the Intervall of Int_NO3.
%        
% Note: Empty CTD data will perform a best fit (FMINSEARCH) for CTD-Values, 
%        a second output variables retrieves these CTD-values: 
%       [ Out , STP ] = CALC_NO3( [] , ... )
%
%-------------------------------------------------------------------------------
%
% Further Readings: Nehir et al. 2021
%                   Sakamoto et al., 2009 and 2017
%                   Zielinski et al., 2011
%                   Fommervault et al., 2015
%
%-------------------------------------------------------------------------------
%
% Initial Sources: Ken Johnson, MBARI, ISUS processing
% 
% Ocean Observatories Initiative CI at GitHub: https://github.com/
% ooici/ion-functions/blob/master/ion_functions/data/matlab_scripts/nutnr/NUTNR_Example_MATLAB_code_20140521_ver_1_00.m
%
%-------------------------------------------------------------------------------

Nin = nargin;

if Nin < 8, int_no3 = []; end
if Nin < 9, int_bsl = []; end
if Nin < 10, scl = [ 100  1000  1  1  ]; end

rms_out = ( Nin > 10 );
 
if isempty(int_no3), int_no3 = [ 217 240 ]; end

excl_int = zeros(0,2);

if size(int_no3,1) >= 2
   excl_int = int_no3(2:end,:);
   int_no3 = int_no3(1,:);
end

   
n1 = size(intens,1);
n2 = size(intens,2);

o1 = ones(n1,1);
o2 = ones(1,n2);

%=======================================================================
if rms_out
%=======================================================================

% Special Input from FMINSEARCH with M as last Input, see further below
% Int-Parameter are here the Index Vectors

   ii_no3 = int_no3;
   ii_bsl = int_bsl;

%=======================================================================
else
%=======================================================================
% Get Wavelengths for no3 and CDOM-Interval
% Form Matrice for linear fits

   if ~( isnumeric(int_no3) & isequal(size(int_no3),[1 2]) )
       error('Intervall for NOX-Fit must be a 2-element numeric.');
   end

   ii_no3 = ( ( int_no3(1) <= wvl ) & ( wvl <= int_no3(2) ) );

   for jj = 1 : size(excl_int,1)
       ii_no3 = ( ii_no3 & ~( ( excl_int(jj,1) < wvl ) & ( wvl < excl_int(jj,2) ) ) );
   end

   ii_no3 = find(ii_no3);

   if isempty(ii_no3)
      error('Can''t find Wavelengths for Data within Intervall for NOx')
   end

   ii_bsl = [];

   if ~isempty(int_bsl)

       if ~( isnumeric(int_bsl) & isequal(size(int_bsl),[1 2]) )
           error('Intervall for Baseline-Fit must be a 2-element numeric.');
       end

       ii_bsl = find( ( int_bsl(1) <= wvl ) & ( wvl <= int_bsl(2) ) );

       if isempty(int_bsl)
          error('Can''t find WavelLengths for Data within Intervall for CDOM')
       end

       M1 = cat( 1 , o2(ii_bsl) , wvl(ii_bsl) );
       M1 = pinv(permute(M1,[2 1]));
       M1 = permute(M1,[2 1]);         % 2 Columns   

       M2 = eno3(:,ii_no3);
       M2 = pinv(permute(M2,[2 1]));
       M2 = permute(M2,[2 1]);         % 1 or 2 Column

       % Length of ii_no3 and ii_bsl may differ
       % Fill Matrices with NaN's

       l1 = prod(size(M1,1));
       l2 = prod(size(M2,1));

       M1 = cat(1,M1,NaN*ones(max(0,l2-l1),size(M1,2)));
       M2 = cat(1,M2,NaN*ones(max(0,l1-l2),size(M2,2)));

       M = cat( 2 , M1 , M2 ); % 1. and 2. Col: Baseline-fit (CDOM); 3. Col: no3-fit

   else

       M = cat( 1 , o2(ii_no3)/scl(1) , wvl(ii_no3)/scl(2) , eno3(:,ii_no3)/scl(3) );
       M = pinv(permute(M,[2 1]));
       M = permute(M,[2 1]);          % 3 Columns 

   end

%=======================================================================
end
%=======================================================================

%-----------------------------------------------------------------------
% Start FMINSEARCH to estimate CTD-Values

if isempty(stp) & ( rms_out == 0 )

   stp_start = [35,20,0];
   opt = optimset('display','off','maxiter',1000);

   stp = NaN * ones(n1,3);

   sc = []; cl = [];
   res = 5; pct = 0;

   if exist('loopdot','file') == 2
      sc = [ 50  5+i ];
      cl = loopdot(sc,n1,'Fit STP records');
   else
      frm = '\rFit STP records %.0f times, please be patient ... %3d%%';
      fprintf(1,['\n' frm],n1,pct);
   end

   for ii = 1 : n1

       stp(ii,:) = fminsearch(@(x) calc_no3(...
             x,wvl,intens(ii,:),ref,esw_temp,esw,eno3,ii_no3,ii_bsl,scl,M), ...
                stp_start, opt );

       if ~isempty(cl)
           loopdot(sc,n1,ii,1,cl);
       else
           p = res * floor( 100 * ii/n1 / res );
           if p > pct
              pct = p;
              fprintf(1,frm,n1,pct);
           end
       end

   end

   if isempty(cl)
      fprintf(1,[frm ' %s\n\n'],n1,pct,'done');
   end

end

%-----------------------------------------------------------------------
% Absorbance at esw-Calibration Temp. and insitu-Temp
% Coefficients from Nehir et al. 2021
A =  1e-7; 
B = -9e-5;
C =  1.6e-3;
D =  2.12e-2;

W =  210; % wavelength offset adjustable parameter within 206 and 212 nm
P =  0.026;   % Pressure Factor from Sakamoto et al. 2017
%P =  0.020;   % Pressure Factor from Fommervault et al. 2015

cf = [ A  B  C  D ];


%-----------------------------------------------------------------------
% Error Check on esw_Temp
if ~( isnumeric(esw_temp) & ( prod(size(esw_temp)) == 1 ) )

   error('Input for esw_temp must be a single numeric or Structure.')

end
  
%-----------------------------------------------------------------------
% Exp-WaveLength-Polynomal Fit

wvr = wvl(1,:) - W;
            
% Correct esw for insitu Temperature 

esw_in_situ = esw(o1,:) .* exp((cf(1) .* wvr.^3 + cf(2)  .* wvr.^2 + cf(3) .* wvr + cf(4) ) .*(stp(:,2*o2)-esw_temp));

% Correct esw for insitu Pressure

esw_in_situ = esw_in_situ .* ( 1 - stp(:,3*o2)./1000 .* P );


%-----------------------------------------------------------------------
% Absorbance 

abs_in_situ                       = intens ./ ref(o1,:);
abs_in_situ(find(abs_in_situ<=0)) = NaN;

abs_in_situ  = -1 * log( abs_in_situ ) / log(10);

% Remove excpected Spectral Component due to sea-salt (bromide)
%  from Measured Absorbance at in-situ Temp.

abs_cor = abs_in_situ - esw_in_situ .* stp(:,1*o2);

%-----------------------------------------------------------------------
% Linear Regressions

if isempty(ii_bsl)
  
   % BaselineCorrection and no3-fit in one Step on same Interval

   c = abs_cor(:,ii_no3) * M;  % Linear Fit

   ii_scl = [ 1  2  3*ones(1,size(c,2)-2) ];

   c = c ./ scl(o1,ii_scl);

   bsl = c(:,[1 2]);     % Offset Slope

   no3 = c(:,3:end);     % Slope only

else

   % CDOM Baseline Correction with Measurements from 240 .. 260nm

   nm = prod(size(ii_bsl));

   bsl = abs_cor(:,ii_bsl) * M(1:nm,[1 2]);  % Linear Fit

   no3 = [];

end


% Baseline Absorbance (CDOM)
abs_bsl_fit = bsl(:,1*o2) + bsl(:,2*o2) .* wvl(o1,:);

% Nitrate Absorbance of Observation
abs_no3 = abs_cor - abs_bsl_fit;

if isempty(no3)

   % no3-fit second

   nm = prod(size(ii_no3));

   no3 = abs_no3(:,ii_no3) * M(1:nm,3:end);  % Linear Fit

end

% Expected NO3 Absorbance from Fit

nx = size(eno3,1);

abs_no3_fit = zeros(n1,n2,nx);

for ii = 1 : nx

    abs_no3_fit(:,:,ii) = no3(:,ii*o2) .* eno3(ii*o1,:);

end

abs_diff_fit = abs_cor - abs_bsl_fit - sum(abs_no3_fit,3);

%---------------------------------------------------------

dft = abs_diff_fit(:,ii_no3);
isn = isnan(dft);

dft(find(isn)) = 0;

diff_rms = sqrt( sum(dft.^2,2) ./ (size(dft,2)-sum(isn,2)) );

%---------------------------------------------------------

if isequal(rms_out,1)
   out = diff_rms;
   return
end

out = struct( 'esw_in_situ' , { esw_in_situ } , ...
              'abs_in_situ' , { abs_in_situ } , ...
              'abs_cor'     , { abs_cor }     , ...
              'int_no3'     , { int_no3 }     , ...
              'ind_no3'     , {  ii_no3 }     , ...
              'wvl_no3'     , { wvl(ii_no3) } , ...
              'int_bsl'     , { int_bsl }     , ...
              'ind_bsl'     , {  ii_bsl }     , ...
              'wvl_bsl'     , { wvl(ii_bsl) } , ...
                  'bsl'     , { bsl }         , ...
                  'no3'     , { no3 }         , ...
              'abs_bsl_fit' , { abs_bsl_fit } , ...
              'abs_bsl_cor' , { abs_in_situ-abs_bsl_fit } , ...
              'abs_no3_fit' , { abs_no3_fit } , ...
             'abs_diff_fit' , { abs_diff_fit } , ...
             'abs_diff_rms' , { diff_rms     } , ...
                   'coeff'  , { cf       } , ...
                 'esw_temp' , { esw_temp } , ...
                 'wvl_offs' , { W        }           );

