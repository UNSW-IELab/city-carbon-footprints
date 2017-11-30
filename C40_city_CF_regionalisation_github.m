%% C40_city_CF_regionalisation_github.m

% This script accompanies the paper "Consumption-Based Carbon Footprints
% of Global Cities" by Wiedmann et al, submitted to Nature Climate Change
% in November 2017.

% The script creates regionalised two-region IO tables for single cities
% and rest-of-nations in GTAP format. It then inserts these two-region IO
% tables into GTAP9 (replacing the nation's 'cross') and calculates carbon
% footprints of city's final demand.

% Data to run the script can be made available upon request.
% Corresponding author email: t.wiedmann@unsw.edu.au


tic; clc; clear;
diary on

%% Define paths 

modelpath    = ('/.../');
citydatapath = ('/.../');
resultspath  = ('/.../');
GTAPdatapath = ('/.../');


% loops for LQ and RAS methods
for LQ = 5:6 ;
for RAS_choice = 0:4 ;


%% Load and prepare data

% disp('.');
% disp('Load GTAP data');
% disp('.');
% fileID = fopen([resultspath 'log.txt']);
% fprintf(fileID,'Load GTAP data');
% fclose(fileID);
% T = evalc(expression)

load([GTAPdatapath 'A.mat']) ;
load([GTAPdatapath 'X.mat']) ;
load([GTAPdatapath 'Y.mat']) ; % This contains five matrices Y_h, Y_g, Y_k, Y_VST and Y_total; all of dimension 7980x140 (one column of total FD for each country) 
load([GTAPdatapath 'VA.mat']) ;
load([GTAPdatapath 'GHGs_2011.mat']) ;
load([GTAPdatapath 'CO2_GTAP_HH_2011.mat']) ;
load([GTAPdatapath 'f.mat']) ; % Vector of CO2e emissions. f=sum(CO2_GTAP_sec,1)+(sum(GHG,1)*1000);


% clear Y_block ;  % Y_block not used further
T = A .* repmat(X',7980,1) ;  % transaction matrix; unit US$m 

% Explanation of GHG emissions
% 
% CO2_GTAP_sec, 5x7980 
% The GTAP file "emissions from domestic production" includes the CO2
% emissions by industry for all countries. The five rows represent the CO2
% from Coal, Oil, Gas, P_C = petroleum and coal products, and gdt = gas
% distribution. 
% Unit is kt CO2e
% 
% Here we sum up the rows to obtain total CO2 from fossil fuels 
%

% GHG, 3x7980, 'Non-CO2 vector by sector and type of gas' 
% The three labels (rows) are the types of gas here, so N2O, CH4 and
% F-gasses 
% Unit is Mt CO2e (according to Anne Owen, 03Mar17)
% Ref: https://www.gtap.agecon.purdue.edu/resources/download/7637.pdf

% CO2_GTAP_HH from CO2_GTAP_HH_2011.mat
% Direct HH emissions data. 5 rows by 140 columns. Columns add up to total

% For this work we use four rows of GHGs, plus one summary row of CO2e (for
% industries; households ignored for now)

% f=sum(CO2_GTAP_sec,1)+(sum(GHG,1)*1000);
E_ind = f ;

% GHG are:
% CO2 from fossil fuels (CO2) 
% Nitrous Oxide (N2O) 
% Methane (CH4) 
% Fluorinated gasses (FGAS) 
% CO2 equivalents (sum of the four rows above) 
% Unit is kt CO2e
clear CO2_GTAP_sec ;
clear GHG ;
% csvwrite([resultspath 'E_ind.csv'], E_ind) ;

% Check balance of GTAP data
Xin  = sum(T,1) + sum(VA,1);
Xout = sum(T,2) + sum(Y_total,2);
GTAPbal = Xout - Xin';

% disp(['Sum of GTAP MRIO imbalance = ' num2str(sum(abs(GTAPbal))./(sum(sum(abs(T)))+sum(sum(abs(VA)))+sum(sum(abs(Y_total))))*100) '% of total production']);
% disp(['Max of GTAP MRIO imbalance = ' num2str(max(200*abs(GTAPbal)./(Xout+Xin'))) '% of total production']); 


%% Total country carbon footprints
% Uses the GTAP data as is.
 
% Calculate DIMs and TIMs
Xin  = Xin  .* (abs(Xin) >1e-3);  % set small values (smaller than 1000$ to zero to prevent large DIMs (turns out to be optimum)
  
DIMsE  = E_ind * diag(1./Xin) ;  % one row of total CO2e DIMs (because Anne's one-row vector f is used)
% five rows of DIMs for five categories of industry Emissions in E_ind block
DIMsE(isinf(DIMsE))=0 ;  % Set inf in DIMs to zero
DIMsE(isnan(DIMsE))=0 ;  % Set nan in DIMs to zero
% csvwrite([resultspath 'DIMsE.csv'], DIMsE) ;
   
% Calculate L 
L = inv(eye(size(A))-A) ;  % Leontief Inverse matrix
clear A ;  % not needed further in the script


%% Preparing data for city-specific analysis

% Define country location in GTAP:
% Example: The first sector of Canada has these coordinates:
% RegionNum = 27
% SectorNum = 1
% RegSecNum = ((RegionNum-1)*57)+SectorNum = 1483
% 1483	27	CAN	Paddy rice	PDR
% 
% Example: The last sector of Germany has these coordinates:
% RegionNum = 62
% SectorNum = 57
% RegSecNum = ((RegionNum-1)*57)+SectorNum = 3534
% 3534	62	DEU	Dwellings	DWE 

% Look up city and country IDs and other data
file =  fopen('City_lookups3.csv'); 
% This look up file has 8 columns:
% (1) row/city number
% (2) city 3 digit letter code (strings)
% (3) city population
% (4) city direct housing emissions
% (5) city direct travel emissions
% (6) country 3 digit letter code (strings)
% (7) country population
% (8) GTAP number = RegionNum

C_lookup = textscan(file, '%u %s %u %u %u %s %u %u','delimiter', ','); % this tells MATLAB if the data is a string %s or a number %u

for n = 1:78 ;

RegionNum = C_lookup{8}(n); % extract the GTAP region number
   
disp('.'); disp('.'); 
disp(['Country/region number: ' num2str(RegionNum) ]) ; 
disp(['City number: ' num2str(n) ]) ; 

RegSecNumStart= ((RegionNum-1)*57)+1 ;
RegSecNumEnd  = ((RegionNum-1)*57)+57 ;

% Load and prepare data for one city at a time
% Order: It is always City first, then RON (rest of nation)!
% So in TOR_GTAP_F.mat and TOR_GTAP_V.mat:
% City:  columns 1:57
% RON: columns 58:114
% And in TOR_GTAP_Y.mat:
% City:  columns 1-4
% RON: columns 5-8
% Ascertain that the first set of 57 values is smaller (i.e. for the City)
% than the second set of 57 values (i.e. for RON).

% Note on city final demand data: 
% The TOR_GTAP_Y file is 7980 x 8 and the column headings are:
% TOR_hhold
% TOR_gov
% TOR_cap
% TOR_vst
% RON_hhold
% RON_gov
% RON_cap
% RON_vst

% Load the n-th city in the list
   %this takes the n-th row of the 2nd col, turns it from an array to a
   %string then uses string concatenate (strcat) to attach to the file name

load([citydatapath strcat(char(C_lookup{2}(n)), '_GTAP_F.mat')]); 
load([citydatapath strcat(char(C_lookup{2}(n)), '_GTAP_V.mat')]); 
load([citydatapath strcat(char(C_lookup{2}(n)), '_GTAP_Y.mat')]); 
load([citydatapath strcat(char(C_lookup{2}(n)), '_GTAP_Yhh.mat')]); 
load([citydatapath strcat(char(C_lookup{2}(n)), '_GTAP_Yoth.mat')]); 

eval(['CITY_GTAP_F = ' strcat(char(C_lookup{2}(n)), '_GTAP_F') ';']) 
eval(['CITY_GTAP_V = ' strcat(char(C_lookup{2}(n)), '_GTAP_V') ';']) 
% eval(['CITY_GTAP_Y = ' strcat(char(C_lookup{2}(n)), '_GTAP_Y') ';']) ; obsolete, i.e. replace with Yhh and Yoth
eval(['CITY_GTAP_Yhh = ' strcat(char(C_lookup{2}(n)), '_GTAP_Yhh') ';']) 
eval(['CITY_GTAP_Yoth = ' strcat(char(C_lookup{2}(n)), '_GTAP_Yoth') ';']) 
CITY_GTAP_Y(:,1) = sum(CITY_GTAP_Yhh,2);
CITY_GTAP_Y(:,2:8) = CITY_GTAP_Yoth;

% combining Yhh and Yoth into one matrix of final demand by COICOP (to
% avoid handling them separately
CITY_GTAP_Y_COICOP = ([CITY_GTAP_Yhh, CITY_GTAP_Yoth]) ;

% CITY_GTAP_Y_COICOP is in two files ? the first, CITY_GTAP_Yhh, is the city hhold spends with 12 columns for 12 COICOP categories:
% �         Food and non-alcoholic beverages          
% �         Alcoholic beverages, tobacco and narcotics         
% �         Clothing and footwear  
% �         Housing, water, electricity, gas and other fuels 
% �         Furnishings, household equipment and routine household maintenance             
% �         Health  
% �         Transport           
% �         Communication               
% �         Recreation and culture 
% �         Education           
% �         Restaurants and hotels
% �         Miscellaneous goods and services
%  
% The second, CITY_GTAP_Yoth, are the other FD categories & 7 columns representing:
% �         City gov fd
% �         City cap fd
% �         City vst fd
% �         RoN hhold fd
% �         RoN gov fd
% �         RoN cap fd
% �         RoN vst fd

disp(['City chosen: ' char(C_lookup{2}(n))]) ; 

disp('.'); 
% disp(['Check whether total CO2e emissions for City     (' num2str(sum(CITY_GTAP_F(:,1:57))) ') are similar to total TRUE CO2e emissions for City (' num2str(sum(CITY_TRUE_F(:,1:57))) ')']) ; 
disp(['Check whether total CO2e emissions for City     (' num2str(sum(CITY_GTAP_F(:,1:57))) ') are smaller than total GTAP CO2e emissions for RON (' num2str(sum(CITY_GTAP_F(:,58:114))) ')']) ; 
% disp(['Check whether total CO2e emissions for City+RON (' num2str(sum(CITY_GTAP_F)) ') are similar to total TRUE CO2e emissions for City+RON (' num2str(sum(CITY_TRUE_F)) ')']) ; 
disp(['Check whether total CO2e emissions for City+RON (' num2str(sum(CITY_GTAP_F)) ') are equal to total national CO2e emissions (' num2str(sum(sum(E_ind(1,RegSecNumStart:RegSecNumEnd)))) ')']) ; 
disp('.'); 
disp(['Check whether total value added for City     (' num2str(sum(CITY_GTAP_V(:,1:57))) ') is smaller than total value added for RON (' num2str(sum(CITY_GTAP_V(:,58:114))) ')']) ; 
disp(['Check whether total value added for City+RON (' num2str(sum(CITY_GTAP_V)) ') is equal to total national value added (' num2str(sum(VA(:,RegSecNumStart:RegSecNumEnd))) ')']) ; 
disp('.'); 
disp(['Check whether FD for City               (' num2str(sum(sum(CITY_GTAP_Y(RegSecNumStart:RegSecNumEnd,1:4)))) ') is smaller than FD for RON (' num2str(sum(sum(CITY_GTAP_Y(RegSecNumStart:RegSecNumEnd,5:8)))) ')']) ; 
disp(['Check whether FD for City by COICOP     (' num2str(sum(sum(CITY_GTAP_Y_COICOP(RegSecNumStart:RegSecNumEnd,1:15)))) ') is smaller than FD for RON (' num2str(sum(sum(CITY_GTAP_Y_COICOP(RegSecNumStart:RegSecNumEnd,16:19)))) ')']) ; 
disp(['Check whether FD for City+RON           (' num2str(sum(sum(CITY_GTAP_Y(RegSecNumStart:RegSecNumEnd,:)))) ') is equal to national FD WITHOUT direct imports to FD (' num2str(sum(Y_total(RegSecNumStart:RegSecNumEnd,RegionNum))) ')']) ; 
disp(['Check whether FD for City+RON by COICOP (' num2str(sum(sum(CITY_GTAP_Y_COICOP(RegSecNumStart:RegSecNumEnd,:)))) ') is equal to national FD WITHOUT direct imports to FD (' num2str(sum(Y_total(RegSecNumStart:RegSecNumEnd,RegionNum))) ')']) ; 


%% Simplified City CF calculated with total city FD only, i.e. without split
% produces results table 

 result(n,1) = DIMsE * L * Y_total(:,RegionNum); % country product footprint in kt CO2e
 result(n,2) = sum(CO2_GTAP_HH(:,RegionNum),1)*1000; % country direct emissions converted to kt CO2e
 result(n,3) = (result(n,1)+result(n,2))*1000/double(C_lookup{7}(n)); % country carbon footprint per capita, incl. HH emissions, in t CO2e/cap
 result(n,4) = double(C_lookup{4}(n)) + double(C_lookup{5}(n)); % city direct emissions in kt CO2e (HH residential and private transport)
 eval(['result(n,5) = DIMsE * L * sum(' strcat(char(C_lookup{2}(n)), '_GTAP_Y(:,1:4)') ',2);'])
   	%this takes the total city demand from the city vector and calculates the simplified
   	%city CF without splitting. We need eval() here to make it work...
 result(n,6) = (result(n,4)+result(n,5))*1000/double(C_lookup{3}(n)); % simple city CF per capita, incl. HH emissions, in t CO2e/cap


%% Splitting the cross
% This section extracts the various parts of the country 'cross' in the
% GTAP-MRIO table and splits them into city (c) and rest-of-nation (r)
% parts
% disp('.'); toc;
% disp('**********************************************************************')
% disp('.'); disp('Splitting the cross'); 

% Defining the parts of the country cross
% disp('.'); disp(['Ascertain that RegionNum is ' num2str(RegionNum)]); 

Tnat = T(RegSecNumStart:RegSecNumEnd,RegSecNumStart:RegSecNumEnd) ; % 57x57 national T matrix

% Load City and RON specific data
VAcit = CITY_GTAP_V(:,1:57) ;
VAron = CITY_GTAP_V(:,58:114) ;
% Define segments of FD in GTAP format (4 columns of FD)
FDcit = CITY_GTAP_Y(RegSecNumStart:RegSecNumEnd,1:4) ;
FDron = CITY_GTAP_Y(RegSecNumStart:RegSecNumEnd,5:8) ;
FDcitIMPupper = CITY_GTAP_Y(1:(RegSecNumStart-1),1:4) ;
FDcitIMPlower = CITY_GTAP_Y((RegSecNumEnd+1):7980,1:4) ;
FDronIMPupper = CITY_GTAP_Y(1:(RegSecNumStart-1),5:8) ;
FDronIMPlower = CITY_GTAP_Y((RegSecNumEnd+1):7980,5:8) ;
% Define segments of FD in COICOP format (12+7 columns of FD)
FDcit_COICOP = CITY_GTAP_Y_COICOP(RegSecNumStart:RegSecNumEnd,1:15) ;
FDron_COICOP = CITY_GTAP_Y_COICOP(RegSecNumStart:RegSecNumEnd,16:19) ;
FDcit_COICOP_IMPupper = CITY_GTAP_Y_COICOP(1:(RegSecNumStart-1),1:15) ;
FDcit_COICOP_IMPlower = CITY_GTAP_Y_COICOP((RegSecNumEnd+1):7980,1:15) ;
FDron_COICOP_IMPupper = CITY_GTAP_Y_COICOP(1:(RegSecNumStart-1),16:19) ;
FDron_COICOP_IMPlower = CITY_GTAP_Y_COICOP((RegSecNumEnd+1):7980,16:19) ;

% Example for Mel FD in COICOP format:
% MEL_hh_GTAP_Y ? this will be an n*12 matrix of Melbourne?s household expenditure by GTAP(down) and COICOP (across)
% MEL_oth_GTAP_Y ? this will be an n*7 matrix where the columns represent Mel gov, Mel cap, Mel vst, RoN hh, Ron gov, Ron cap and Ron vst
% So Melbourne?s footprint will be shown broken down into 15 segments: 12 household COICOP products, gov, cap and vst.

% FDcit  = FDcit  .* (FDcit>0);  % set negative values to zero % commented out because there are no negative values in GTAP FD
% FDron  = FDron  .* (FDron>0);  % set negative values to zero % commented out because there are no negative values in GTAP FD
FDnat = FDcit + FDron; VAnat = VAcit + VAron;

IMPnatIDupper = T(1:RegSecNumStart-1,RegSecNumStart:RegSecNumEnd) ; % upper part of national imports to ID; 57 columns
IMPnatIDlower = T(RegSecNumEnd+1:7980,RegSecNumStart:RegSecNumEnd); % lower part of national imports to ID; 57 columns
EXPnatIDleft = T(RegSecNumStart:RegSecNumEnd,1:RegSecNumStart-1) ; % left part of national exports to ID; 57 rows
EXPnatIDright = T(RegSecNumStart:RegSecNumEnd,RegSecNumEnd+1:7980) ; % right part of national exports to ID; 57 rows

IMPnatFDupper = [Y_g(1:RegSecNumStart-1,RegionNum), Y_h(1:RegSecNumStart-1,RegionNum), Y_k(1:RegSecNumStart-1,RegionNum), Y_VST(1:RegSecNumStart-1,RegionNum)] ;    % upper part of national imports to FD; 4 columns
IMPnatFDlower = [Y_g(RegSecNumEnd+1:7980,RegionNum), Y_h(RegSecNumEnd+1:7980,RegionNum), Y_k(RegSecNumEnd+1:7980,RegionNum), Y_VST(RegSecNumEnd+1:7980,RegionNum)]; % lower part of national imports to FD; 4 columns
EXPnatFDleft = Y_total(RegSecNumStart:RegSecNumEnd,1:RegionNum-1) ; % left part of national exports to FD; 57 rows; total FD sufficient
EXPnatFDright = Y_total(RegSecNumStart:RegSecNumEnd,RegionNum+1:140) ; % right part of national exports to FD; 57 rows; total FD sufficient

% Balance check for country data cross before it is split
Xoutnat= sum(EXPnatIDleft,2) + sum(Tnat,2) + sum(EXPnatIDright,2) + sum(EXPnatFDleft,2) + sum(FDnat,2) + sum(EXPnatFDright,2) ; % total output of country cross; parts in MRIO order
Xinnat = sum(IMPnatIDupper,1) + sum(Tnat,1) + sum(IMPnatIDlower,1) + sum(VAnat,1); % total input of country cross; parts in MRIO order
IObalnat = Xoutnat - Xinnat';
disp('.'); 
disp(['Sum of country cross IO imbalance BEFORE splitting = ' num2str(sum(abs(IObalnat))./(sum(sum(abs(Tnat)))+sum(sum(abs(VAnat)))+sum(sum(abs(FDnat))))*100) '% of total production.']);
disp(['Max of country cross IO imbalance BEFORE splitting = ' num2str(max(200*abs(IObalnat)./(Xoutnat+Xinnat'))) '% of total production.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Treg_choice 0: Splitting Options T0 and Y0 (i.e. separting city and rest-of-nation tables as well as interregional trade matrices)
% In Treg0, a location quotient approach is employed.
% LQs are interpreted in a way that if LQ_j^s>=1 then region s is self-sufficient in its production of commodity j. 
% columns are split by LQ (tau), rows are split by VA.
 
% Choose LQ method here:
LQ=6;
 
switch LQ
% Trading coefficients (tau)   
    case 1
    % A) Simple LQ
    tau_T_cc = min(1, repmat(  ( (VAcit./sum(VAcit)) ./ (VAnat./sum(VAnat)) )' ,[1 size(Tnat,1)] ) ); 
    tau_T_rr = min(1, repmat(  ( (VAron./sum(VAron)) ./ (VAnat./sum(VAnat)) )' ,[1 size(Tnat,1)] ) );
    case 2
    % B) Cross-industry LQ (doesn't need repmat-ing with [1 size(Tnat,1)]; only for T, not for y)
    for i=1:size(Tnat,1); for j=1:size(Tnat,1);
        tau_T_cc(i,j) = min(1, (VAcit(1,i)./VAcit(1,j)) ./ (VAnat(1,i)./VAnat(1,j)) );
        tau_T_rr(i,j) = min(1, (VAron(1,i)./VAron(1,j)) ./ (VAnat(1,i)./VAnat(1,j)) );
    end; end;
    case 3
    % C) Symmetric cross-industry LQ (doesn't need repmat-ing with [1 size(Tnat,1)]; only for T, not for y)
    for i=1:size(Tnat,1); for j=1:size(Tnat,1);
        tau_T_cc(i,j) = 2-2./(min(1, (VAcit(1,i)./VAcit(1,j)) ./ (VAnat(1,i)./VAnat(1,j)) ) +1);
        tau_T_rr(i,j) = 2-2./(min(1, (VAron(1,i)./VAron(1,j)) ./ (VAnat(1,i)./VAnat(1,j)) ) +1);
    end; end;
    case 4
    % D) Round?s semi-logarithmic LQ (doesn't need repmat-ing with [1 size(Tnat,1)]; only for T, not for y)
    for i=1:size(Tnat,1); for j=1:size(Tnat,1);
        tau_T_cc(i,j) = min(1, (VAcit(1,i)./sum(VAcit)) ./ (VAnat(1,i)./sum(VAnat)) ./ log2(1+ (VAcit(1,j)./sum(VAcit)) ./ (VAnat(1,j)./sum(VAnat))) ); 
        tau_T_rr(i,j) = min(1, (VAron(1,i)./sum(VAron)) ./ (VAnat(1,i)./sum(VAnat)) ./ log2(1+ (VAron(1,j)./sum(VAron)) ./ (VAnat(1,j)./sum(VAnat))) ); 
    end; end;
    case 5
    % E) Flegg's LQ
    for i=1:size(Tnat,1); for j=1:size(Tnat,1);
        if i==j;
            tau_T_cc(i,j) = min(1, (log2(1+sum(VAcit)/sum(VAnat)))^0.3 .* (VAcit(1,i)./sum(VAcit)) ./ (VAnat(1,i)./sum(VAnat)) ); 
            tau_T_rr(i,j) = min(1, (log2(1+sum(VAron)/sum(VAnat)))^0.3 .* (VAron(1,i)./sum(VAron)) ./ (VAnat(1,i)./sum(VAnat)) );
        else
            tau_T_cc(i,j) = min(1, (log2(1+sum(VAcit)/sum(VAnat)))^0.3 .* (VAcit(1,i)./VAcit(1,j)) ./ (VAnat(1,i)./VAnat(1,j)) );
            tau_T_rr(i,j) = min(1, (log2(1+sum(VAron)/sum(VAnat)))^0.3 .* (VAron(1,i)./VAron(1,j)) ./ (VAnat(1,i)./VAnat(1,j)) );
        end;
    end; end;
    case 6
    % F) Flegg's Adjusted LQ
    for i=1:size(Tnat,1); for j=1:size(Tnat,1);
        if i==j;
            tau_T_cc(i,j) = (log2(1+sum(VAcit)/sum(VAnat)))^0.3 .* (VAcit(1,i)./sum(VAcit)) ./ (VAnat(1,i)./sum(VAnat)); 
            tau_T_rr(i,j) = (log2(1+sum(VAron)/sum(VAnat)))^0.3 .* (VAron(1,i)./sum(VAron)) ./ (VAnat(1,i)./sum(VAnat));
        else
            tau_T_cc(i,j) = (log2(1+sum(VAcit)/sum(VAnat)))^0.3 .* (VAcit(1,i)./VAcit(1,j)) ./ (VAnat(1,i)./VAnat(1,j));
            tau_T_rr(i,j) = (log2(1+sum(VAron)/sum(VAnat)))^0.3 .* (VAron(1,i)./VAron(1,j)) ./ (VAnat(1,i)./VAnat(1,j));
        end;
        tau_SLQ_cc = repmat(  ( (VAcit./sum(VAcit)) ./ (VAnat./sum(VAnat)) )' ,[1 size(Tnat,1)] ); 
        tau_SLQ_rr = repmat(  ( (VAron./sum(VAron)) ./ (VAnat./sum(VAnat)) )' ,[1 size(Tnat,1)] );
        tau_T_cc = min(1, tau_T_cc .* (1.*double(tau_SLQ_cc<=1) + log2(1+tau_SLQ_cc).*double(tau_SLQ_cc>1)) );
        tau_T_rr = min(1, tau_T_rr .* (1.*double(tau_SLQ_rr<=1) + log2(1+tau_SLQ_rr).*double(tau_SLQ_rr>1)) );
    end; end;
 
end;
 
% Off-diagonal coefficients
tau_T_rc = (1 - tau_T_cc); % if tau_T_cc = 1, then c is self-sufficient, and does not need imports from r
tau_T_cr = (1 - tau_T_rr); % if tau_T_rr = 1, then r is self-sufficient, and does not need imports from c
 
% regionalising Tnat with Option T0
Treg0 = [ [ Tnat .* tau_T_cc .* repmat((VAcit./VAnat)',[1 size(Tnat,2)])     , ...
            Tnat .* tau_T_cr .* repmat((VAron./VAnat)',[1 size(Tnat,2)]) ]   ; ...
          [ Tnat .* tau_T_rc .* repmat((VAcit./VAnat)',[1 size(Tnat,2)])     , ...
            Tnat .* tau_T_rr .* repmat((VAron./VAnat)',[1 size(Tnat,2)]) ] ] ;
VAreg0 = [VAcit , VAron];
 
% regionalising FD with LQ: VA (rows) and trading coefficients (columns) (Option Y0)
tau_T_cc = min(1, (VAcit./sum(VAcit)) ./ (VAnat./sum(VAnat)) ); 
tau_T_rr = min(1, (VAron./sum(VAron)) ./ (VAnat./sum(VAnat)) );
tau_T_rc = (1 - tau_T_cc); % if tau_T_cc = 1, then c is self-sufficient, and does not need imports from r
tau_T_cr = (1 - tau_T_rr); % if tau_T_rr = 1, then r is self-sufficient, and does not need imports from c
FDreg0 = [ FDcit .* repmat(tau_T_cc',[1 size(FDcit,2)]),...
           FDron .* repmat(tau_T_cr',[1 size(FDron,2)]) ; ...
           FDcit .* repmat(tau_T_rc',[1 size(FDcit,2)]),...
           FDron .* repmat(tau_T_rr',[1 size(FDron,2)]) ] ;
      
% regionalising FD by COICOP with VA (Option Y1):
FDreg0_COICOP = [ FDcit_COICOP .* repmat(tau_T_cc',1,15),...
                  FDron_COICOP .* repmat(tau_T_cr',1,4) ; ...
                  FDcit_COICOP .* repmat(tau_T_rc',1,15),...
                  FDron_COICOP .* repmat(tau_T_rr',1,4) ] ;
       
% splitting national exports to final demand in other regions (Option Y0):
% first left part
EXPtoFDleftsplit0 = [ EXPnatFDleft .* repmat((VAcit./VAnat)',[1 size(EXPnatFDleft,2)]) ;
                      EXPnatFDleft .* repmat((VAron./VAnat)',[1 size(EXPnatFDleft,2)]) ] ;

% then right part
EXPtoFDrightsplit0 = [ EXPnatFDright .* repmat((VAcit./VAnat)',[1 size(EXPnatFDright,2)]) ;
                       EXPnatFDright .* repmat((VAron./VAnat)',[1 size(EXPnatFDright,2)]) ] ;
 
% splitting national exports to intermediate demand in other regions (Option T0)
% first left part
EXPtoIDleftsplit0 = [ EXPnatIDleft .* repmat((VAcit./VAnat)',[1 size(EXPnatIDleft,2)]) ;
                      EXPnatIDleft .* repmat((VAron./VAnat)',[1 size(EXPnatIDleft,2)]) ] ;

% then right part
EXPtoIDrightsplit0 = [ EXPnatIDright .* repmat((VAcit./VAnat)',[1 size(EXPnatIDright,2)]) ;
                       EXPnatIDright .* repmat((VAron./VAnat)',[1 size(EXPnatIDright,2)]) ] ;          

% splitting imports from other regions to national intermediate demand (Option T0)  
% first upper part
IMPtoIDuppersplit0 = [ IMPnatIDupper .* repmat(VAcit./VAnat,[size(IMPnatIDupper,1) 1]), IMPnatIDupper .* repmat(VAron./VAnat,[size(IMPnatIDupper,1) 1]) ] ;
              
% then lower part
IMPtoIDlowersplit0 = [ IMPnatIDlower .* repmat(VAcit./VAnat,[size(IMPnatIDlower,1) 1]), IMPnatIDlower .* repmat(VAron./VAnat,[size(IMPnatIDlower,1) 1]) ] ;
              
% Balance check for country data cross after it has been split with Option1
disp('.'); 
disp(['Total transaction volume of T before splitting: ' num2str(sum(sum(Tnat))) '; after splitting: ' num2str(sum(sum(Treg0))) ' with Option 0.']);

Xoutreg0  = sum(EXPtoIDleftsplit0,2) + sum(Treg0,2) + sum(EXPtoIDrightsplit0,2) + sum(EXPtoFDleftsplit0,2) + sum(FDreg0,2) + sum(EXPtoFDrightsplit0,2) ; % total output of split country cross; parts in MRIO order
Xinreg0   = sum(IMPtoIDuppersplit0,1) + sum(Treg0,1) + sum(IMPtoIDlowersplit0,1) + sum(VAreg0,1); % total input of split country cross; parts in MRIO order
IObalreg0 = Xoutreg0 - Xinreg0';

Treg_choice = 0;

    Treg  = Treg0 ;
    VAreg = VAreg0 ;
    FDreg = FDreg0 ;
    FDreg_COICOP = FDreg0_COICOP ;
    Xoutreg = Xoutreg0 ;
    Xinreg  = Xinreg0 ;
    EXPtoIDleftsplit  = EXPtoIDleftsplit0 ;
    EXPtoIDrightsplit = EXPtoIDrightsplit0 ;
    EXPtoFDleftsplit  = EXPtoFDleftsplit0 ;
    EXPtoFDrightsplit = EXPtoFDrightsplit0 ;
    IMPtoIDuppersplit = IMPtoIDuppersplit0 ;
    IMPtoIDlowersplit = IMPtoIDlowersplit0 ;
    disp('.'); 
    disp(['Splitting Option ' num2str(Treg_choice) ' (Treg_choice 0)']);
    disp(['Splitting Option T0 and Y0 (LQ method: rows are split by VA and quadrants by trading coefficients)']) ; 
    disp(['LQ method chosen: ' num2str(LQ)]);

clear Treg0 ;
clear VAreg0 ;
clear FDreg0 ;
clear FDreg0_COICOP ;
clear Xoutreg0 ;
clear Xinreg0 ;
clear EXPtoIDleftsplit0 ;
clear EXPtoIDrightsplit0 ;
clear EXPtoFDleftsplit0 ;
clear EXPtoFDrightsplit0 ;
clear IMPtoIDuppersplit0 ;
clear IMPtoIDlowersplit0 ;


%% RAS Balancing 
% Now choose RAS balancing with FD fixed (RAS_choice == 1), with VA fixed
% (RAS_choice == 2) or with both FD and VA fixed (RAS_choice == 3).
% RAS == 4 means neither FD nor VA are fixed.
% To switch off (RAS0) just write RAS_choice = 0 ;
 
RAS_choice = 0 ;

if RAS_choice == 1 ;    
    %% RAS Balancing with FD fixed
    disp(['Balancing Option 1: RAS balancing with FD fixed (Option RAS1)']); disp('.'); 
    for loop = 1:200 ;  % loop counter; number of balancing steps

        % Capital X stands for global total in/output, 
        % Small x stands for regional (local) total in/output

        % First, determine what row totals in must be to balance with Xinreg
        % when FD is kept fixed!
        xoutreg = Xinreg' - sum(FDreg,2);

        % Then, determine the scaling factor between actual and target xoutreg
        scalar_out = xoutreg ./ ( sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(EXPtoFDrightsplit,2) ) ;
        scalar_out(isnan(scalar_out))=1 ;  % Set NaN to one
        scalar_out(isinf(scalar_out))=1 ;  % Set NaN to one
        scalar_out = scalar_out  .* (scalar_out>0) ;  % set negative values to zero
        scalar_out0= ones(114,1) .* (scalar_out<1e-10) ;  % set negative values to one
        scalar_out = scalar_out + scalar_out0 ;  % add ones instead of zeros
        % scalar_out
%     	disp(['sum of scalar_out in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_out))) ' ']) ;
    	
        % The actual rows are now scaled so that they adhere to xoutreg
        EXPtoIDleftsplit  = EXPtoIDleftsplit  .* repmat(scalar_out,[1 size(EXPtoIDleftsplit,2)]) ;
        Treg              = Treg              .* repmat(scalar_out,[1 size(Treg,2)]) ;
        EXPtoIDrightsplit = EXPtoIDrightsplit .* repmat(scalar_out,[1 size(EXPtoIDrightsplit,2)]) ;
        EXPtoFDleftsplit  = EXPtoFDleftsplit  .* repmat(scalar_out,[1 size(EXPtoFDleftsplit,2)]) ;
        EXPtoFDrightsplit = EXPtoFDrightsplit .* repmat(scalar_out,[1 size(EXPtoFDrightsplit,2)]) ;
                
        % For columns, the procedure is different, because we allow VA to change (not fixed)!
        % The only condition is that Xinreg must become equal to Xoutreg
        % Determine what global column totals must be to balance with Xoutreg
        Xinreg = ( sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) )' ;
        
        % Then, determine the scaling factor between actual and target Xinreg
        scalar_in = Xinreg ./ (sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1)) ;        
        scalar_in(isnan(scalar_in))=1 ;  % Set NaN to one
        scalar_in(isinf(scalar_in))=1 ;  % Set NaN to one
        scalar_in = scalar_in  .* (scalar_in>0) ;  % set negative values to zero
        scalar_in0= ones(1,114).* (scalar_in<1e-10) ;  % set negative values to one
        scalar_in = scalar_in + scalar_in0 ;  % add ones instead of zeros
        % scalar_in
%     	disp(['sum of scalar_in in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_in))) ' ']);
    	
        % The actual columns of elements of the vertical part of the cross are now scaled so that they adhere to Xinreg
        IMPtoIDuppersplit = IMPtoIDuppersplit .* repmat(scalar_in,[size(IMPtoIDuppersplit,1) 1]) ;
        Treg              = Treg              .* repmat(scalar_in,[size(Treg,1) 1]) ;
        IMPtoIDlowersplit = IMPtoIDlowersplit .* repmat(scalar_in,[size(IMPtoIDlowersplit,1) 1]) ;
        VAreg             = VAreg             .* repmat(scalar_in,1,1) ;
        
    end  
    
    % Balance check for country data cross after it has been balance with FD fixed 2
%     disp('.'); 
%     disp(['Total transaction volume of T after RAS balancing with FD fixed: ' num2str(sum(sum(Treg))) ]);

    Xoutreg  = sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) ; % total output of balanced country cross; parts in MRIO order
    Xinreg   = sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1); % total input of balanced country cross; parts in MRIO order
    IObalreg = Xoutreg - Xinreg';

%     disp('.'); 
    disp(['Total national input: ' num2str(sum(sum(Xinnat))) '; total national output: ' num2str(sum(sum(Xoutnat))) '.']);
    disp(['Total regional input: ' num2str(sum(sum(Xinreg))) '; total regional output: ' num2str(sum(sum(Xoutreg))) '.']);
    disp('.'); 
    disp(['Sum of country cross IO imbalance AFTER balancing with FD fixed = ' num2str(sum(abs(IObalreg))./(sum(sum(Xoutreg+Xinreg'))/2)*100) '% of total output.']);
    disp(['Max of country cross IO imbalance AFTER balancing with FD fixed = ' num2str(max(200*abs(IObalreg)./(Xoutreg+Xinreg'))) '% of total output.']); disp('.');
        
%     csvwrite([resultspath 'Treg_balanced.csv'], Treg);
    
    disp('.'); toc;
    disp('**********************************************************************')

else if RAS_choice == 2 ;
	%% RAS Balancing with VA fixed
    disp(['Balancing Option 2: RAS balancing with VA fixed (Option RAS2)']); disp('.'); 
	for loop = 1:200 ;  % loop counter; number of balancing steps

        % Capital X stands for global total in/output, 
        % Small x stands for regional (local) total in/output

        % First, determine what column totals in Treg must be to balance with Xoutreg
        % This means we keep VA fixed!
        xinreg = Xoutreg' - sum(VAreg,1);

        % Then, determine the scaling factor between actual and target xinreg
        scalar_in = xinreg ./ (sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1)) ;
        scalar_in(isnan(scalar_in))=1 ;  % Set NaN to one
        scalar_in(isinf(scalar_in))=1 ;  % Set NaN to one
        scalar_in = scalar_in  .* (scalar_in>0) ;  % set negative values to zero
        scalar_in0= ones(1,114).* (scalar_in<1e-10) ;  % set negative values to one
        scalar_in = scalar_in + scalar_in0 ;  % add ones instead of zeros
        % scalar_in
    	disp(['sum of scalar_in in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_in))) ' ']) ;
       
        % The actual colums are now scaled so that they adhere to xinreg
        Treg              = Treg              .* repmat(scalar_in,[size(Treg,1) 1]) ;  % check whether Treg = Treg .* repmat(scalar_in,sizeTreg(1:1),1) is the same!
        IMPtoIDuppersplit = IMPtoIDuppersplit .* repmat(scalar_in,[size(IMPtoIDuppersplit,1) 1]) ;
        IMPtoIDlowersplit = IMPtoIDlowersplit .* repmat(scalar_in,[size(IMPtoIDlowersplit,1) 1]) ;
        % Treg

        % For rows, the procedure is different, because we allow FD to change (not fixed)!
        % The only condition is that Xoutreg must become equal to Xinreg
        % Determine what global row totals must be to balance with Xinreg
        Xoutreg = (sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1))' ;

        % Then, determine the scaling factor between actual and target Xoutreg
        scalar_out = Xoutreg ./ (sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2)) ;
        scalar_out(isnan(scalar_out))=1 ;  % Set NaN to one
        scalar_out(isinf(scalar_out))=1 ;  % Set NaN to one
        scalar_out = scalar_out  .* (scalar_out>0) ;  % set negative values to zero
        scalar_out0= ones(114,1) .* (scalar_out<1e-10) ;  % set negative values to one
        scalar_out = scalar_out + scalar_out0 ;  % add ones instead of zeros
        % scalar_out
    	disp(['sum of scalar_out in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_out))) ' ']) ;
    	
        % The actual rows of horizontal part of the cross are now scaled so that they adhere to Xoutreg      
        EXPtoIDleftsplit  = EXPtoIDleftsplit  .* repmat(scalar_out,[1 size(EXPtoIDleftsplit,2)]) ;
        Treg              = Treg              .* repmat(scalar_out,[1 size(Treg,2)]) ;
        EXPtoIDrightsplit = EXPtoIDrightsplit .* repmat(scalar_out,[1 size(EXPtoIDrightsplit,2)]) ;
        EXPtoFDleftsplit  = EXPtoFDleftsplit  .* repmat(scalar_out,[1 size(EXPtoFDleftsplit,2)]) ;
        FDreg             = FDreg             .* repmat(scalar_out,[1 size(FDreg,2)]) ; 
        FDreg_COICOP	  = FDreg_COICOP      .* repmat(scalar_out,[1 size(FDreg_COICOP,2)]) ;  % FDreg_COICOP also needs to be scaled
        EXPtoFDrightsplit = EXPtoFDrightsplit .* repmat(scalar_out,[1 size(EXPtoFDrightsplit,2)]) ;
        % Treg
    end  
    
    % Balance check for country data cross after it has been balance with
    % VA fixed (RAS2)
%     disp('.'); 
%     disp(['Total transaction volume of T after RAS balancing with VA fixed: ' num2str(sum(sum(Treg))) ]);

    Xoutreg  = sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) ; % total output of balanced country cross; parts in MRIO order
    Xinreg   = sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1); % total input of balanced country cross; parts in MRIO order
    IObalreg = Xoutreg - Xinreg';

    disp('.'); 
    disp(['Total national input: ' num2str(sum(sum(Xinnat))) '; total national output: ' num2str(sum(sum(Xoutnat))) '.']);
    disp(['Total regional input: ' num2str(sum(sum(Xinreg))) '; total regional output: ' num2str(sum(sum(Xoutreg))) '.']);
    disp('.'); 
    disp(['Sum of country cross IO imbalance AFTER balancing with VA fixed = ' num2str(sum(abs(IObalreg))./(sum(sum(Xoutreg+Xinreg'))/2)*100) '% of total output.']);
    disp(['Max of country cross IO imbalance AFTER balancing with VA fixed = ' num2str(max(200*abs(IObalreg)./(Xoutreg+Xinreg'))) '% of total output.']); disp('.');
        
%     csvwrite([resultspath 'Treg_balanced.csv'], Treg);
    
    disp('.'); toc;
    disp('**********************************************************************')
    
else if RAS_choice == 3 ;
	%% RAS Balancing with FD and VA fixed
    disp(['Balancing Option 3: RAS balancing with both FD and VA fixed (Option RAS3)']); disp('.');  
	for loop = 1:200 ;  % loop counter; number of balancing steps

        % Capital X stands for global total in/output, 
        % Small x stands for regional (local) total in/output

        % First, determine what row totals in must be to balance with Xinreg
        % For this, we first need to calculate the  Xinreg (it changes in
        % every loop)
        Xinreg  = sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1); % total input of  country cross; parts in MRIO order
        
        % when FD is kept fixed!
        xoutreg = Xinreg' - sum(FDreg,2);

        % Then, determine the scaling factor between actual and target xoutreg
        scalar_out = xoutreg ./ ( sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(EXPtoFDrightsplit,2) ) ;
        scalar_out(isnan(scalar_out))=1 ;  % Set NaN to one
        scalar_out(isinf(scalar_out))=1 ;  % Set NaN to one
        scalar_out = scalar_out  .* (scalar_out>0) ;  % set negative values to zero
        scalar_out0= ones(114,1) .* (scalar_out<1e-10) ;  % set negative values to one
        scalar_out = scalar_out + scalar_out0 ;  % add ones instead of zeros
        % scalar_out
    	disp(['sum of scalar_out in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_out))) ' ']) ;
    	
        % The actual rows are now scaled so that they adhere to xoutreg
        EXPtoIDleftsplit  = EXPtoIDleftsplit  .* repmat(scalar_out,[1 size(EXPtoIDleftsplit,2)]) ;
        Treg              = Treg              .* repmat(scalar_out,[1 size(Treg,2)]) ;
        EXPtoIDrightsplit = EXPtoIDrightsplit .* repmat(scalar_out,[1 size(EXPtoIDrightsplit,2)]) ;
        EXPtoFDleftsplit  = EXPtoFDleftsplit  .* repmat(scalar_out,[1 size(EXPtoFDleftsplit,2)]) ;
        EXPtoFDrightsplit = EXPtoFDrightsplit .* repmat(scalar_out,[1 size(EXPtoFDrightsplit,2)]) ;
         
        % Then we determine what column totals in Treg must be to balance with Xoutreg
        % For this, we first need to calculate the new Xoutreg:
        Xoutreg  = sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) ; % total output of country cross; parts in MRIO order

        % Now we can keep VA fixed by:
        xinreg = Xoutreg' - sum(VAreg,1);

        % Then, determine the scaling factor between actual and target xinreg
        scalar_in = xinreg ./ (sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1)) ;
        scalar_in(isnan(scalar_in))=1 ;  % Set NaN to one
        scalar_in(isinf(scalar_in))=1 ;  % Set NaN to one
        scalar_in = scalar_in  .* (scalar_in>0) ;  % set negative values to zero
        scalar_in0= ones(1,114).* (scalar_in<1e-10) ;  % set negative values to one
        scalar_in = scalar_in + scalar_in0 ;  % add ones instead of zeros
        % scalar_in
    	disp(['sum of scalar_in in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_in))) ' ']) ;
       
        % The actual colums are now scaled so that they adhere to xinreg
        Treg              = Treg              .* repmat(scalar_in,[size(Treg,1) 1]) ;  % check whether Treg = Treg .* repmat(scalar_in,sizeTreg(1:1),1) is the same!
        IMPtoIDuppersplit = IMPtoIDuppersplit .* repmat(scalar_in,[size(IMPtoIDuppersplit,1) 1]) ;
        IMPtoIDlowersplit = IMPtoIDlowersplit .* repmat(scalar_in,[size(IMPtoIDlowersplit,1) 1]) ;
        % Treg
   end
    
    % Balance check for country data cross after it has been balanced with
    % VA fixed (RAS2)
%     disp('.'); 
%     disp(['Total transaction volume of T after RAS balancing both FD and VA fixed: ' num2str(sum(sum(Treg))) ]);

    Xoutreg  = sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) ; % total output of balanced country cross; parts in MRIO order
    Xinreg   = sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1); % total input of balanced country cross; parts in MRIO order
    IObalreg = Xoutreg - Xinreg';

    disp('.'); 
    disp(['Total national input: ' num2str(sum(sum(Xinnat))) '; total national output: ' num2str(sum(sum(Xoutnat))) '.']);
    disp(['Total regional input: ' num2str(sum(sum(Xinreg))) '; total regional output: ' num2str(sum(sum(Xoutreg))) '.']);
    disp('.'); 
    disp(['Sum of country cross IO imbalance AFTER balancing with FD and VA fixed = ' num2str(sum(abs(IObalreg))./(sum(sum(Xoutreg+Xinreg'))/2)*100) '% of total output.']);
    disp(['Max of country cross IO imbalance AFTER balancing with FD and VA fixed = ' num2str(max(200*abs(IObalreg)./(Xoutreg+Xinreg'))) '% of total output.']); disp('.');
        
%     csvwrite([resultspath 'Treg_balanced.csv'], Treg);
    
    disp('.'); toc;
    disp('**********************************************************************')
    
else if RAS_choice == 4 ;
	%% RAS Balancing with nothing fixed
    disp(['Balancing Option 4: RAS balancing without constraints (FD and VA variable; Option RAS4)']); disp('.');  
	for loop = 1:200 ;  % loop counter; number of balancing steps

        % Capital X stands for global total in/output,     
    
        % The first condition is that Xoutreg must become equal to Xinreg
        % Determine what global row totals must be to balance with Xinreg
        Xoutreg = (sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1))' ;

        % Then, determine the scaling factor between actual and target Xoutreg
        scalar_out = Xoutreg ./ (sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2)) ;
        scalar_out(isnan(scalar_out))=1 ;  % Set NaN to one
        scalar_out(isinf(scalar_out))=1 ;  % Set NaN to one
        scalar_out = scalar_out  .* (scalar_out>0) ;  % set negative values to zero
        scalar_out0= ones(114,1) .* (scalar_out<1e-10) ;  % set negative values to one
        scalar_out = scalar_out + scalar_out0 ;  % add ones instead of zeros
        % scalar_out
    	disp(['sum of scalar_out in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_out))) ' ']) ;
    	
        % The actual rows of horizontal part of the cross are now scaled so that they adhere to the new Xoutreg      
        EXPtoIDleftsplit  = EXPtoIDleftsplit  .* repmat(scalar_out,[1 size(EXPtoIDleftsplit,2)]) ;
        Treg              = Treg              .* repmat(scalar_out,[1 size(Treg,2)]) ;
        EXPtoIDrightsplit = EXPtoIDrightsplit .* repmat(scalar_out,[1 size(EXPtoIDrightsplit,2)]) ;
        EXPtoFDleftsplit  = EXPtoFDleftsplit  .* repmat(scalar_out,[1 size(EXPtoFDleftsplit,2)]) ;
        FDreg             = FDreg             .* repmat(scalar_out,[1 size(FDreg,2)]) ;  
        FDreg_COICOP	  = FDreg_COICOP      .* repmat(scalar_out,[1 size(FDreg_COICOP,2)]) ;  % FDreg_COICOP also needs to be scaled
        EXPtoFDrightsplit = EXPtoFDrightsplit .* repmat(scalar_out,[1 size(EXPtoFDrightsplit,2)]) ;
        % Treg
    
        % The second condition is that Xinreg must become equal to Xoutreg
        % Determine what global column totals must be to balance with Xoutreg
        Xinreg = ( sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) )' ;
        
        % Then, determine the scaling factor between actual and target Xinreg
        scalar_in = Xinreg ./ (sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1)) ;        
        scalar_in(isnan(scalar_in))=1 ;  % Set NaN to one
        scalar_in(isinf(scalar_in))=1 ;  % Set NaN to one
        scalar_in = scalar_in  .* (scalar_in>0) ;  % set negative values to zero
        scalar_in0= ones(1,114).* (scalar_in<1e-10) ;  % set negative values to one
        scalar_in = scalar_in + scalar_in0 ;  % add ones instead of zeros
        % scalar_in
    	disp(['sum of scalar_in in loop ' num2str(loop) ' is ' num2str(sum(sum(scalar_in))) ' ']);
    	
        % The actual columns of elements of the vertical part of the cross are now scaled so that they adhere to Xinreg
        IMPtoIDuppersplit = IMPtoIDuppersplit .* repmat(scalar_in,[size(IMPtoIDuppersplit,1) 1]) ;
        Treg              = Treg              .* repmat(scalar_in,[size(Treg,1) 1]) ;
        IMPtoIDlowersplit = IMPtoIDlowersplit .* repmat(scalar_in,[size(IMPtoIDlowersplit,1) 1]) ;
        VAreg             = VAreg             .* repmat(scalar_in,1,1) ;
    
    end
    
        % Balance check for country data cross after it has been balanced 
        Xoutreg  = sum(EXPtoIDleftsplit,2) + sum(Treg,2) + sum(EXPtoIDrightsplit,2) + sum(EXPtoFDleftsplit,2) + sum(FDreg,2) + sum(EXPtoFDrightsplit,2) ; % total output of balanced country cross; parts in MRIO order
        Xinreg   = sum(IMPtoIDuppersplit,1) + sum(Treg,1) + sum(IMPtoIDlowersplit,1) + sum(VAreg,1); % total input of balanced country cross; parts in MRIO order
        IObalreg = Xoutreg - Xinreg';

        disp('.');
        disp(['Total national input: ' num2str(sum(sum(Xinnat))) '; total national output: ' num2str(sum(sum(Xoutnat))) '.']);
        disp(['Total regional input: ' num2str(sum(sum(Xinreg))) '; total regional output: ' num2str(sum(sum(Xoutreg))) '.']);
        disp('.'); 
        disp(['Sum of country cross IO imbalance AFTER balancing with FD and VA variable = ' num2str(sum(abs(IObalreg))./(sum(sum(Xoutreg+Xinreg'))/2)*100) '% of total output.']);
        disp(['Max of country cross IO imbalance AFTER balancing with FD and VA variable = ' num2str(max(200*abs(IObalreg)./(Xoutreg+Xinreg'))) '% of total output.']);
        
%     csvwrite([resultspath 'Treg_balanced.csv'], Treg);
    
    disp('.'); toc;
    disp('**********************************************************************')    
    
   end
  end
 end
end


%% Reassembling enlarged GTAP table (with City and RON cross)

% Extracting GTAP table quadrants outside of cross:
GTAPupperleft  = T(1:(RegSecNumStart-1) , 1:(RegSecNumStart-1)) ;
GTAPupperright = T(1:(RegSecNumStart-1) , (RegSecNumEnd+1):7980) ;
GTAPlowerleft  = T((RegSecNumEnd+1):7980, 1:(RegSecNumStart-1)) ;
GTAPlowerright = T((RegSecNumEnd+1):7980, (RegSecNumEnd+1):7980) ;

% Assembling enlarged GTAP MRIO
T_large = [ [ GTAPupperleft, IMPtoIDuppersplit, GTAPupperright ] ;...
            [ EXPtoIDleftsplit, Treg, EXPtoIDrightsplit ] ;...
            [ GTAPlowerleft, IMPtoIDlowersplit, GTAPlowerright ] ] ;

VA_large = [ VA(:,1:(RegSecNumStart-1)), VAreg, VA(:,(RegSecNumEnd+1):7980) ] ;

FD_ROW_upperleft  = Y_total(1:(RegSecNumStart-1) , 1:(RegionNum-1)) ;
FD_ROW_upperright = Y_total(1:(RegSecNumStart-1) , (RegionNum+1):140) ;
FD_ROW_lowerleft  = Y_total((RegSecNumEnd+1):7980, 1:(RegionNum-1)) ;
FD_ROW_lowerright = Y_total((RegSecNumEnd+1):7980, (RegionNum+1):140) ;

FD_city_large = [ FDcitIMPupper; FDreg(:,1:4); FDcitIMPlower ] ;
FD_ron_large  = [ FDronIMPupper; FDreg(:,5:8); FDronIMPlower ] ;

FD_total_large = [ [ FD_ROW_upperleft, FDcitIMPupper, FDronIMPupper, FD_ROW_upperright ] ;...
                   [ EXPtoFDleftsplit, FDreg, EXPtoFDrightsplit ] ;...
                   [ FD_ROW_lowerleft, FDcitIMPlower, FDronIMPlower, FD_ROW_lowerright ] ] ;
               
FD_COICOP_large = [ [ FDcit_COICOP_IMPupper, FDron_COICOP_IMPupper ] ;...
                   [ FDreg_COICOP ] ;...
                   [ FDcit_COICOP_IMPlower, FDron_COICOP_IMPlower ] ] ;
                   
Xin_large = sum(T_large,1) + VA_large ;
Xout_large = sum(T_large,2) + sum(FD_total_large,2) ;
GTAPlargeIObal = Xout_large - Xin_large' ;

% disp('.'); 
% disp(['Check whether total X of original GTAP MRIO (' num2str(sum(sum(X))) ') is total Xin_large (' num2str(sum(Xin_large,2)) ') and total Xout_large (' num2str(sum(Xout_large,1)) ') of enlarged GTAP MRIO']) ;
% disp('.'); 
% disp(['Sum of large MRIO imbalance = ' num2str(sum(abs(GTAPlargeIObal))./(sum(sum(Xout_large+Xin_large'))/2)*100) '% of total output.']);
% disp(['Max of large MRIO imbalance = ' num2str(max(200*abs(GTAPlargeIObal)./(Xout_large+Xin_large'))) '% of total output.']); 
% disp('.'); toc ;
% disp('**********************************************************************')
% 


%% Carbon footprints of specific cities and rest of nations

% set small values (smaller than 1000$) in total input to zero to prevent large DIMs 
Xin_large = Xin_large  .* (abs(Xin_large) >1e-3) ;  

% Coefficient matrix of enlarged GTAP MRIO
A_large = T_large * diag(1./Xin_large); % Coefficient matrix A for MRIO. 
A_large(isinf(A_large))=0 ;  % Set inf in matrix A to zero
A_large(isnan(A_large))=0 ;  % Set nan in matrix A to zero
clear T_large ;  % not needed any longer

% assemble large emissions vector (currently only one row of CO2?
% emissions), row1 for CO2 and row5 for CO2e
E_ind_large = [ E_ind(1,1:(RegSecNumStart-1)), CITY_GTAP_F, E_ind(1,(RegSecNumEnd+1):7980) ] ; % this is the default case when no 'true' emissions are known
% E_ind_large = [ E_ind(1,1:(RegSecNumStart-1)), CITY_TRUE_F, E_ind(1,(RegSecNumEnd+1):7980) ] ; % use this if 'true' emissions are known!
    
% Calculate DIMs of enlarged GTAP MRIO
DIMsE_large  = E_ind_large * diag(1./Xin_large) ;  

% five rows of DIMs for five categories of industry Emissions in E_ind block
DIMsE_large(isinf(DIMsE_large))=0 ;  % Set inf in DIMs to zero
DIMsE_large(isnan(DIMsE_large))=0 ;  % Set nan in DIMs to zero
% csvwrite([resultspath 'DIMsE_large.csv'], DIMsE_large) ;
% csvwrite([resultspath 'DIMsE_CITY_RON.csv'], DIMsE_large(:,RegSecNumStart:(RegSecNumEnd+57))) ;

% Calculate L_large
L_large = inv(eye(size(A_large))-A_large) ;  % enlarged Leontief Inverse matrix
clear A_large ;  % not needed further in the script

% Calculate TIMs of enlarged GTAP MRIO
TIMsE_large = DIMsE_large * L_large ;
% csvwrite([resultspath 'TIMsE_large.csv'], TIMsE_large) ;
% csvwrite([resultspath 'TIMsE_CITY_RON.csv'], TIMsE_large(:,RegSecNumStart:(RegSecNumEnd+57))) ;

% CFs of all final demand
CF_CO2e_large  = repmat(TIMsE_large(1,:)',1,(140-1+8)) .* FD_total_large ;

% CFs of City by COICOP
CF_CO2e_COICOP_large  = repmat(TIMsE_large(1,:)',1,19) .* FD_COICOP_large ;

% Position of electricity sector in CCM fields (sector 43)
ELY_pos = zeros(57,57) ;
ELY_pos(43,43) = 1 ;

% Full City Carbon Map (CCM) and Scopes 1, 2 and 3
CCM_large = diag(DIMsE_large) * L_large * diag(sum(FD_city_large,2)) ;
CCM_large_ELY = CCM_large .* repmat(ELY_pos,141,141) ;

Scope1       = sum(CITY_GTAP_F(:,1:57)) ;
Scope1_ELY   = CITY_GTAP_F(:,43) ;
Scope2       = sum(CCM_large(RegSecNumStart+42,RegSecNumStart:RegSecNumEnd)) - CCM_large(RegSecNumStart+42,RegSecNumStart+42) + sum(sum(CCM_large_ELY)) ;
Scope1_in_CF = sum(sum(CCM_large(RegSecNumStart:RegSecNumEnd,:))) - sum(CCM_large(RegSecNumStart+42,RegSecNumStart:RegSecNumEnd)) + CCM_large(RegSecNumStart+42,RegSecNumStart+42) - sum(CCM_large_ELY(RegSecNumStart+42,:)) ;
Scope1_EXP   = Scope1_in_CF - Scope1 ;
Scope3       = sum(sum(CCM_large)) - Scope1_in_CF - Scope2 ;

% Explanation of COICOP format for final demand:
% MEL_hh_GTAP_Y ? this will be an n*12 matrix of Melbourne?s household expenditure by GTAP(down) and COICOP (across)
% MEL_oth_GTAP_Y ? this will be an n*7 matrix where the columns represent Mel gov, Mel cap, Mel vst, RoN hh, Ron gov, Ron cap, and Ron vst
% So Melbourne?s footprint will be shown broken down into 15 segments: 12 household COICOP products, gov, cap and vst.

% The dot product matrix sum of e*L .* SYD_GTAP_hhY to get a results matrix that is n*12. 
% This can then be summed down the columns to show impact by COICOP category (or sum part of the column to get domestic and imports).

% City and RON CF aggregated to 57x4
City_CF_CO2e_57 = repmat(eye(57),1,141) * CF_CO2e_large(:,RegionNum:(RegionNum+3)) ;
RON_CF_CO2e_57  = repmat(eye(57),1,141) * CF_CO2e_large(:,(RegionNum+4):(RegionNum+7)) ;

% City CF aggregated to 57x4, only City FD for city products (Stream 1+3a)
City_CF_CO2e_57_CITY_prod = CF_CO2e_large(RegSecNumStart:RegSecNumEnd,RegionNum:(RegionNum+3)) ;
% City CF aggregated to 57x4, only City FD for RON products
% (regional/national imports) (part of Stream 4a)
City_CF_CO2e_57_RON_prod = CF_CO2e_large((RegSecNumEnd+1):(RegSecNumEnd+57),RegionNum:(RegionNum+3)) ;
% City CF aggregated to 57x4, only City FD for ROW products (international imports) (part of Stream 4a)
City_CF_CO2e_57_ROW_prod = CF_CO2e_large(:,RegionNum:(RegionNum+3)) ;
City_CF_CO2e_57_ROW_prod((RegSecNumStart):(RegSecNumEnd+57),:) = zeros(114,4) ;
City_CF_CO2e_57_ROW_prod = repmat(eye(57),1,141) * City_CF_CO2e_57_ROW_prod ;

% City CF and RON by COICOP aggregated to 57x4
City_CF_CO2e_COICOP_57 = repmat(eye(57),1,141) * CF_CO2e_COICOP_large(:,1:15) ;
% RON_CF_CO2e_COICOP_57  = repmat(eye(57),1,141) * CF_CO2e_COICOP_large(:,16:19) ;  % This is not a COICOP breakdown (only 1+3 columns)
 
disp('.'); 
disp(['Check whether total CO2e footprint of City (' num2str(sum(sum(CF_CO2e_large(:,RegionNum:(RegionNum+3))))) ') plus RON (' num2str(sum(sum(CF_CO2e_large(:,(RegionNum+4):(RegionNum+7))))) ') is equal to national CF (' num2str(sum(sum(result(n,1)))) ')']) ; 
disp(['Check whether total CO2e footprint of City (' num2str(sum(sum(CF_CO2e_large(:,RegionNum:(RegionNum+3))))) ') is equal to CCM (' num2str(sum(sum(CCM_large))) ')']) ; 
disp(['Total CO2e footprint of City by COICOP = ' num2str(sum(sum(CF_CO2e_COICOP_large(:,1:15)))) ' or ' num2str(sum(sum(City_CF_CO2e_57))) ' kt CO2e']) ; 
disp(['Total Scope 1 of City                  = ' num2str(Scope1) ' kt CO2e']) ; 
disp(['Total Scope 1 in CF of City            = ' num2str(Scope1_in_CF) ' kt CO2e']) ; 
disp(['Total Scope 1 exported from City       = ' num2str(Scope1_EXP) ' kt CO2e']) ; 
disp(['City electricity sector emissions      = ' num2str(Scope1_ELY) ' kt CO2e']) ; 
disp(['Total Scope 2 of City                  = ' num2str(Scope2) ' kt CO2e']) ; 
disp(['Total Scope 3 of City                  = ' num2str(Scope3) ' kt CO2e']) ; 
disp(['Check whether Scope 3 (incl. electricity) (' num2str(sum(sum(CCM_large))-sum(sum(CCM_large(RegSecNumStart:RegSecNumEnd,:)))) ') is larger than total EEI to city (' num2str(sum(sum(City_CF_CO2e_57_RON_prod)) + sum(sum(City_CF_CO2e_57_ROW_prod))) '), because Scope 3 is larger than Stream 4a (by 3a)!']) ; 
disp(['Check whether Scope1(CF)+2+3 (' num2str(Scope1_in_CF + Scope2 + Scope3) ') is equal to City CF (' num2str(sum(sum(CF_CO2e_large(:, RegionNum   :(RegionNum+3))))) ')']) ; 
disp('.'); 
disp(['Total CO2e footprint of RoN            = ' num2str(sum(sum(CF_CO2e_large(:,(RegionNum+4):(RegionNum+7))))) ' kt CO2e']) ; 
% disp(['Total CO2e footprint of RoN  by COICOP = ' num2str(sum(sum(CF_CO2e_COICOP_large(:,16:19)))) ' kt CO2e']) ;
disp(['Total CO2e footprint of Country        = ' num2str(sum(sum(CF_CO2e_large(:, RegionNum   :(RegionNum+7))))) ' kt CO2e']) ; 

%csvwrite([resultspath 'City_CF_CO2e_large.csv'], CF_CO2e_large(:,RegionNum:(RegionNum+3))) ;
csvwrite([resultspath strcat(char(C_lookup{2}(n)), '_CF_CO2e57_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '.csv')], [City_CF_CO2e_57, zeros(57,1), City_CF_CO2e_57_CITY_prod, zeros(57,1), City_CF_CO2e_57_RON_prod, zeros(57,1), City_CF_CO2e_57_ROW_prod]) ;
csvwrite([resultspath strcat(char(C_lookup{2}(n)), '_CF_CO2e57_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '_COICOP.csv')], City_CF_CO2e_COICOP_57) ;

% Legend for City_CF_CO2e_57.csv is:
% first four columns are Total City CF aggregated to 57x4
% then one column of zeros
% next four columns are City CF, only City FD for city products (local goods and services)
% then one column of zeros
% next four columns are City CF, only City FD for RON products (regional/national imports)
% then one column of zeros
% next four columns are City CF, only City FD for ROW products (international imports)

% Generating summary results for table

 result(n,7) = sum(sum(City_CF_CO2e_57)) ; % city CBA after split, excl. HH direct emissions, kt CO2e
 result(n,8) = sum(sum(City_CF_CO2e_COICOP_57)) ; % city CF by COICOP, excl. HH direct emissions, kt CO2e
 result(n,9) = (result(n,4)+result(n,7))*1000/double(C_lookup{3}(n)); % city CF per capita after split and including HH direct emissions, in t CO2e/cap
 result(n,10)= (result(n,4)+result(n,8))*1000/double(C_lookup{3}(n)); % city CF per capita after split and including HH direct emissions, in t CO2e/cap
 result(n,11) = Scope1 ;       % Scope 1, excl. HH direct emissions
 result(n,12) = Scope1_EXP ;   % Scope 1 exported
 result(n,13) = Scope1_in_CF ; % Scope 1 remaining in CF, excl. HH direct emissions
 result(n,14) = Scope2 ;       % Scope 2
 result(n,15) = Scope3 ;       % Scope 3
 result(n,16) = Scope1_in_CF + Scope2 + Scope3 ; % total CF from scopes but without HH direct emissions
 result(n,17) = double(C_lookup{4}(n)) + double(C_lookup{5}(n)); % city direct emissions in kt CO2e (HH residential and private transport)
 result(n,18) = result(n,4)*1000/double(C_lookup{3}(n));  %  HH direct emissions per cap in t CO2e/cap
 result(n,19) = result(n,11)*1000/double(C_lookup{3}(n)); % Scope 1, excl. HH direct emissions, per cap
 result(n,20) = result(n,12)*1000/double(C_lookup{3}(n)); % Scope 1 exported, excl. HH direct emissions, per cap
 result(n,21) = result(n,13)*1000/double(C_lookup{3}(n)); % Scope 1 remaining in CF, excl. HH direct emissions, per cap
 result(n,22) = result(n,14)*1000/double(C_lookup{3}(n)); % Scope 2, excl. HH direct emissions, per cap
 result(n,23) = result(n,15)*1000/double(C_lookup{3}(n)); % Scope 3, excl. HH direct emissions, per cap
 result(n,24) = result(n,18)+result(n,21)+result(n,22)+result(n,23);  % Total city CF incl. HH direct emissions, per cap
 result(n,25) = sum(sum(RON_CF_CO2e_57)) ; % RON CBA after split, excl. HH direct emissions, kt CO2e
 result(n,26) = result(n,2) - result(n,4) ; % RON HH emissions are Country HH minus City HH emission, kt CO2e
 result(n,27) = (result(n,25)+result(n,26))*1000/double(C_lookup{7}(n)-C_lookup{3}(n)); % RON CF per capita after split and including HH direct emissions, in t CO2e/cap
   
%  result(n,11) = (strcat('T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'RAS', num2str(RAS_choice))) ;
%  result(n,12) = strcat(num2str(date))' ; % this and the next line lead to the same error message
%  result(n,12) = strcat(char(num2str(date)))' ; % this and the next line lead to the same error message
%  eval(['result(n,12) = strcat(num2str(date))']) ;

% fig = figure('Name',num2str(n),'NumberTitle','off');
% mesh(Treg) ;
%axis([0 114 0 114 0 1e5]) ;
% saveas(fig, [resultspath strcat(char(C_lookup{2}(n)), '_Treg_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date))]) ;
% saveas(fig, [resultspath strcat(char(C_lookup{2}(n)), '_Treg_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '.png')]) ;
% saveas(fig, [resultspath strcat(char(C_lookup{2}(n)), '_Treg_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '.jpg')]) ;
% saveas(fig, [resultspath strcat(char(C_lookup{2}(n)), '_Treg_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '.tif')]) ;
csvwrite([resultspath 'City_CF_results_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '.csv'], result) ; % if emissions from GTAP are used (this appears twice as needed here if validation is not done)
% csvwrite([resultspath 'City_CF_results_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'LQ', num2str(LQ), 'RAS', num2str(RAS_choice), '_', num2str(date), '_GTAPEm_CC', num2str(CurrencyConversion), '.csv'], result) ; % if emissions from GTAP are used (this appears twice as needed here if validation is not done)
% csvwrite([resultspath 'City_CF_results_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'RAS', num2str(RAS_choice), '_', num2str(date), '_TrueEm_CC', num2str(CurrencyConversion), '.csv'], result) ; % if true emissions are known (this appears twice as needed here if validation is not done)

SplitCrossSum = [Treg, (sum(EXPtoIDleftsplit,2)+sum(EXPtoIDrightsplit,2)), (sum(EXPtoFDleftsplit,2)+sum(EXPtoFDrightsplit,2)), FDreg;
                 (sum(IMPtoIDuppersplit,1)+sum(IMPtoIDlowersplit,1)), 0, 0, (sum(FDcitIMPupper,1)+sum(FDcitIMPlower,1)), (sum(FDronIMPupper,1)+sum(FDronIMPlower,1));
                 sum(VAreg,1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ;
                                 
% csvwrite([resultspath strcat(char(C_lookup{2}(n)), '_SplitCrossSum_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'RAS', num2str(RAS_choice), '_', num2str(date), '_GTAPEm_CC', num2str(CurrencyConversion), '.csv')], SplitCrossSum) ;  % if emissions from GTAP are used 
% csvwrite([resultspath strcat(char(C_lookup{2}(n)), '_SplitCrossSum_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'RAS', num2str(RAS_choice), '_', num2str(date), '_TrueEm_CC', num2str(CurrencyConversion), '.csv')], SplitCrossSum) ;  % if true emissions are known

% disp('loop time: '); toc ;
disp('**********************************************************************')

end

end

end
 
% csvwrite([resultspath 'City_CF_results_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'RAS', num2str(RAS_choice), '_', num2str(date), '_GTAPEm_CC', num2str(CurrencyConversion), '.csv'], result) ; % if emissions from GTAP are used
% csvwrite([resultspath 'City_CF_results_T', num2str(Treg_choice), 'Y', num2str(Treg_choice), 'RAS', num2str(RAS_choice), '_', num2str(date), '_TrueEm_CC', num2str(CurrencyConversion), '.csv'], result) ; % if true emissions are known
% 

disp('.'); toc ;
disp('**********************************************************************')

diary off

