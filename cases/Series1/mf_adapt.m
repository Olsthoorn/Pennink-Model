%% Pennink1 - density flow in Pennnink's (1905) sand box model
% Experiments series 1
%
% The computation takes about 50 sec on a 2 GHz mac
%
%see http://www.citg.tudelft.nl/live/pagina.jsp?id=68e12562-a4d2-489a-b82e-deca5dd32c42&lang=en
%
% in his sand-box experiment, Pennink studies freshwater flow between from
% rainwater recharge to aa extraction canal at the left.
% He then injects ink at 2 points at the top of the model
% and shows in a series of photos the movement of the ink.
%
% This model can be run by mt3dms as well as by seawat.
% To easier match the actual times of the tests, we will use time in minutes
% instead of days.
%
% TO 090312 100523 100722

clear variables; close all;

basename='PenninkSeries1';
save name basename

fprintf('Tests Pennink (1915), series 1, 2 Sep 1904\n');

GREP = 'STRESS PERIOD';

% Countour of sand mass obtained form photo by ginput
xSand =[
   -2.8568
   68.3082
   68.6354
   55.7112
   21.8465
   20.3741
   13.8302
    7.6135
   -1.2208
];

zSand =[
   -2.6501
   -3.0078
   65.3106
   58.6934
   57.6204
   56.3684
   50.8243
   44.7436
   44.2071
   ];

xCanL =[
   -2.3660
    5.6503
    6.7955
    8.4314
    9.2494
    9.2494
   -2.8568
   ];

zCanL =[
   37.4110
   37.5898
   43.1340
   53.1493
   59.4088
   63.8799
   64.5953
   ];

xCanR=[
   68.6354
   60.7827
   61.6007
   62.2551
   62.7459
   63.4003
   68.4718
];

zCanR =[
   65.3106
   62.4491
   51.1820
   46.1743
   43.3128
   41.5244
   41.3456
];


%% The grid, the box is 96x96 cm and 1.8 cm thick
FRESHWATER=0;    % Relative minimum concentration
SEAWATER  =1;    % Relative maximum concentration
k=86500/(24*60); % [cm/min] calibrated
peff=0.38;       % [-] effective porosity, calibrated

Scale=65/113;    % [cm/mm] foto p32
MW   =65;        % [cm] Width of model. Pennink p6
MH   =65;        % [cm]
zCap =52;        % [cm] Top of full capillary zone (see description)
hCanL=45.2;      % [cm] foto p32

%% Grid is always 3D
dx=1.0;          % [cm] grid cell width
dy=1.8;          % [cm] grid cell length = depth of model
dz=1.0;          % [cm] grid cell size vertical

xGr=0:dx:MW;     % [cm] grid
yGr=[0 dy];      % [cm] grid
zGr=MH:-dz:0;    % [cm] grid

gr = gridObj(xGr,yGr,zGr);

%% Point sources locations
xyzInk=[
    54.5 0 55
    38.5 0 55 ];

idx = xyzindex(xyzInk,gr);

%% Times of the photos in Pennink (1915) 2 sept 1904
%  these times are used in workbook, sheet PER
year=1904; month=5;  % day unknown, first day in known month assumed
TIMES=datenum([
    year,month,1, 9,22,0;   % ink added    year,month,1,10,19;   % opname p12
    year,month,1,10,39,0;   % p14
    year,month,1,11,09,0;   % p16
    year,month,1,11,12,0;   % some ink drops added to show flow path
    year,month,1,11,39,0;   % p18
]);

%% Model arrays

IBOUND = gr.const(0);

% use IBOUND as zone array
IBOUND(inpolyz(gr.XM,gr.ZM,xSand,zSand)) =  1;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanL,zCanL)) = -2;
IBOUND(inpolyz(gr.XM,gr.ZM,xCanR,zCanR)) =  3;  % not a fixed head

% dry above canal
IBOUND(IBOUND == -2  & gr.ZM>hCanL)=0;
IBOUND(IBOUND ==  3  & gr.ZM>hCanL)=0;

HK =  gr.const(k);
VK =  gr.const(k);

HK(HK>0 & gr.ZM>zCap)=k/10; % HK=0 above full capillary zone

% peff in saturated and unsaturation zone
PEFF = gr.const(peff);
PEFF( gr.ZM>zCap)=peff/3;           % unsaturated

ICBUND = IBOUND;

STRTHD = gr.const(hCanL);
STCONC = gr.const(FRESHWATER);

ink1 = 5; IBOUND(idx(1)) = ink1;
ink2 = 6; IBOUND(idx(2)) = ink2; 


zoneVals = {ink1 0 0;ink2 0 0};
zoneConc = {'C_Ink1';'C_ink2'};

[~,PNTSRC] = bcnZone(basename,'CCC',IBOUND,zoneVals,zoneConc);

%% RCH is 4.2 L/h, converted to cm3/min

W=45; % [cm] width of rain added to top of model
M=37; % [cm] center of rain added to top of model (mid between 2 screws see photo)
N=4.2*1e3/60/W/dy;  % [cm/min]

NPER = getPeriods(basename);

RECH=zeros(gr.Ny,gr.Nx);                  % specify recharge for all stress periods
RECH(:,gr.Xm>=M-0.5*W & gr.Xm<=M+0.5*W,:)=N;   % only where rain is applied in model
RECH = bsxfun(@times,RECH,YS(1:NPER));

save underneath xSand zSand xCanL zCanL % needed to for PSIMASK in mf_analyze