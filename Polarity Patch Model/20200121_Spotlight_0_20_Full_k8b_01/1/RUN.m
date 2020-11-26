close all;
clear all;
tic
tot_start=tic;
global k1a
global k1b
global k2a
global k2b
global k3
global k4a
global k4b
global k5a
global k5b
global k6a
global k6b
global k7
global eta
global dx
global dx2
global k8a                      %rate of pheromone binding to receptor complex
global k8b                      %rate of pheromone unbinding to receptor complex
global k9                       %rate of nucleotide exchange in Cdc42D->T assisted by liganded receptor RecA

%Use a random seed

load('seed.mat');
rng(seed);

%% variables for time-stepping
dt    = 0.05;                               % diffusion time step in seconds
sim   = 1.0*60*60.0;                            % simulation time in (s)
display(['Simulation time is ' num2str(sim/60) '(min), with a time step ' num2str(dt) ' (s)']); %Masha debug
simdt = sim/dt;                             % number of dt time steps
tsave = 5;                                 % data is going to be saved every tsave seconds
display(['Data is going to be saved every ' num2str(tsave) ' (s)']);

%% variables for spatial discretization
N = 100;                                    % number of spatial grid points
Cellsize = 5;                               % sphere diameter (microns)
Internalcompartment = 0.7;                  % relative size (area) of internal compartment vs plasma membrane
dx  = Cellsize*sqrt(pi)/N;                  % spatial meshsize (1D)
dx2 = sqrt(Internalcompartment)*dx;         % spatial meshsize for internal compartment (1D)

%% initalising diffusion mechanism
nRsteps = 100;                              % take nRsteps reaction steps per diffusion step
dt2 = dt/nRsteps;                           % reaction time-step
totcnt = 1;

%% initial concentrations
if(exist('Cdc42T_end1.mat','file')==2)      %check if we have run pre-equilibration and can load equilibrated concentraitons
    version='1';
else version='0';
end

%% Pheromone constant concentration profile

Phermin = 0.000; Phermax = 0.020; %Pheromone conceentration profile max and min in uMolar

Pher = mean([Phermax,Phermin])*ones(100,100); %uniform pheromone concentration; also helps to initialize Pher matrix if you are using a gradient

%Initialize parameters for pheromone gradient

Pherxdir = linspace(Phermin,Phermax,N);
[PherCol , PherRows] = meshgrid([1:N]);
Pherslope = (Phermax - Phermin)/(N - 1);
Pheroffset = Pherxdir(10) - Pherslope*10; %model pheromone gradient as Pher(x) = mx + c


display(['Average pheromone concentration is ' num2str(mean(mean(Pher))) ' (uM)']);
%save(['Pher_end_' version '.mat'],'Pher');

%% reaction constants
mult  = 16;

k1a = 10;                       %s-1,       BemGEFc -> BemGEF
k1b = 10;                       %s-1,       BemGEF -> BemGEFc

% GEF
k2a = 0.16;                     %uM-1.s-1,  BemGEF + Cdc42D -> Cdc42T
k3  = 0.35;                     %uM-1.s-1,  BemGEF42 + Cdc42D -> Cdc42T
%GAP
k2b = 0.63;                     %s-1,       Cdc42T -> Cdc42D
k2bvesc = zeros(N,N);           %s-1,       Cdc42T -> Cdc42D due to vesicles 
k2bdeliv = 0.3;                 %s-1,       Cdc42T -> Cdc42D GAP activity delivered by exocytic vesicle
k2bhalflife = 10;             %s, half life

k4a = 10;                       %uM-1.s-1,  BemGEF + Cdc42T -> BemGEF42
k4b = 10;                       %s-1,       BemGEF42 -> BemGEF + Cdc42T
k7  = 10;                       %uM-1.s-1,  BemGEFc + Cdc42T -> BemGEF42

% Cdc42D on
k5a = mult*9;                   %s-1,       RDI42c -> RDI42
%k6b = mult*5;                   %s-1,       RDI42 -> Cdc42D + RDIc
% Cdc42D off
%k6a = mult*15;                  %uM-1.s-1,  Cdc42D + RDIc -> RDI42
k5b = mult*1.3;                 %s-1,       RDI42 -> RDI42c

%Rec + alpha <-> RecA %Masha dev oct2012
KD8=6*10^(-3);                  % KD of pheromone binding reaction Raths et al. 1988 PMID 2846561 gives 20 nM, Jenness, D.D., Burkholder, A.C., and Hartwell, L.H. (1986) Mol. Cell. Biol. 6, 318-320 give 6 nM. Take an average btwn the two, 10 nM.
k8b=0.01;                      %Raths et al. 1988 PMID 2846561
k8a=(k8b/KD8)*Pher;             %Keep KD constant in case individual rates need to be varied

%Cdc42D + RecA ->(k9)-> Cdc42T + RecA
k9=0.08;
if(mean(mean(Pher))==0)         %fix for cases when we simulate 0 pheromone. the interpolation does not work well with zero profiles
    Pher=ones(N).*10^(-4);
    k9=0;
    k8a=(k8b/KD8)*Pher; 
end
display(['Rate of Cdc42D->Cdc42T, facilitated by RecA, k9 = ' num2str(k9)]);

%% define parameters for actin cable attachment by gaussian around actin cable centroid

actcolcoord = repmat(1:N,[ N , 1 ]);
actrowcoord = repmat([1:N]',[ 1 , N ]);
actinsd = 0.24/dx; % 0.24um is the standard deviation of 
%Spa2 patch -- 20190113_SETMAWP_Treadmill\IN\MeasureGaussianAroundActinCable
%but is subject to change

k_gauss = 1.5; %strength of Gaussian
%% conversion between 3D cytosol and 2D membrane
eta0=0.01;                      %eta: Vm/Vc, membrane/cytoplasm volume correction
eta =eta0;                      % notes in '08_endocytic_vesicle_30Dec2010.doc'
mem_depth  = (Cellsize/2)*(1-(1/((eta0+1)^(1/3))));
cyt_mult   = ((Cellsize/2)^3)/(eta0+1);
Rnew_mult  = N/2/(pi^0.5);

%% compute diffusion-coeff times dt/dx^2
Dconst = 0.0045;                % diffusion coefficient
Diff1  = Dconst * dt/dx^2;      % diffusion matrix, includes dt/dx^2 factor
e        = ones(N,1);
Dxx      = spdiags([e, -2*e, e], -1:1, N, N);
Dxx(1,N) = 1;                   % periodic boundary conditions
Dxx(N,1) = 1;
Dxx      = Dxx*Diff1;

Mxx = kron(speye(N),Dxx);
Myy = kron(Dxx,speye(N));
Hop = speye(N^2) - Mxx - Myy;   % heat operator, I-Dxx

%% compute separate diffusion-coeff time dt/dx^2 for receptors

Dconstrec = 0.0005;                % diffusion coefficient
Diff1rec  = Dconstrec * dt/dx^2;      % diffusion matrix, includes dt/dx^2 factor
erec        = ones(N,1);
Dxxrec      = spdiags([erec, -2*erec, erec], -1:1, N, N);
Dxxrec(1,N) = 1;                   % periodic boundary conditions
Dxxrec(N,1) = 1;
Dxxrec      = Dxxrec*Diff1rec;

Mxxrec = kron(speye(N),Dxxrec);
Myyrec = kron(Dxxrec,speye(N));
Hoprec = speye(N^2) - Mxxrec - Myyrec;   % heat operator, I-Dxx

%% computer separate diffusion-coeff time dt/dx^2 for vsnares

Dconstvsnare = 0.0025;                % diffusion coefficient
Diff1vsnare  = Dconstvsnare * dt/dx^2;      % diffusion matrix, includes dt/dx^2 factor
evsnare        = ones(N,1);
Dxxvsnare      = spdiags([evsnare, -2*evsnare, evsnare], -1:1, N, N);
Dxxvsnare(1,N) = 1;                   % periodic boundary conditions
Dxxvsnare(N,1) = 1;
Dxxvsnare      = Dxxvsnare*Diff1vsnare;

Mxxvsnare = kron(speye(N),Dxxvsnare);
Myyvsnare = kron(Dxxvsnare,speye(N));
Hopvsnare = speye(N^2) - Mxxvsnare - Myyvsnare;   % heat operator, I-Dxx

%% initialize vesicle mechanism
% Parameters related to vesicle fusion/fission
cable       = [];               % hold upto MAXcable cable co-ordinates
MAXcable    = 10;               % maximum number of cables
detatchPROB = dt/(60.);         % (time step)/(residence time),	if time step < residence time
attachPROB  = dt;               % time step, if time step < 1, (prob 1 if time step = 1sec)

nVInlen  = N/50;                % no. of bins (1D) used to represent an exocytic
nVOutlen = N/100;               % no. of bins (1D) used to represent an endocytic  vesicle

LotHiV = 100;                % time to switch from low to high vesicle activity
HitLoV = 20;                % time to switch from high to low vesicle activity
avgvrate = 2/0.6;                % mean period for endocytosis, sec (= 15 sec actin patch/25 patches in cell)
lowvrate = avgvrate*0.5;
highvrate = ( avgvrate*(LotHiV + HitLoV) - lowvrate*LotHiV )/HitLoV; % calculate high vesicle rate

highperiodOut = 1/lowvrate;            % high period, low rate of endocytosis
lowperiodOut = 1/highvrate;           % low period, high rate of endocytosis

periodOut = lowperiodOut;       % period for endocytosis, mean is 0.3
periodIn = 4*periodOut;         % mean period for exocytosis, sec
patch_fill = 20;                % fill level for endocytosis, used for sink_case =2
patch_fill_give_up = 24;        % even if the patch doesn't fill up with cargo,

endo_red_time = 15;             % 15 s before endocytosis when spots are "red"
exo_cnt = 0;                    % count exocytic events

%% Defining hill coefficients for Cdc42 probability and ~RecA for cable attachement
n1=4;                   %exponent and threshold for RecGEF
n2=4;                   %exponent and threshold for Cdc42T
nv=3;                   %softer power dependence for vSNARE hill function
%display([' vSNARE nv = ' num2str(nv) ', Cdc42 n2 = ' num2str(n2)]);

%% Loop over restarts
Nruns=1;
display(['This node will run ' num2str(Nruns) ' ' num2str(sim/60/60) ' (hr) long simulations.']);
%tau=-1*ones([1 Nruns]);
tot_Nruns_start=tic;

%%
for runi=1:Nruns
    tot_runi_start=tic;
    
    initialize_system_storage;
    
    disp('Entering simulation loop');
    %% start the i run
    for loop=1:simdt                % looping over dt time steps
        
        
        %%
        curr_t = loop*dt;
        
           %% decide how quickly vesicle activity will happen
        
        die = rand();
        if die < dt/LotHiV && periodOut == highperiodOut %switch from low to high rate
            periodOut = lowperiodOut; %higher rate
            periodIn = 4*periodOut;
        end
        die = rand();
        if die < dt/HitLoV && periodOut == lowperiodOut %switch from high to low rate
            periodOut = highperiodOut; %lower rate
            periodIn = 4*periodOut;
        end
        
        %% Add GAP by vesicles to the mix
        
        k2bvesc = k2bvesc - dt*log(2)/k2bhalflife*(k2bvesc - 0.01 ); %setting rate constant = log(2)/halflife. The 0.01 is so that this never becomes zero or it will mess with interpolation.
        k2btotal = k2b + k2bvesc; %add GAP activity and vesicles
        
        %%
        Reaction_step;
        
        Cargo_and_diffusion_step;
        
        if(mod(curr_t,900)==0) % progress report
            display(['    Run #' num2str(runi) ': Current time step: ',num2str(curr_t/60) ' (min)']);
        end
        
        if(mod(curr_t , tsave*6)==0 && curr_t>900)  % take a snap shot every tsave seconds
        
            %% write down some observables
            
            if ~exist('Rec_store','var'), Rec_store = []; end
            Rec_store = cat(3, Rec_store ,  Rec ); 
             
            if ~exist('RecA_store','var'), RecA_store = []; end
            RecA_store = cat(3, RecA_store ,  RecA ); 
             
            if ~exist('Total_42T_store','var'), Total_42T_store = []; end
            Total_42T_store = cat(3, Total_42T_store ,  (Cdc42T+BemGEF42) ); %I am taking total Cdc42T 
            
            if ~exist('Pher_store','var'), Pher_store = []; end
            Pher_store = cat(3, Pher_store,  (Pher) ); %I am taking total Cdc42T 
            
        end
        
        if(mod(curr_t,tsave)==0 && curr_t>900) % take a snap shot every tsave seconds
            
            %% calculate and write down observables
           
            record_cdc42t_reca_centers;
            
            %recenter gradient to where the patch is
            
            % spotlight patch
            
            if ~exist('ploc','var'), ploc = 0; end %generate variable ploc
            
            if numel(c42cen(:,1)) > 1 %make sure you have at least two points recorded before attempting to recreate the gradient
                
                cen_diff = c42cen(end,2) - c42cen(end-1,2); %difference in distance travelend by patch along x axis
                
                %account for when patch goes over the edge
                
                if cen_diff < -50, cen_diff = 100 + cen_diff; %when patch goes over 100th pixel
                elseif cen_diff > 50 , cen_diff = - (100 - cen_diff); %when patch goes below 1st pixel
                end
                
                ploc = ploc + cen_diff; %ploc keeps track of how much the patch has moved (this gets around periodic boundary conditions)
                
            end
            
            Colmov = PherCol - ( c42cen(end,2)  - 50.5 ); Colmov = mod(  Colmov  , N  ) ; % the centroid-finding code is set up to assign values from 0 to 100
            Pher = Pherslope*Colmov + Pheroffset + ploc*Pherslope; % the ploc*Pherslope is what makes this the spottlight paradigm
            Pher(Pher<0) = 0; % anything that becomes negative is zero
            
            %Rec + alpha <-> RecA %Masha dev oct2012
            KD8=6*10^(-3);                  % KD of pheromone binding reaction Raths et al. 1988 PMID 2846561 gives 20 nM, Jenness, D.D., Burkholder, A.C., and Hartwell, L.H. (1986) Mol. Cell. Biol. 6, 318-320 give 6 nM. Take an average btwn the two, 10 nM.
            k8b=0.01;                      %Raths et al. 1988 PMID 2846561
            k8a=(k8b/KD8)*Pher;             %Keep KD constant in case individual rates need to be varied
            
            %Cdc42D + RecA ->(k9)-> Cdc42T + RecA
            k9=0.08;
            if(mean(mean(Pher))==0)         %fix for cases when we simulate 0 pheromone. the interpolation does not work well with zero profiles
                Pher=ones(N).*10^(-4);
                k9=0;
                k8a=(k8b/KD8)*Pher; 
            end
            
        end
    end
   
 
end
save('AllProteins_Pher_10nM');
%save('CheckProfiles.mat','Pher','Cdc42T','Cdc42D','BemGEF42','BemGEF','BemGEFc','Cdc42Dc','Cdc42Dic','Rec','Recic','RecA' , 'Rec_store' , 'RecA_store' , 'Total_42T_store' , 'Pher_store'); %save data for verification
display(['--> Total time for ' num2str(Nruns) ' runs = ' num2str(toc(tot_Nruns_start)/60/60) ' (hr).']);
