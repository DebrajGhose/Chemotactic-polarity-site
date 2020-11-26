
%% Initialize storage for the output.
clear 'c42cen' 'recacen' 'k42_store' 'kra_store' 'endo_cell' 'cable' 'probv_endo';
totcnt = 1;
curr_t = 0;
fpt_flag=0; 
tsi=0;
Fwarp = zeros(N+1); % storage for use in interpolation that needs wrap around func values
c42cen=[];recacen=[]; 
endo_cell        = struct('endo_time',[],'pick_time',[],'vSNARE',[],'Rec',[],'RecA',[],'Cdc42T',[], ...
    'Cdc42D',[], 'BemGEF42',[], 'BemGEF',[], 'posx',[],'posy',[],'dx',[],'num',0);
% Cdc42_store   = fopen(['data_Cdc42T_time_course_run' num2str(runi) '.txt'],    'w'); 
% RecA_store    = fopen(['data_RecGEF_time_course_run' num2str(runi) '.txt'],    'w');
% Rec_store     = fopen(['data_Rec_time_course_run' num2str(runi) '.txt'],    'w');
% neg_store     = fopen(['data_negative_run' num2str(runi) '.txt'],              'w'); % going to hold the time points and protein when it goes negative
mass_store    = fopen(['data_mass_run' num2str(runi) '.txt'],                  'w'); % going to hold total mass
% endo_loc      = fopen(['data_endopatches_run' num2str(runi) '.txt'],           'w');
% vsnare_store  = fopen(['data_vsnare_time_course_run' num2str(runi) '.txt'],    'w');
%cable_loc     = fopen(['data_cable_location_run' num2str(runi) '.txt'],        'w');
Internalcompartment = 0.7;                  % relative size (area) of internal compartment vs plasma membrane 
dx  = Cellsize*sqrt(pi)/N;                  % spatial meshsize (1D)
dx2 = sqrt(Internalcompartment)*dx;
eta0=0.01;                     %eta: Vm/Vc, membrane/cytoplasm volume correction
eta =eta0;                     % notes in '08_endocytic_vesicle_30Dec2010.doc'
mem_depth  = (Cellsize/2)*(1-(1/((eta0+1)^(1/3))));
cyt_mult   = ((Cellsize/2)^3)/(eta0+1);
Rnew_mult  = N/2/(pi^0.5);

%% initial concentrations

load('Init_AllProteins_Pher_10nM.mat');

%% Load the desired initial position and translate all membrane concentration profiles (with pbc wrap)
% % ee=1:N;
% % [Y, X]=meshgrid(ee);
% % Rdist=sqrt((X-N/2).^2 +(N/2-Y).^2);
% % [indxs,indys]=find((Rdist>N/2-1)&(Rdist<N/2+1));        %N/2 away, but in a random spot
% % i0=randi([1 length(indxs)]);
% % xnew0=indxs(i0);ynew0=indys(i0);
% % x=importdata('input.dat'); 
% % xnew0=x(1);
% % ynew0=x(2);
% xnew0=N/2;ynew0=N/2;
% %display(['Randomizing the initial peak position: (x0,y0) = (' num2str(xnew0) ',' num2str(ynew0) ')']);
% display(['   Shifting the peak position to (x0,y0) = (' num2str(xnew0) ',' num2str(ynew0) ')']);
% [indx,indy]=find(Cdc42T==max(max(Cdc42T)));
% if (length(indx)>1) indx=indx(1); end
% if (length(indy)>1) indy=indy(1); end
% delx=(xnew0-indx);
% dely=(ynew0-indy);
% cdct=Cdc42T;cdcd=Cdc42D; gdicdc=GDI42;
% bem=BemGEF; bem42=BemGEF42;vsn=vSNARE;reca=RecA; rec=Rec;newindx=0;newindy=0;
% for i=1:N
%     for j=1:N
%         if (i+delx<0)
%             newindx=N+i+delx;
%         elseif (i+delx==0 | i+delx==N)
%             newindx=N;
%         else
%             newindx=mod(i+delx,N);
%         end
%         if (j+dely<0)
%             newindy=N+j+dely;
%         elseif (j+dely==0 | j+dely==N)
%             newindy=N;
%         else
%             newindy=mod(j+dely,N);
%         end
%         Cdc42T(  newindx,newindy)=cdct(  i,j);
%         Cdc42D(  newindx,newindy)=cdcd(  i,j);
%         GDI42(   newindx,newindy)=gdicdc(i,j);
%         BemGEF(  newindx,newindy)=bem(   i,j);
%         BemGEF42(newindx,newindy)=bem42( i,j);
%         vSNARE(  newindx,newindy)=vsn(   i,j);
%         RecA(    newindx,newindy)=reca(  i,j);
%         Rec(     newindx,newindy)=rec(   i,j);
%     end
% end

%% Calclate the global mass __________________________________________
Cdc42s = Cdc42Dc ...                                    % Cdc42D cytoplasmic content 
+ eta*(mean(mean(Cdc42D+Cdc42T+BemGEF42))) ...   % Cdc42 membrane content 
+ eta*(Cdc42Dic*(dx2^2)/(dx^2)) ...                    % Cdc42D on inner membrane           
+ eta*sum(endo_cell.Cdc42T(1:endo_cell.num)   .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2) ...     % Cdc42T in vesicle                            
+ eta*sum(endo_cell.Cdc42D(1:endo_cell.num)   .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2) ...     % Cdc42D in vesicle 
+ eta*sum(endo_cell.BemGEF42(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);        % BemGEF42 in vesicle 

Bem1s = BemGEFc ...                                    % cytoplasmic content scaled to outer membrane
+ eta*(mean(mean(BemGEF+BemGEF42))) ...                % Bem1 membrane content  
+ eta*sum(endo_cell.BemGEF(1:endo_cell.num)   .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2) ...     % Bem1 vesicle content 
+ eta*sum(endo_cell.BemGEF42(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);        % BemGEF42 vesicle content                                                        

vSNAREs        = eta*mean(mean(vSNARE)) ...              % v-SNARE membrane content
+ eta*(vSNAREic*(dx2^2)/(dx^2)) ...                   % v-SNARE on inner membrane content 
+ eta*sum(endo_cell.vSNARE(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);          % v-SNARE vesicle content 

Recs        = eta*mean(mean(Rec)) ...              % Receptor complex membrane content %Masha dev oct2012
+ eta*(Recic*(dx2^2)/(dx^2)) ...                   % Receptor complex on inner membrane content 
+ eta*sum(endo_cell.Rec(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);          % Receptor complex vesicle content 

RecAs        = eta*mean(mean(RecA)) ...              % Receptor complex membrane content %Masha dev oct2012 
+ eta*sum(endo_cell.RecA(1:endo_cell.num) .* (endo_cell.dx(1:endo_cell.num).^2))/(dx^2)/(N^2);          % Receptor complex vesicle content 
fprintf(mass_store,'%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n',0,Cdc42s, Bem1s,vSNAREs,Recs,RecAs);

%% Clean vesicle containers endo and exo
cable       = [];       % hold upto MAXcable cable co-ordinates
% record_vSNARE = [];     % record [v-SNARE] at time of endocytosis
% record_endot = [];      % record time of endocytosis
% record_lifetime = [];   % record lifetime of each sink
% record_Rec = [];     % record [Rec] at time of endocytosis
% record_RecA = [];     % record pheromone bound receptor [RecA] at time of endocytosis

%exocytic thresholds for hill distribution
if(max(max(RecA))~=0)     k1=0.5*max(max(RecA)); else    k1=1; end %exponent and threshold for RecA, threshold shosen as a half of <Rec> in the absense of pheromone 1.5 microMolar
k2=0.5*max(max(Cdc42T+BemGEF42)); %about a 1/2 of the maximum Cdc42T concentration in a well formed peak in our simulations
k42_store=zeros([1 ceil(sim/tsave)+1]);
k42_store(1)=k2;
% kra_store=zeros([1 ceil(sim/tsave)+1]);
% kra_store(1)=k1;
if(max(max(vSNARE))~=0)     kv=0.5*max(max(vSNARE)); else    kv=1; end 
kvsn_store=zeros([1 ceil(sim/tsave)+1]);
kvsn_store(1)=kv;

%% report at the end
display(['--> Initialization lasted ' num2str(toc(tot_runi_start)) ' (s). Starting simulations...']); %Masha debug
