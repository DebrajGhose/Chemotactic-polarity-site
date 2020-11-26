% Euler_step



% membrane
Cdc42Th    = Cdc42T    + dts_i.* ((k2a*BemGEF+k3*BemGEF42+k9*RecA).*Cdc42D - k2btotal.*Cdc42T - k4a*BemGEF.*Cdc42T + k4b*BemGEF42 - k7*BemGEFc.*Cdc42T);   % Added GDP->GTP  exchange by Reca (k9) %Masha dev sept2012

BemGEF42h     = BemGEF42     + dts_i.*(k4a*BemGEF.*Cdc42T - k4b*BemGEF42 + k7*BemGEFc.*Cdc42T);             

BemGEFh    = BemGEF    + dts_i.*(k1a*BemGEFc - k1b*BemGEF - k4a*BemGEF.*Cdc42T + k4b*BemGEF42);

Cdc42Dh    = Cdc42D    + dts_i.*(k5a*Cdc42Dc - k5b*Cdc42D + k2btotal.*Cdc42T - (k2a*BemGEF+k3*BemGEF42+k9*RecA).*Cdc42D ); % Added GDP->GTP  exchange by Reca (k9) %Masha dev sept2012

% cytoplasm

BemGEFch    = BemGEFc    + dts_i.*mean(mean((eta*(k1b*BemGEF - (k1a+k7*Cdc42T).*BemGEFc))));

Cdc42Dch  = Cdc42Dc  + dts_i.*mean(mean((eta*(k5b*Cdc42D  - k5a*Cdc42Dc) + ...
                     eta*(k5b*Cdc42Dic - k5a*Cdc42Dc   )*dx2*dx2/dx/dx )));

%internal compartment

Cdc42Dich = Cdc42Dic + dts_i.*(k5a*mean(mean(Cdc42Dc)) - k5b*Cdc42Dic);

%Masha dev sept2012
Recich = Recic;
Rech  = Rec   + dts_i.*(-k8a.*Rec + k8b*RecA); %free receptor on the membrane convertion to being bound to pheromone
RecAh  = RecA       + dts_i.*( k8a.*Rec - k8b*RecA); %receptor bound to pheromone time evolution

