function [t_values,m_values,g_values,j_values,Rcom_values,Rsid_values,Riron_values,dmdt]...
    = simulation_community_interaction4(num_strain,num_sid,a,m0,Riron_supply,t_end,h,barcode_syn,barcode_rec)
% a general function to simulate the dynamics of a community
% interaction based on bar code
% change ep k1 k2 to alternate the steady state sid value
% alternate the step size
% from simulation_community_interaction3_tune.m
% reduce the memory usage

arguments
    num_strain = 2 % number of strains
    num_sid = 2 % number of sid
    a = 0.5*ones(2,num_strain) % strategies of strains
    m0=ones(1,num_strain) % initial mass of strains
    Riron_supply=1 % iron supplied
    t_end=200 % simulation range
    h = 0.05 % step size
    barcode_syn=[1,0;...
            0,1]
    barcode_rec=[1,0;...
            0,1]
end

t0 = 0;             % Initial time
d = 0.01;            % dilution rate
vm = 0.1;             % rate cofficient for complex influxes
Km = 1;             % affinity cofficient
k1 = 1e-2;           % iron-sid association
k2 = 1e-7;             % iron-sid dissociation
gm = 20;            % growth coefficient: gamma
ep = 0.2;             % production coefficient: epsilon
r = 1;              % biomass per cell volume

% Initialize arrays to store the results
t_values = (t0:h:t_end).';
t_record = (t0:t_end).';
t_record_index = 0;
m_values = zeros(size(t_record,1)-1,num_strain); % biomass, each column: one strain
g_values = zeros(size(t_record,1)-1,num_strain); % growth rate
j_values = zeros(size(t_record,1)-1,num_strain);          % fluxes of sid
Rcom_values = zeros(size(t_record,1)-1,num_sid); % complex
Rsid_values = zeros(size(t_record,1)-1,num_sid); % siderophore
Riron_values = zeros(size(t_record,1)-1,1);      % iron

m_values(1,:) = m0;
Riron_values(1) = Riron_supply;

m = m_values(1,:); % biomass, each column: one strain
g = g_values(1,:); % growth rate
j = j_values(1,:);          % fluxes of sid
Rcom = Rcom_values(1,:); % complex
Rsid = Rsid_values(1,:); % siderophore
Riron = Riron_values(1,:);      % iron

for i = 1:length(t_values) - 1

%     if i == 458
%         qzy=1
%     end

    % iterate biomass m
    j = vm.*Rcom*barcode_rec; % fluxes of sid
    g = gm.*a(1,:).*j; % growth rate
    
    dmdt = m.*(g-d).*h;
    m = m + dmdt;

    % iterate siderophore Rsid
    dRsid_dt = (- d.*Rsid + ep.*a(2,:).*m*(barcode_syn.')./r ...
        - k1.*Riron.*Rsid + k2.*Rcom) .*h;
    Rsid = Rsid + dRsid_dt;
    

    % iterate complex Rcom
    dRcom_dt = (- d.*Rcom - m./r.*vm*(barcode_rec.').*Rcom ...
        + k1.*Rsid.*Riron - k2.*Rcom).*h;
    Rcom = Rcom + dRcom_dt;
    

    % iterate iron Riron
    dRiron_dt = (d*(Riron_supply-Riron) ...
        +sum(-k1*Rsid*Riron + k2*Rcom,2)).*h;
    Riron = Riron + dRiron_dt;
    

    if any(Rcom<0)
        qzy =1
    end

    if t_values(i)>=t_record(t_record_index+1)
        j_values(t_record(t_record_index+1)+1,:) = j;
        g_values(t_record(t_record_index+1)+1,:) = g;
        m_values(t_record(t_record_index+1)+1,:) = m;
        Rsid_values(t_record(t_record_index+1)+1,:) = Rsid;
        Rcom_values(t_record(t_record_index+1)+1,:) = Rcom;
        Riron_values(t_record(t_record_index+1)+1) = Riron;
        t_record_index = t_record_index+1;
    end

end

end