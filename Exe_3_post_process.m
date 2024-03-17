clear all; close all; clc;

% structures initialization
RHO = struct("REG_out", zeros(3,2), "REG_in", zeros(3,2), "REG_in_350K", zeros(3,2), "REG_in_350K_Sab", zeros(3,2));
BETA_eff = struct("REG_out", zeros(3,2), "REG_in", zeros(3,2), "REG_in_350K", zeros(3,2), "REG_in_350K_Sab", zeros(3,2));
K_eff = struct("REG_out", zeros(3,2), "REG_in", zeros(3,2), "REG_in_350K", zeros(3,2), "REG_in_350K_Sab", zeros(3,2));

% files and scenario to analyze
filename = {"REG_out\Triga_res.m"; "REG_in\Triga_res.m"};
scenario = {'REG_out'; 'REG_in'};

for i = 1:length(filename)

    run(filename{i})

    % rho, k, beta = 3x2 --> 1° row: (-), 2° row: (pcm), 3° row: $
    rho = zeros(3,2);
    beta_eff = zeros(3,2);
    k_eff = zeros(3,2);
    
    % beta effective
    beta_eff(1,:) = ADJ_MEULEKAMP_BETA_EFF(1,1:2); 
    beta_eff(1,2) = beta_eff(1,2) * beta_eff(1,1);
    beta_eff(2,:) = beta_eff(1,:) * 1e5; % pcm
    
    disp(['In the case of ' scenario{i} ':'])
    disp(['Beta effective is: ' num2str(beta_eff(2,1)) ' +/- ' num2str(beta_eff(2,2)) ' pcm'])
    
    % k effective 
    k_eff(1,:) = IMP_KEFF(1,:);
    k_eff(1,2) = k_eff(1,2) * k_eff(1,1);
    
    % reactivity
    rho(1,1) = (k_eff(1,1)-1) / k_eff(1,1);
    rho(1,2) = k_eff(1,2) / k_eff(1,1)^2;
    
    rho(2,:) = rho(1,:) * 1e5; % pcm
    
    rho(3,1) = rho(1,1) / beta_eff(1,1); % $
    rho(3,2) = sqrt((rho(1,2)/beta_eff(1,1)).^2 + (-rho(1,1)/beta_eff(1,1).^2 * beta_eff(1,2)).^2); % $
    
    disp(['k effective is: ' num2str(k_eff(1,1)) ' +/- ' num2str(k_eff(1,2))])
    disp(['Reactivity is: ' num2str(rho(2,1)) ' +/- ' num2str(rho(2,2)) ' pcm'])
    disp(['Reactivity is: ' num2str(rho(3,1)) ' +/- ' num2str(rho(3,2)) ' $'])
    
    if scenario{i} == "REG_out"
        BETA_eff.REG_out = beta_eff;
        RHO.REG_out = rho;
        K_eff.REG_out = k_eff;
    elseif scenario{i} == "REG_in"
        BETA_eff.REG_in = beta_eff;
        RHO.REG_in = rho;
        K_eff.REG_in = k_eff;
    end

    clearvars -except filename scenario BETA_eff RHO K_eff
    fprintf('\n')
end

fprintf('\n')




%%
filename = {"REG_in\Triga_res.m"; "Doppler\Triga_res.m"; "Doppler_Sab\Triga_res.m"};
scenario = {'Environmental temperature (300 K)'; 'Higher temperature (350 K) without S(a,b)'; 'Higher temperature (350 K) considering S(a,b) variation'};

for i = 1:length(filename)

    run(filename{i})

    % rho, k, beta = 3x2 --> 1° row: (-), 2° row: (pcm), 3° row: $
    rho = zeros(3,2);
    k_eff = zeros(3,2);
    
    disp(['In the case of ' scenario{i} ':'])
    
    % k effective 
    k_eff(1,:) = IMP_KEFF(1,:);
    k_eff(1,2) = k_eff(1,2) * k_eff(1,1);
    
    % reactivity
    rho(1,1) = (k_eff(1,1)-1) / k_eff(1,1);
    rho(1,2) = k_eff(1,2) / k_eff(1,1)^2;
    
    rho(2,:) = rho(1,:) * 1e5; % pcm
    
    disp(['k effective is: ' num2str(k_eff(1,1)) ' +/- ' num2str(k_eff(1,2))])
    disp(['Reactivity is: ' num2str(rho(2,1)) ' +/- ' num2str(rho(2,2)) ' pcm'])

    fprintf('\n')

    if i==2
        RHO.REG_in_350K = rho;
        K_eff.REG_in_350K = k_eff;
    elseif i==3
        RHO.REG_in_350K_Sab = rho;
        K_eff.REG_in_350K_Sab = k_eff;
    end

    clearvars -except filename scenario BETA_eff RHO K_eff
    fprintf('\n\n')

end

CRW_REG(1) = RHO.REG_out(3,1) - RHO.REG_in(3,1);              % $
CRW_REG(2) = sqrt( RHO.REG_out(3,2)^2 + RHO.REG_in(3,2)^2 );  % $

delta_temp = 350 - 300; % K

alpha_fuel_350(1) = (RHO.REG_in_350K(2,1) - RHO.REG_in(2,1)) / delta_temp;               % pcm/°C
alpha_fuel_350(2) = sqrt( RHO.REG_in_350K(2,2)^2 + RHO.REG_in(2,2)^2 ) / delta_temp;     % pcm/°C

alpha_fuel_350_Sab(1) = (RHO.REG_in_350K_Sab(2,1) - RHO.REG_in(2,1)) / delta_temp;               % pcm/°C
alpha_fuel_350_Sab(2) = sqrt( RHO.REG_in_350K_Sab(2,2)^2 + RHO.REG_in(2,2)^2 ) / delta_temp;     % pcm/°C

disp(['The prompt fuel feedback coefficient, measured from 300K to 350K, is: ' num2str(alpha_fuel_350(1)) ' +/- ' num2str(alpha_fuel_350(2)) ' pcm/°C'])
disp(['The prompt fuel feedback coefficient, measured from 300K to 350K considering Sab, is: ' num2str(alpha_fuel_350_Sab(1)) ' +/- ' num2str(alpha_fuel_350_Sab(2)) ' pcm/°C'])



%% REG calibration

% #1 = estimate value
% #2 = tolerance

k_cal(1,1) = K_eff.REG_in(1,1); % -4
k_cal(1,2) = K_eff.REG_in(1,2);

k_cal(2,1) = 1.00269; % 1
k_cal(2,2) = 0.00067 * k_cal(2,1);

k_cal(3,1) = 1.00343; % 6
k_cal(3,2) = 0.00069 * k_cal(3,1);

k_cal(4,1) = 1.00542; % 11
k_cal(4,2) = 0.00068 * k_cal(4,1);

k_cal(5,1) = 1.00632; % 12.5
k_cal(5,2) = 0.00097 * k_cal(5,1) ;

k_cal(6,1) = 1.00820; % 16
k_cal(6,2) = 0.00063 * k_cal(6,1);

k_cal(7,1) = 1.00966; % 21
k_cal(7,2) = 0.00078 * k_cal(7,1);

k_cal(8,1) = 1.01148; % 26
k_cal(8,2) = 0.00062 * k_cal(8,1);

k_cal(9,1) = 1.01401; % 31
k_cal(9,2) = 0.00064 * k_cal(9,1);

k_cal(10,1) = K_eff.REG_out(1,1); % 34
k_cal(10,2) = K_eff.REG_out(1,2);

k_cal(11,1) = 1.01533; % 39
k_cal(11,2) = 0.00091 * k_cal(11,1);

rho_cal(:,1) = (k_cal(:,1)-1) ./ k_cal(:,1); % (-)
rho_cal(:,2) = k_cal(:,2) ./ k_cal(:,1).^2;

rho_cal(:,1) = rho_cal(:,1) / BETA_eff.REG_in(1,1); % ($)
rho_cal(:,2) = rho_cal(:,1) .* sqrt( (rho_cal(:,2)./rho_cal(:,1)).^2 + (BETA_eff.REG_in(1,2)/BETA_eff.REG_in(1,1)).^2 );

reg_pos_23 = [0 7.2 12.6 15.7 18.1 20.1 24.6 29.0 34.7 38.1];
reg_exp_23 = [0 0.1427 0.3446 0.4819 0.5792 0.6691 0.8805 1.0401 1.1684 1.1909];
reg_exp_23_err = [0 0.0039 0.0076 0.0083 0.0086 0.0097 0.0114 0.0124 0.0129 0.0129];



reg_pos_65 = [116 260 344 417 480 559 645 821];
reg_exp_65 = [0 23.5 47.5 74.5 97 122 145.5 164.5];
reg_exp_65_err = sqrt( reg_exp_65/4 );

reg_pos_65 = 38 * (reg_pos_65 - reg_pos_65(1))/(reg_pos_65(end) - reg_pos_65(1));
reg_exp_65 = reg_exp_65/100;
reg_exp_65_err = reg_exp_65_err/100;



heights = [-4 1 6 11 12.5 16 21 26 31 34 39]-1;
figure
errorbar(heights(2:end), rho_cal(2:end,1)-rho_cal(2,1), 1.96*rho_cal(2:end,2), "r", 'LineWidth',1)
hold on
errorbar(reg_pos_65, reg_exp_65, 1.96*reg_exp_65_err, "g", 'LineWidth',1)
errorbar(reg_pos_23, reg_exp_23, 1.96*reg_exp_23_err, "b", 'LineWidth',1)
grid on
title('REG calibration curve', 'FontSize',12)
ylabel('Reactivity ($)', 'FontSize',12)
xlabel('Height from fuel element bottom (cm)', 'FontSize',12)
legend(["Serpent", "1965 experiment", "2023 experiment"], 'Location','northwest')

% Experimental CRW
CRW_65 = reg_exp_65(end)-reg_exp_65(1);
CRW_65_err = sqrt( reg_exp_65_err(end)^2 + reg_exp_65_err(1)^2 );
CRW_23 = reg_exp_23(end)-reg_exp_23(1);
CRW_23_err = sqrt( reg_exp_23_err(end)^2 + reg_exp_23_err(1)^2 );

disp(['The CRW of REG, measured at 300 K, is: ' num2str(rho_cal(end,1)-rho_cal(2,1)) ' +/- ' num2str(2*rho_cal(end,2)) ' $'])
disp(['The CRW of REG, measured in 1965, is: ' num2str(CRW_65) ' +/- ' num2str(CRW_65_err) ' $'])
disp(['The CRW of REG, measured in 2023, is: ' num2str(CRW_23) ' +/- ' num2str(CRW_23_err) ' $'])

%% Void coefficient
% In case of water filled central channel:
k_water_cc(1) = 1.00448;
k_water_cc(2) = 3.2e-4 * k_water_cc(1);

rho_water_cc(1) = (k_water_cc(1)-1) / k_water_cc(1); % (-)
rho_water_cc(2) = k_water_cc(2) / k_water_cc(1)^2;

rho_water_cc = rho_water_cc * 1e5; % pcm

%rho_water_cc(1) = rho_water_cc(1) / BETA_eff.REG_in(1,1); % ($)
%rho_water_cc(2) = rho_water_cc(1) .* sqrt( (rho_water_cc(2)./rho_water_cc(1)).^2 + (BETA_eff.REG_in(1,2)/BETA_eff.REG_in(1,1)).^2 );

volume_cc = pi*1.79^2 * 35.56; % cm3

alpha_void(1) = -(rho_water_cc(1) - RHO.REG_in(2,1)) / volume_cc;
alpha_void(2) = sqrt(rho_water_cc(2)^2 + RHO.REG_in(2,2)^2) / volume_cc;

fprintf('\n')
disp(['The void coefficient is ' num2str(alpha_void(1)) ' +/- ' num2str(alpha_void(2)) ' pcm/cm^3'])

%% Neutron flux spectra inside universes
%run("REG_in\Triga_res.m")
run("REG_in\Triga_res.m")

figure
for i=1:6
    % energie = (MICRO_E(i,1:end-1) + MICRO_E(i,2:end)) / 2;
    % errorbar(log10(energie), INF_MICRO_FLX(i,1:2:end-1), 2*INF_MICRO_FLX(i,2:2:end).*INF_MICRO_FLX(i,1:2:end-1))
    energie = (MICRO_E(i,1:end-1) + MICRO_E(i,2:end)) / 2;
    errorbar(energie, INF_MICRO_FLX(i,1:2:end-1), 2*INF_MICRO_FLX(i,2:2:end).*INF_MICRO_FLX(i,1:2:end-1))
    hold on
end
set(gca,'XScale','Log')
grid on
legend(["Fuel"; "Water"; "SHIM"; "REG"; "Graphite"; "Cladding"])
hold off



% energy index that separates the thermal and fast intervals
index = find(MICRO_E(1,:) == 6.25e-7);

% Logarithmic energy intervals
th_logE_intervals = log10(MICRO_E(1,2:index)) - log10(MICRO_E(1,1:index-1));
f_logE_intervals = log10(MICRO_E(1,index+1:end)) - log10(MICRO_E(1,index:end-1));

% Total thermal and fast fluxes in the universes
for i=1:5
    tot_th_flux(i) = sum( INF_MICRO_FLX(1,1:2:2*(index-1)-1) .* th_logE_intervals ) ;
    tot_f_flux(i) =  sum( INF_MICRO_FLX(1,index*2-1:2:end-1) .* f_logE_intervals ) ;
end




% 2 groups cross sections: 1°,3° rows = th and f XS; 2°,4° rows =
% corresponding errors
XS = struct("Total", zeros(4,5), ...   % (1/cm)
    "Capture", zeros(4,5), ...
    "Fission", zeros(4,5), ...
    "Absorption", zeros(4,5), ... % = capture + fission
    "Removal", zeros(4,5), ...    % = group removal + absorption
    "Scattering_gg", zeros(4,5),...  % = total - removal
    "OutScattering", zeros(4,5),...  % = removal - absorption
    "Scattering", zeros(4,5)...  % = total - absorption
    );

for i=1:5
    XS.Total(:,i) = INF_TOT(i,:);
    XS.Capture(:,i) = INF_CAPT(i,:);
    XS.Fission(:,i) = INF_FISS(i,:);
    XS.Absorption(:,i) = INF_ABS(i,:);
    XS.Removal(:,i) = INF_REMXS(i,:);
    XS.Scattering_gg(:,i) = XS.Total(:,i) - XS.Removal(:,i);
    XS.OutScattering(:,i) = XS.Removal(:,i) - XS.Absorption(:,i);
    XS.Scattering(:,i) = XS.Total(:,i) - XS.Absorption(:,i);
end

fields = fieldnames(XS);

% Now energy groups are in descending order, I want them in ascending
% order:
invert_mat = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]; % I switch the rows 1,2 with 3,4 and viceversa

for i=1:length(fields)
    XS.(fields{i}) = invert_mat*XS.(fields{i});
end

% Calculate errors
for i=1:length(fields)
    if fields{i}=="Scattering_gg"
        XS.(fields{i})(2,:) = sqrt( XS.Total(2,:).^2 + XS.Removal(2,:).^2 );
        XS.(fields{i})(4,:) = sqrt( XS.Total(4,:).^2 + XS.Removal(4,:).^2 );
    elseif fields{i}=="OutScattering"
        XS.(fields{i})(2,:) = sqrt( XS.Absorption(2,:).^2 + XS.Removal(2,:).^2 );
        XS.(fields{i})(4,:) = sqrt( XS.Absorption(4,:).^2 + XS.Removal(4,:).^2 );
    elseif fields{i}=="Scattering"
        XS.(fields{i})(2,:) = sqrt( XS.Total(2,:).^2 + XS.Absorption(2,:).^2 );
        XS.(fields{i})(4,:) = sqrt( XS.Total(4,:).^2 + XS.Absorption(4,:).^2 );
    else
        XS.(fields{i})(2,:) = XS.(fields{i})(2,:) .* XS.(fields{i})(1,:);
        XS.(fields{i})(4,:) = XS.(fields{i})(4,:) .* XS.(fields{i})(3,:);
    end
end

Diffusion_coeff = zeros(4,5);
for i=1:5
    Diffusion_coeff(:,i) = INF_DIFFCOEF(i,:);
end
Diffusion_coeff = invert_mat*Diffusion_coeff;
Diffusion_coeff(2,:) = Diffusion_coeff(2,:).*Diffusion_coeff(1,:);
Diffusion_coeff(4,:) = Diffusion_coeff(4,:).*Diffusion_coeff(3,:);




%% Microscopic cross sections 
%clear all
% run("REG_in\Triga_xs0.m");

%clearvars -except E -regexp \w*xs

%fields = fieldnames(mic_xs);

%plot(log10(E), [i1001_03c_xs])%; i1001_03s_xs; i5010_03c_xs; i5011_03c_xs; i6012_03c_xs; i6012_03s_xs; i7014_03c_xs])

%legend(fields{2:end})
% %%
% figure
% loglog(E, [mCladding_xs(:,1)'; mControlRods_xs(:,1)'; mFuel_xs(:,1)'; mGraphite_xs(:,1)'; mVacuum_xs(:,1)'; mWater_xs(:,1)'], 'LineWidth',2)
% legend(["Cladding", "CRs", "Fuel", "Graphite", "Vacuum", "Water"])
% grid on
% %%
% figure
% loglog(E, [mWater_xs(:,1)'; i1001_03c_xs(:,1)'; i1001_03s_xs(:,1)'])
% legend(["water", "03c", "03s"])
% grid on

