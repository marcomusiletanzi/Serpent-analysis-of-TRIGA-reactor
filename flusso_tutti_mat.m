%#########################################
%## Initial checking and pre-processing ##
%#########################################

%## Check that the detector file exists

if (exist("..\prova ese 3\Triga_det0.m", "file") ~= 2)
  disp("Could not find infinite_det0.m from current folder! Cannot do analysis.")
end

%## Run the detector output file to bring the results to workspace

run ..\'prova ese 3'\Triga_det0.m;

%## Check that the detector output exist

%#####################################
%## Plot the energy-integrated flux ##
%#####################################

%## Scale the energy integrated flux to a maximum of 1.0


% DETEnergyDetFuel(:,11) = DETEnergyDetFuel(:,11)/max(DETEnergyDetFuel(:,11));
% DETEnergyDetWater(:,11) = DETEnergyDetWater(:,11)/max(DETEnergyDetWater(:,11));
% DETEnergyDetREG(:,11) = DETEnergyDetREG(:,11)/max(DETEnergyDetREG(:,11));
% DETEnergyDetSHIM(:,11) = DETEnergyDetSHIM(:,11)/max(DETEnergyDetSHIM(:,11));
% DETEnergyDetRef(:,11) = DETEnergyDetRef(:,11)/max(DETEnergyDetRef(:,11));

total = max(max([DETEnergyDetFuel(:,11),DETEnergyDetWater(:,11),DETEnergyDetREG(:,11),DETEnergyDetSHIM(:,11),DETEnergyDetRef(:,11)]));
DETEnergyDetFuel(:,11) = DETEnergyDetFuel(:,11)/total;
DETEnergyDetWater(:,11) = DETEnergyDetWater(:,11)/total;
DETEnergyDetREG(:,11) = DETEnergyDetREG(:,11)/total;
DETEnergyDetSHIM(:,11) = DETEnergyDetSHIM(:,11)/total;
DETEnergyDetRef(:,11) = DETEnergyDetRef(:,11)/total;




%## Plot

figure(1);
errorbar(DETEnergyDetFuelE(:,3), DETEnergyDetFuel(:,11), 2*DETEnergyDetFuel(:,11).*DETEnergyDetFuel(:,12),'b.');
hold on
errorbar(DETEnergyDetWaterE(:,3), DETEnergyDetWater(:,11), 2*DETEnergyDetWater(:,11).*DETEnergyDetWater(:,12),'r.');
errorbar(DETEnergyDetREGE(:,3), 1e2*DETEnergyDetREG(:,11), 2*DETEnergyDetREG(:,11).*DETEnergyDetREG(:,12),'g.');
errorbar(DETEnergyDetSHIME(:,3), 1e2*DETEnergyDetSHIM(:,11), 2*DETEnergyDetSHIM(:,11).*DETEnergyDetSHIM(:,12),'y.');
errorbar(DETEnergyDetRefE(:,3), DETEnergyDetRef(:,11), 2*DETEnergyDetRef(:,11).*DETEnergyDetRef(:,12),'k.');

%## Set axes

set(gca,'XScale','log');
set(gca,'YScale','linear');
set(gca,'XTick',[1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e0,1e2]);
%set(gca,'YTick',[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1]);
set(gca,'FontSize',16);
set(gca,'XLim',[1e-12 1e2]);
set(gca,'YLim',[0 1.1]);

%## Make the plot a bit nicer

xlabel('Energy (MeV)')
ylabel('Neutron flux (a.u.)')
legend(["Fuel", "Water", "REG x 100", "SHIM x 100", "Reflector"])
grid on

box on
