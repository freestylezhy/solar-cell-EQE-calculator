function phi = photonFlux_AM15G_eV(absorptivity)

q = 1.602176565e-19;      %unit chage
k = 1.3806505e-23;  %Bolzmann constant
h = 6.6260693e-34;  %Plank constant
c = 299792458;      %speek of light

global AM15G
if isempty(AM15G)
    AM15G = load('AM15G.txt');
end

AM15G_lambda = AM15G(:,1);  %nm
AM15G_power = AM15G(:,2);   %W/m^2/nm
AM15G_eV = h*c./(AM15G_lambda*1e-9)/q;  %eV
AM15G_photonFluxDensity = AM15G_power*1e9/1e4*c*h/q^2./(AM15G_eV.^3); %#/cm^2/eV
AM15G_eV = flipud(AM15G_eV);
AM15G_photonFluxDensity = flipud(AM15G_photonFluxDensity);

% %%% debug
%     figure(1)
%     plot(AM15G_eV, AM15G_photonFluxDensity)
% %%% end of debug

AM15G_absorptivity = interp1(absorptivity(:,1), absorptivity(:,2), AM15G_eV);
AM15G_absorptivity(isnan(AM15G_absorptivity)) = 0;
phi = trapz(AM15G_eV, AM15G_absorptivity.*AM15G_photonFluxDensity);


