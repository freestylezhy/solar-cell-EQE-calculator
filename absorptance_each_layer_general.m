close all; clear all; clc

addpath ./material_parameters
addpath ./measurement_results

data = importdata('CdTe.txt');
eV_CdTe = data.data(:,1);
n_CdTe = data.data(:,2);
k_CdTe = data.data(:,3);
lambda_CdTe = 1.24./eV_CdTe;

data = importdata('InSb.txt');
lambda_InSb = data.data(:,1)/1e3;
n_InSb = data.data(:,2);
k_InSb = data.data(:,3);

data = importdata('p-type_a-Si_Holman.txt');
lambda_aSi = data.data(:,1)/1e3;
n_aSi = data.data(:,2);
k_aSi = data.data(:,3);

data = importdata('Mathieu_measured_ITO.txt');
lambda_ITO = data.data(:,1)/1e3;
n_ITO = data.data(:,2);
k_ITO = data.data(:,3);

data = importdata('Mathieu_measured_SiC.txt');
lambda_aSiC = data.data(:,1)/1e3;
n_aSiC = data.data(:,2);
k_aSiC = data.data(:,3);


%%
lambda_set = 0.3:0.001:1.6;
eV_set = 1.24./lambda_set;
N0 = 1;     %complex refractive index of air

%   layer structure
%
%   index       material  thickness (um)
%   0           air         inf
%   1           ITO         0.07
%   2           aSi         0.01
%   3           MgCdTe      0.01
%   4           CdTe        1
%   5           MgCdTe      0.05
%   6           CdTe        0.5
%   7           InSb        inf

%thickness of each layer
d1 = 0.07;
d2 = 0.01;
d3 = 0.01;
d4 = 1;
d5 = 0.05;
d6 = 0.5;
d = [d1 d2 d3 d4 d5 d6];

sampleNo = 4;

%n&k of each layer
Eg_windows = 2.04;
Eg_BSF = 1.93; 
N1_set = interp1(lambda_ITO, n_ITO, lambda_set,'linear','extrap') - 1i*interp1(lambda_ITO, k_ITO, lambda_set,'linear','extrap');
N2_set = interp1(lambda_aSi, n_aSi, lambda_set,'linear','extrap') - 1i*interp1(lambda_aSi, k_aSi, lambda_set,'linear','extrap');
%N2_set = interp1(lambda_aSiC, n_aSiC, lambda_set,'linear','extrap') - 1i*interp1(lambda_aSiC, k_aSiC, lambda_set,'linear','extrap');
N3_set = interp1(eV_CdTe, n_CdTe, eV_set-(Eg_windows-1.5),'linear','extrap') - 1i*interp1(eV_CdTe, k_CdTe, eV_set-(Eg_windows-1.5),'linear','extrap');
N4_set = interp1(lambda_CdTe, n_CdTe, lambda_set,'linear','extrap') - 1i*interp1(lambda_CdTe, k_CdTe, lambda_set,'linear','extrap');
N5_set = interp1(eV_CdTe, n_CdTe, eV_set-(Eg_BSF-1.5),'linear','extrap') - 1i*interp1(eV_CdTe, k_CdTe, eV_set-(Eg_BSF-1.5),'linear','extrap');
N6_set = N4_set;
N7_set = interp1(lambda_InSb, n_InSb, lambda_set,'linear','extrap') - 1i*interp1(lambda_InSb, k_InSb, lambda_set,'linear','extrap');

for i = 1:length(lambda_set)
    lambda = lambda_set(i);
    N1 = N1_set(i);
    N2 = N2_set(i);
    N3 = N3_set(i);
    N4 = N4_set(i);
    N5 = N5_set(i);
    N6 = N6_set(i);
    N7 = N7_set(i);
    
    %normal incidence
    nta = [N1 N2 N3 N4 N5 N6 N7];
    num_layers = length(nta);
    nta0 = N0;
    for j = 1:num_layers-1
        delta(j) = 2*pi*nta(j)*d(j)/lambda;
        cMat{j} = [cos(delta(j)) 1i*sin(delta(j))/nta(j);
            1i*nta(j)*sin(delta(j))  cos(delta(j))];
    end
    
    cMat_all = eye(2);
    for j = num_layers-1:-1:1
        cMat_all = cMat{j}*cMat_all;
    end
    BC = cMat_all*[1;nta(end)];
    Y = BC(2)/BC(1);
    rou = (nta0-Y)/(nta0+Y);
    R = rou*conj(rou);
    R_set(i) = R;
    
    B = BC(2);
    C = BC(1);
    T_set(i) = 4*nta0*real(nta(end))/((nta0*B+C)*conj(nta0*B+C));
    
    EH{num_layers} = [1;nta(end)];
    for j = num_layers-1:-1:1
        EH{j} = cMat{j}*EH{j+1};
    end
    
    for j = 1:num_layers
        intensity_raw(j) = 1/2*real(EH{j}(1)*conj(EH{j}(2)));
    end
    normalization_factor = (1-R)/intensity_raw(1);
    intensity = normalization_factor*intensity_raw;
    for j = 1:num_layers-1
        abs_set(j,i) = intensity(j)-intensity(j+1);
    end
    abs_set(num_layers,i) = intensity(num_layers);
    
    
    
end

figure(1)
plot(eV_set, R_set)
figure(2)
plot(lambda_set, R_set)

figure(3)
area(lambda_set, [abs_set; R_set]');
xlim([0.3 1.5])
ylim([0,1.1])

figure(1003)
abs_set_new = abs_set;
abs_set_new(1, :) = abs_set(4, :);
abs_set_new([2 3 4], :) = abs_set([1 2 3], :);
area(lambda_set, [abs_set_new; R_set]');
xlim([0.3 0.9])
ylim([0,1.1])
colormap Jet

%area(lambda_set, [abs_set; R_set]');
X1 = lambda_set*1e3;
ymatrix1 = [abs_set_new; R_set]'*100;
createfigure_Jsc_loss(X1, ymatrix1)
xlim([300 900])
ylim([0 100])
xlabel('Wavelength (nm)')
ylabel('Reflectence, Transmittance, Absorptance (%)')

for i = 1:num_layers
    figure(i+100)
    plot(eV_set, abs_set(i,:))
end

for i = 1:num_layers
    figure(i+200)
    plot(lambda_set, abs_set(i,:))
end

%% Jsc in each layer
%addpath D:\research\solar_cell_simulation_matlab\spectrum
q=1.60217733e-19;       % C     electron change
index = eV_set>1.5;
for i = 1:num_layers
    absorptivity = [eV_set(index)' abs_set(i,index)'];
    Jsc_set(i) = q*photonFlux_AM15G_eV(absorptivity);
end

index = eV_set>1.5;
absorptivity = [eV_set(index)' R_set(index)'];
Jsc_Reflection = q*photonFlux_AM15G_eV(absorptivity);

Jsc_set
Jsc_Reflection
sum(Jsc_set) + Jsc_Reflection

%%
%compared measured EQE and absorptance
figure(1001)
plot(lambda_set*1e3, abs_set(4,:)*100)
hold all
data = importdata('measured_EQE_A1746_all.txt');
lambda_EQE = data.data(:,1);
EQE = data.data(:,sampleNo);
plot(lambda_EQE, EQE*100)

figure(1001)
plot(lambda_set*1e3, (1-R_set)*100)
hold all
data = importdata('measured_1-R_A1746_all.txt');
lambda_one_minus_R = data.data(:,1);
one_minus_R = data.data(:,sampleNo);
plot(lambda_one_minus_R, one_minus_R*100)

xlim([300, 1500])


xlabel('Wavelength (nm)')
ylabel('1-R and EQE (%)')

legend('Calculated absorptance', 'EQE', 'Calculated 1-R', 'Measured 1-R')

