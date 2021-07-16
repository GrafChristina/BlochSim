%%% Simulate Bloch Mc-Connell equations using various examples
% See
% C. Graf, A. Rund, C.S. Aigner, R. Stollberger,
% Accuracy and Performance Analysis for Bloch and Bloch-McConnell
% simulation methods
% Journal of Magnetic Resonance 329(3):107011
% doi: 10.1016/j.jmr.2021.107011
%%%

clc
clear
 

example = 3;
if example == 3                     % 2 Pool CEST, water & solute, no MT
    load('init_gauss.mat');         % Gaussian RF pulse
%     load('init_block.mat');       % Block RF pulse
elseif example == 4                 % 2 Pool CEST, water & MT, using basic Bloch-McConnell
    load('Example4.mat');
elseif example == 5                 % 3 Pool CEST, water, solute & MT, using basic Bloch-McConnell
    load('Example5.mat');
elseif example == 6                 % 2 Pool Pulsed MT, i.e. water & MT
    load('Example6.mat');
elseif example == 7                 % 3 Pool Pulsed MT, i.e. water, solute & MT
    load('Example7.mat');
end




M_BMcC_s=bmcc_symmetric_splitting(d);

M_BMcC_as=bmcc_asymmetric_splitting(d);


figure
plot(-d.xZspec,M_BMcC_s(3,:,end),'LineWidth',3)
hold on
plot(-d.xZspec,M_BMcC_as(3,:,end),'LineWidth',3)
xlabel('z-spectrum in ppm')
ylabel('water magnetization in a.u.')
legend('SY','ASY')

ax=gca;
set(gca, 'FontSize', 25)
set(gca,'XTick',[-5 -2.5 0 2.5 5])
set(gca,'YTick',[-1 -0.5 0 0.5 1])