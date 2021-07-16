%%% Simulate Bloch equations using various examples
% See
% C. Graf, A. Rund, C.S. Aigner, R. Stollberger,
% Accuracy and Performance Analysis for Bloch and Bloch-McConnell
% simulation methods
% Journal of Magnetic Resonance 329(3):107011
% doi: 10.1016/j.jmr.2021.107011
%%%

clear
clc

%%%% Function calls Bloch simulation with asymmetric and symmetric operator
%%%% splitting

example = 2; % 1: sinc RF and trapezoidal Gs, 2: SMS pulse

if example==1
    load('Example1.mat');
elseif example==2
    load('Example2.mat');
end

T1=[10^-9 1331 400 832 1420]; % without relax, grey matter, tendons, white matter, muscle  at 3T
T2=[10^-9 110 5 79.6 31.7];
relax_type = {'Without Relaxation', 'Grey Matter', 'Tendons', 'White Matter', 'Muscle'};
relax=[0 1 1 1 1]; % flag, set to 1 if relaxation should be included in the simulation

    
for i=1:5
    d.T1=T1(i);
    d.T2=T2(i);
    d.relax=relax(i);
    M_sy=bloch_symmetric_splitting(u,v,w,d);
    
    M_asy=bloch_asymmetric_splitting(u,v,w,d);
    
        
    figure
    subplot(1,2,1)
    plot(d.xdis*10^3,M_asy(1,:,end),'LineWidth',3)
    hold on
    plot(d.xdis*10^3,M_asy(2,:,end),'LineWidth',3)
    plot(d.xdis*10^3,M_asy(3,:,end),'LineWidth',3)
    xlabel('distance in mm')
    ylabel('magnetization in a.u.')
    ax=gca;
    set(gca, 'FontSize', 25)
    legend('Mx','My','Mz')
    title('Asymmetric Splitting')
    if example ==1
        set(gca,'XTick',[-5 -2.5 0 2.5 5])
    else
        set(gca,'XTick',[-60 -30 0 30 60])
    end
    set(gca,'YTick',[-1 -0.5 0 0.5 1])
    
    subplot(1,2,2)
    plot(d.xdis*10^3,M_sy(1,:,end),'LineWidth',3)
    hold on
    plot(d.xdis*10^3,M_sy(2,:,end),'LineWidth',3)
    plot(d.xdis*10^3,M_sy(3,:,end),'LineWidth',3)
    xlabel('distance in mm')
    ylabel('magnetization in a.u.')
    ax=gca;
    set(gca, 'FontSize', 25)
    legend('Mx','My','Mz')
    title('Symmetric Splitting')
    if example ==1
        set(gca,'XTick',[-5 -2.5 0 2.5 5])
    else
        set(gca,'XTick',[-60 -30 0 30 60])
    end
    set(gca,'YTick',[-1 -0.5 0 0.5 1])
    
    t1 = suptitle(relax_type(i));
    t1.FontSize = 40;
end
