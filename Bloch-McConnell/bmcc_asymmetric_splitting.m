function[M]=bmcc_asymmetric_splitting(d)
%%% Asymmetric Operator Splitting using a different number of pools and 2
% different models
% See
% C. Graf, A. Rund, C.S. Aigner, R. Stollberger,
% Accuracy and Performance Analysis for Bloch and Bloch-McConnell
% simulation methods
% Journal of Magnetic Resonance 329(3):107011
% doi: 10.1016/j.jmr.2021.107011
%%%

np=d.np;                                                            % Number of Pools

if d.MT==1                                                          % Henkelman model
    
    % Calculate magnetization vector M using an asymmetric operator splitting method
    if np==2
        M=bmcc_asymmetric_splitting_2_large_scale_MT_2(d);    % 2 Pools
    elseif np==3
        M=bmcc_asymmetric_splitting_3_vectorised_large_scale_MT_2(d);    % 3 Pools
    else
        disp('Not implemented');
    end
    
else                                                                % classical Cest model
    
    if np==2
        M=bmcc_asymmetric_splitting_2_vectorised_large_scale(d);         % 2 Pools
        
    elseif np==3
        M=bmcc_asymmetric_splitting_3_vectorised_large_scale(d);         % 3 Pools
    else
        disp('Not implemented');
    end
    
end

end