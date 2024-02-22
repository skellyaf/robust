function [inverse_stm] = invert_stm(stm, simparams)




if strcmp(simparams.dynSys,'2bp')
    % Symplectic unit matrix
    J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];
    % Symplectic
    inverse_stm = -J * stm' * J;
% elseif strcmp(simparams.dynSys,'cr3bp')
else
    % Not symplectic
    inverse_stm = inv(stm);
end



end