x = simparams.x0;
% x = x_opt(:);

% % % 
%% Objective function gradient

if 1

%     dx = 1e-6;
    dx = sqrt(eps);
%     dx = 1e-10;
    testGradJ = zeros(length(x(:)),1);

    [obj_x, ~] = obj_min_tcm(x,simparams);

    for i = 1:length(x(:))
%     parfor i = 1:length(x(:))
        i
        if i == 59
            stopp = 1;
        end
%     for i = length(x(:))
%             dxrand = dx*randn;
            dxrand = dx*1;
%             dxrand = dx*x(i);
            xdx = x;
            xdx(i) = x(i) + dxrand;
%             testGradJ(i) = (  obj_min_tcm(xdx,simparams) - obj_min_tcm(x,simparams)  )/dxrand; 
            [obj_xdx, ~] = obj_min_tcm(xdx,simparams);
            testGradJ(i) = (  obj_xdx - obj_x  )/dxrand; 

    end
    testGradJ
    [~,analyticalGradJ] = obj_min_tcm(x,simparams);
%     analyticalGradJ = analyticalGradJ'

    [testGradJ, analyticalGradJ]
    p=2;
    
end



%%%%%%%%%%%%%%%%% VERIFYING EQUALITY CONSTRAINT GRADIENTS
if 1
    dx = sqrt(eps);
    outputGradients = 1;
    [~,nonpCeq] = constraint_circ_x0xf(x, simparams);
    testGradceq = [];
    
    for j = 1:length(x(:)) 
        j
        xdx = x;
        dxd = dx;
%         dxd = dx*x(j);
        xdx(j) = x(j) + dxd;
        [~,pertCeq] = constraint_circ_x0xf(xdx, simparams);
    
        testGradceq(:,end+1) = (  pertCeq - nonpCeq  )/dxd; 

        if j == 127
            pp=1;
        end
    
    end
    testGradceq
    [~, ~, ~, ceqGrad_an] = constraint_circ_x0xf(x, simparams)
    ceqGrad_an = ceqGrad_an';
    diffGradceq = testGradceq - ceqGrad_an

    (ceqGrad_an - testGradceq) ./ testGradceq

    

end


%%%%%%%%%%%%%%%%%% VERIFYING INEQUALITY CONSTRAINT GRADIENTS

dx = sqrt(eps);
[nonpCin] = constraint_circ_x0xf(x, simparams);

testGradcin = [];
    for j = 1:length(x(:)) %3
        j

        if j == 134
            ppp=1;
        end

        xdx = x;
        xdx(j) = x(j) + dx;
        [pertCin] = constraint_circ_x0xf(xdx, simparams);        
        testGradcin(:,end+1) = (  pertCin - nonpCin  )/dx; 

        
        
    end

[~, ~, cinGrad_an, ~] = constraint_circ_x0xf(x, simparams);

cinGrad_an = cinGrad_an';


pp=1;
