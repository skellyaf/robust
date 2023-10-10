function [ret] = existsAndTrue(inputName, simparams)
%existsAndTrue Checks if the inputName string exists and is true. Returns true
%if so, false if either aren't true.

if isfield(simparams, inputName)
    if simparams.(inputName)
        ret = true;
    else
        ret = false;
    end
else
    ret = false;
end

end