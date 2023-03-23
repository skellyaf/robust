function [X] = subdivide_segment(x_t, t, n)
%subdivide_segment splits a single trajectory segment up into n coast
%segments, separated evenly by number of indices across the time vector t.

X = [];
X = [x_t(1,:)'];

idx_sep = floor( length(t) / n );

iter = 1;
i = idx_sep;
X(7,1) = t(i);

while iter    
    
    if i > length(t) - idx_sep
        if i - idx_sep > 0
            X(7,end) = t(end) - t(i - idx_sep);
        end
        iter = 0;
    else
        X(1:6,end+1) = x_t(i,:)';
        X(7,end) = t(i + idx_sep) - t(i);
    end
    i = i + idx_sep;   
    

end



end