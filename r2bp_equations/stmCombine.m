function stmkj = stmCombine(stm_i, j, k)
% Input - tensor of STMs that are sequential in time
% Output - a combination of multiple STMs into a single longer time STM
% based on the desired indices (from index j to index k)

assert(j <= k,"Error: Requesting STM to go backwards");

stmkj = eye(size(stm_i,1));

for i = j:k
    stmkj = stm_i(:,:,i) * stmkj;
end

end