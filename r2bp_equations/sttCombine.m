function sttkj = sttCombine(stm_i, stt_i, j, k, stmkj, sttkj)
% Input - tensor of STTs that are sequential in time
% input - tensor of STMs that are sequential for the same time periods
% Optional inputs - initial STM and STT to combine with
% Output - a combination of multiple STTs into a single longer time STT
% based on the desired indices (from index j to index k)

assert(j <= k,"Error: Requesting STM to go backwards");
assert(size(stm_i,3) == size(stt_i,4),'Error: STM and STT not the same number of segments!')

if nargin == 4
    stmkj = eye(6);
    sttkj = zeros(6,6,6);
end

for i = j:k

    sttkj = tensorCombine(stmkj, sttkj, stm_i(:,:,i), stt_i(:,:,:,i));
    stmkj = stm_i(:,:,i) * stmkj;

end

end