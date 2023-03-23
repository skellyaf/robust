stt = randn(6,6,6);

v1 = randn(6,1);

% v2 = randn(6,1);
v2 = v1;

mult = zeros(6,1);

for i = 1:6
    for k1 = 1:6
        for k2 = 1:6
            mult(i, k1, k2) = mult(i, k1, k2) + stt(i, k1, k2) * v1(k1) * v1(k2);
        end
    end
end

mult2 = tmult(stt, v1*v1');

mult - mult2;


mult3 = zeros(6,6,6);
for i = 1:6
    mult3(i,:,:) = squeeze(stt(i,:,:)) * v1*v1';
end




vmult = zeros(6,6);
% svmult = sym('vm', [6 6]);
sv1 = sym('a', [6 1], 'real');


for k1 = 1:6
    for k2 = 1:6
        svmult(k1,k2) =  sv1(k1) * sv1(k2);
    end
end


sstt = sym('s',[ 6 6 6], 'real')
smult2 = sym('z', [6 1], 'real')
% smult2 = zeros(6,6,6);

for i = 1:6
    for k1 = 1:6
        for k2 = 1:6
%             smult2(i, k1, k2) = smult2(i, k1, k2) + sstt(i, k1, k2) * sv1(k1) * sv1(k2);
            smult2(i) = smult2(i) + sstt(i, k1, k2) * sv1(k1) * sv1(k2);
        end
    end
end






smult3 = sym('z', [6 1], 'real')
% smult2 = zeros(6,6,6);

for i = 1:6
    for k1 = 1:6
        for k2 = 1:6
%             smult2(i, k1, k2) = smult2(i, k1, k2) + sstt(i, k1, k2) * sv1(k1) * sv1(k2);
            smult3(i) = smult3(i) + sv1(k1) * sv1(k2);
        end
    end
end








sa = sym('a',[3 1],'real')
sv = sym('v',[3 1],'real')


t1=(sa+sv) * (sa+sv)'
t2=sa*sa' + sv*sv' + sa*sv' + sv*sa'
simplify(t1-t2)
