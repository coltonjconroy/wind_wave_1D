
function [alph,bet,crk] = SSPRK(p,nrk)

%--------------------------------------------------------------------------

order = nrk;

if nrk == 1 & order == 1
    cfl = 1/(2*p+1);
    alph = ones([1 1]);
    bet = ones([1 1]);
    crk = [1];

end

if nrk == 2 & order == 2
    cfl = 1/(2*p+1);
    alph = [1,0;.5,.5];
    bet = [1,0;0,.5];
    crk = [ 0 1 ];
end

if nrk == 3 & order == 3
    alph = [1,0,0;3/4,1/4,0;1/3,0,2/3];
    bet = [1,0,0;0,1/4,0;0,0,2/3];
    crk = [ 0 1/2 1 ];
end

if order > 3

    alph = [1,0,0;3/4,1/4,0;1/3,0,2/3];
    bet = [1,0,0;0,1/4,0;0,0,2/3];
    crk = [ 0 1/2 1 ];
end