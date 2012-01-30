function res = diffDensity(obj, istate)
% diffDensity(indo, state_number)
% Calculate change in density matrix (in atomic orbitals) for 
% state istate of an indo object (1= ground state) 

norb = obj.norb;
nfilled = obj.nfilled;
nempty = obj.norb-obj.nfilled;
% convert wavefunction into format wf(a,r)
wf = zeros(nfilled,nempty);
for isci = 1:obj.nscibasis
   a = obj.ehsci(isci,1);          % hole label
   r = obj.ehsci(isci,2)-nfilled;  % elec label (LUMO = 1)
   wf(a,r) = obj.wfsci(isci,istate);
end

% First get the change in density in molecular orbitals
rhomo = zeros(norb,norb);
% the filled portion is: -sum_s wf(a,s) wf(b,s) = -wf*wf'
rhomo(1:nfilled,1:nfilled) = -wf*(wf');
rempty = (nfilled+1):norb;
% the empty portion is: +sum_a wf(a,r) wf(a,s) = +wf'*wf
rhomo(rempty,rempty) = +wf'*wf

% now transfrom to atomic orbitals  orb(atomic, molecular)
res = orb * (rhomo*(orb'));

end