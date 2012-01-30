field = 0:0.01:0.05;
phi1 = 0;
phi2 = 180;
p1i = 1;
p2i = 1;

dp = cell(1,length(field));
dpgs = cell(1,length(field));
hlgap = cell(1,length(field));
eexc = cell(1,length(field));
igs = cell(1,length(field));
opint = cell(1,length(field));
tpint = cell(1,length(field));
dp2exc = cell(1,length(field));

for i = 1:numel(p1i)
    for j = 1:numel(p2i)
        %field = [0:0.1:(floor(10*maxdata(p1i(i),p2i(j)) - 0.1) / 10) secondthreshold(p1i(i),p2i(j)):0.1:(secondthreshold(p1i(i),p2i(j))+1)];
        dptemp = zeros(1,length(field));
        dpgstemp = zeros(1,length(field));
        dp2exctemp = zeros(1,length(field));
        hltemp = zeros(1,length(field));
        eexctemp = zeros(length(field),25);
        igstemp = zeros(1,length(field));
        opinttemp = zeros(length(field),25);
        tpinttemp = zeros(length(field),25);
        
        for k = 1:numel(field)
            myece = ECEParams('AM1', {phi1(i); phi2(j)}, true, field(k), 100, 25);
            myexp = EnergyCalcExp(myece,...
                's:\NoAngles\12-merPPV.dat',...
                [6 94],...
                's:\NoAngles\scratch\12-merPPV-scratch.mat',...
                'o:\',...
                'p:\',...
                false,...
                ['12-merPPV MaxFieldTest - ',num2str(field(k)),'VA']);
            myexp.run('quiet','noprogress');
            
            myexp.data(1).load_to_memory('indo','load');
            
            tmp = myexp.get_field('indo.dipole',1,1,phi1(i),phi2(j));
            dpgstemp(k) = sum(tmp .^ 2, 1) .^ (0.5);
            
            tmp = myexp.get_field('indo.dipole',2,2,phi1(i),phi2(j));
            dptemp(k) = sum(tmp .^ 2, 1) .^ (0.5);
            
            tmp = myexp.get_field('indo.dipole',3,3,phi1(i),phi2(j));
            dp2exctemp(k) = sum(tmp .^ 2, 1) .^ (0.5);
            
            nfill = myexp.get_field('indo.nfilled',phi1(i),phi2(j));
            tmp = myexp.get_field('indo.orbE',[nfill nfill+1],phi1(i),phi2(j));
            hltemp(k) = tmp(2) - tmp(1);
            
            eexctemp(k,1:25) = myexp.get_field('Eexc',:,phi1(i),phi2(j));
            igstemp(k) = myexp.get_field('indo.esci',1,phi1(i),phi2(j));
            opinttemp(k,1:25) = myexp.get_field('Tint',:,phi1(i),phi2(j));
            tmp = squeeze(myexp.get_field('indo.r',2,:,:,phi1(i),phi2(j)));
            for l = 1:25
                tpinttemp(k,l) = (eexctemp(k,l) - eexctemp(k,2)) * sum(tmp(l,:) .^ 2);
            end
        end
        
        dp{i,j} = dptemp;
        hlgap{i,j} = hltemp;
        eexc{i,j} = eexctemp;
        igs{i,j} = igstemp;
        opint{i,j} = opinttemp;
        tpint{i,j} = tpinttemp;
        dp2exc{i,j} = dp2exctemp;
        dpgs{i,j} = dpgstemp;
    end
end
% 
% %% 
% stilenergy = zeros(6,25);
% stilgstint = zeros(6,25);
% stilfestint = zeros(6,25,3);
% sfi = zeros(6,25);
% 
% for i = 1:6
%     % amp = myexp(i).get_field('ampac.Hf') * 1/23.06;
%     indE = myexp(i).get_field('indo.esci',1);
%     stilenergy(i,:) = myexp(i).get_field('Eexc',:) + indE;
%     stilgstint(i,:) = myexp(i).get_field('Tint',:);
%     stilfestint(i,:,:) = myexp(i).get_field('indo.r',2,:,:);
%     for j = 1:25;
%         sfi(i,j) = (stilenergy(i,j) - stilenergy(i,2)) * sum(stilfestint(i,j,:) .^ 2, 3);
%     end
% end

%% 
figure(8)
hold on

% xaxis = 0:0.01:0.05;
i = 1; j = 1;
%xaxis = [0:0.1:(floor(10*maxdata(p1i(i),p2i(j)) - 0.1) / 10) secondthreshold(p1i(i),p2i(j)):0.1:(secondthreshold(p1i(i),p2i(j))+1)];
% xaxis = xaxis(1:8);
xaxis = field;

for k = 1:length(xaxis)
    plot(xaxis(k), eexc{i,j}(k,:) + igs{i,j}(k), 'g^');
end

maxopint = max(opint{i,j}(:));
for k = 1:length(xaxis)
    for l = 1:25
        plot(xaxis(k), eexc{i,j}(k,l) + igs{i,j}(k), 'bo', 'MarkerSize', (opint{i,j}(k,l)*30 / maxopint) + 1e-3);
    end
end

maxtpint = max(tpint{i,j}(:));

for k = 1:length(xaxis)
    for l = 1:25
        if (tpint{i,j}(k,l) > 0)
            plot(xaxis(k), eexc{i,j}(k,l) + igs{i,j}(k), 'rs', 'MarkerSize', (tpint{i,j}(k,l)*30 / maxtpint) + 1e-3);
        end
    end
end