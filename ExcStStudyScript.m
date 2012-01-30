i = 1;
myexp = repmat(EnergyCalcExp(), [1 8]);
for field = 0:0.1:0.7
    myece = ECEParams('AM1', {0; [90 93 135 180]}, true, field, 100, 25);
    myexp(i) = EnergyCalcExp(myece,...
        's:\2 angles\stilbene.dat',...
        [6 14],...
        ['s:\2 angles\stilbene-ExcStStudy-',num2str(field),'VA.mat'],...
        'o:\',...
        'p:\',...
        false,...
        ['Stilbene Excited State Study - ',num2str(field),'VA']);
    myexp(i).run(true, 60, true);
    i = i + 1;
end
%%
i = 1;
myexp = repmat(EnergyCalcExp(), [1 6]);
for field = 0:0.01:0.05
    myece = ECEParams('AM1', {}, true, field, 500, 25);
    myexp(i) = EnergyCalcExp(myece,...
        'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\Exp\NoAngles\12-merPPV.dat',...
        [6 94],...
        ['C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\Exp\NoAngles\stilbene-ExcStStudy-',num2str(field),'VA.mat'],...
        'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\INDOLib\',...
        'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\GSLib\',...
        false,...
        ['12-merPPV Excited State Study - ',num2str(field),'VA']);
    myexp(i).run(true, 60, true);
    i = i + 1;
end
%%
stilenergy = zeros(6,25);
stilgstint = zeros(6,25);
stilfestint = zeros(6,25,3);
sfi = zeros(6,25);

for i = 1:6
    % amp = myexp(i).get_field('ampac.Hf') * 1/23.06;
    indE = myexp(i).get_field('indo.esci',1);
    stilenergy(i,:) = myexp(i).get_field('Eexc',:) + indE;
    stilgstint(i,:) = myexp(i).get_field('Tint',:);
    stilfestint(i,:,:) = myexp(i).get_field('indo.r',2,:,:);
    for j = 1:25;
        sfi(i,j) = (stilenergy(i,j) - stilenergy(i,2)) * sum(stilfestint(i,j,:) .^ 2, 3);
    end
end

%% 
figure(5)
hold on

xaxis = 0:0.01:0.05;

for i = 1:25
    plot(xaxis, stilenergy(:,i), 'g^');
end

maxgsint = max(max(stilgstint));
for i = 1:6
    for j = 1:25
        plot(xaxis(i), stilenergy(i,j), 'bo', 'MarkerSize', (stilgstint(i,j)*30 / maxgsint) + 1e-3);
    end
end

maxfsi = max(max(sfi));

for i = 1:6
    for j = 1:25
        if (sfi(i,j) > 0)
            plot(xaxis(i), stilenergy(i,j), 'rs', 'MarkerSize', (sfi(i,j)*30 / maxfsi) + 1e-3);
        end
    end
end
    