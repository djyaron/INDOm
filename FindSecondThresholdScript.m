secondthreshold = zeros(31,31);
i = 1;
pow = [10 5 1 0.1 0.01];
for phi1 = 0:3:90
    j = 1;
    for phi2 = 90:3:180
        field = 5;
        powidx = 1;
        getOut = false;
        while (~getOut)
            myece = ECEParams('AM1', {phi1; phi2}, true, field, 100, 25);
            myexp = EnergyCalcExp(myece,...
                's:\2 angles\stilbene.dat',...
                [6 14],...
                ['s:\2 angles\SecondThreshold\stilbene-SecondThreshold-',num2str(field),'VA.mat'],...
                'o:\',...
                'p:\',...
                false,...
                ['Stilbene Second Threshold Test - ',num2str(field),'VA']);
            myexp.run('quiet');
            if (myexp.fail_iterator ~= 1)
                field = field + (pow(powidx));
            elseif (myexp.fail_iterator == 1 && powidx < 5)
                field = field - (pow(powidx));
                powidx = powidx + 1;
                field = field + (pow(powidx));
            elseif (myexp.fail_iterator == 1 && powidx == 5)
                getOut = true;
            end
        end
        secondthreshold(i,j) = field;
        j = j + 1;
    end
    i = i + 1;
end
phi1 = 0:3:90;
phi2 = 90:3:180;
desc = 'secondthreshold is (phi1,phi2)'; 
save('s:\2 angles\SecondThreshold\secondthreshold.mat', 'secondthreshold', 'phi1', 'phi2', 'desc');
        