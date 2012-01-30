maxdata = zeros(31,31);
i = 1;
for phi1 = 0:3:90
    j = 1;
    for phi2 = 90:3:180
        field = 0;
        pow = 0;
        getOut = false;
        while (~getOut)
            myece = ECEParams('AM1', {phi1; phi2}, true, field, 100, 25);
            myexp = EnergyCalcExp(myece,...
                's:\2 angles\stilbene.dat',...
                [6 14],...
                ['s:\2 angles\MaxFieldTest\stilbene-MaxFieldTest-',num2str(field),'VA.mat'],...
                'o:\',...
                'p:\',...
                false,...
                ['Stilbene MaxFieldTest - ',num2str(field),'VA']);
            myexp.run(true, 60, true);
            if (myexp.fail_iterator == 1)
                field = field + (1 * 10^(pow));
            elseif (myexp.fail_iterator > 1 && pow > -2)
                field = field - (1 * 10^(pow));
                pow = pow - 1;
                field = field + (1 * 10^(pow));
            elseif (myexp.fail_iterator > 1 && pow == -2)
                getOut = true;
            end
        end
        maxdata(i,j) = field - 0.01;
        j = j + 1;
    end
    i = i + 1;
end
phi1 = 0:3:90;
phi2 = 90:3:180;
desc = 'Maxdata is (phi1,phi2)'; 
save('s:\2 angles\MaxFieldTest\maxdata.mat', 'maxdata', 'phi1', 'phi2', 'desc');
        