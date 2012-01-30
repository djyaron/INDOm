% for fs = 0.138:0.001:0.2
%     myece = ECEParams('AM1', {}, true, fs, 800, 25);
%     myexp = EnergyCalcExp(myece,...
%         's:\NoAngle\8-merPPV.dat',...
%         [6 62],...
%         ['s:\NoAngle\ExploreDipoleNewAlg\8-merPPV-NewAlg-',num2str(fs),'VA.mat'],...
%         'o:\',...
%         'p:\',...
%         false,...
%         ['Explore Dipole with new HF Algorithm - ',num2str(fs),'VA']);
% %     for phi2 = [90 93 135 180]
% %         thiskey = myexp.keygen('indo',[0 phi2]);
% %         if (exist(['o:\',DataHash(thiskey),'.mat'],'file'))
% %             load(['o:\',DataHash(thiskey),'.mat'],'indo','key');
% %             if (isequal(key,thiskey))
% %                 % if (indo.indo_succeed == 0)
% %                     delete(['o:\',DataHash(thiskey),'.mat']);
% %                 % end
% %             else
% %                 load(['o:\',DataHash(DataHash(thiskey)),'.mat'],'indo','key');
% %                 if (isequal(key,thiskey))
% %                     % if (indo.indo_succeed == 0)
% %                         delete(['o:\',DataHash(DataHash(thiskey)),'.mat']);
% %                     % end
% %                 else
% %                     keyboard;
% %                 end
% %             end
% %         end
% %     end
%     myexp.run('quiet');
% end
% 
% %% Stilbene
% xaxis = [0:0.001:0.166];
% num_ang = 1;
% 
% dp = zeros(num_ang,length(xaxis));
% dpgs = zeros(num_ang,length(xaxis));
% hlgap = zeros(num_ang,length(xaxis));
% eexc = zeros(num_ang,length(xaxis),25);
% igs = zeros(num_ang,length(xaxis));
% opint = zeros(num_ang,length(xaxis),25);
% tpint = zeros(num_ang,length(xaxis),25);
% dp2exc = zeros(num_ang,length(xaxis));
% dp3exc = zeros(num_ang,length(xaxis));
% dp4exc = zeros(num_ang,length(xaxis));
% 
% fidx = 1;
% for fs = xaxis
%     myexp = EnergyCalcExp.ReadMATFile(['s:\noangle\ExploreDipoleNewAlg\8-merPPV-NewAlg-',num2str(fs),'VA.mat'], false);
%    
%     for i = 1:num_ang
%         myexp.data(i).load_to_memory('indo','load');
%     end
%     
%     pidx = 1;
%     for phi2 = [90 93 135 180]
%             
%         tmp = myexp.get_field('indo.dipole',1,1,0,phi2);
%         dpgs(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
% 
%         tmp = myexp.get_field('indo.dipole',2,2,0,phi2);
%         dp(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
% 
%         tmp = myexp.get_field('indo.dipole',3,3,0,phi2);
%         dp2exc(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
%         
%         tmp = myexp.get_field('indo.dipole',4,4,0,phi2);
%         dp3exc(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
%         
%         tmp = myexp.get_field('indo.dipole',5,5,0,phi2);
%         dp4exc(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
% 
%         nfill = myexp.get_field('indo.nfilled',0,phi2);
%         tmp = myexp.get_field('indo.orbE',[nfill nfill+1],0,phi2);
%         hlgap(pidx,fidx) = tmp(2) - tmp(1);
% 
%         eexc(pidx,fidx,1:25) = myexp.get_field('Eexc',:,0,phi2);
%         igs(pidx,fidx) = myexp.get_field('indo.esci',1,0,phi2);
%         opint(pidx,fidx,1:25) = myexp.get_field('Tint',:,0,phi2);
%         tmp = squeeze(myexp.get_field('indo.r',2,:,:,0,phi2));
%         for l = 1:25
%             tpint(pidx,fidx,l) = (eexc(pidx,fidx,l) - eexc(pidx,fidx,2)) * sum(tmp(l,:) .^ 2);
%         end
%         
%         pidx = pidx + 1;
%     end
%     
%     fidx = fidx + 1;
% end

%% 8-merPPV - No Twist

xaxis = [0:0.001:0.166];
num_ang = 1;

dp = zeros(num_ang,length(xaxis));
dpgs = zeros(num_ang,length(xaxis));
hlgap = zeros(num_ang,length(xaxis));
eexc = zeros(num_ang,length(xaxis),25);
igs = zeros(num_ang,length(xaxis));
opint = zeros(num_ang,length(xaxis),25);
tpint = zeros(num_ang,length(xaxis),25);
dp2exc = zeros(num_ang,length(xaxis));
dp3exc = zeros(num_ang,length(xaxis));
dp4exc = zeros(num_ang,length(xaxis));
dpmag = zeros(num_ang,length(xaxis));

fidx = 1;
for fs = xaxis
    myexp = EnergyCalcExp.ReadMATFile(['s:\noangle\ExploreDipoleNewAlg\8-merPPV-NewAlg-',num2str(fs),'VA.mat'], false);
   
    for i = 1:num_ang
        myexp.data(i).load_to_memory('indo','load');
    end
    
    pidx = 1;
    % for phi2 = [90 93 135 180]
            
        tmp = myexp.get_field('indo.dipole',1,1);
        dpgs(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

        tmp = myexp.get_field('indo.dipole',2,2);
        dp(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

        tmp = myexp.get_field('indo.dipole',3,3);
        dp2exc(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
        
        tmp = myexp.get_field('indo.dipole',4,4);
        dp3exc(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
        
        tmp = myexp.get_field('indo.dipole',5,5);
        dp4exc(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);
        
        tmp = myexp.get_field('indo.dipole',25,25);
        dpmag(pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

        nfill = myexp.get_field('indo.nfilled');
        tmp = myexp.get_field('indo.orbE',[nfill nfill+1],:);
        hlgap(pidx,fidx) = tmp(2) - tmp(1);

        eexc(pidx,fidx,1:25) = myexp.get_field('Eexc',:,:,:);
        igs(pidx,fidx) = myexp.get_field('indo.esci',1,:);
        opint(pidx,fidx,1:25) = myexp.get_field('Tint',:,:);
        [~,opmax] = max(opint(pidx,fidx,:));
        tmp = squeeze(myexp.get_field('indo.r',opmax,:,:,:));
        for l = 1:25
            tpint(pidx,fidx,l) = (eexc(pidx,fidx,l) - eexc(pidx,fidx,2)) * sum(tmp(l,:) .^ 2);
        end
        
        % pidx = pidx + 1;
    % end
    
    fidx = fidx + 1;
end

%%
figure(9)
hold on

j = 1;

for k = 1:length(xaxis)
    plot(xaxis(k), reshape(eexc(j,k,:),1,[]) + igs(j,k), 'g^');
end

maxopint = max(opint(j,:));
for k = 1:length(xaxis)
    for l = 1:25
        plot(xaxis(k), eexc(j,k,l) + igs(j,k), 'bo', 'MarkerSize', (opint(j,k,l)*30 / maxopint) + 1e-3);
    end
end

maxtpint = max(tpint(j,:));

for k = 1:length(xaxis)
    for l = 1:25
        if (tpint(j,k,l) > 0)
            plot(xaxis(k), eexc(j,k,l) + igs(j,k), 'rs', 'MarkerSize', (tpint(j,k,l)*30 / maxtpint) + 1e-3);
        end
    end
end

%%
figure(17);
hold on;
% xaxis = 0:0.1:4;
k = 1;
%scatter(xaxis,reshape(dpgs(k,1:length(xaxis)),1,[]),'MarkerEdgeColor','g','Marker','*');
scatter(xaxis,reshape(dp(k,1:length(xaxis)),1,[]),'MarkerEdgeColor','b','Marker','o');
%scatter(xaxis,reshape(dp2exc(k,1:length(xaxis)),1,[]),'MarkerEdgeColor','r','Marker','+');
%scatter(xaxis,reshape(dp3exc(k,1:length(xaxis)),1,[]),'MarkerEdgeColor','k','Marker','d');
%scatter(xaxis,reshape(dp4exc(k,1:length(xaxis)),1,[]),'MarkerEdgeColor','c','Marker','^');
scatter(xaxis,reshape(dpmag(k,1:length(xaxis)),1,[]),'MarkerEdgeColor','m','Marker','x');