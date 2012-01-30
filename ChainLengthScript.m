cl_params = [3, 22; 4, 30; 5, 38; 6, 46; 7, 54];

for chain_id = 4:5
%     fid = fopen('C:\Documents and Settings\Linda Group\My Documents\Dropbox\dyes\christian\remote_job_progress.txt','a+');
%     fprintf(fid,'%s',[num2str(cl_params(chain_id)),'-merPPV job started at ',datestr(now)]);
%     fclose(fid);
%     fid = fopen('C:\Documents and Settings\Linda Group\My Documents\Dropbox\dyes\christian\remote_job_progress.txt','a+');
    system('subst s: /d');
    system('subst o: /d');
    system('subst p: /d');
    
    rootdir = ['C:\Users\Christian\Documents\Research\Yaron\dyes2\data\',...
        num2str(cl_params(chain_id,1)),'-merPPV\'];
    SFolder = 'Exp';
    OFolder = 'INDOLib';
    PFolder = 'GSLib';

    system(['subst s: "',rootdir,SFolder,'"']);
    system(['subst o: "',rootdir,OFolder,'"']);
    system(['subst p: "',rootdir,PFolder,'"']);
    for fs = 0:0.001:0.2
        myece = ECEParams('AM1', {}, true, fs, 2000, 25);
        myexp = EnergyCalcExp(myece,...
            ['s:\NoAngle\',num2str(cl_params(chain_id,1)),'-merPPV.dat'],...
            [3 cl_params(chain_id,2)],...
            ['s:\NoAngle\ExploreDipoleNewAlg\',num2str(cl_params(chain_id,1)),'-merPPV-NewAlg-',num2str(fs),'VA.mat'],...
            'o:\',...
            'p:\',...
            false,...
            ['Explore Dipole with new HF Algorithm - ',num2str(fs),'VA']);
        myexp.run('quiet');
    end
%     fprintf(fid,'%s',[num2str(cl_params(chain_id)),'-merPPV job completed at ',datestr(now)]);
%     fclose(fid);
end

%% Collect data

xaxis = [0:0.001:0.2];
num_ang = 1;
num_chains = 6;

dp = NaN(num_chains,num_ang,length(xaxis));
dpgs = NaN(num_chains,num_ang,length(xaxis));
hlgap = NaN(num_chains,num_ang,length(xaxis));
eexc = NaN(num_chains,num_ang,length(xaxis),25);
igs = NaN(num_chains,num_ang,length(xaxis));
opint = NaN(num_chains,num_ang,length(xaxis),25);
tpint = NaN(num_chains,num_ang,length(xaxis),25);
dp2exc = NaN(num_chains,num_ang,length(xaxis));
dp3exc = NaN(num_chains,num_ang,length(xaxis));
dp4exc = NaN(num_chains,num_ang,length(xaxis));
dpmag = NaN(num_chains,num_ang,length(xaxis));

cidx = 1;
for chain_len = 3:8
    fidx = 1;
    system('subst s: /d');
    system('subst o: /d');
    system('subst p: /d');

    rootdir = ['C:\Users\Christian\Documents\Research\Yaron\dyes2\data\',...
        num2str(chain_len),'-merPPV\'];
    SFolder = 'Exp';
    OFolder = 'INDOLib';
    PFolder = 'GSLib';

    system(['subst s: "',rootdir,SFolder,'"']);
    system(['subst o: "',rootdir,OFolder,'"']);
    system(['subst p: "',rootdir,PFolder,'"']);
    
    for fs = xaxis
        
        myexp = EnergyCalcExp.ReadMATFile(['s:\noangle\ExploreDipoleNewAlg\',...
            num2str(chain_len),'-merPPV-NewAlg-',num2str(fs),'VA.mat'], false);

        for i = 1:num_ang
            myexp.data(i).load_to_memory('indo','load');
        end

        pidx = 1;
        % for phi2 = [90 93 135 180]
        if (myexp.data(i).indo_succeed)
            tmp = myexp.get_field('indo.dipole',1,1);
            dpgs(cidx,pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',2,2);
            dp(cidx,pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',3,3);
            dp2exc(cidx,pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',4,4);
            dp3exc(cidx,pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',5,5);
            dp4exc(cidx,pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

            tmp = myexp.get_field('indo.dipole',25,25);
            dpmag(cidx,pidx,fidx) = sum(tmp .^ 2, 1) .^ (0.5);

            nfill = myexp.get_field('indo.nfilled');
            tmp = myexp.get_field('indo.orbE',[nfill nfill+1],:);
            hlgap(cidx,pidx,fidx) = tmp(2) - tmp(1);

            tmp = myexp.get_field('Eexc',:,:,:);
            if (~any(isnan(tmp)))
                eexc(cidx,pidx,fidx,1:25) = tmp;
            end
            tmp = myexp.get_field('indo.esci',1,:);
            if (~any(isnan(tmp)))
                igs(cidx,pidx,fidx) = tmp;
            end
            tmp = myexp.get_field('Tint',:,:);
            if (~any(isnan(tmp)))
                opint(cidx,pidx,fidx,1:25) = tmp;
                [~,opmax] = max(opint(cidx,pidx,fidx,:));
                tmp = squeeze(myexp.get_field('indo.r',opmax,:,:,:));
                for l = 1:25
                    tpint(cidx,pidx,fidx,l) = (eexc(cidx,pidx,fidx,l) - eexc(cidx,pidx,fidx,2)) * sum(tmp(l,:) .^ 2);
                end
            end
        end
            

            % pidx = pidx + 1;
        % end

        fidx = fidx + 1;
    end
    cidx = cidx + 1;
end
%%
figure(9)
hold on

i = 6;      % Chain length index
j = 1;      % Phi index (always 1)

xaxis = 0:0.001:0.2;

for k = 1:length(xaxis)
    plot(xaxis(k), reshape(eexc(i,j,k,:),1,[]) + igs(i,j,k), 'g^');
end

% maxopint = max(max(opint(i,j,1:length(xaxis),:)));
% for k = 1:length(xaxis)
%     for l = 1:25
%         if (~isnan(opint(i,j,k,l)) && opint(i,j,k,l) >= 0)
%             plot(xaxis(k), eexc(i,j,k,l) + igs(i,j,k), 'bo', 'MarkerSize', (opint(i,j,k,l)*30 / maxopint) + 1e-3);
%         end
%     end
% end
% 
% maxtpint = max(max(tpint(i,j,1:length(xaxis),:)));
% 
% for k = 1:length(xaxis)
%     for l = 1:25
%         if (tpint(i,j,k,l) >= 0 && ~isnan(tpint(i,j,k,l)))
%             plot(xaxis(k), eexc(i,j,k,l) + igs(i,j,k), 'rs', 'MarkerSize', (tpint(i,j,k,l)*30 / maxtpint) + 1e-3);
%         end
%     end
% end

topint = log10(opint+1e-3);
minopint = min(min(topint(i,j,1:length(xaxis),:)));

topint = topint + abs(minopint);

maxopint = max(max(topint(i,j,1:length(xaxis),:)));


for k = 1:length(xaxis)
    for l = 1:25
        if (~isnan(topint(i,j,k,l)) && topint(i,j,k,l) >= 0)
            plot(xaxis(k), eexc(i,j,k,l) + igs(i,j,k), 'bo', 'MarkerSize', (topint(i,j,k,l)*20 / maxopint) + 1e-3);
        end
    end
end


ttpint = log10(tpint+1e-3);
ttpint = ttpint .* ~imag(ttpint);
mintpint = max(-1*min(ttpint(i,j,1:length(xaxis),:)));
ttpint = ttpint + abs(mintpint);
maxtpint = max(max(ttpint(i,j,1:length(xaxis),:)));

try
for k = 1:length(xaxis)
    for l = 1:25
        if (ttpint(i,j,k,l) >= 0 && ~isnan(ttpint(i,j,k,l)))
            plot(xaxis(k), eexc(i,j,k,l) + igs(i,j,k), 'rs', 'MarkerSize', (ttpint(i,j,k,l)*20 / maxtpint) + 1e-3);
        end
    end
end
catch exc
    keyboard
end
%%
figure(17);
hold on;
% xaxis = 0:0.1:4;
i = 1;
k = 1;

scatter(xaxis,reshape(dpgs(i,k,1:length(xaxis)),1,[]),'MarkerEdgeColor','g','Marker','*');
scatter(xaxis,reshape(dp(i,k,1:length(xaxis)),1,[]),'MarkerEdgeColor','b','Marker','o');
scatter(xaxis,reshape(dp2exc(i,k,1:length(xaxis)),1,[]),'MarkerEdgeColor','r','Marker','+');
scatter(xaxis,reshape(dp3exc(i,k,1:length(xaxis)),1,[]),'MarkerEdgeColor','k','Marker','d');
scatter(xaxis,reshape(dp4exc(i,k,1:length(xaxis)),1,[]),'MarkerEdgeColor','c','Marker','^');
scatter(xaxis,reshape(dpmag(i,k,1:length(xaxis)),1,[]),'MarkerEdgeColor','m','Marker','x');