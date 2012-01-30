% root = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\stilbenetest\';
% myece = ECEParams('AM1', {0; 90:5:180}, true, 0, 100, 25);
% myexp = EnergyCalcExp(myece,...
%     [root, 'stilbene.dat'],...
%     [6 14],...
%     [root, 'exp\stilbenefinertestrun.mat'],...
%     [root, 'INDO'],...
%     [root, 'GS']);
% myexp.run(true, 50, true);

%% Ampac Data Conversion
% root = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\twistedphi3\data\';
% newroot = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\GSLib\';
% datfile = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\twistedphi3\data\3-merPPV.dat';
% 
% if (~exist(newroot, 'dir'))
%     mkdir(newroot);
% end
% 
% thiszmat = lower(fileread(datfile));
% 
% % lib = containers.Map('KeyType','char','ValueType','any');
% 
% S = load([newroot, 'catalog.mat'], 'lib');
% lib = S.lib;
% 
% molecule = '3-merPPV';
% 
% key = struct('method',[],'zmat',[],'phi',[]);
% key.method = 'am1';
% 
% files = dir([root, '*.out']);
% fn = cell(1,numel(files));
% fn = {files(1:end).name};
% 
% for i = 1:numel(fn)
%     dataout = struct();
%     phi = [];
%     curfn = fn{i};
%     [tok, ~] = regexpi(curfn, [molecule, '_', key.method, '_([0-9]+).*_([0-9]+).*_([0-9]+).*\.out'], 'tokens');
%     if (~isempty(tok))
%         for j = 1:length(tok{1})
%             phi(end+1) = str2num(cell2mat(tok{1}(j)));
%         end
%     else
%         phi = [];
%     end
%     
%     key.zmat = thiszmat;
% 
%     if (isempty(phi))
%         key.phi = [];
%     else
%         key.phi = reshape(phi, 1, []);
%         for j = 1:length(phi)
%             key.zmat = strrep(key.zmat, ['phi', num2str(j)], num2str(phi(j), '%3.6f'));
%         end
%     end
%     key.zmat = textscan(key.zmat,'%s','delimiter','\n');
%     key.zmat = {key.zmat{1}{4:end-1}};
%     key.zmat = cellfun(@(x)textscan(x, '%s'), key.zmat);
% 
%     for j = 1:length(key.zmat)
%         key.zmat{j}{2} = num2str(str2double(key.zmat{j}{2}),'%3.6f');
%         key.zmat{j}{4} = num2str(str2double(key.zmat{j}{4}),'%3.6f');
%         key.zmat{j}{6} = num2str(str2double(key.zmat{j}{6}),'%3.6f');
%     end
%     
%     hash = DataHash(key);
%     getOut = false;
%     
%     while (~getOut)
%         if (lib.isKey(hash))
%             hash = DataHash(hash);
%         else
%             getOut = true;
%         end
%     end
%     
%     lib(hash) = key;
%     
%     dataout.molecule = molecule;
%     dataout.key = key;
%     dataout.hash = hash;
%     dataout.out = fileread([root, curfn]);
%     delete([root, curfn]);
%     if (exist([root, curfn(1:end-4), '.vis'], 'file'))
%         dataout.vis = fileread([root, curfn(1:end-4), '.vis']);
%         delete([root, curfn(1:end-4), '.vis']);
%     else
%         dataout.vis = [];
%     end
%     if (exist([root, curfn(1:end-4), '.arc'], 'file'))
%         dataout.arc = fileread([root, curfn(1:end-4), '.arc']);
%         delete([root, curfn(1:end-4), '.arc']);
%     else
%         data.arc = [];
%     end
%     if (exist([root, curfn(1:end-4), '.dat'], 'file'))
%         dataout.dat = fileread([root, curfn(1:end-4), '.dat']);
%         delete([root, curfn(1:end-4), '.dat']);
%     else
%         dataout.dat = [];
%     end
%     
%     save([newroot, hash, '.mat'], '-struct', 'dataout');
% end
% 
% save([newroot, 'catalog.mat'], 'lib');
% 
% clear classes;

%% INDO Data Conversion
ampacroot = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\GSLib\';
newroot = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\INDOLib\';
datfile = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\twistedphi3\data\3-merPPV.dat';

if (~exist(newroot, 'dir'))
    mkdir(newroot);
end

zmatin = fileread(datfile);

S = load([ampacroot, 'catalog.mat']);
Alib = S.lib;

% S = load([newroot, 'catalog.mat']);
% Ilib = S.lib;

% Ilib = containers.Map('KeyType','char','ValueType','any');

field = 0;
root = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\3-merPPV\twistedphi3\data\';

files = dir([root, '*.ido']);
fn = cell(1,numel(files));
fn = {files(1:end).name};

for i = 1:numel(fn)
    dataout = struct();
    key = struct();
    dataout.molecule = '3-merPPV';
    key.method = 'am1';
    key.zmat = lower(zmatin);

    phi = [];
    curfn = fn{i};
    [tok, ~] = regexpi(curfn, [dataout.molecule, '_', key.method, '_([0-9]+).*_([0-9]+).*_([0-9]+).*\.ido'], 'tokens');
    
    if (~isempty(tok))
        for j = 1:length(tok{1})
            phi(end+1) = str2num(cell2mat(tok{1}(j)));
        end
    else
        phi = [];
    end

    if (isempty(phi))
        key.phi = [];
    else
        key.phi = reshape(phi, 1, []);
        for j = 1:length(phi)
            key.zmat = strrep(key.zmat, ['phi', num2str(j)], num2str(phi(j), '%3.6f'));
        end
    end

    key.zmat = textscan(key.zmat,'%s','delimiter','\n');
    key.zmat = {key.zmat{1}{4:end-1}};
    key.zmat = cellfun(@(x)textscan(x, '%s'), key.zmat);

    for j = 1:length(key.zmat)
        key.zmat{j}{2} = num2str(str2double(key.zmat{j}{2}),'%3.6f');
        key.zmat{j}{4} = num2str(str2double(key.zmat{j}{4}),'%3.6f');
        key.zmat{j}{6} = num2str(str2double(key.zmat{j}{6}),'%3.6f');
    end

    Ahash = DataHash(key);
    getOut = false;

    while (~getOut)
        if (~isequal(Alib(Ahash),key))
            Ahash = DataHash(Ahash);
        else
            getOut = true;
        end
    end

    key.Efield_mag = field;
    if (field == 0)
        key.Efield_type = 'c';
        key.Efield_params = [0 0 0];
    else
        key.Efield_type = 'p';
        key.Efield_params = [6 37];
    end

    Ihash = DataHash(key);
    getOut = false;
    while (~getOut)
        % if (Ilib.isKey(Ihash))
        if (exist([newroot, Ihash, '.mat'], 'file'))
            Ihash = DataHash(Ihash);
        else
            getOut = true;
        end
    end

    dataout.key = key;
    dataout.hash = Ihash;

    

    if (field ~= 0)
        S = load([ampacroot, Ahash, '.mat'], 'arc','out');
        ampac = parseAmpacFromText(S.arc, S.out);
        myece = ECEParams();
        myece.setEFVectorFromAtomNums(ampac, key.Efield_params);
        configIn = struct('charge',0,'norbs',100,'nstates',25,'field', myece.Efield_vector * field);
    else
        configIn = struct('charge',0,'norbs',100,'nstates',25,'field', [0 0 0]);
    end

    dataout.indo = Indo.LoadExistingData([root, curfn], configIn, ampacroot, Ahash);

    
    % Ilib(Ihash) = dataout.key;
    save([newroot, Ihash, '.mat'], '-struct', 'dataout');
    delete([root, curfn]);
end

% lib = Ilib;
% save([newroot, 'catalog.mat'], 'lib');

clear classes;

%% Fix INDO
% root = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\INDOLib\newlib\';
% Aroot = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\GSLib\';
% files = dir([root, '*.mat']);
% Afiles = dir([Aroot, '*.mat']);
% 
% for i = 1:numel(files)
%     if (~strcmp(files(i).name, 'catalog.mat'))
%         I = load([root, files(i).name]);
%         Ahash = I.indo.jobName;
%         runFullCheck = false;
%         
%         if (exist([Aroot, Ahash, '.mat'], 'file'))
%             A = load([Aroot, Ahash, '.mat']);
%             if (~isequal(struct('molecule',A.key.molecule,'method',A.key.method,'phi',A.key.phi),...
%                     struct('molecule',I.key.molecule,'method',I.key.method,'phi',I.key.phi)))
%                 runFullCheck = true;
%             end
%         else
%             runFullCheck = true;
%         end
%         
%         if (runFullCheck)
%             k = numel(Afiles);
%             j = 1;
%             while (j<k)
%                 if (~strcmp(Afiles(j).name, 'catalog.mat'))
%                     A = load([Aroot, Afiles(j).name]);
%                     if (isequal(struct('molecule',A.key.molecule,'method',A.key.method,'phi',A.key.phi),...
%                             struct('molecule',I.key.molecule,'method',I.key.method,'phi',I.key.phi)))
%                         I.indo.jobName = A.hash;
%                         I.key.zmat = A.key.zmat;
%                         Ihash = DataHash(I.key);
%                         getOut = false;
%                         
%                         while (~getOut)
%                             if (exist([root, '\newnewlib\', Ihash, '.mat'], 'file'))
%                                 Ihash = DataHash(Ihash);
%                             else
%                                 getOut = true;
%                             end
%                         end
%                         
%                         I.hash = Ihash;
%                         
%                         save([root, '\newnewlib\', Ihash, '.mat'], '-struct', 'I');
%                         delete([root, files(i).name]);
%                                 
%                         j=k;
%                     end
%                 end
%                 j=j+1;
%             end
%         end
%     end
% end
%%
% root = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\';
% files = dir([root, 'Exp\*.mat']);
% 
% for i = 1:numel(files)
%     S = load([root,'\Exp\', files(i).name], 'obj');
%     S.obj.indodatapath = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\INDOLib\';
%     S.obj.ampacdatapath = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\GSLib\';
%     S.obj.MATFilename = ['C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\Exp\',S.obj.MATFilename,'.mat'];
%     S.obj.molecule = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\MiddleTwist\12-merPPV.dat';
%     
%     Alib = load([root, 'GSLib\catalog.mat']);
%     S.obj.ampac_lib = Alib.lib;
%     Ilib = load([root, 'INDOLib\catalog.mat']);
%     S.obj.indo_lib = Ilib.lib;
%     
%     S.obj.run(true, 50, true);
% end
%%
% root = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\INDOLib\';
% newroot = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\12-merPPV\INDOLib\newlib\';
% 
% files = dir([root, '*.mat']);
% catalog = containers.Map('KeyType','char','ValueType','any');
% 
% for i = 1:numel(files)
%     S = load([root, files(i).name]);
%     S.key = rmfield(S.key, 'molecule');
%     S.molecule = '12-merPPV';
%     
%     hash = DataHash(S.key);
%     getOut = false;
%     
%     while (~getOut)
%         if (exist([newroot, hash, '.mat'], 'file'))
%             hash = DataHash(hash);
%         else
%             getOut = true;
%         end
%     end
%     
%     S.hash = hash;
%     S.indo.jobName = alookup(S.indo.jobName);
%     save([newroot, hash, '.mat'], '-struct', 'S');
%     catalog(hash) = S.key;
% end
% 
% lib = catalog;
% save([newroot, 'catalog.mat'], 'lib');
% 
% clear classes;