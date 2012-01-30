boandrdata = cell(109,5);
idx = 1;

% for fs = 0:0.001:0.2
while false
    myece = ECEParams('AM1', {}, true, fs, 800, 25,...
        's:\NoAngle\params.txt');
    
    if (fs ~= 0)
        myece.dm_guess = ['o:\', myexp.data(1).indo_hash, '-dm.bin'];
    end
    
    myexp = EnergyCalcExp(myece,...
        's:\NoAngle\8-merPPV.dat',...
        [6 62],...
        ['s:\NoAngle\8-merPPV-',num2str(fs),'VA.mat'],...
        'o:\',...
        'p:\',...
        false,...
        ['Bond Order calcs - ',num2str(fs),'VA']);
    myexp.run('quiet');
    
    dat = myexp.data(1).indoOutput;
    dat = textscan(dat,'%s','delimiter','\n');
    dat = dat{1};
    dat = textscan(dat{2},'%s');
    dat = dat{1};
    
    efield = [str2double(dat{1}), str2double(dat{2}), str2double(dat{3})];
    
    myexp.data(1).load_to_memory('indo','load');
    indo = myexp.data(1).raw_indo;
    myexp.data(1).generate_ampac_file('out', 'p:\');
    
    [rgs, rex, deltar, ~, numruns] = OptExcStateStructure('p:\', myexp.data(1).ampac_hash, 'o:\', myexp.data(1).indo_hash, efield, indo);
    
    boandrdata{idx, 1} = fs;
    boandrdata{idx, 2} = rgs;
    boandrdata{idx, 3} = rex;
    boandrdata{idx, 4} = deltar;
    boandrdata{idx, 5} = numruns;
    idx = idx + 1;
end

% keyboard;
%%
nidx = 1;
AmpacEXE = 'C:\Program Files\Semichem, Inc\Ampac-9.2\ampac.exe';
    OpenBabelEXE = 'C:\Program Files (x86)\OpenBabel-2.3.1\babel.exe';

    load(['s:\NoAngle\8-merPPV-0VA.mat'],'obj');
    
    obj.data(1).generate_ampac_file('out','p:\');
    obj.data(1).load_to_memory('indo','load');
    indo = obj.data(1).raw_indo;
    
for fs = 0:0.001:0.175
    
    
    load(['s:\NoAngle\8-merPPV-',num2str(fs),'VA.mat'],'obj');
    
    obj.data(1).generate_ampac_file('out','p:\');
    
    ampac_pathonly = 'p:\';
    indo_pathonly = 'o:\';
    ampac_nameonly = obj.data(1).ampac_hash;
    indo_nameonly = obj.data(1).indo_hash;
    
    ampac_filepath = [ampac_pathonly, ampac_nameonly];
    indo_filepath = [indo_pathonly, indo_nameonly];
    % copyfile([ampac_filepath, '.out'],[ampac_filepath, '-new.out']);
    % copyfile([indo_filepath, '-dm.bin'],[indo_filepath, '-new-dm.bin']);
    
%     if (length(varargin) > 1)
%         indo = varargin{2};
%         indo_nameonly = [indo_nameonly, '-new'];
%         indo_filepath = [indo_pathonly, indo_nameonly];
%     else
%         copyfile([indo_filepath, '.ido'],[indo_filepath, '-new.ido']);
%         indo_nameonly = [indo_nameonly, '-new'];
%         indo_filepath = [indo_pathonly, indo_nameonly];
%         indo = Indo.LoadExistingData([indo_filepath,'.ido'],[],[],[]);
%     end
%        
%     ampac_nameonly = [ampac_nameonly, '-new'];
%     ampac_filepath = [ampac_pathonly, ampac_nameonly];
    
    
    getOut = false;
    oldbls = [];
    numruns = 0;

%    while (~getOut)
    
    fdat = fileread([ampac_filepath, '.out']);

    zmat = textscan(fdat,'%s','delimiter','\n');
    zmat = zmat{1};

    idx = 0;
    i = 1;
    foundonce = false;

    while (idx == 0)
        %if (strcmp(strtrim(zmat{i}), 'GEOMETRY OPTIMISED : ENERGY MINIMISED'))
            %while (idx == 0)
                if (strcmp(strtrim(zmat{i}), '(I)                   NA:I          NB:NA:I      NC:NB:NA:I     NA    NB    NC'))
                    if (foundonce)
                        idx = i+1;
                        eidx = 0;
                        i = i+2;
                        while (eidx == 0)
                            if (strcmp(zmat{i},''))
                                eidx = i;
                            else
                                i = i + 1;
                            end
                        end
                    else
                        foundonce = true;
                        i = i+1;
                    end
                else
                    i = i+1;
                end
            %end
        %else
            %i = i + 1;
        %end
    end

    natoms = eidx - idx;
    zmat = {zmat{idx:eidx-1}};
    zmatrix = cell(natoms, 8);

    temp = textscan(zmat{1}, '%s');
    temp = temp{1};
    zmatrix{1,1} = temp{1};
    zmatrix{1,2} = temp{2};
    zmatrix{1,3} = '';
    zmatrix{1,4} = '';
    zmatrix{1,5} = '';
    zmatrix{1,6} = '';
    zmatrix{1,7} = '';
    zmatrix{1,8} = '';

    temp = textscan(zmat{2}, '%s');
    temp = temp{1};
    zmatrix{2,1} = temp{1};
    zmatrix{2,2} = temp{2};
    zmatrix{2,3} = temp{3};
    zmatrix{2,4} = '';
    zmatrix{2,5} = '';
    zmatrix{2,6} = temp{5};
    zmatrix{2,7} = '';
    zmatrix{2,8} = '';

    temp = textscan(zmat{3}, '%s');
    temp = temp{1};
    zmatrix{3,1} = temp{1};
    zmatrix{3,2} = temp{2};
    zmatrix{3,3} = temp{3};
    zmatrix{3,4} = temp{5};
    zmatrix{3,5} = '';
    zmatrix{3,6} = temp{7};
    zmatrix{3,7} = temp{8};
    zmatrix{3,8} = '';

    for i = 4:natoms
        temp = textscan(zmat{i}, '%s');
        temp = temp{1};
        zmatrix{i,1} = temp{1};
        zmatrix{i,2} = temp{2};
        zmatrix{i,3} = temp{3};
        zmatrix{i,4} = temp{5};
        zmatrix{i,5} = temp{7};
        zmatrix{i,6} = temp{9};
        zmatrix{i,7} = temp{10};
        zmatrix{i,8} = temp{11};
    end
    
    newzmatrix = zmatrix;
    
    %% Write MOPAC input file and convert to SYBYL Mol2 Format

    fid = fopen([ampac_filepath, '.mopin'],'w');

    fprintf(fid, '%s\r\n\r\n\r\n', 'INSERT KEYWORDS HERE');
    % fprintf(fid, '%s\r\n\r\n', 'Title');
    % '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n'
    fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{1,2}, '0.000000', '1', '0.000000', '1', '0.000000', '1', '0','0','0');

    if (natoms > 1)
        fprintf(fid,'%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{2,2}, zmatrix{2,3}, '1',...
            '0.000000', '1', '0.000000', '1', zmatrix{2,6},'0','0');
    end
    if (natoms > 2)
        fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{3,2}, zmatrix{3,3}, '1',...
            zmatrix{3,4}, '1', '0.000000', '1', zmatrix{3,6}, zmatrix{3,7},'0');
    end
    if (natoms > 3)
        for i = 4:natoms
            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', zmatrix{i,2}, zmatrix{i,3}, '1',...
                zmatrix{i,4}, '1', zmatrix{i,5}, '1', zmatrix{i,6}, zmatrix{i,7}, zmatrix{i,8});
        end
    end

    fclose(fid);

    [~, ~] = system(['"', OpenBabelEXE, '" -imopin "', ampac_filepath, '.mopin" -omol2 "', ampac_filepath, '.mol2"']);

    %% Get all bonds

    mf = fileread([ampac_filepath, '.mol2']);

    mf = textscan(mf,'%s','delimiter','\n');
    mf = mf{1};

    bidx = find(cellfun(@(x)strcmp(x, '@<TRIPOS>BOND'), mf)) + 1;

    mf = {mf{bidx:end}};

    mf = cellfun(@(x)textscan(x,'%s'), mf);

    atomsbonds = cell(natoms,1);

    for i = 1:length(mf)
        atomsbonds{str2double(mf{i}(2))}(end+1) = str2double(mf{i}(3));
        atomsbonds{str2double(mf{i}(3))}(end+1) = str2double(mf{i}(2));
    end

    %% Get rings

    rings = {};     % known rings
    ringset = [];   % current path we're following
    atomids = [0];  % atom numbers of known ring atoms
    gotaring = false;   % bool to let us get out of our current path if we've found a ring

    for i = 1:length(atomsbonds)    % loop over all atom numbers
        if (~ismember(i, atomids))  % check to see if atom is one of the known ring atoms
            ringset(end+1) = i;     % First atom in path
            for j = 1:length(atomsbonds{i}) % Loop over all paths from atom i
                if (length(atomsbonds{atomsbonds{i}(j)}) > 1)   % Check for backtracking
                    ringset(end+1) = atomsbonds{i}(j);  % Select a path from i
                    if (~gotaring)  % Keep going if we're not trying to get out from a found ring
                        for k = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom j
                            if (atomsbonds{ringset(end)}(k) ~= i && ~gotaring)  % Check for backtracking & getting out
                                ringset(end+1) = atomsbonds{ringset(end)}(k); % Select path from j
                                for l = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom k
                                    if (atomsbonds{ringset(end)}(l) ~= ringset(end-1) && ~gotaring) % Check for backtracking
                                        ringset(end+1) = atomsbonds{ringset(end)}(l);   % Select path from k
                                        for m = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom l
                                            if (atomsbonds{ringset(end)}(m) ~= ringset(end-1) && ~gotaring) % Check for backtracking
                                                ringset(end+1) = atomsbonds{ringset(end)}(m);   % Select path from l
                                                for n = 1:length(atomsbonds{ringset(end)})  % Loop over all paths from atom m
                                                    if (atomsbonds{ringset(end)}(n) ~= ringset(end-1))  % Check for backtracking
                                                        if (~ismember(i, atomsbonds{ringset(end)})) % Check to see if any path fom m is to atom i
                                                            if (~gotaring)  % If not, check to get out and continue
                                                                    ringset(end+1) = atomsbonds{ringset(end)}(n); % Select a path from m
                                                                    if (ismember(i, atomsbonds{ringset(end)}) && atomsbonds{ringset(end)}(n) ~= ringset(end-1)) % Check if any path from atom n is to atom i
                                                                        rings{end+1} = ringset; % If so, we found a ring, record it
                                                                        atomids = unique([atomids ringset]);
                                                                        gotaring = true; % and get out
                                                                    end
                                                                    ringset = ringset(1:end-1); % remove n from current path
                                                            end
                                                        elseif (~gotaring)  % if we're not trying to get out and we found a 5-membered ring, record it and get out
                                                            rings{end+1} = ringset;
                                                            atomids = unique([atomids ringset]);
                                                            gotaring = true;
                                                        end
                                                    end
                                                end 
                                                ringset = ringset(1:end-1); % remove m from current path
                                            end
                                        end
                                        ringset = ringset(1:end-1); % remove l from current path
                                    end
                                end
                                ringset = ringset(1:end-1); % remove k from current path
                            end
                        end

                    end
                    ringset = ringset(1:end-1); % remove j from current path
                end
            end
        end
        ringset = ringset(1:end-1); % remove i from current path
        gotaring = false; % start looking again, excluding atoms we marked as in a ring already
    end

    %% Get z-matrix angles

    ringonly = [];
    ringinc = [];
    noring = [];

    for i = 1:natoms-2
        aset = [str2double(zmatrix{i+2,1}) str2double(zmatrix{i+2,6}) str2double(zmatrix{i+2,7})];
        switch (length(intersect(aset, atomids)))
            case 3
                ringonly(end+1,:) = aset;
            case 2
                ringinc(end+1,:) = aset;
            otherwise
                noring(end+1,:) = aset;          
        end
    end

    allbanums = [ringonly; ringinc; noring];

    %% Get all C-C and C-H bonds

    CCbonds = {};
    CHbonds = {};

    for i = 1:length(mf)
        if (strcmpi(zmatrix{str2double(mf{i}(2)),2}, 'c') && strcmpi(zmatrix{str2double(mf{i}(3)),2}, 'c'))
            myset = sort([str2double(mf{i}(2)) str2double(mf{i}(3))]);
            if (~any(cellfun(@(x)isequal(x, myset), CCbonds)))
                CCbonds{end+1} = myset;
            end
        elseif ((strcmpi(zmatrix{str2double(mf{i}(2)),2}, 'c') && strcmpi(zmatrix{str2double(mf{i}(3)),2}, 'h')) || ...
                (strcmpi(zmatrix{str2double(mf{i}(2)),2}, 'h') && strcmpi(zmatrix{str2double(mf{i}(3)),2}, 'c')))
            myset = sort([str2double(mf{i}(2)) str2double(mf{i}(3))]);
            if (~any(cellfun(@(x)isequal(x, myset), CHbonds)))
                CHbonds{end+1} = myset;
            end
        end
    end

    temp = reshape(cell2mat(CCbonds),2, length(CCbonds))';
    [b,ix] = sort(temp,1);
    CCbonds = [b(:,1), temp(ix(:,1), 2)];

    temp = reshape(cell2mat(CHbonds),2, length(CHbonds))';
    [b,ix] = sort(temp,1);
    CHbonds = [b(:,1), temp(ix(:,1), 2)];
    
    
    oldbls = cellfun(@(x)str2double(x),{zmatrix{2:end,3}});
    rgs = oldbls;
    rex = rgs;
    gsbo = [];
    gsbobyatom = [];

    % while (~getOut)
        numruns = numruns + 1;
        % newzmatrix = zmatrix;
        
        %% Read in ground state density matrix file

        fid = fopen([indo_filepath,'-dm.bin']);

        matdim = fread(fid,1,'int');

        dm = fread(fid,[matdim,matdim],'double');

        fclose(fid);
        
        fid = fopen([indo_filepath,'-new-dm.bin']);

        matdim = fread(fid,1,'int');

        esdm = fread(fid,[matdim,matdim],'double');

        fclose(fid);

        %% Calculate bond order of C-C bonds

        gsbo = size(CCbonds,2);
        gsbobyatom = zeros(natoms);

        for i = 1:size(CCbonds,1)
            blmat = dm(indo.aorbAtom == CCbonds(i,1), indo.aorbAtom == CCbonds(i,2)) .^ 2;
            gsbo(i) = sum(blmat(:));
            gsbobyatom(CCbonds(i,1), CCbonds(i,2)) = gsbo(i);
            gsbobyatom(CCbonds(i,2), CCbonds(i,1)) = gsbo(i);
        end
        
        esbo = size(CCbonds,2);
        esbobyatom = zeros(natoms);

        for i = 1:size(CCbonds,1)
            blmat = esdm(indo.aorbAtom == CCbonds(i,1), indo.aorbAtom == CCbonds(i,2)) .^ 2;
            esbo(i) = sum(blmat(:));
            esbobyatom(CCbonds(i,1), CCbonds(i,2)) = esbo(i);
            esbobyatom(CCbonds(i,2), CCbonds(i,1)) = esbo(i);
        end
        
        newblsbyatom = 1.54-0.3*log(esbobyatom);  % Pauling bond order
        newzmatrix = zmatrix;
        for i = 2:natoms
            if (strcmpi(zmatrix{i,2},'c'))
                newzmatrix{i,3} = num2str(newblsbyatom(i,str2double(newzmatrix{i,6})));
            end
        end
        
        newbls = cellfun(@(x)str2double(x),{newzmatrix{2:end,3}});
        
        
        %% Get current bond lengths and check them
        
        getOut = true;
            % deltar = rex - rgs;
            deltar = zeros(natoms-1,1);
            atom_nums_in_bonds = zeros(natoms-1,2);
            
            for (i = 1:natoms-1)
                trgs = gsbobyatom(str2double(zmatrix{i+1,1}),str2double(zmatrix{i+1,6}));
                trex = esbobyatom(str2double(zmatrix{i+1,1}),str2double(zmatrix{i+1,6}));
                atom_nums_in_bonds(i,1:2) = [str2double(zmatrix{i+1,1}), str2double(zmatrix{i+1,6})];
                deltar(i) = 0.3*log(trgs / trex);
            end
            
            rgs = rgs';
            rex = rex';
            
            boandrdata{nidx, 1} = fs;
            boandrdata{nidx, 2} = rgs;
            boandrdata{nidx, 3} = rex;
            boandrdata{nidx, 4} = deltar;
            boandrdata{nidx, 5} = 0;
            nidx = nidx + 1;
    % end
end