% %% Calculates which carbons are bonded by finding carbons that are <160 pm apart
% bonds = [];
% for i = 1:ampac.natom - 1
%     for j = i+1:ampac.natom
%         if (strcmpi(ampac.element{i},'c') && strcmpi(ampac.element{j},'c'))
%             bl = sum((ampac.r(:,i) - ampac.r(:,j)) .^ 2) ^ 0.5;
%             if (bl < 1.6)
%                 bonds(1,end+1) = i;
%                 bonds(2,end) = j;
%                 bonds(3,end) = bl;
%             end
%         end
%     end
% end
% %% Draws structure to verify all bonds are accounted for
% 
% figure(7);
% hold on;
% axis equal;
% axis([min(ampac.r(1,:)) max(ampac.r(1,:)) min(ampac.r(2,:)) max(ampac.r(2,:)) min(ampac.r(3,:))-1 max(ampac.r(3,:))+1]);
% for i = 1:ampac.natom
%     if (strcmpi(ampac.element{i},'c'))
%         scatter3(ampac.r(1,i),ampac.r(2,i),ampac.r(3,i),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b');
%     end
% end
% 
% for i = 1:size(bonds,2)
%     dat = [ampac.r(:,bonds(1,i)) ampac.r(:,bonds(2,i))];
%     line(dat(1,:),dat(2,:),dat(3,:),'color','r')
% end

%% Read in optimized ground state z-matrix
function [rgs, rex, deltar, atom_nums_in_bonds, numruns] = OptExcStateStructure(ampac_pathonly, ampac_nameonly, indo_pathonly, indo_nameonly, varargin)
    
% ampac_pathonly = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\DMG-8mer\';
% ampac_nameonly = '8-merPPVampac';
    
    % EXECUTABLE PATHWAYS, NO QUOTES
    AmpacEXE = 'C:\Program Files\Semichem, Inc\Ampac-9.2\ampac.exe';
    OpenBabelEXE = 'C:\Program Files (x86)\OpenBabel-2.3.1\babel.exe';
    
    
    ampac_filepath = [ampac_pathonly, ampac_nameonly];
    indo_filepath = [indo_pathonly, indo_nameonly];
    copyfile([ampac_filepath, '.out'],[ampac_filepath, '-new.out']);
    copyfile([indo_filepath, '-dm.bin'],[indo_filepath, '-new-dm.bin']);
    
    if (length(varargin) > 1)
        indo = varargin{2};
        indo_nameonly = [indo_nameonly, '-new'];
        indo_filepath = [indo_pathonly, indo_nameonly];
    else
        copyfile([indo_filepath, '.ido'],[indo_filepath, '-new.ido']);
        indo_nameonly = [indo_nameonly, '-new'];
        indo_filepath = [indo_pathonly, indo_nameonly];
        indo = Indo.LoadExistingData([indo_filepath,'.ido'],[],[],[]);
    end
       
    ampac_nameonly = [ampac_nameonly, '-new'];
    ampac_filepath = [ampac_pathonly, ampac_nameonly];
    
    
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

    while (~getOut)
        numruns = numruns + 1;
        % newzmatrix = zmatrix;
        
        %% Read in ground state density matrix file

        fid = fopen([indo_filepath,'-dm.bin']);

        matdim = fread(fid,1,'int');

        dm = fread(fid,[matdim,matdim],'double');

        fclose(fid);

        %% Calculate bond order of C-C bonds

        bo = size(CCbonds,2);
        bobyatom = zeros(natoms);

        for i = 1:size(CCbonds,1)
            blmat = dm(indo.aorbAtom == CCbonds(i,1), indo.aorbAtom == CCbonds(i,2)) .^ 2;
            bo(i) = sum(blmat(:));
            bobyatom(CCbonds(i,1), CCbonds(i,2)) = bo(i);
            bobyatom(CCbonds(i,2), CCbonds(i,1)) = bo(i);
        end
        
        if (isempty(gsbo))
            gsbo = bo;
            gsbobyatom = bobyatom;
        end

        % Columns 1-3: Atom 1, Atom 2, Calculated GS Bond Order
        atomsandbo = [CCbonds(:,1) CCbonds(:,2) bo'];

        %% Calculate new bond lengths and ring angles

        newbls = [atomsandbo(:,1:2) (1.54-0.3*log(bo'))];
        newblsbyatom = 1.54-0.3*log(bobyatom);  % Pauling bond order

        newringang = cell(length(rings),1);

        for i = 1:length(newringang)
            blset = zeros(length(rings{i}),1);
            for j = 1:length(rings{i})-1
                blset(j) = newblsbyatom(rings{i}(j), rings{i}(j+1));
            end
            blset(end) = newblsbyatom(rings{i}(end), rings{i}(1));

            ropt = RingOptimizer(blset);
            [in,~,~,~] = ropt.run('deg');
            newringang{i} = in;
        end

        internalringangles = cell(sum(cellfun(@(x)length(x), rings)), 2);

        i = 1;

        for j = 1:length(rings)
            ratoms = [rings{j} rings{j}(1:2)];
            for k = 1:length(ratoms)-2
                internalringangles{i,1} = ratoms(k:(k+2));
                internalringangles{i,2} = newringang{j}(k);
                i = i + 1;
            end
        end


        %% Assign new bond lengths and angles

        allbas = zeros(natoms,natoms,natoms);
        for i = 3:natoms
            allbas(i, str2double(zmatrix{i,6}), str2double(zmatrix{i,7})) = str2double(zmatrix{i,4});
            allbas(str2double(zmatrix{i,7}), str2double(zmatrix{i,6}), i) = str2double(zmatrix{i,4});
        end

        % Write new bond lengths to new z-matrix
        for i = 2:natoms
            if (strcmpi(zmatrix{i,2},'c'))
                newzmatrix{i,3} = num2str(newblsbyatom(i,str2double(newzmatrix{i,6})));
            end
        end

        % Get internal ring angles which match the atoms that connect to
        % rings
        results = cell(length(ringinc),1);
        ira = reshape(cell2mat({internalringangles{:,1}}),3,[])';
        for i = 1:length(ringinc)
            centeratom = ringinc(i,2);

        %     for j = 1:length(allbanums)
        %         if (~all(allbanums(j,:) == ringinc(i,:)))
        %             if (allbanums(j,2) == centeratom)
        %                 results{i}(end+1) = j;
        %             end
        %         end
        %     end

            for j = 1:length(ira)
                if (~all(ira(j,:) == ringinc(i,:)))
                    if (ira(j,2) == centeratom)
                        results{i}(end+1) = j;
                    end
                end
            end
        end

        % Split the internal angle in half and make that the bond angle. If
        % the center atom is part of two rings, subtract the two internal
        % angles from 360 to get the bond angle necessary
        for i = 1:length(results)
            if (length(results{i}) == 1)
                nang = internalringangles{results{i}(1),2};
                nang = 180 - nang / 2;
            else
                nang = internalringangles{results{i}(1),2};
                nang = 360 - nang - internalringangles{results{i}(2),2};
            end
            allbas(ringinc(i,1),ringinc(i,2),ringinc(i,3)) = nang;
            allbas(ringinc(i,3),ringinc(i,2),ringinc(i,1)) = nang;
        end

        % Assign only internal ring angles their new values
        for i = 1:length(ringonly)
            loc = find(cellfun(@(x)isequal(x,ringonly(i,:)), internalringangles(:,1)));
            if (isempty(loc))
                loc = find(cellfun(@(x)isequal(x,fliplr(ringonly(i,:))), internalringangles(:,1)));
            end

            allbas(ringonly(i,1),ringonly(i,2),ringonly(i,3)) = internalringangles{loc,2};
            allbas(ringonly(i,3),ringonly(i,2),ringonly(i,1)) = internalringangles{loc,2};
        end

        % Assign new bond angles into the new z-matrix
        for i = 3:natoms
            nang = allbas(i, str2double(newzmatrix{i,6}), str2double(newzmatrix{i,7}));
            newzmatrix{i,4} = num2str(nang);
        end
        
        %% Get current bond lengths and check them
        
        temp = cellfun(@(x)str2double(x),{newzmatrix{2:end,3}});
        diff = max(abs(temp - rex));
        rex = temp;
        
        disp(['OptExcStStruct: run ', num2str(numruns), ', diff ', num2str(diff)]);
        
        if (diff < 1e-4)
            getOut = true;
            % deltar = rex - rgs;
            deltar = zeros(natoms-1,1);
            atom_nums_in_bonds = zeros(natoms-1,2);
            
            for (i = 1:natoms-1)
                trgs = gsbobyatom(str2double(zmatrix{i+1,1}),str2double(zmatrix{i+1,6}));
                trex = bobyatom(str2double(zmatrix{i+1,1}),str2double(zmatrix{i+1,6}));
                atom_nums_in_bonds(i,1:2) = [str2double(zmatrix{i+1,1}), str2double(zmatrix{i+1,6})];
                deltar(i) = 0.3*log(trgs / trex);
            end
            
            rgs = rgs';
            rex = rex';
        end
        
        if (~getOut)
        %% Write out new z-matrix as Ampac DAT for single-point calculation


            fid = fopen([ampac_filepath, '.dat'],'w');

            fprintf(fid, '%s\r\n%s\r\n%s\r\n', 'AM1 rhf singlet 1scf t=auto geom=ok', 'Title', 'Comment');
            % fprintf(fid, '%s\r\n\r\n\r\n', 'INSERT KEYWORDS HERE');
            % fprintf(fid, '%s\r\n\r\n', 'Title');
            % '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n'
            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{1,2}, '0.000000', '1', '0.000000', '1', '0.000000', '1', '0','0','0');

            if (natoms > 1)
                fprintf(fid,'%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{2,2}, newzmatrix{2,3}, '1',...
                    '0.000000', '1', '0.000000', '1', newzmatrix{2,6},'0','0');
            end
            if (natoms > 2)
                fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{3,2}, newzmatrix{3,3}, '1',...
                    newzmatrix{3,4}, '1', '0.000000', '1', newzmatrix{3,6}, newzmatrix{3,7},'0');
            end
            if (natoms > 3)
                for i = 4:natoms
                    fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s\r\n', newzmatrix{i,2}, newzmatrix{i,3}, '1',...
                        newzmatrix{i,4}, '1', newzmatrix{i,5}, '1', newzmatrix{i,6}, newzmatrix{i,7}, newzmatrix{i,8});
                end
            end

            fprintf(fid, '%s %s  %s  %s  %s  %s  %s  %s  %s  %s', '0', '0.000000', '0', '0.000000', '0', '0.000000', '0', '0','0','0');

            fclose(fid);


            %% Run single point calculation to generate OUT file and then send it to INDO

            [~,~] = system(['"', AmpacEXE, '" "', ampac_filepath, '.dat"']);
            
            if (~isempty(varargin) && ~isempty(varargin{1}))
                efield = varargin{1};
            else
                efield = [0 0 0];
            end
            
            res = [];
            res.charge = 0;
            res.norbs = 500;
            res.nstates = 25;
            res.field = efield;
            res.initial_shiftc = 80.0;
            res.initial_shift_step = 1.0;
            res.min_shift_step = 0.1;
            res.max_shift_step = 10.0;
            res.initial_second_shift_step = 0.5;
            res.min_second_shift_step = 0.01;
            res.max_second_shift_step = 0.5;
            res.initial_eeint = 1.0;
            res.initial_eestep = 0.0;
            res.min_eestep = 0.0;
            res.max_eestep = 0.0;
            res.initial_conv = 1e-3;
            res.min_conv = 1e-10;
            res.max_inner_iter = 1500;
            res.max_iter = 150000;
            res.dm_guess = [indo_filepath, '-dm.bin'];
            res.try_default_first = true;
            res.output_dm = true;
            res.pot_file = [];

            indo = Indo(res, ampac_pathonly, ampac_nameonly);
            if (~strcmp(ampac_filepath, indo_filepath))
                movefile([ampac_pathonly, ampac_nameonly, '.ido'], [indo_pathonly, indo_nameonly, '.ido']);
                movefile([ampac_pathonly, ampac_nameonly, '-dm.bin'], [indo_pathonly, indo_nameonly, '-dm.bin']);
            end
            % delete([ampac_pathonly, ampac_nameonly, '.ipf']);
        end
    end
end


    