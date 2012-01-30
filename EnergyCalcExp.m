classdef EnergyCalcExp < handle
%--------------------------------------------------------------------------------------------------%
%   Class Name: EnergyCalcExp
%   Description: Runs AMPAC and INDO calculations from template files
%   Written by: Christian Legaspi, Carnegie Mellon University
%   Date: June 10, 2011
%   Comments: If you don't know how to operate this...read the function headers
%       and comments on properties and hopefully that will help you out. Related classes are
%       ECEParams (written by CML), ECEDataStuct (written by CML), Indo (written by Dave Yaron),
%       DataHash (written by Jan Simon), waitbar (written by Andrew, modified by CML).
%       Inspiration was taken from pausebutton (written by Murphy O'Brien) but no code was directly
%       used from the original work.
%   Revisions:
%       7/6/11: Completion of hash file name implementation. Replace the command window updates
%           of progress and with GUI progress bar with option to pause, save and quit on the fly.
%           Periodic save still retained in case of power outage or crash or what not. (CML)
%--------------------------------------------------------------------------------------------------%
    
    properties (SetAccess = private, Transient)
        ampac_lib;              % containers.Map. Library of Ampac hashes and their parameters
        indo_lib;               % containers.Map. Library of INDO hashes and their parameters
    end
    
    properties (SetAccess = private)
        zmatrix_template;       % char[]: stores the template from the DAT input file
        mol_charge;             % int: molecular charge as dictated by the input file (assumed 0 if not present in file)
        progress;               % unsigned int[1,n]: 1xn matrix indicating the indices of calculation progress
                                % where each n represents PHI(n).
                                % progress = finishpt when complete.
        finishpt;               % unsigned int[1,n]: 1xn matrix indicating the size of each PHI(n) parameter matrix
        bigO;                   % unsigned int: total number of configurations to be calculated
        data;                   % ECEDataStruct[]: n-dimensional matrix containing the calculated data stored in ECEDataStruct objects
        fail_bucket;            % struct[1,m]: 1xm struct matrix which stores information about failed calculations. Struct contains:
                                % phiID: 1xn matrix with absolute angle values where error occurred.
                                % type: char array containing either 'indo' or 'ampac' indicating the type of calc failure.
                                %       If the type is 'ampac' and INDO calcs were to be run, there will be no indo calcs
                                %       available.
                                % msg: char array description of the error that occurred
        fail_iterator;          % unsigned int: keeps track of the next index available in fail_bucket.
                                % length(fail_bucket) = m = fail_iterator - 1
        
        ampacexe = '"c:\Program Files\Semichem, Inc\Ampac-9.2\ampac.exe"';    % Path to the AMPAC executable, including " marks
    end
    
    properties (SetAccess = public)
        params;                 % ECEParams: class object containing information about the experiment to be run.
        indodatapath;           % char[]: path where INDO calc files will be stored/retrieved (libary directory)
        ampacdatapath;          % char[]: path where AMPAC files will be stored/retrieved (library/directory)
        MATFilename;            % char[]: path to loaction where MAT file will be saved for this experiment
        molecule;               % char[]: path to template file for the experiment
        axis_params;            % unsigned int[1,2|3]: contains atom numbers to create E field vector. If empty, a custom
                                % vector must be specified in params. If length = 2, vector will be created parallel to the
                                % vector between the two atom numbers. If length = 3, vector will be created normal to the
                                % plane formed by the three atom numbers. Order of these arguments does matter.
        description;            % char[]. A short description about the experiment for bookkeeping purposes
    end
    methods (Access = private)  % Functions only accessible by functions of the class object
%--------------------------------------------------------------------------------------------------%
% Function:     run_Ampac
% Description:  Runs the Ampac software package to produce optimized geometry and ground state energy
%               information. If a filename exists which matches the experiment parameters, data is
%               read in from these files.
% Input Args:
%       obj:    EnergyCalcExp calling object. For the sake of clarity, all functions contained
%               within this class will have 'obj' as their EnergyCalcExp calling object if they
%               aren't static functions. This parameter does not need to be included in the argument
%               list but is assigned when the function is called. 'obj' is passed by reference and
%               functions similar to the 'this' pointer in C++. All other input args are passed by
%               value, unless they are specified in the output args as well. This explanation
%               will be omitted in other functions listed below.
%       phi:    int/double[1,n]. Contains the angle values for each PHI(n) angle being varied on
%               this calculation.
%       suppress_output: bool. TRUE suppresses output to the command window.
% Output:  
%       ampac:  struct. Contains all fields returned by parseAmpac in addition to the job name
%               (char[]: jobname), calculation success (bool: ampac_succeed) and any error message
%               (char[]: ampac_err_msg)
%--------------------------------------------------------------------------------------------------%
        function ampac = run_Ampac(obj, phi, suppress_output)
            key = obj.keygen('ampac', phi);
            hash = DataHash(key);
            data_exist = false;
            ampac.hash = [];
            
            while (isempty(ampac.hash))
                if (~isempty(obj.ampac_lib))
                    if (obj.ampac_lib.isKey(hash))
                        if (~isequal(obj.ampac_lib(hash), key))
                            hash = DataHash(hash);
                        else
                            ampac.hash = hash;
                            data_exist = true;
                        end
                    else
                        ampac.hash = hash;
                    end
                else
                    if (exist([obj.ampacdatapath, hash, '.mat'], 'file'))
                        S = load([obj.ampacdatapath, hash, '.mat'], 'key');
                        if (~isequal(S.key, key))
                            hash = DataHash(hash);
                        else
                            ampac.hash = hash;
                            data_exist = true;
                        end
                    else
                        ampac.hash = hash;
                    end
                end
            end
            
            if (~data_exist)
                data_out = struct();
                data_out.hash = hash;
                data_out.key = key;

                filetext = strrep(obj.zmatrix_template, 'METHOD', obj.params.method);

                for i = 1:length(phi)
                    filetext = strrep(filetext, ['PHI', num2str(i)], num2str(phi(i)));
                end 

                dat_file = [obj.ampacdatapath, hash, '.dat'];
                fid1 = fopen(dat_file,'w');
                fwrite(fid1, filetext, 'char');
                fclose(fid1);
                
                data_out.dat = filetext;

                if (~suppress_output)
                    disp(['Running AMPAC calculation on ',jobname]);
                end
                [~, result] = system([obj.ampacexe,' ', dat_file]); 
                
                if (~exist([obj.ampacdatapath, hash, '.arc'], 'file') && ...
                        exist([obj.ampacdatapath, hash, '.out'], 'file'))
                    
                    outfiletext = fileread([obj.ampacdatapath, hash, '.out']);
                    data_out.out = outfiletext;
                    data_out.vis = [];
                    data_out.arc = []; %#ok<STRNU>
                    ampac.ampac_succeed = false;
                    if (~isempty(regexpi(outfiletext, '.*TO CONTINUE CALCULATION SPECIFY "GEO-OK".*', 'start')))
                        ampac.ampac_err_msg = 'Geometry failure.';
                    else
                        ampac.ampac_err_msg = 'An unknown error occurred. See *.out file.';
                    end

                    if (~suppress_output)
                        disp(['AMPAC failed. ',ampac.ampac_err_msg]);
                    end
                else
                    data_out.out = fileread([obj.ampacdatapath, hash, '.out']);
                    data_out.vis = fileread([obj.ampacdatapath, hash, '.vis']);
                    data_out.arc = fileread([obj.ampacdatapath, hash, '.arc']);
                    
                    if (~suppress_output)
                        disp('AMPAC finished');
                    end
                    % Parse AMPAC data
                    if (~suppress_output)
                        disp('Parsing AMPAC data');
                    end
                    ampac = parseAmpacFromText(data_out.arc, data_out.out);
                    ampac.ampac_succeed = true;
                    ampac.hash = hash;
                    
                    if (~suppress_output)
                        disp('Done parsing AMPAC');
                    end
                end
                
                if (~isempty(obj.ampac_lib))
                    obj.ampac_lib(hash) = key;
                    lib = obj.ampac_lib; 
                    save([obj.ampacdatapath, 'catalog.mat'], 'lib', '-append');
                end
                
                save([obj.ampacdatapath, hash, '.mat'], '-struct', 'data_out');
            else
                S = load([obj.ampacdatapath, hash, '.mat'], 'out', 'arc');
                if (~isempty(S.arc))
                    ampac = parseAmpacFromText(S.arc, S.out);
                    ampac.hash = hash;
                    ampac.ampac_succeed = true;
                    
                    fid = fopen([obj.ampacdatapath, hash, '.out'], 'w');
                    fwrite(fid, S.out, 'char');
                    fclose(fid);
                else
                    ampac.ampac_succeed = false;
                    if (~isempty(regexpi(out, '.*TO CONTINUE CALCULATION SPECIFY "GEO-OK".*', 'start')))
                        ampac.ampac_err_msg = 'Geometry failure.';
                    else
                        ampac.ampac_err_msg = 'An unknown error occurred. See *.out file.';
                    end

                    if (~suppress_output)
                        disp(['AMPAC failed. ',ampac.ampac_err_msg]);
                    end
                end
            end
        end
    end
    
    methods (Static)    % Functions that can be used without creating an instance of the class
%--------------------------------------------------------------------------------------------------%
% Function:     get_subscripts
% Description:  Takes a linear index and the size matrix of the array and returns the subscripts.
%               Same as the Matlab function ind2sub() but the output type is different.
% Input Args:
%       size_of_matrix: unsigned int[1,n]. Contains the size of the array the linear index is from
%       lin_idx: unsigned int. Linear index to get subscripts for.
% Output:  
%       pos_mat: unsigned int[1,n]. Matrix containing the subscripts represented by the linear index
%--------------------------------------------------------------------------------------------------%
        function pos_mat = get_subscripts(size_of_matrix, lin_idx)
            pos_mat = zeros(1,length(size_of_matrix));
            
            for i = 1:length(size_of_matrix)
                pos_mat(i) = mod(lin_idx - 1, size_of_matrix(i)) + 1;
                lin_idx = ((lin_idx - pos_mat(i)) / size_of_matrix(i)) + 1;
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     get_linear_index
% Description:  Takes a matrix of subscripts and the size matrix of the array and returns the
%               linear index represented by the subscripts.
% Input Args:
%       size_of_matrix: unsigned int[1,n]. Contains the size of the array the linear index is from
%       pos_mat: unsigned int[1,n]. Matrix containing the subscripts of the location.
% Output:  
%       idx: unsigned int. Linear index of the subscript set.
%--------------------------------------------------------------------------------------------------%
        function idx = get_linear_index(size_of_matrix, pos_mat, reverse)
            if (length(pos_mat) == 1)
                idx = pos_mat(1);
                return;
            end
            
            if (nargin == 3 && reverse)     % Calculates linear index in reverse
                idx = 0;

                for i = 1:(length(size_of_matrix) - 1)
                    innerloop = pos_mat(i) - 1;
                    for j = (i+1) : length(size_of_matrix)
                        innerloop = innerloop * size_of_matrix(j);
                    end
                    idx = idx + innerloop;
                end 

                idx = idx + pos_mat(end);
            else
                % Calculates the actual linear index MATLAB understands
                idx = 0;

                for i = (length(size_of_matrix)):-1:2
                    innerloop = pos_mat(i) - 1;
                    for j = (i-1) : -1 : 1
                        innerloop = innerloop * size_of_matrix(j);
                    end
                    idx = idx + innerloop;
                end 

                idx = idx + pos_mat(1);
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     ReadMATFile
% Description:  Loads experiment data saved by this class to a MAT file.
% Input Args:
%       matfile: char[]. Path of the MAT file.
% Output:  
%       obj: EnergyCalcExp object. Contents of the file.
%--------------------------------------------------------------------------------------------------%
        function obj = ReadMATFile(matfile, varargin)
            if (nargin > 1)
                loadcats = varargin{1};
            else
                loadcats = true;
            end
            
            if (~exist(matfile,'file'))
                throw(MException('EnergyCalcExp:ReadMATFile:File_DNE',...
                    'EnergyCalcExp: MAT File does not exist or cannot be accessed'));
            end
            
            mydata = load(matfile);
            obj = mydata.obj;
            
            if (loadcats)
                obj.load_catalogs();
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     generate_catalog
% Description:  Generates a library catalog from the files in a directory
% Input Args:
%       path: char[]. Path to the directory to be analyzed for cataloging
%       varargin: bool. TRUE means to silence queries about libraries >10,000 files. Optional.
% Output:  
%       successful: bool. TRUE means it was successful
%--------------------------------------------------------------------------------------------------%
        function successful = generate_catalog(path, varargin)
            if (isempty(varargin))
                silent = false;
            else
                silent = varargin{1};
            end
            
            if (path(end) ~= '\')
                path(end+1) = '\';
            end
            
            files = dir([path, '*.mat']);
            fn = cell(1,numel(files)); 
            fn = {files(1:end).name};
            fn = cellfun(@(x)regexpi(x, '^([A-Fa-f0-9]){32}\.mat', 'match'), fn, 'UniformOutput', false);
            fn = {fn{~(cellfun(@(x)isempty(x), fn))}}; 
            fn = cat(2, fn{1:end});
            
            if (length(fn) > 10000 && ~silent)
                res = input(['This library contains >10,000 eligible files.\nIt is ',...
                    'not recommended to create libraries larger than 10,000 files due to ',...
                    'memory usage concerns and speed issues.\nAre you sure you want to ',...
                    'proceed?\n(y/n)'], 's');
                if (~isempty(res) && strcmpi(res(1), 'n'))
                    successful = 0;
                    return;
                end
            end              
            
            lib = containers.Map('KeyType','char','ValueType','any');
            
            if (~isempty(fn))
                for i = 1:length(fn)
                    S = load([path, fn{i}], 'hash', 'key');
                    if (length(fieldnames(S)) == 2)
                        lib(S.hash) = S.key;
                    end
                end
            end
            
            save([path, 'catalog.mat'], 'lib');
            
            successful = 1;       
        end
            
    end
    
    methods             % Public functions that can be called with an instance of a class
%--------------------------------------------------------------------------------------------------%
% Function:     Constructor
% Description:  Builds an instance of the EnergyCalcExp object. Can be built without any arguments.
% Input Args:
%       parameters: ECEParams object. Contains the experiment parameters.
%       mol: char[]. Full path to the molecule template file.
%       axis_params: unsigned int[1,2|3]. contains atom numbers to create E field vector. If empty, a custom
%               vector must be specified in params. If length = 2, vector will be created parallel to the
%               vector between the two atom numbers. If length = 3, vector will be created normal to the
%               plane formed by the three atom numbers.
%       MAT: char[]. Full path to the save location and file name for the experiment file.
%       indo: char[]. Path to the INDO library
%       ampac: char[]. Path where Ampac library
%       load_libs: bool. TRUE attempts to load the library catalogs for Ampac and INDO. If they do not exist,
%               new libraries are created.
%       varargin: char[]. Description about the experiment. Optional.
% Output:  
%       obj: EnergyCalcExp object. The constructed object.
%--------------------------------------------------------------------------------------------------%
        function obj = EnergyCalcExp(parameters, mol, axis_params, MAT, indo, ampac, load_libs, varargin)
            if (nargin > 4)
                if (~exist(mol,'file'))
                    throw(MException('EnergyCalcExp:Constructor:Bad_Template_File',...
                        'EnergyCalcExp: Bad input template file'));
                end
                
                obj.params = parameters;
                obj.MATFilename = MAT;
                obj.molecule = mol;
                obj.axis_params = axis_params;
                
                if (~isempty(varargin))
                    obj.description = varargin{1};
                end

                if (indo(end) ~= '\')
                    indo = [indo, '\'];
                end
                obj.indodatapath = indo;
                
                if (load_libs)
                    if (~exist([obj.indodatapath, 'catalog.mat'], 'file'))
                        lib = containers.Map('KeyType','char','ValueType','any');
                        obj.indo_lib = lib;
                        if (~exist(obj.indodatapath, 'dir'))
                            mkdir(obj.indodatapath);
                        end
                        save([obj.indodatapath, 'catalog.mat'], 'lib');
                    else
                        S = load([obj.indodatapath, 'catalog.mat'], 'lib');
                        obj.indo_lib = S.lib;
                    end
                else
                    obj.indo_lib = [];
                end

                if (ampac(end) ~= '\')
                    ampac = [ampac, '\'];
                end
                obj.ampacdatapath = ampac;
                
                if (load_libs)
                    if (~exist([obj.ampacdatapath, 'catalog.mat'], 'file'))
                        lib = containers.Map('KeyType','char','ValueType','any');
                        obj.ampac_lib = lib;
                        if (~exist(obj.ampacdatapath, 'dir'))
                            mkdir(obj.ampacdatapath);
                        end
                        save([obj.ampacdatapath, 'catalog.mat'], 'lib');
                    else
                        S = load([obj.ampacdatapath, 'catalog.mat'], 'lib');
                        obj.ampac_lib = S.lib;
                    end
                else
                    obj.ampac_lib = [];
                end

                obj.zmatrix_template = fileread(obj.molecule);

                idx1 = regexpi(obj.zmatrix_template, 'charge\s*=\s*', 'once', 'end');

                if (isempty(idx1))
                    obj.mol_charge = 0;
                else
                    idx2 = regexpi(obj.zmatrix_template(idx1 + 1 : end), '\s', 'once');
                    obj.mol_charge = str2num(obj.zmatrix_template(idx1 + 1 : idx1 + idx2 - 1)); %#ok<ST2NM>
                end
                
                obj.fail_iterator = 1;
                
                if (~load_libs)
                    warning('off','backtrace');
                    warning('EnergyCalcExp:LibraryNotLoaded',['The library catalogs were not loaded. This is recommended for libraries',...
                        ' with >10,000 entries due to memory concerns. However this means seek time will ',...
                        'increase for previous entries.']); 
                    warning('on','backtrace');
                end
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     run
% Description:  Runs the calculations specified by the parameters in the class object.
% Input Args:
%       obj: EnergyCalcExp calling object.
%       varargin: Specify the following optional properties as strings (case-indifferent):
%           'quiet': Prevents display of excess text from being output to the command window and
%                   shows only the progress indicator (unless used in conjunction with 'NoProgress').
%           'SaveAfter': Next argument must be an integer representing the number of calculations to run
%                   before saving the progress to the MAT file.
%           'overwrite': Has no effect if the calculation has never been run. If the calculation
%                   has been loaded from a MAT file, this argument starts the calculation from the beginning and
%                   reads in the data files anew.
%           'NoProgress': Prevents the loading of the progress indicator window. Useful if the calculations are
%                   being run in an outside loop and the run() function only does a few calculations. However,
%                   this also means that the pause button will be unavailable.
% Output:  
%       successful: bool. TRUE means that everything went OK. FALSE means there were calculation errors.
%--------------------------------------------------------------------------------------------------%
        function successful = run(obj, varargin)
            quiet_mode = any(cellfun(@(x)isequal(lower(x),'quiet'), varargin)); %#ok<*EFIND>
            force_overwrite = any(cellfun(@(x)isequal(lower(x),'overwrite'), varargin));
            no_prog_window = any(cellfun(@(x)isequal(lower(x),'noprogress'), varargin));
            
            save_count = 0;
            if (any(cellfun(@(x)isequal(lower(x),'saveafter'), varargin)))
                try
                    save_val = varargin{find(cellfun(@(x)isequal(lower(x),'saveafter'), varargin)) + 1};
                catch exception
                    if (strcmp(exception.identifier, 'MATLAB:badsubscript'))
                        throw(MException('EnergyCalcExp:run:BadSaveAfterArg',...
                            'You must specify an integer value following the argument ''SaveAfter'''));
                    else
                        rethrow(exception);
                    end
                end
                if (isnumeric(save_val) && ~isinf(save_val) && ~isnan(save_val))
                    save_count = floor(real(save_val));
                else
                    throw(MException('EnergyCalcExp:run:BadSaveAfterArg',...
                        'You must specify an integer value following the argument ''SaveAfter'''));
                end
            end
                
            if (isequal(obj.progress, obj.finishpt) && ~isempty(obj.progress) && ~isempty(obj.finishpt) &&...
                    ~force_overwrite)
                warning('off','backtrace');
                warning('EnergyCalcExp:CalcsAlreadyComplete',['Calculations have already been marked as completed. ',...
                    'To force reloading of the data files, please specify TRUE ',...
                    'for the third argument in run(). To force recalculation, ',...
                    'please manually remove the files or change the library directory.']); %#ok<*WNTAG>
                warning('on','backtrace');
                successful = true;
                return;
            end
            
            save_iterator = 0;  % Counter to determine when to save checkpoint files
            
            if (isempty(obj.params.phi))
                obj.run_opt_only(quiet_mode);
                return;
            end
            
            if (isempty(obj.progress) || force_overwrite)
                obj.progress = ones(1,length(obj.params.phi));                
                obj.finishpt = [];
                
                for i = 1:length(obj.params.phi)
                    obj.finishpt(i) = length(obj.params.phi{i});
                end
                
                obj.bigO = prod(obj.finishpt);
                obj.data = repmat(ECEDataStruct(), obj.finishpt);
                obj.fail_bucket = repmat(struct('phiID', [], 'type', [], 'msg', []), [numel(obj.data) 1]);
            end
            
            indoparams = [];
            
            if (obj.params.doIndo)
                warning('off','MATLAB:structOnObject');
                indoparams = struct(obj.params);
                warning('on','MATLAB:structOnObject');
                indoparams.charge = obj.mol_charge;
                indoparams.norbs = obj.params.norbs_basis_set;
                indoparams.nstates = obj.params.nexcstates;
            end
            
            getOut = false;
            % mytic = tic();
            
            if (~no_prog_window)
                waithdl = waitbar(0, 'Calculations running...');
            end
            
            while (~getOut)
                phi_params = zeros(1,length(obj.params.phi));
                for i = 1:length(obj.params.phi)
                    phi_params(i) = obj.params.phi{i}(obj.progress(i));
                end
                
                ampac = obj.run_Ampac(phi_params, quiet_mode);
                
                if (ampac.ampac_succeed)
                    if (obj.params.doIndo)
                        if (~no_prog_window)
                            waithdldata = get(waithdl,'userdata');
                            if (waithdldata.bail == 0)
                                waitbar((prod(obj.progress) - 1) / obj.bigO, waithdl);
                            else
                                if (waithdldata.bail == 2)
                                    obj.save_data();
                                end

                                if (exist([obj.ampacdatapath, ampac.hash, '.vis'], 'file'))
                                    delete([obj.ampacdatapath, ampac.hash, '.vis']);
                                end
                                if (exist([obj.ampacdatapath, ampac.hash, '.out'], 'file'))
                                    delete([obj.ampacdatapath, ampac.hash, '.out']);
                                end
                                if (exist([obj.ampacdatapath, ampac.hash, '.arc'], 'file'))
                                    delete([obj.ampacdatapath, ampac.hash, '.arc']);
                                end
                                if (exist([obj.ampacdatapath, ampac.hash, '.dat'], 'file'))
                                    delete([obj.ampacdatapath, ampac.hash, '.dat']);
                                end

                                disp('Calculation ended by user.');
                                close(waithdl,'force');
                                successful = 0;
                                return;
                            end
                        end
                        
                        if (isempty(obj.params.Efield_vector) && ...
                            ~(obj.params.Efield_mag == 0 || isempty(obj.params.Efield_mag)) && ...
                            (isempty(obj.axis_params)))

                            throw(MException('EnergyCalcExp:run:No_Efield_Vector',...
                                'EnergyCalcExp: No data provided for E field vector'));
                        elseif (obj.params.Efield_mag == 0 || isempty(obj.params.Efield_mag))
                            obj.params.setCustomEFVector([0 0 0]);
                        end

                        if (isempty(obj.params.Efield_type) || (obj.params.Efield_type ~= 'c'))
                            obj.params.setEFVectorFromAtomNums(ampac, obj.axis_params);
                        end

                        indoparams.field = obj.params.Efield_vector * obj.params.Efield_mag;
                        
                        key = obj.keygen('indo', phi_params);
                        hash = DataHash(key);
                        data_exist = false;
                        indo_out = struct();
                        indo_out.key = key;
                        indo_out.hash = [];

                        while (isempty(indo_out.hash))
                            if (~isempty(obj.indo_lib))
                                if (obj.indo_lib.isKey(hash))
                                    if (~isequal(obj.indo_lib(hash), key))
                                        hash = DataHash(hash);
                                    else
                                        indo_out.hash = hash;
                                        data_exist = true;
                                    end
                                else
                                    indo_out.hash = hash;
                                end
                            else
                                if (exist([obj.indodatapath, hash, '.mat'], 'file'))
                                    S = load([obj.indodatapath, hash, '.mat'], 'key');
                                    if (~isequal(S.key, key))
                                        hash = DataHash(hash);
                                    else
                                        indo_out.hash = hash;
                                        data_exist = true;
                                    end
                                else
                                    indo_out.hash = hash;
                                end
                            end
                        end

                        if (~data_exist)
                            if (~quiet_mode)
                                disp('Running INDO calculation');
                            end
%                             if (obj.ampacdatapath(end-1) == ':')
%                                 indo = Indo(indoparams, obj.ampacdatapath(1:end), ampac.hash);
%                             else
                                indo = Indo(indoparams, obj.ampacdatapath(1:end-1), ampac.hash);
%                             end
                            if (indoparams.output_dm)
                                movefile([obj.ampacdatapath, ampac.hash, '-dm.bin'], ...
                                    [obj.indodatapath, hash, '-dm.bin']);
                            end
                                
                            indo_out.indo = indo;
                            
                            if (~isempty(obj.indo_lib))
                                obj.indo_lib(hash) = key;
                                lib = obj.indo_lib; 
                                save([obj.indodatapath, 'catalog.mat'], 'lib', '-append');
                            end
                            
                            save([obj.indodatapath, hash, '.mat'], '-struct', 'indo_out');
                        else
                            S = load([obj.indodatapath, hash, '.mat'], 'indo');
                            indo = S.indo;
                        end
                        
                        if (~indo.indo_succeed)
                            if (~quiet_mode)
                                disp(['INDO Failed. ', indo.indo_err_msg]);
                            end
                            obj.fail_bucket(obj.fail_iterator) = ...
                                struct('phiID', phi_params, 'type', 'indo', 'msg', indo.indo_err_msg);
                            obj.fail_iterator = obj.fail_iterator + 1;
                        elseif (~quiet_mode)
                            disp(['INDO finished. Done parsing INDO.', indo.indo_err_msg]);
                        end
                        
                        obj.data(obj.get_linear_index(obj.finishpt, obj.progress)) = ...
                            ECEDataStruct(ampac, indo, phi_params, obj.ampacdatapath, ampac.hash,...
                            obj.indodatapath, indo_out.hash);

                    else
                        obj.data(obj.get_linear_index(obj.finishpt, obj.progress)) = ...
                            ECEDataStruct(ampac, [], phi_params, obj.ampacdatapath, ampac.hash,...
                            [], []);
                    end
                else
                    obj.fail_bucket(obj.fail_iterator) = ...
                        struct('phiID', phi_params, 'type', 'ampac', 'msg', ampac.ampac_err_msg);
                    obj.fail_iterator = obj.fail_iterator + 1;
                    obj.data(obj.get_linear_index(obj.finishpt, obj.progress)) = ...
                        ECEDataStruct(ampac, [], phi_params, obj.ampacdatapath, ampac.hash, [], []);
                end
                
                if (exist([obj.ampacdatapath, ampac.hash, '.vis'], 'file'))
                    delete([obj.ampacdatapath, ampac.hash, '.vis']);
                end
                if (exist([obj.ampacdatapath, ampac.hash, '.out'], 'file'))
                    delete([obj.ampacdatapath, ampac.hash, '.out']);
                end
                if (exist([obj.ampacdatapath, ampac.hash, '.arc'], 'file'))
                    delete([obj.ampacdatapath, ampac.hash, '.arc']);
                end
                if (exist([obj.ampacdatapath, ampac.hash, '.dat'], 'file'))
                    delete([obj.ampacdatapath, ampac.hash, '.dat']);
                end
                if (exist([obj.ampacdatapath, ampac.hash, '.ido'], 'file'))
                    delete([obj.ampacdatapath, ampac.hash, '.ido']);
                end
                
                obj.progress(end) = obj.progress(end) + 1;
                
                if (~no_prog_window)
                    waithdldata = get(waithdl,'userdata');
                    if (waithdldata.bail == 0)
                        waitbar((prod(obj.progress) - 1) / obj.bigO, waithdl);
                    else
                        if (waithdldata.bail == 2)
                            obj.save_data();
                        end

                        disp('Calculation ended by user.');
                        close(waithdl,'force');
                        successful = 0;
                        return;
                    end
                end
                
                for i = length(obj.progress):-1:1
                    if (obj.progress(i) > obj.finishpt(i) && i ~= 1)
                        obj.progress(i) = 1;
                        obj.progress(i-1) = obj.progress(i-1) + 1;
                    end
                    
                    if (obj.progress(1) > obj.finishpt(1))
                        obj.progress = obj.finishpt;
                        getOut = true;
                    end
                end
                
                if (save_iterator == save_count - 1 && ~getOut)
                    if (~quiet_mode)
                        disp('Saving Checkpoint MAT file...');
                    end
                    [tok, ~] = regexpi(obj.MATFilename, '(.+)\\.+\.mat', 'tokens');
                    if (~exist(cell2mat(tok{1}(1)), 'dir'))
                        mkdir(cell2mat(tok{1}(1)));
                    end
                    matfile = obj.MATFilename;
                    save(matfile, 'obj');
                    
                    if (~quiet_mode)
                        disp(['Data saved to ',matfile]);
                    end
                    save_iterator = 0;

%                     elapsed = toc(mytic);
%                     completed_iter = obj.get_linear_index(obj.finishpt, obj.progress, true);
%                     time_remain = (elapsed / save_count) * (obj.bigO - completed_iter);
%                     dhms_time.s = round(mod(time_remain, 60));
%                     time_remain = floor(time_remain / 60);
%                     dhms_time.m = mod(time_remain, 60);
%                     time_remain = floor(time_remain / 60);
%                     dhms_time.h = mod(time_remain, 24);
%                     time_remain = floor(time_remain / 24);
%                     dhms_time.d = time_remain;

%                     disp(['Appx. time remaining: ', num2str(dhms_time.d), 'd ',...
%                         num2str(dhms_time.h), 'h ', num2str(dhms_time.m), 'm ', ...
%                         num2str(dhms_time.s), 's'])
%                     disp([sprintf('%0.2f',(completed_iter * 100 / obj.bigO)), '% complete']);
%                     mytic = tic();
                else
                    save_iterator = save_iterator + 1;
                end
                    
            end
            
            mindata = min(cell2mat({obj.data(:).Ehf}));

            for i = 1:obj.bigO
                if (~isempty(obj.data(i).Ehf))
                    obj.data(i).Ehf = (obj.data(i).Ehf - mindata) * (1/23.04);
                end
            end
                 
            if (obj.fail_iterator ~= 1)
                warning('off','backtrace');
                warning('EnergyCalcExp:CalcFailuresDetected',[num2str(obj.fail_iterator - 1),' calculation failure(s) were detected. ',...
                    'See the fail_bucket for further details.']);
                warning('on','backtrace');
                obj.fail_bucket = obj.fail_bucket(1:obj.fail_iterator - 1);
            else
                obj.fail_bucket = [];
            end
            
            disp('Saving Final MAT file...');
            [tok, ~] = regexpi(obj.MATFilename, '(.+)\\.+\.mat', 'tokens');
            if (~exist(cell2mat(tok{1}(1)), 'dir'))
                mkdir(cell2mat(tok{1}(1)));
            end
            matfile = obj.MATFilename;
            save(matfile, 'obj');
            disp(['Data saved to ',matfile]);
            
            if (~isempty(obj.ampac_lib))
                disp('Saving AMPAC library catalog...');
                lib = obj.ampac_lib; 
                save([obj.ampacdatapath, 'catalog.mat'], 'lib');
                disp('AMPAC library catalog saved');
            end
            
            if (~isempty(obj.indo_lib))
                disp('Saving INDO library catalog...');
                lib = obj.indo_lib; 
                save([obj.indodatapath, 'catalog.mat'], 'lib');
                disp('INDO library catalog saved');
            end
            
            if (~no_prog_window)
                close(waithdl,'force');
            end
            successful = (obj.fail_iterator == 1);
        end
%--------------------------------------------------------------------------------------------------%
% Function:     Plot3DData
% Description:  Plots 3-D graph of the data
% Input Args:
%       obj: EnergyCalcExp calling object
%       figure_num: unsigned int. Handle for the plot window
%       surface_not_contour: bool. TRUE plots a surface/contour plot by surfc(). FALSE plots a contour plot.
%       n_state: unsigned int. Value for the state to plot, 0 being the ground state.
%       varargin: cell array. Contains values of all angles desired, 'x' for x-axis, 'y' for y-axis.
%               A note about varargin: varargin contains any excess arguments specified after the
%               defined arguments. So, to specify the angle values, just include them as separate aruments.
%               Ex. For an experiment with 3 PHI angles: myECE.Plot3DData(1, true, 1, 'x', 'y', 90)
%               PHI1 is on the x-axis, PHI2 is on the y-axis and PHI3 = 90 degrees.
% Output:  
%       successful: bool. Basically a worthless output. TRUE just means it got to the end of the function.
%--------------------------------------------------------------------------------------------------%
        function successful = Plot3DData(obj, figure_num, surface_not_contour, n_state, varargin)
            if (length(obj.params.phi) == 1)
                throw(MException('EnergyCalcExp:Plot3DData:CannotMake3D',...
                    'EnergyCalcExp: There are not enough PHI parameters to make a 3-D plot.'));
            end
            
            x_axis = 0;
            y_axis = 0;
            
            for i = 1:length(varargin)
                if (isequal(varargin{i}, 'x') && x_axis == 0)
                    x_axis = i;
                    varargin{i} = ':';
                elseif (isequal(varargin{i}, 'y') && y_axis == 0)
                    y_axis = i;
                    varargin{i} = ':';
                elseif (~isnumeric(varargin{i}) || ~isequal(size(varargin{i}),[1 1]))
                    throw(MException('EnergyCalcExp:Plot3DData:Invalid_Input',...
                         'EnergyCalcExp: You have specified an invalid input for plotting.'));
                end
            end
            
            if (x_axis == 0 || y_axis == 0)
                throw(MException('EnergyCalcExp:Plot3DData:No_XY',...
                        'EnergyCalcExp: No X and/or Y specified for plotting.'));
            end
            
            if (n_state + 1 > obj.params.nexcstates || n_state < 0)
                throw(MException('EnercyCalcExp:Plot3DData:Invalid_State',...
                    'EnergyCalcExp: Invalid state specified for plotting.'));
            end
            
            figure(figure_num);
            hold on;
            
            GrdStEnergy = obj.get_field('Ehf', varargin{1:end});
            ExcEnergyData = GrdStEnergy + obj.get_field('Eexc', n_state + 1, varargin{1:end});
            
            if (x_axis > y_axis)
                GrdStEnergy = GrdStEnergy';
                ExcEnergyData = ExcEnergyData';
            end
            
            if (n_state == 0)
                if (surface_not_contour)
                    surfc(obj.params.phi{x_axis}, obj.params.phi{y_axis}, GrdStEnergy');
                    zlabel('Ground State Energy (eV)')
                else
                    contour(obj.params.phi{x_axis}, obj.params.phi{y_axis}, GrdStEnergy');
                    title('Ground State Energy (eV)')
                end
            else
                if (surface_not_contour)
                    surfc(obj.params.phi{x_axis}, obj.params.phi{y_axis}, ExcEnergyData');
                    zlabel('First Excited State Energy (eV)')
                else
                    contour(obj.params.phi{x_axis}, obj.params.phi{y_axis}, ExcEnergyData');
                    title('First Excited State Energy (eV)')
                end
            end

            xlabel(['PHI',num2str(x_axis)]);
            ylabel(['PHI',num2str(y_axis)]);

            figure(figure_num);

            successful = 1;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     Plot2DData
% Description:  Plots 2-D graph of the data
% Input Args:
%       obj: EnergyCalcExp calling object
%       figure_num: unsigned int. Handle for the plot window
%       plot_intensity: bool. TRUE plots intensity as circles around the excited state data points
%       n_states: int[]. Array of states to plot. 0 is the ground state.
%       varargin: cell array. Contains values of all angles desired, 'x' for x-axis. See note about varargin above.
% Output:  
%       successful: bool. Basically a worthless output. TRUE just means it got to the end of the function.
%--------------------------------------------------------------------------------------------------%
        function successful = Plot2DData(obj, figure_num, plot_intensity, n_states, varargin)

            x_axis = 0;
            
            for i = 1:length(varargin)
                if (isequal(varargin{i}, 'x') && x_axis == 0)
                    x_axis = i;
                    varargin{i} = ':';
                elseif (~isnumeric(varargin{i}) || ~isequal(size(varargin{i}),[1 1]))
                    throw(MException('EnergyCalcExp:Plot2DData:Invalid_Input',...
                         'EnergyCalcExp: You have specified an invalid input for plotting.'));
                end
            end
            
            if (x_axis == 0)
                throw(MException('EnergyCalcExp:Plot2DData:No_X',...
                        'EnergyCalcExp: No X-axis specified for plotting.'));
            end
            
            if (max(n_states) + 1 > obj.params.nexcstates || min(n_states) < 0)
                throw(MException('EnercyCalcExp:Plot2DData:Invalid_State',...
                    'EnergyCalcExp: Invalid state specified for plotting.'));
            end

            figure(figure_num);
            hold on;
            
            EnergyData = repmat(obj.get_field('Ehf', varargin{1:end}), 1, length(n_states)) + ...
                obj.get_field('Eexc', n_states + 1, varargin{1:end});
            
            if (obj.fail_iterator ~= 1 && ismember(0, n_states))
                EnergyData(:, n_states == 0) = obj.get_field('Ehf', varargin{1:end});
            end
            
            if (plot_intensity)
                IntensityData = obj.get_field('Tint', n_states + 1, varargin{1:end});
                maxInt = max(max(IntensityData));
                IntensityData = 30*IntensityData / maxInt + 1e-3;
            end
                           
            for is = 1:length(n_states)
                if (n_states(is) == 0)
                    plot(obj.params.phi{x_axis}, EnergyData(:,is), 'g^');
                else
                    plot(obj.params.phi{x_axis}, EnergyData(:,is), 'b.');
                end
                
                if (plot_intensity && (n_states(is) ~= 0) && (length(IntensityData) ~= 1))
                    for iangles = 1:length(obj.params.phi{x_axis})
                        if (~isnan(IntensityData(iangles, is)))
                            plot(obj.params.phi{x_axis}(iangles), EnergyData(iangles,is),'ro' ...
                                 ,'MarkerSize',IntensityData(iangles,is));
                        end
                    end
                end
            end
            
            xlabel(['PHI', num2str(x_axis)]);
            ylabel('Energy (eV)');
            successful = 1;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     get_data
% Description:  Retrieves stored data based on the angle values specified.
% Input Args:
%       obj: EnergyCalcExp calling object
%       varargin: cell array. Contains values of all angles desired in the output.
% Output:  
%       res: ECEDataStruct[]. n-dimensional ECEDataStruct matrix containing the values specified. All dims with a
%               size of 1 are removed. The dims of the matrix are in numerical order of the PHI
%               values which were larger than 1 in size.
%               Ex. For 3 PHI values, myECE.get_data([90 92 106 111 120 151], 180, 0:2:20)
%               will return a 6x11 matrix (PHI1, PHI3)
%--------------------------------------------------------------------------------------------------%
        function res = get_data(obj, varargin)
            if (length(varargin) ~= length(obj.params.phi))
                throw(MException('EnergyCalcExp:get_data:Not_Enough_Args',...
                    'EnergyCalcExp: Not enough arguments specified for PHI'));
            end
            
            index = cell(1,length(varargin));
            
            for i = 1:length(varargin)
                if (varargin{i} == ':')
                    index{i} = 1:length(obj.params.phi{i});
                else
                    [~, loc] = ismember(obj.params.phi{i}, varargin{i});
                    [~, c] = find(loc ~= 0);
                    
                    if (length(c) ~= length(varargin{i}))
                        throw(MException('EnergyCalcExp:get_data:Not_Member_of_PHI',...
                            ['EnergyCalcExp: One of the members specified for PHI',...
                            num2str(i),' does not exist in the data set.']));
                    end
                    
                    index{i} = c;
                end                
            end
            
            mystruct = struct('type',{'()'},'subs',{index});
            subset = subsref(obj.data, mystruct);
            newsize = size(subset);
            newsize = newsize(newsize ~= 1);
            
            if (isempty(newsize))
                res = reshape(subset,1,1);
            elseif (length(newsize) == 1)
                res = reshape(subset,[],newsize);
            else
                res = reshape(subset,newsize);
            end                
        end
%--------------------------------------------------------------------------------------------------%
% Function:     get_field
% Description:  Retrieves data stored in fields based on the angle values specified.
% Input Args:
%       obj: EnergyCalcExp calling object
%       field: char[]. Name of the field. If the field is listed in the EnergyCalcExp class definition,
%               form is 'fieldname'. If it is in the raw data, form is 'ampac' or 'indo' followed by the field name
%               or function separated by a period. Ex. 'ampac.element'. A list of available fields not listed in
%               this class definition is listed within the function as mem_loadable_field and mem_loadable_functions.
%               NOTE: For large data retrievals from fields/functions which are not in the class definition,
%               as the number of data points increases, the retrieval time drastically increases because the data
%               must be loaded from the harddisk instead of from memory.
%       varargin: cell array. The first section contains any arguments/subscripts required for the fields/functions.
%               The last arguments should contain values of all angles desired in the output.
% Output:  
%       res: double[]. n-dimensional matrix containing the values specified. Data returned from function calls
%               will be loaded into 'res' as a cell array instead of a matrix. All dims with a
%               size of 1 are removed. The dims of the matrix are in numerical order of the PHI
%               values which were larger than 1 in size. The first dimensions are for the PHI values
%               and the last dims are for the arguments/subscripts in the fields.
%               Ex. For 3 PHI values and one subscript:
%               myECE.get_field('Eexc', 1:5, 90:100, 180, 0:2:28)
%               will return a 11x15x5 matrix (PHI1, PHI3, subscript)
%--------------------------------------------------------------------------------------------------%
        function res = get_field(obj, field, varargin)
            mem_loadable_fields = {'ampac.r'; 'ampac.element'; 'ampac.natom';...
                'ampac.Hf'; 'indo.norb'; 'indo.aorbAtom'; 'indo.aorbType';...
                'indo.nfilled'; 'indo.orbE'; 'indo.orb'; 'indo.nsci';...
                'indo.nscibasis'; 'indo.esci'; 'indo.r'; 'indo.wfsci';...
                'indo.ehsci'; 'indo.hfE'};
            mem_loadable_functions = {'indo.get_osc'; 'indo.dipole'};
            
            if (ismember(field, fieldnames(obj.data(1))))
                type_of_field = 1;
            elseif (ismember(field, mem_loadable_fields))
                type_of_field = 2;
            elseif (ismember(field, mem_loadable_functions))
                type_of_field = 3;
            else
                throw(MException('EnergyCalcExp:get_field:Invalid_Field',...
                    'EnergyCalcExp: No field of that name available.'));
            end
            
            if (length(varargin) < length(obj.params.phi))
                throw(MException('EnergyCalcExp:get_field:Not_Enough_Args',...
                    'EnergyCalcExp: Not enough arguments specified for PHI'));
            end
            
            mydat = obj.get_data(varargin{end-length(obj.params.phi)+1:end});
            
            if (length(varargin) == length(obj.params.phi))
                varargin = {};
            else
                varargin = {varargin{1:(end-length(obj.params.phi))}}; 
            end
            
            switch (type_of_field)
                case 1
                    res = reshape({mydat(:).(deblank(field))}, size(mydat));
                    ind = cellfun(@isempty, res);
                    checkforempty = find(ind, 1);
                    
                    if (~isempty(checkforempty))
                        ind = find(ind == 0);
                        if (isempty(ind))
                            res = NaN;
                            return;
                        end
                        sizefield = size(cell2mat(res(ind(1))));
                        sqzsizefield = sizefield(sizefield ~= 1);
                        argsize = zeros(1, length(varargin));
                        if (~isempty(varargin))
                            for i = 1:length(varargin)
                                if (ischar(varargin{i}))
                                    if (length(sizefield) == 2)
                                        varargin{i} = 1:sqzsizefield; %#ok<AGROW>
                                    end
                                end
                                argsize(i) = numel(varargin{i});
                            end
                            for i = 1:numel(res)
                                if (~isempty(res{i}))
                                    res{i} = res{i}(varargin{1:end});
                                else
                                    res{i} = ones([1 argsize]) * NaN;
                                    res{i} = squeeze(res{i});
                                end
                            end
                        else
                            for i = 1:numel(res)
                                if (isempty(res{i}))
                                    res{i} = ones(sizefield) * NaN;
                                end
                            end
                        end
                        
                        res = cell2mat({res{:}}); 
                        res = reshape(res, [argsize size(mydat)]);
                        res = squeeze(res);
                        res = shiftdim(res, length(argsize(argsize ~= 1)));
                    else
                        sizefield = size(mydat(1).(deblank(field)));
                        res = cell2mat({res{:}}); 
                        res = reshape(res, [sizefield size(mydat)]);
                        res = squeeze(res);
                        
                        if (length(find((size(res) ~= 1) == 1)) > 1)
                            res = shiftdim(res, length(sizefield(sizefield ~= 1)));
                        end

                        subs = cell(1,length(find((size(res) ~= 1) == 1)));

                        for i = 1:length(find((size(res) ~= 1) == 1))
                            if (i <= length(find(size(mydat) ~= 1)))
                                subs{i} = 1:size(res,i);
                            else
                                if (~isempty(varargin))
                                    subs{i} = varargin{i-length(find(size(mydat) ~= 1))};
                                else
                                    subs{i} = [1]; %#ok<NBRAK>
                                end
                            end
                        end

                        mystruct = struct('type',{'()'},'subs',{subs});
                        res = subsref(res, mystruct);
                    end
                case 2
                    [tok, ~] = regexp(field, '(\w*)\.(.*)', 'tokens');
                    func = str2func(cell2mat(['get_',tok{1}(1)]));
                    field = cell2mat(tok{1}(2));
                    % sizeoffield = size(func(mydat(1), field));
                    % sizeoffield = sizeoffield(sizeoffield ~= 1);
                    
                    res = cell(size(mydat));
                    
                    for i = 1:numel(mydat)
                        % myloc = EnergyCalcExp.get_subscripts(size(mydat), i);
                        res{i} = func(mydat(i), field);
                    end
                    
                    ind = cellfun(@isempty, res);
                    ind = find(ind == 0);
                    if (isempty(ind))
                        sizefield = 1;
                    else
                        sizefield = size(res{ind(1)});
                    end
                    
                    for i = 1:numel(res)
                        if (isempty(res{i}))
                            res{i} = ones(sizefield) * NaN;
                        end
                    end
                    
                    try
                        checkres = cell2mat({res{1:end}}); %#ok<*CCAT1>
                        res = reshape(checkres, [sizefield(sizefield ~= 1) size(mydat)]);
                        res = squeeze(res);
                        if (length(find((size(res) ~= 1) == 1)) > 1)
                            res = shiftdim(res, length(sizefield(sizefield ~= 1)));
                        end
                        res = subsref(res, struct('type',{'()'},'subs',...
                            {cat(2,num2cell(repmat(':',1,length(find((size(mydat) ~= 1) == 1)))),varargin{1:end})}));
                    catch exception
                        if (~strcmp('MATLAB:getReshapeDims:notSameNumel', exception.identifier) && ...
                            ~strcmp('MATLAB:catenate:dimensionMismatch', exception.identifier))
                        
                            rethrow(exception);
                        end
                        for i = 1:numel(res)
                            if (ismember(i, ind))
                                sizeofres = size(res{i});
                                if (length(sizeofres(sizeofres ~= 1)) == 1)
                                    if (isequal(find(sizeofres ~= 1), find(size(res{i}(varargin{1:end})) ~= 1)))
                                        res{i} = res{i}(varargin{1:end});
                                    else
                                        res{i} = res{i}(varargin{1:end});
                                        res{i} = res{i}';
                                    end
                                end                                   
                            else
                                res{i} = NaN;
                            end
                        end
                    end
                case 3
                    [tok, ~] = regexp(field, '(\w*)\.(.*)', 'tokens');
                    datatype = cell2mat(tok{1}(1));
                    fieldname = cat(2,'raw_',datatype);
                    func = str2func(cell2mat(tok{1}(2)));
                    
                    res = cell(size(mydat));
                    for i = 1:numel(mydat)
                        if (mydat(i).(deblank(cat(2, datatype, '_succeed'))) == 1)
                            mydat(i).load_to_memory(datatype,'load');
                            if (isempty(varargin))
                                res{i} = func(mydat(i).(deblank(fieldname)));
                            else
                                res{i} = func(mydat(i).(deblank(fieldname)), varargin{1:end});
                            end
                            mydat(i).load_to_memory(datatype,'unload');
                        else
                            res{i} = NaN;
                        end
                    end
                    
                    try
                        ind = cellfun(@(x)isequal(x, NaN), res);
                        ind = find(ind == 0);
                        if (isempty(ind))
                            sizefield = 1;
                        else
                            sizefield = size(res{ind(1)});
                        end
                        checkres = cell2mat({res{1:end}});
                        res = reshape(checkres, [sizefield(sizefield ~= 1) size(mydat)]);
                        if (length(find((size(res) ~= 1) == 1)) > 1)
                            res = shiftdim(res, length(sizefield(sizefield ~= 1)));
                        end
                    catch exception
                        if (~strcmp('MATLAB:getReshapeDims:notSameNumel', exception.identifier) && ...
                            ~strcmp('MATLAB:catenate:dimensionMismatch', exception.identifier))
                        
                            rethrow(exception);
                        end
                    end
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     update_paths
% Description:  Changes the Ampac and/or INDO paths if the files have been moved.
% Input Args:
%       obj: EnergyCalcExp calling object
%       newampac: char[]. Full path to the new AMPAC directory. If [], the AMPAC path will not be changed.
%       newindo: char[]. Full path to the new INDO directory. If [], the INDO path will not be changed.
%       save_mat_file: bool. TRUE saves the MAT file to reflect the changes.
% Output:  
%       successful: bool. Basically useless output. TRUE means we got to the end of the function.
%--------------------------------------------------------------------------------------------------%
        function successful = update_paths(obj, newampac, newindo, save_mat_file)
            if (~isempty(newampac))
                if (isempty(regexpi(newampac, '(\w+:\\.*)|(\\\\.+\\.+)', 'start')) || ...
                    (regexpi(newampac, '(\w+:\\.*)|(\\\\.+\\.+)', 'start') ~= 1))
                    
                    throw(MException('EnergyCalcExp:update_paths:Invalid_Ampac_Path',...
                        'EnergyCalcExp: Invalid Ampac path. Path must be complete and either a network or local disk path.'));
                end
                
                if (newampac(end) ~= '\')
                    newampac = [newampac, '\'];
                end
                
                if (exist([newampac, 'catalog.mat'], 'file'))
                    S = load([newampac, 'catalog.mat'], 'lib');
                    alib = S.lib;
                    ampachashes = {obj.data(:).ampac_hash};

                    if (~ismember(ampachashes, alib.keys()))
                        throw(MExcpetion('EnergyCalcExp:update_paths:NewAmpacLibMissingData',...
                            'EnergyCalcExp: New Ampac library does not contain necessary data.'));
                    end
                
                    
                else
                    files = dir([newampac, '*.mat']);
                    fn = cell(1,numel(files)); 
                    fn = {files(1:end).name};
                    fn = cellfun(@(x)regexpi(x, '^([A-Fa-f0-9]){32}\.mat', 'once', 'tokens'), fn, 'UniformOutput', false);
                    fn = {fn{~(cellfun(@(x)isempty(x), fn))}};
                    fn = cat(2, fn{1:end});
                    
                    ampachashes = {obj.data(:).ampac_hash};
                    
                    if (length(fn) >= numel(obj.data))
                        nonmember = ismember(ampachashes, fn);
                        nonmember = nonmember(nonmember == 0);
                        
                        if (~isempty(nonmember))
                            throw(MExcpetion('EnergyCalcExp:update_paths:NewAmpacLibMissingData',...
                                'EnergyCalcExp: New Ampac library does not contain necessary data.'));
                        end
                    else
                        throw(MExcpetion('EnergyCalcExp:update_paths:NewAmpacLibMissingData',...
                            'EnergyCalcExp: New Ampac library does not contain necessary data.'));
                    end
                end
                
                obj.ampacdatapath = newampac; 
            end
            
            if (~isempty(newindo))
                if (isempty(regexpi(newindo, '(\w+:\\.*)|(\\\\.+\\.+)', 'start')) || ...
                    (regexpi(newindo, '(\w+:\\.*)|(\\\\.+\\.+)', 'start') ~= 1))
                    
                    throw(MException('EnergyCalcExp:update_paths:Invalid_Indo_Path',...
                        'EnergyCalcExp: Invalid Indo path. Path must be complete and either a network or local disk path.'));
                end
                
                if (newindo(end) ~= '\')
                    newindo = [newindo, '\'];
                end
                
                if (exist([obj.newindo, 'catalog.mat'], 'file'))
                    S = load([newindo, 'catalog.mat'], 'lib');
                    ilib = S.lib;
                    indohashes = {obj.data(:).indo_hash};
                    if (~ismember(indohashes, ilib.keys()))
                        throw(MExcpetion('EnergyCalcExp:update_paths:NewINDOLibMissingData',...
                            'EnergyCalcExp: New INDO library does not contain necessary data.'));
                    end
                else
                    files = dir([newindo, '*.mat']);
                    fn = cell(1,numel(files)); %#ok<*NASGU>
                    fn = {files(1:end).name(1:end-4)};
                    fn = cellfun(@(x)regexpi(x, '^([A-Fa-f0-9]){32}\.mat', 'match'), fn, 'UniformOutput', false);
                    fn = {fn{~(cellfun(@(x)isempty(x), fn))}};
                    fn = cat(2, fn{1:end});
                    
                    indohashes = {obj.data(:).indo_hash};
                    
                    if (length(fn) >= numel(obj.data))
                        nonmember = ismember(indohashes, fn);
                        nonmember = nonmember(nonmember == 0);
                        
                        if (~isempty(nonmember))
                            throw(MExcpetion('EnergyCalcExp:update_paths:NewAmpacLibMissingData',...
                                'EnergyCalcExp: New Ampac library does not contain necessary data.'));
                        end
                    else
                        throw(MExcpetion('EnergyCalcExp:update_paths:NewAmpacLibMissingData',...
                            'EnergyCalcExp: New Ampac library does not contain necessary data.'));
                    end
                end
                
                obj.indodatapath = newindo;
            end

            for i = 1:numel(obj.data)
                obj.data(i).update_paths(newampac, newindo);
            end
            
            if (save_mat_file)
                matfile = obj.MATFilename;
                save(matfile, 'obj');
            end
            
            successful = 1;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     run_opt_only
% Description:  Runs an optimization calculation. Runs if no PHI angles are specified.
% Input Args:
%       obj: EnergyCalcExp calling object
%       quiet_mode: bool. TRUE causes output to be suppressed to the command window.
% Output:  
%       successful: bool. TRUE means there were no errors in the calculation.
%--------------------------------------------------------------------------------------------------%
        function successful = run_opt_only(obj, quiet_mode)            
            obj.data = ECEDataStruct();
            obj.fail_bucket = struct('phiID', [], 'type', [], 'msg', []);
            obj.bigO = 1;
            
            ampac = obj.run_Ampac([], quiet_mode);
            
            if (ampac.ampac_succeed)
                if (obj.params.doIndo)
                    warning('off','MATLAB:structOnObject');
                    indoparams = struct(obj.params);
                    warning('on','MATLAB:structOnObject');
                    indoparams.charge = obj.mol_charge;
                    indoparams.norbs = obj.params.norbs_basis_set;
                    indoparams.nstates = obj.params.nexcstates;
                    
                    if (isempty(obj.params.Efield_vector) && ...
                        ~(obj.params.Efield_mag == 0 || isempty(obj.params.Efield_mag)) && ...
                        (isempty(obj.axis_params)))

                        throw(MException('EnergyCalcExp:run:No_Efield_Vector',...
                            'EnergyCalcExp: No data provided for E field vector'));
                    elseif (obj.params.Efield_mag == 0 || isempty(obj.params.Efield_mag))
                        obj.params.setCustomEFVector([0 0 0]);
                    end

                    if (isempty(obj.params.Efield_type) || (obj.params.Efield_type ~= 'c'))
                        obj.params.setEFVectorFromAtomNums(ampac, obj.axis_params);
                    end

                    indoparams.field = obj.params.Efield_vector * obj.params.Efield_mag;
                    
                    key = obj.keygen('indo', []);
                    hash = DataHash(key);
                    data_exist = false;
                    indo_out = struct();
                    indo_out.key = key;
                    indo_out.hash = [];

                    while (isempty(indo_out.hash))
                        if (~isempty(obj.indo_lib))
                            if (obj.indo_lib.isKey(hash))
                                if (~isequal(obj.indo_lib(hash), key))
                                    hash = DataHash(hash);
                                else
                                    indo_out.hash = hash;
                                    data_exist = true;
                                end
                            else
                                indo_out.hash = hash;
                            end
                        else
                            if (exist([obj.indodatapath, hash, '.mat'], 'file'))
                                S = load([obj.indodatapath, hash, '.mat'], 'key');
                                if (~isequal(S.key, key))
                                    hash = DataHash(hash);
                                else
                                    indo_out.hash = hash;
                                    data_exist = true;
                                end
                            else
                                indo_out.hash = hash;
                            end
                        end
                    end

                    if (~data_exist)
                        if (~quiet_mode)
                            disp('Running INDO calculation');
                        end
                        if (obj.ampacdatapath(end-1) == ':')
                            indo = Indo(indoparams, obj.ampacdatapath(1:end), ampac.hash);
                        else
                            indo = Indo(indoparams, obj.ampacdatapath(1:end-1), ampac.hash);
                        end
                        indo_out.indo = indo;
                        
                        if (indoparams.output_dm)
                            movefile([obj.ampacdatapath, ampac.hash, '-dm.bin'], ...
                                [obj.indodatapath, hash, '-dm.bin']);
                        end
                        
                        if (~isempty(obj.indo_lib))
                            obj.indo_lib(hash) = key;
                            lib = obj.indo_lib;
                            save([obj.indodatapath, 'catalog.mat'], 'lib', '-append');
                        end
                        
                        save([obj.indodatapath, hash, '.mat'], '-struct', 'indo_out');
                    else
                        S = load([obj.indodatapath, hash, '.mat'], 'indo');
                        indo = S.indo;
                    end

                    if (~indo.indo_succeed)
                        if (~quiet_mode)
                            disp(['INDO Failed. ', indo.indo_err_msg]);
                        end
                        obj.fail_bucket(obj.fail_iterator) = ...
                            struct('phiID', [], 'type', 'indo', 'msg', indo.indo_err_msg);
                        obj.fail_iterator = obj.fail_iterator + 1;
                    elseif (~quiet_mode)
                        disp(['INDO finished. Done parsing INDO.', indo.indo_err_msg]);
                    end

                    obj.data(1) = ...
                        ECEDataStruct(ampac, indo, [], obj.ampacdatapath, ampac.hash,...
                        obj.indodatapath, indo_out.hash);

                else
                    obj.data(1) = ...
                        ECEDataStruct(ampac, [], [], obj.ampacdatapath, ampac.hash, [], []);
                end
            else
                obj.fail_bucket(obj.fail_iterator) = ...
                    struct('phiID', [], 'type', 'ampac', 'msg', ampac.ampac_err_msg);
                obj.fail_iterator = obj.fail_iterator + 1;
                obj.data(1) = ...
                    ECEDataStruct(ampac, [], [], obj.ampacdatapath, ampac.hash, [], []);
            end
            
            if (exist([obj.ampacdatapath, ampac.hash, '.vis'], 'file'))
                delete([obj.ampacdatapath, ampac.hash, '.vis']);
            end
            if (exist([obj.ampacdatapath, ampac.hash, '.out'], 'file'))
                delete([obj.ampacdatapath, ampac.hash, '.out']);
            end
            if (exist([obj.ampacdatapath, ampac.hash, '.arc'], 'file'))
                delete([obj.ampacdatapath, ampac.hash, '.arc']);
            end
            if (exist([obj.ampacdatapath, ampac.hash, '.dat'], 'file'))
                delete([obj.ampacdatapath, ampac.hash, '.dat']);
            end
            if (exist([obj.ampacdatapath, ampac.hash, '.ido'], 'file'))
                delete([obj.ampacdatapath, ampac.hash, '.ido']);
            end
            
            if (obj.fail_iterator == 1 || ~strcmp(obj.fail_bucket(1).type, 'ampac'))
                obj.data(1).Ehf = obj.data(1).Ehf * (1/23.06);
            end
            if (obj.fail_iterator ~= 1)
                warning('off','backtrace');
                warning('EnergyCalcExp:CalcFailuresDetected',['A calculation failure was detected. ',...
                    'The data should still be accessible through the ',...
                    'built-in functions, but plotting may fail. See the fail_bucket for further details.']);
                warning('on','backtrace');
            else
                obj.fail_bucket = [];
            end
            
            disp('Saving Final MAT file...');
            [tok, ~] = regexpi(obj.MATFilename, '(.+)\\.+\.mat', 'tokens');
            if (~exist(cell2mat(tok{1}(1)), 'dir'))
                mkdir(cell2mat(tok{1}(1)));
            end
            matfile = obj.MATFilename;
            save(matfile, 'obj');
            disp(['Data saved to ',matfile]);
            
            if (~isempty(obj.ampac_lib))
                disp('Saving AMPAC library catalog...');
                lib = obj.ampac_lib;
                save([obj.ampacdatapath, 'catalog.mat'], 'lib');
                disp('AMPAC library catalog saved');
            end
            
            if (~isempty(obj.indo_lib))
                disp('Saving INDO library catalog...');
                lib = obj.indo_lib;
                save([obj.indodatapath, 'catalog.mat'], 'lib');
                disp('INDO library catalog saved');
            end
            
            successful = (obj.fail_iterator == 1);
        end
%--------------------------------------------------------------------------------------------------%
% Function:     save_data
% Description:  Saves any changes made to the MAT file output by the EnergyCalcExp object.
% Input Args:
%       obj: EnergyCalcExp calling object
% Output:  
%       successful: bool. TRUE means we got to the end of the function.
%--------------------------------------------------------------------------------------------------%
        function successful = save_data(obj)
            disp('Saving MAT file...');
            matfile = obj.MATFilename;
            save(matfile, 'obj');
            disp(['Data saved to ',matfile]);
            
            successful = 1;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     keygen
% Description:  Generates a struct containing the data necessary for generating the hash key
%               for data storage in the library
% Input Args:
%       obj: EnergyCalcExp calling object
%       key_type: char[]. Specifies what type of library key to generate. 'ampac' or 'indo'
%       phi: int/double[]: array containing the angle values for the key to be generated
% Output:  
%       key: struct. Contains the fields for the key.
%--------------------------------------------------------------------------------------------------%
        function key = keygen(obj, key_type, phi)
            key.method = lower(obj.params.method);
            key.zmat = lower(obj.zmatrix_template);
            
            if (isempty(phi))
                key.phi = [];
            else
                key.phi = reshape(phi, 1, []);
                for i = 1:length(phi)
                    key.zmat = strrep(key.zmat, ['phi', num2str(i)], num2str(phi(i), '%3.6f'));
                end
            end
            
            key.zmat = textscan(key.zmat,'%s','delimiter','\n');
            key.zmat = {key.zmat{1}{4:end-1}};
            key.zmat = cellfun(@(x)textscan(x, '%s'), key.zmat);
            
            for i = 1:length(key.zmat)
                key.zmat{i}{2} = num2str(str2double(key.zmat{i}{2}),'%3.6f');
                key.zmat{i}{4} = num2str(str2double(key.zmat{i}{4}),'%3.6f');
                key.zmat{i}{6} = num2str(str2double(key.zmat{i}{6}),'%3.6f');
            end            
            
            if (strcmp(key_type, 'indo'))
                key.Efield_mag = obj.params.Efield_mag;
                
                if (key.Efield_mag ~= 0)
                    if (obj.params.Efield_type == 'c')
                        key.Efield_type = 'c';
                        key.Efield_params = reshape(obj.params.Efield_vector, [1 3]);
                    elseif (length(obj.axis_params) == 2)
                        key.Efield_type = 'p';
                        key.Efield_params = reshape(obj.axis_params, [1 2]);
                    else
                        key.Efield_type = 'n';
                        key.Efield_params = reshape(obj.axis_params, [1 3]);
                    end
                else
                    key.Efield_type = 'c';
                    key.Efield_params = [0 0 0];
                end
            elseif (~strcmp(key_type, 'ampac'))
                throw(MException('EnergyCalcExp:keygen:BadType',...
                    'EnergyCalcExp: Please specify "ampac" or "indo" for the type'));
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     load_catalogs
% Description:  Loads libraries from saved path names. Useful if the experiment is loaded without
%               using the ReadMATFile() function.
% Input Args:
%       obj: EnergyCalcExp calling object
% Output:  
%       successful: bool. TRUE if both libraries were loaded successfully.
%--------------------------------------------------------------------------------------------------%
        function successful = load_catalogs(obj)
            successful = true;
            
            if (exist([obj.ampacdatapath, 'catalog.mat'], 'file'))
                S = load([obj.ampacdatapath, 'catalog.mat'], 'lib');
                obj.ampac_lib = S.lib;
            else
                warning('off','backtrace');
                warning('EnergyCalcExp:NoAmpacLibrary',['WARNING: Ampac library catalog not found at saved location. ',...
                    'Please run update_paths() if this is not intentional.']);
                warning('on','backtrace');
                obj.ampac_lib = [];
                successful = false;
            end
            
            if (exist([obj.indodatapath, 'catalog.mat'], 'file'))
                S = load([obj.indodatapath, 'catalog.mat'], 'lib');
                obj.indo_lib = S.lib;
            else
                warning('off','backtrace');
                warning('EnergyCalcExp:NoIndoLibrary',['INDO library catalog not found at saved location. ',...
                    'Please run update_paths() if this is not intentional.']);
                warning('on','backtrace');
                obj.indo_lib = [];
                successful = false;
            end
        end     
    end
end

