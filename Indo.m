classdef Indo < handle
    properties (SetAccess = private)
        % Input parameters
        config     % see defaultConfig() for contents
        dataPath   % directory for the data (do not end with \)
        jobName    % will read structure from ampac out file  jobname.out
        % and store indo results in jobname.ido
        % Output success/fail
        indo_succeed;   % TRUE means Indo calcs succeeded, FALSE means an exception was thrown
        indo_err_msg;   % Error msg from output
        % Atomic basis set:
        norb       % number of atomic basis functions, and hf orbitals
        aorbAtom   % (1,i) atom on which the ith atomic orbital resides
        aorbType   % (1,i) type of ith atomic orbital
        %  {s1=0,s2=1,p2x=2,p2y=3,p2z=4,s3=5,p3x=6,p3y=7,p3z=8}
        % Hartree Fock Results:
        nfilled    % number of filled molecular orbitals
        hfE        % HF ground state energy
        orbE       % (1,i) energy of ith orbital
        orb        % (i,j) ith component of jth orbital
        indoOutput % output from indo
        % SCI results
        nsci        % number of sci states, with first being ground state
        nscibasis   % number of basis functions (first being ground state)
        esci        % (1,i) energies of ith sci state
        r           %( i,j, icomp) transition (position) operator icomp = x,y,z
        wfsci       % (i,j)  ith component of jth state
        ehsci       % (i,1) hole of the ith SCI basis function (0 if GS)
        % (i,2) = elec of the ith SCI basis function (0 if GS)

        % indoExe = '"c:\mscpp\demo-dci-backup.exe"';
        indoExe = '"c:\mscpp\demo-dci\Release\demo-dci.exe"';
    end % properties
    properties (Transient)
        osc         % (1,i) oscillator strength from gs to state i
        rx          % r(:,:,1) for backwards compatibility, don't use now
        ry          % r(:,:,2)
        rz          % r(:,:,3)
    end 
    methods (Access = private)
        function readIndo(obj, inputfile)
            if (nargin > 1)
                filename = inputfile;
            else
                filename = [obj.dataPath,'\',obj.jobName,'.ido'];
            end
            
            fid1 = fopen(filename);
            if (fid1 == -1)
               error(['in Indo.readIndo, could not find file: ',filename]);
            end
            obj.norb = fread(fid1,1,'integer*4');
            obj.aorbAtom = fread(fid1,[1,obj.norb],'integer*4');
            obj.aorbAtom = obj.aorbAtom +1; % C++ starts count at 0 instead of 1
            obj.aorbType = fread(fid1,[1,obj.norb],'integer*4');

            ntest = fread(fid1,1,'integer*4');
            if (ntest ~= obj.norb)
               error('atomic and fock basis sizes differ');
            end
            obj.nfilled = fread(fid1,1,'integer*4');
            obj.hfE  = fread(fid1,1,'real*8');
            obj.orbE = fread(fid1,[1,obj.norb],'real*8');
            obj.orb = fread(fid1,[obj.norb,obj.norb],'real*8');

            obj.nsci = fread(fid1,1,'integer*4');
            obj.nscibasis = fread(fid1,1,'integer*4');
            obj.esci = fread(fid1,[1,obj.nsci],'real*8');
            ntest = fread(fid1,[1,2],'integer*4');
            obj.r = zeros(obj.nsci,obj.nsci,3);
            obj.r(:,:,1) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
            ntest = fread(fid1,[1,2],'integer*4');
            obj.r(:,:,2) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
            ntest = fread(fid1,[1,2],'integer*4');
            obj.r(:,:,3) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
            temp = fread(fid1,[2,obj.nscibasis],'integer*4');
            obj.ehsci = temp' +1; % +1 fixes the counting from 0 issue
            obj.wfsci = fread(fid1,[obj.nscibasis,obj.nsci],'real*8');
            fclose(fid1);
        end
    end
    methods (Static)
        function res = defaultConfig()
            res.charge = 1;
            res.norbs = 100;
            res.nstates = 25;
            res.field = [0,0,0];
            res.initial_shiftc = 0.0;
            res.initial_shift_step = 0.0;
            res.min_shift_step = 0.0;
            res.max_shift_step = 0.0;
            res.initial_second_shift_step = 0.0;
            res.min_second_shift_step = 0.0;
            res.max_second_shift_step = 0.0;
            res.initial_eeint = 1.0;
            res.initial_eestep = 0.0;
            res.min_eestep = 0.0;
            res.max_eestep = 0.0;
            res.initial_conv = 1e-10;
            res.min_conv = 1e-10;
            res.max_inner_iter = 300;
            res.max_iter = 10000;
            res.dm_guess = 'default';
            res.try_default_first = false;
            res.output_dm = false;
            res.pot_file = [];
        end
        
        function obj = LoadExistingData(filename, ConfigIn, dataPathIn, JobNameIn)
            try
            if (~exist(filename,'file'))
                throw(MException('Indo:LoadExistingData:FileDNE',...
                    'Indo: Cannot locate INDO file'));
            end
            catch exception
                disp(exception);
            end
            
            obj = Indo();
            
            listing = dir(filename);
            if (listing.bytes == 0)
                obj.indo_succeed = false;
                obj.indo_err_msg = 'File size is 0 bytes. Calc failed on previous run.';
            else
                obj.readIndo(filename);
                obj.indo_succeed = true;
                obj.indo_err_msg = [];
            end         
            
            obj.config = ConfigIn;
            obj.jobName = JobNameIn;
            obj.dataPath = dataPathIn;
        end
    end
    methods       
        function res = Indo(ConfigIn, dataPathIn, jobNameIn)
            if (nargin < 1)
                res.config = Indo.defaultConfig();
            else
                res.config = ConfigIn;
            end
            if (nargin < 2)
                res.dataPath = 'data';
            else
                res.dataPath = dataPathIn;
            end
            if (nargin < 3)
                res.jobName = 'jobname';
            else
                res.jobName = jobNameIn;
            end
            if (nargin > 0)
%                 jobstring = [res.indoExe,' "',res.dataPath,'\',res.jobName,'"', ...
%                     ' ',num2str(res.config.charge), ...
%                     ' ',num2str(res.config.norbs), ...
%                     ' ',num2str(res.config.nstates)];
%                 if (norm(res.config.field) > 0.0)
%                     jobstring = [jobstring, ...
%                         ' ', num2str(res.config.field(1)), ...
%                         ' ', num2str(res.config.field(2)), ...
%                         ' ', num2str(res.config.field(3))];
%                 end
                
                res.write_param_file();
                jobstring = [res.indoExe,' "',res.dataPath,'\',res.jobName,'.ipf"'];
                                
                % disp(['about to do: ',jobstring]);
                [~, result] = system(jobstring);
                delete([res.dataPath,'\',res.jobName,'.ipf']);
                
                res.indoOutput = result;
                if (~isempty(regexpi(result, '.*Repulse exception.*', 'start')))
                    [tok, ~] = regexpi(result, '.*Repulse exception(.*)', 'tokens');
                    res.indo_err_msg = num2str(cell2mat(tok{1}(1)));
                    res.indo_succeed = false;
                    res.norb = [];
                    res.aorbAtom = [];
                    res.aorbType = [];
                    res.nfilled = [];
                    res.hfE = [];
                    res.orbE = [];
                    res.orb = [];
                    res.nsci = [];
                    res.nscibasis = [];
                    res.esci = [];
                    res.r = [];
                    res.wfsci = [];
                    res.ehsci = [];
                else
                    res.indo_succeed = true;
                    res.readIndo();
                end
            end
        end % INDO constructor
        function res = get_osc(obj)
            res = zeros(1,obj.nsci);
            for i=1:obj.nsci
                res(1,i) = (obj.esci(1,i)-obj.esci(1,1)) * ...
                   ( obj.r(1,i,1)^2 + obj.r(1,i,2)^2 + obj.r(1,i,3)^2 );
            end
        end
        function res = get_rx(obj)
            res = obj.r(:,:,1);
        end
        function res = get_ry(obj)
            res = obj.r(:,:,2);
        end
        function res = get_rz(obj)
            res = obj.r(:,:,3);
        end
        function res = dipole(obj,istate,jstate)
            % returns a vector that is the dipole(istate=jstate)
            % or transition moment (istate ~= jstate)
            res = reshape(obj.r(istate,jstate,:),[3,1]);
        end
        function successful = write_param_file(obj)
            fid = fopen([obj.dataPath, obj.jobName, '.ipf'], 'w');
            
            if (fid == -1)
                throw(MException('Indo:Write_Param_File:IOError','Indo: Cannot create param file.'));
            end
            
            fprintf(fid, '%s\r\n', ['jobname = "', obj.dataPath, obj.jobName, '"']);

            if (strcmpi(obj.config.dm_guess, 'default'))
                fprintf(fid, '%s\r\n', 'dm_guess = default');
            else
                fprintf(fid, '%s\r\n', ['dm_guess = "', obj.config.dm_guess, '"']);
            end
            
            if (~isempty(obj.config.pot_file))
                fprintf(fid, '%s\r\n', ['pot_file = "', obj.config.pot_file, '"']);
            end
            
            fprintf(fid, '%s%i\r\n', 'charge = ', obj.config.charge);
            fprintf(fid, '%s%u\r\n', 'norbs = ', obj.config.norbs);
            fprintf(fid, '%s%u\r\n', 'nstates = ', obj.config.nstates);
            fprintf(fid, '%s%u\r\n', 'max_inner_iter = ', obj.config.max_inner_iter);
            fprintf(fid, '%s%u\r\n', 'max_iter = ', obj.config.max_iter);
            fprintf(fid, '%s%u\r\n', 'try_default_first = ', obj.config.try_default_first);
            fprintf(fid, '%s%u\r\n', 'output_dm = ', obj.config.output_dm);
            
            fprintf(fid, '%s%6.10f\r\n', 'efieldx = ', obj.config.field(1));
            fprintf(fid, '%s%6.10f\r\n', 'efieldy = ', obj.config.field(2));
            fprintf(fid, '%s%6.10f\r\n', 'efieldz = ', obj.config.field(3));
            
            fprintf(fid, '%s%3.10f\r\n', 'initial_shiftC = ', obj.config.initial_shiftc);
            fprintf(fid, '%s%3.10f\r\n', 'initial_shift_step = ', obj.config.initial_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'min_shift_step = ', obj.config.min_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'max_shift_step = ', obj.config.max_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'initial_second_shift_step = ', obj.config.initial_second_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'min_second_shift_step = ', obj.config.min_second_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'max_second_shift_step = ', obj.config.max_second_shift_step);
            fprintf(fid, '%s%3.10f\r\n', 'initial_eeint = ', obj.config.initial_eeint);
            fprintf(fid, '%s%3.10f\r\n', 'initial_eestep = ', obj.config.initial_eestep);
            fprintf(fid, '%s%3.10f\r\n', 'max_eestep = ', obj.config.max_eestep);
            fprintf(fid, '%s%3.10f\r\n', 'min_eestep = ', obj.config.min_eestep);
            
            fprintf(fid, '%s%3.10e\r\n', 'initial_conv = ', obj.config.initial_conv);
            fprintf(fid, '%s%3.10e\r\n', 'min_conv = ', obj.config.min_conv);
            
            fclose(fid);
            
            successful = 1;
        end
            
    end % methods
end % class