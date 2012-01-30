classdef ECEDataStruct < handle
%--------------------------------------------------------------------------------------------------%
%   Class Name: ECEDataStruct
%   Description: Holds experimental data for a single calculation within a EnergyCalcExp object
%   Written by: Christian Legaspi, Carnegie Mellon University
%   Date: June 10, 2011
%   Comments: If you don't know how to operate this...read the function headers
%       and comments on properties and hopefully that will help you out. Related classes are
%       EnergyCalcExp and ECEParams (written by CML) and Indo (written by Dave Yaron).
%   Revisions: None
%--------------------------------------------------------------------------------------------------%
    properties (SetAccess = private)
        ampac_file_path;    % char[]. Full path to the Ampac library
        indo_file_path;     % char[]. Full path to the INDO library
        ampac_hash;         % char[]. Hash key for the Ampac data in the Ampac library
        indo_hash;          % char[]. Hash key for the INDO data in the INDO library
        raw_ampac;          % struct. Contains the full Ampac data when it is loaded into memory
                            %       using the load_to_memory() function. Otherwise, it is empty.
        raw_indo;           % Indo object. Contains the full INDO data when it is loaded into memory
                            %       using the load_to_memory() function. Otherwise, it is empty.
    end
    
    properties
        % Identifiers and verifying success
        phiID;              % int/double[1,n]: Matrix containing the PHI(n) angle values represented by this calculation
        indo_succeed;       % bool. TRUE if the INDO calculation ran without any errors
        ampac_succeed;      % bool. TRUE if the Ampac calculation ran without any errors
        
        % Energy and transition data
        Ehf;                % double: Hartree-Fock energy from Ampac
        Eexc;               % double[1,n]: Normalized excited state energies (n states) from INDO/SCI
        Tint;               % double[1,n]: Normalized transition intensity for each excited state
        Tang;               % double[1,n]: Transition angle in degrees relative to the vector of the most intense transition
        indoOutput;         % char[]. Output from INDO program. Can be read to determine error if ~indo_succeed
    end
    
    methods           % Public functions which can be run by an instance of the class
%--------------------------------------------------------------------------------------------------%
% Function:     load_to_memory
% Description:  Loads or unloads raw Ampac or INDO data to/from the ECEDataStruct object. Useful
%               if you wish to access fields and data which are not stored natively.
% Input Args:
%       obj:    ECEDataStruct calling object. For the sake of clarity, all functions contained
%               within this class will have 'obj' as their ECEDataStruct calling object if they
%               aren't static functions. This parameter does not need to be included in the argument
%               list but is assigned when the function is called. 'obj' is passed by reference and
%               functions similar to the 'this' pointer in C++. All other input args are passed by
%               value, unless they are specified in the output args as well. This explanation
%               will be omitted in other functions listed below.
%       aori: char[]. 'ampac' specifies Ampac data and 'indo' specified INDO data for loading/unloading.
%       loadunload: char[]. 'load' specifies loading to memory. Anything else unloads the data.
% Output:  
%       successful: unsigned int. 1 means we loaded successfully. 2 means we unloaded successfully.
%--------------------------------------------------------------------------------------------------%
        function successful = load_to_memory(obj, aori, loadunload)
            if (strcmp(aori,'ampac'))
                if (isempty(obj.raw_ampac) && strcmp(loadunload, 'load'))
                    if (exist([obj.ampac_file_path, obj.ampac_hash, '.mat'], 'file'))
                        hashdat = load([obj.ampac_file_path, obj.ampac_hash, '.mat']);
                        obj.raw_ampac = parseAmpacFromText(hashdat.arc, hashdat.out);
                        successful = 1;
                    else
                        throw(MException('ECEDataStruct:load_to_memory:BadAmpacLibFile',...
                            'ECEDataStruct: Ampac library file does not exist.'));
                    end
                elseif (~strcmp(loadunload, 'load'))
                    obj.raw_ampac = [];
                    successful = 2;
                end
            elseif (strcmp(aori,'indo'))
                if (isempty(obj.indo_file_path))
                    throw(MException('ECEDataStruct:load_to_memory:No_INDO_Data',...
                        'ECEDataStruct: There is no INDO data in this experiment'));
                end
                if (isempty(obj.raw_indo) && strcmp(loadunload, 'load'))
                    if (exist([obj.indo_file_path, obj.indo_hash, '.mat'], 'file'))
                        hashdat = load([obj.indo_file_path, obj.indo_hash, '.mat']);
                        obj.raw_indo = hashdat.indo;
                        successful = 1;
                    else
                        throw(MException('ECEDataStruct:load_to_memory:BadINDOLibFile',...
                            'ECEDataStruct: INDO library file does not exist.'));
                    end
                elseif (~strcmp(loadunload, 'load'))
                    obj.raw_indo = [];
                    successful = 2;
                end
            else
                throw(MException('ECEDataStruct:load_to_memory:Bad_Argument',...
                    'ECEDataStruct: Please specify "ampac" or "indo"'));
            end            
        end
%--------------------------------------------------------------------------------------------------%
% Function:     get_ampac
% Description:  Retrieves fields from raw Ampac data. Loads it to memory and unloads it if
%               load_to_memory() has not been run.
% Input Args:
%       obj:    ECEDataStruct calling object.
%       property: char[]. Name of the field to load
%       varargin: cell array. The subscripts to retrieve from the specified field. See note in 
%               EnergyCalcExp about varargin if unsure how to proceed.
% Output:  
%       res: int/double[]. Multidimensional matrix containing the data requested.
%--------------------------------------------------------------------------------------------------%
        function res = get_ampac(obj, property, varargin)
            if (~obj.ampac_succeed)
                res = [];
                return;
            end
            if (isempty(obj.raw_ampac))
                obj.load_to_memory('ampac','load');
                
                if (isempty(find(strcmp(fieldnames(obj.raw_ampac),property), 1)))
                    throw(MException('ECEDataStruct:get_ampac:Invalid_Property',...
                        'ECEDataStruct: The property specified is not contained in the AMPAC data.'));
                end
                    
                mystruct = struct('type',{'.';'()'},'subs',{property; varargin});
                res = subsref(obj.raw_ampac, mystruct);

                obj.load_to_memory('ampac','unload');
            else
                if (isempty(find(strcmp(fieldnames(obj.raw_ampac),property), 1)))
                    throw(MException('ECEDataStruct:get_ampac:Invalid_Property',...
                        'ECEDataStruct: The property specified is not contained in the AMPAC data.'));
                end
                
                mystruct = struct('type',{'.';'()'},'subs',{property; varargin});
                res = subsref(obj.raw_ampac, mystruct);
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     get_indo
% Description:  Retrieves fields from raw INDO data. Loads it to memory and unloads it if
%               load_to_memory() has not been run.
% Input Args:
%       obj:    ECEDataStruct calling object.
%       property: char[]. Name of the field to load
%       varargin: cell array. The subscripts to retrieve from the specified field. See note in
%               EnergyCalcExp about varargin if unsure how to proceed.
% Output:  
%       res: int/double[]. Multidimensional matrix containing the data requested.
%--------------------------------------------------------------------------------------------------%
        function res = get_indo(obj, property, varargin)
            if (~obj.indo_succeed)
                res = [];
                return;
            end
            if (isempty(obj.raw_indo))
                obj.load_to_memory('indo','load');
                
                if (isempty(find(strcmp(fieldnames(obj.raw_indo),property), 1)))
                    throw(MException('ECEDataStruct:get_indo:Invalid_Property',...
                        'ECEDataStruct: The property specified is not contained in the INDO data.'));
                end
                
                mystruct = struct('type',{'.';'()'},'subs',{property; varargin});
                res = subsref(obj.raw_indo, mystruct);

                obj.load_to_memory('indo','unload');
            else
                if (isempty(find(strcmp(fieldnames(obj.raw_indo),property), 1)))
                    throw(MException('ECEDataStruct:get_indo:Invalid_Property',...
                        'ECEDataStruct: The property specified is not contained in the INDO data.'));
                end
                
                mystruct = struct('type',{'.';'()'},'subs',{property; varargin});
                res = subsref(obj.raw_indo, mystruct);
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     Constructor
% Description:  Creates an instance of the ECEDataStruct class. Can be built empty.
% Input Args:
%       ampac: struct. Output from EnergyCalcExp.run_Ampac(). parseAmpac() is missing necessary fields.
%       indo: Indo object. Output from Indo calculations.
%       phi_info: int/double[]. 1xn matrix containing PHI(n) angles run for this calculation
%       ampac_path: char[]. Path to the folder where the Ampac data is stored
%       a_hash: char[]. Hash key for the Ampac library
%       indo_path: char[]. Path to the folder where the INDO data is stored
%       i_hash: char[]. Hash key for the INDO library
% Output:  
%       obj: ECEDataStruct object. Constructed object.
%--------------------------------------------------------------------------------------------------%
        function obj = ECEDataStruct(ampac, indo, phi_info, ampac_path, a_hash, indo_path, i_hash)
            if (nargin > 0)
                obj.phiID = phi_info;
                obj.ampac_file_path = ampac_path;
                obj.ampac_hash = a_hash;
                obj.indo_file_path = indo_path;
                obj.indo_hash = i_hash;
                obj.raw_ampac = [];
                obj.raw_indo = [];
                
                if (ampac.ampac_succeed)
                    obj.Ehf = ampac.Hf;
                    obj.ampac_succeed = true;
                else
                    obj.ampac_succeed = false;
                    obj.Ehf = [];
                end

                if (~isempty(indo))
                    obj.indoOutput = indo.indoOutput;
                    
                    if (indo.indo_succeed)
                        obj.Eexc = indo.esci(1:end) - indo.esci(1);

                        for i = 1:length(obj.Eexc)
                            obj.Tint(i) = obj.Eexc(i) * (indo.r(1,i,1)^2 + indo.r(1,i,2)^2 + indo.r(1,i,3)^2);
                        end

                        [~, maxInt] = max(obj.Tint);

                        for i = 1:length(obj.Eexc)
                            obj.Tang(i) = acosd( dot(indo.r(:,i), indo.r(:,maxInt)) / (norm(indo.r(:,maxInt)) * norm(indo.r(:,i))));
                        end
                        
                        obj.indo_succeed = true;
                    else
                        obj.indo_succeed = false;
                        obj.Tint = [];
                        obj.Eexc = [];
                        obj.Tang = [];
                    end
                else
                    obj.Tint = [];
                    obj.Eexc = [];
                    obj.Tang = [];
                end             
            end
        end
%--------------------------------------------------------------------------------------------------%
% Function:     update_paths
% Description:  Updates the locations of the files which contain the raw data for these calculations
% Input Args:
%       obj:    ECEDataStruct calling object.
%       newampac: char[]. New path to the Ampac data. If [], the current path will not change
%       newindo: char[]. New parth to the INDO data. If [], current path will not change.
% Output:  
%       successful: bool. Useless output. TRUE means we got to the end.
%--------------------------------------------------------------------------------------------------%
        function successful = update_paths(obj, newampac, newindo)          
            if (~isempty(newampac))                            
                if (newampac(end) ~= '\')
                    newampac = [newampac, '\'];
                end
                
                if (exist([newampac, obj.ampac_hash, '.mat'], 'file'))
                    obj.ampac_file_path = newampac;
                else
                    throw(MException('ECEDataStruct:update_paths:BadAmpacLibFile',...
                        'ECEDataStruct: Ampac library file does not exist.'));
                end                
            end
            
            if (~isempty(newindo))                
                if (newindo(end) ~= '\')
                    newindo = [newindo, '\'];
                end
                
                if (exist([newindo, obj.indo_hash, '.mat'], 'file'))
                    obj.indo_file_path = newindo;
                else
                    throw(MException('ECEDataStruct:update_paths:BadINDOLibFile',...
                        'ECEDataStruct: INDO library file does not exist.'));
                end 
            end
            
            successful = 1;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     generate_ampac_file
% Description:  Generates out, vis, dat or arc files from library file
% Input Args:
%       obj:    ECEDataStruct calling object.
%       file_type: char[]. String specifying file to load. 'vis', 'out', 'dat' or 'arc'
%       varargin: cell array (char[]). Optional argument to specify output path. Default is to output
%               to the library directory.
% Output:  
%       successful: bool. Useless output. TRUE means we got to the end.
%--------------------------------------------------------------------------------------------------%
        function successful = generate_ampac_file(obj, file_type, varargin)
            if (isempty(varargin))
                output_path = obj.ampac_file_path;
            else
                output_path = varargin{1};
                if (output_path(end) ~= '\')
                    output_path(end+1) = '\';
                end
            end
            
            if (strcmp(file_type, 'dat') || strcmp(file_type, 'out') || strcmp(file_type, 'vis') || ...
                    strcmp(file_type, 'arc'))
                if (exist([obj.ampac_file_path, obj.ampac_hash, '.mat'], 'file'))
                    if (~exist(output_path, 'dir'))
                        mkdir(output_path);
                    end
                    hashdat = load([obj.ampac_file_path, obj.ampac_hash, '.mat']);
                    fid = fopen([output_path, obj.ampac_hash, '.', file_type], 'w');
                    fwrite(fid, hashdat.(deblank(file_type)), 'char');
                    fclose(fid);
                else
                    throw(MException('ECEDataStruct:generate_ampac_file:BadAmpacLibFile',...
                        'ECEDataStruct: Ampac library file does not exist.'));
                end
            else
                throw(MException('ECEDataStruct:generate_ampac_file:InvalidFileType',...
                    'ECEDataStuct: Bad file type specified for load.'));
            end
            
            successful = 1;
        end
    end
end

