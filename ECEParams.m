classdef ECEParams < handle
%--------------------------------------------------------------------------------------------------%
%   Class Name: ECEParams
%   Description: Holds experiment parameters for EnergyCalcExp objects
%   Written by: Christian Legaspi, Carnegie Mellon University
%   Date: June 10, 2011
%   Comments: If you don't know how to operate this...read the function headers
%       and comments on properties and hopefully that will help you out. Related classes are
%       EnergyCalcExp and ECEDataStuct (written by CML) and Indo (written by Dave Yaron).
%   Revisions: None
%--------------------------------------------------------------------------------------------------%
    
    properties (SetAccess = public)
        method;                 % char[]: method to be used for calculations, ex. 'AM1'
        phi;                    % double{[]}: mx1 cell array where each m cell contains a 1xn matrix of the
                                %       angle values to be substituted for PHI(m) in the template
        doIndo;                 % bool: TRUE runs INDO calculations
        Efield_mag;             % int/double: magnitude of the electric field applied in INDO in volts per Angstrom.
                                %       0 for no field. Negative for a field in the opposite direction.
        norbs_basis_set;        % unsigned int: number of orbitals to be included in the basis set
        nexcstates;             % unsigned int: number of states to be calculated by INDO
        
        % INDO HF Convergence Parameters
        initial_shiftc;         % Starting value for shift operator constant. Double.
        initial_shift_step;     % Starting value for the step by which to change the shift operator constant. Double.
        min_shift_step;         % Smallest shift step. Double.
        max_shift_step;         % Largest shift step. Double.
        initial_second_shift_step;  % After the shift operator constant reaches a value of 1, this is the step size that will be used to reach zero. Double.
        min_second_shift_step;  % Smallest shift step after 1 is reached. Double.
        max_second_shift_step;  % Largest shift step after 1 is reached. Double.
        initial_eeint;          % Starting value for the percentage of electron-electron interactions to include. Double.
        initial_eestep;         % Staring step value for the e-e interaction. Double.
        min_eestep;             % Smallest step size for e-e interaction percentage. Double.
        max_eestep;             % Largest step size for e-e interaction percentage. Double.
        initial_conv;           % Starting convergence value threshold. When the difference between two successive density matrices is below this value 10 times in a row without increasing, the matrix is considered converged for that set of parameters. Double.
        min_conv;               % When the shift operator constant is zero and e-e interaction is 1.0, this is the convergence threshold that must be met. Double.
        max_inner_iter;         % Maximum number of iterations to run for a given set of parameters. This number will be exceeded if the difference between the successive density matrices continues to decrease. Integer.
        max_iter;               % Maximum number of total iterations for the whole process. Integer.
        dm_guess;               % Full path to the binary file containing the initial guess at the density matrix we wish to use. If we want to use the default guess, use "default" as the value for this parameter. Paths must be enclosed in quotes if they contain a space.
        try_default_first;      % If you want to try the default parameters (shift operator constant = 0, e-e interaction = 1) to see if it can succeed without the use of the specialized parameters, specify 1 for this argument. Otherwise specify 0.
        output_dm;               % TRUE if you want to output the density matrix as a binary file. Boolean
        pot_file;               % Full path to the binary file containing the atom potentials and fields. If the path contains spaces, enclose it in quotes. Omit this parameter if you do not wish to use a potential file. This parameter takes precedence over (but does not overwrite) efieldx, efieldy and efieldz.
    end
    
    properties (SetAccess = private)
        Efield_vector;          % int/double[1,3]: Vector along which the electric field in INDO will be applied
        Efield_type;            % char: 'c' for custom, 'p' for parallel and 'n' for normal (normal to plane formed by atoms)
        Efield_description;     % char[]: Description of the vector.
    end
    
    methods                % Public functions that can be accessed by an instance of the class
%--------------------------------------------------------------------------------------------------%
% Function:     setEFVectorFromAtomNums
% Description:  Takes atom numbers and ampac data and creates an electric field vector along the
%               vector formed by two atoms or normal to the plane formed by three atoms. Sets the
%               normalized vector in the parameters.
% Input Args:
%       obj:    ECEParams calling object. For the sake of clarity, all functions contained
%               within this class will have 'obj' as their ECEParams calling object if they
%               aren't static functions. This parameter does not need to be included in the argument
%               list but is assigned when the function is called. 'obj' is passed by reference and
%               functions similar to the 'this' pointer in C++. All other input args are passed by
%               value, unless they are specified in the output args as well. This explanation
%               will be omitted in other functions listed below.
%       ampac: struct. Contains the data output by EnergyCalcExp.run_Ampac() or parseAmpac()
%       axis_params: unsigned int[1,2|3]. contains atom numbers to create E field vector. If length = 2,
%               vector will be created parallel to the vector between the two atom numbers. If length = 3,
%               vector will be created normal to the plane formed by the three atom numbers.
% Output:  
%       fvect: double/int[1,3]. Matrix containing the normalized electric field vector [x,y,z]
%--------------------------------------------------------------------------------------------------%
        function fvect = setEFVectorFromAtomNums(obj, ampac, axis_params)
            if (length(axis_params) == 2)
                EFMN1 = axis_params(1);
                EFMN2 = axis_params(2);
                EFMN3 = [];            
            elseif (length(axis_params) == 3)
                EFMN1 = axis_params(1);
                EFMN2 = axis_params(2);
                EFMN3 = axis_params(3);
            else
                throw(MException('ECEParams:setEFVectorFromAtomNums:Not_Enough_Params',...
                    'ECEParams: There are not enough atom numbers to generate a vector.'));
            end
            
            if (EFMN1 > size(ampac.r, 2) || EFMN1 < 1 || EFMN2 > size(ampac.r, 2) || EFMN2 < 1)
                throw(MException('ECEParams:setEFVectorFromAtomNums:Bad_Atom_Num',...
                    'ECEParams: Atom number specified does not exist in the specified Z-matrix'));
            end
            
            AB = ampac.r(:, EFMN1) - ampac.r(:, EFMN2);
            fvect = AB;

            if (~isempty(EFMN3))
                if (EFMN3 > size(ampac.r, 2) || EFMN3 < 1)
                throw(MException('ECEParams:setEFVectorFromAtomNums:Bad_Atom_Num',...
                    'ECEParams: Atom number specified does not exist in the specified Z-matrix'));
                end
                
                AC = ampac.r(:, EFMN1) - ampac.r(:, EFMN3);
                fvect = cross(AB,AC);
                obj.Efield_description = ['Normal to plane formed by atoms ',...
                    num2str(EFMN1), ', ' , num2str(EFMN2), ' and ', num2str(EFMN3)];
                obj.Efield_type = 'n';
            else
                obj.Efield_description = ['Vector formed by atoms ',...
                    num2str(EFMN1), ' (head) and ' , num2str(EFMN2), ' (tail)'];
                obj.Efield_type = 'p';
            end
            
            fvect = fvect' / norm(fvect);
            
            obj.Efield_vector = fvect;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     setCustomEFVector
% Description:  Takes a custom electric field vector, normalizes and stores it.
% Input Args:
%       obj: ECEParams calling object.
%       vectin: int/double[1,3]. Electric field vector (x,y,z)
% Output:  
%       fvect: int/double[1,3]. Normalized electric field vector (x,y,z)
%--------------------------------------------------------------------------------------------------%
        function fvect = setCustomEFVector(obj, vectin)           
            if (length(vectin) ~= 3)
                throw(MException('ECEParams:setCustomEFVector:Bad_Vector_Size',...
                    'ECEParams: Custom vector size is not three-dimensional'));
            end
            
            if (norm(vectin) ~= 0)
                vectin = vectin / norm(vectin);
            end
            
            vectin = reshape(vectin, [1 3]);
            
            fvect = vectin;
            obj.Efield_vector = vectin;
            
            obj.Efield_description = ['Custom E field: (',...
                num2str(vectin(1)), ', ', num2str(vectin(2)), ...
                ', ', num2str(vectin(3)), ')'];
            obj.Efield_type = 'c';
        end
        
        function successful = read_param_template(obj, path)
            obj.set_default_hf_params();
            if (~exist(path,'file'))
                successful = 0;
                return;
            end
            
            param_text = fileread(path);
            param_text = textscan(param_text,'%s','delimiter','\n');
            param_text = param_text{1};
            param_text = cellfun(@(x)textscan(x, '%s', 'delimiter', '='), param_text);
            param_text = cellfun(@(x)strtrim(x), param_text, 'UniformOutput', false);
            
            successful = 1;
            fn = fieldnames(obj);
            do_not_read = {'jobname', 'charge', 'norbs', 'nstates', 'efieldx', 'efieldy', 'efieldz'};
            treat_as_str = {'dm_guess', 'pot_file'};
            
            for i = 1:length(param_text)
                param_text{i}{1} = lower(param_text{i}{1});
                if (~ismember(param_text{i}{1}, do_not_read) && ...
                        ismember(param_text{i}{1}, fn) && ...
                        ~ismember(param_text{i}{1}, treat_as_str))
                    eval(['obj.', param_text{i}{1}, ' = ', param_text{i}{2},';']);
                elseif (~ismember(param_text{i}{1}, do_not_read) && ...
                        ismember(param_text{i}{1}, fn) && ...
                        ismember(param_text{i}{1}, treat_as_str))
                    eval(['obj.', param_text{i}{1}, ' = ''', param_text{i}{2},''';']);
                elseif (~ismember(param_text{i}{1}, do_not_read) && ~ismember(param_text{i}{1}, fn))
                    successful = -1;
                end
            end
        end
        
        function successful = set_default_hf_params(obj)
            obj.initial_shiftc = 0.0;
            obj.initial_shift_step = 0.0;
            obj.min_shift_step = 0.0;
            obj.max_shift_step = 0.0;
            obj.initial_second_shift_step = 0.0;
            obj.min_second_shift_step = 0.0;
            obj.max_second_shift_step = 0.0;
            obj.initial_eeint = 1.0;
            obj.initial_eestep = 0.0;
            obj.min_eestep = 0.0;
            obj.max_eestep = 0.0;
            obj.initial_conv = 1e-10;
            obj.min_conv = 1e-10;
            obj.max_inner_iter = 300;
            obj.max_iter = 10000;
            obj.dm_guess = 'default';
            obj.try_default_first = false;
            obj.output_dm = false;
            obj.pot_file = [];
            
            successful = true;
        end
%--------------------------------------------------------------------------------------------------%
% Function:     Constructor
% Description:  Builds an ECEParams object. Can build empty object if no args specified.
% Input Args:
%       method: char[]. method to be used for calculations, ex. 'AM1'
%       phi: double{[]}}. mx1 cell array where each m cell contains a 1xn matrix of the
%               angle values to be substituted for PHI(m) in the template
%       doIndo: bool. TRUE runs INDO calculations
%       Efield_mag: int/double. magnitude of the electric field applied in INDO in volts per Angstrom.
%               0 for no field. Negative for a field in the opposite direction.
%       norbs: unsigned int. number of orbitals to be included in the basis set
%       nstates: unsigned int. number of states to be calculated by INDO
%       set_vector_type: char. 'c' for custom, 'p' for parallel and 'n' for normal. This and the following two arguments
%               can be omitted and the electric field vector can be set during EnergyCalcExp.run(), but if you
%               want to specify a vector now, you may do so. 'c' takes the contents of 'components' to be x,y,z components
%               of the electric field vector you wish to have. 'p' takes the two items in 'components' to be the
%               atom numbers between which to draw a directional vector for the electric field. 'n' takes three
%               items in 'components', makes a plane from them and then finds a vector normal to that plane.
%       components: unsigned int[1,2|3] or int/double[1,3]. Matrix containing the necessary parameters listed above.
%       ampac: struct. Ampac data struct provided from EnergyCalcExp.run_Ampac() or parseAmpac(). Only necessary if
%               set_vector_type is 'p' or 'n'.
% Output:  
%       obj: ECEParams. Newly constructed object.
%--------------------------------------------------------------------------------------------------%
        function obj = ECEParams(method, phi, doIndo, Efield_mag, norbs, nstates, param_file_template, set_vector_type, components, ampac)
            if (nargin > 1)
                if (~isempty(phi))
                    if (~iscell(phi))
                        phi = mat2cell(phi, ones(1, length(phi)));
                    end
                else
                    phi = {};
                end
    
                if (nargin > 3 && doIndo)
                    if (norbs < 1)
                        throw(MException('ECEParams:Constructor:Bad_Basis_Set',...
                            'ECEParams: You must have a positive integer for the basis set orbitals'));
                    end

                    if (nstates < 1)
                        throw(MException('ECEParams:Constructor:Bad_Num_Exc_States',...
                            'ECEParams: You must have a positive integer for the number of excited states'));
                    end
                    
                    if (nargin > 6)
                        param_res = obj.read_param_template(param_file_template);
                        if (param_res == 0)
                            throw(MException('ECEParams:Constructor:Bad_Param_Template',...
                                'ECEParams: Paramter template file could not be read'));
                        elseif (param_res == -1)
                            throw(MException('ECEParams:Constructor:Unrecognized_Fields',...
                                ['ECEParams: Your template parameter file contained unrecognized',...
                                ' fields. Any other matching fields were read in.']));
                        end
                    else
                        obj.set_default_hf_params();
                    end

                    if (nargin > 9)
                        if ((length(components) ~= 2 && strcmp(set_vector_type, 'p')) || ...
                            (length(components) ~= 3 && strcmp(set_vector_type, 'n')))

                            throw(MException('ECEParams:Constructor:Bad_Atom_Number',...
                                'ECEParams: Incorrect number of atom number parameters'));
                        end                

                        if (nargin == 8 && (strcmp(set_vector_type, 'p') || strcmp(set_vector_type, 'n')))
                            throw(MException('ECEParams:Constructor:No_Ampac_Data',...
                                'ECEParams: No Ampac data provided for position data'));
                        end

                        if (strcmp(set_vector_type, 'c')) % Custom vector
                            obj.setCustomEFVector(components);
                        elseif (strcmp(set_vector_type, 'p')) % Parallel vector to atom numbers
                            obj.setEFVectorFromAtomNums(ampac, components(1), components(2));
                        elseif (strcmp(set_vector_type, 'n')) % Normal vector to plane formed by atom numbers
                            obj.setEFVectorFromAtomNums(ampac, components(1), components(2), components(3));
                        end
                    else
                        obj.Efield_vector = [];
                        obj.Efield_type = [];
                        obj.Efield_description = 'Vector not set';
                    end

                    obj.Efield_mag = Efield_mag;
                    obj.norbs_basis_set = norbs;
                    obj.nexcstates = nstates;
                end

                obj.method = method;
                obj.phi = phi;
                obj.doIndo = doIndo;
            end
        end
    end
    
end

