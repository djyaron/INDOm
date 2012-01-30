classdef RingOptimizer < handle
%--------------------------------------------------------------------------------------------------%
%   Class Name: RingOptimizer
%   Description: Calculates the optimum internal angles of any polygon
%       given the lengths of the sides
%   Written by: Dave Yaron & Christian Legaspi, Carnegie Mellon University
%   Date: November 3, 2011
%   Comments: For inputs r(1) to r(n) where r(i) is connected with r(i+1) and r(i-1) (and r(n) connects to r(1))...
%       Inner angles are defined as:
%       theta(i) = angle between r(i) & r(i+1) inside the polygon and theta(n) is
%       the angle between r(n) & r(1).
%
%       Outer angles output from the run() function are defined as:
%       out(i) = positive angle (0 < out(i) < 360 deg) between r(i) and
%       the x axis where r(1) is defined to be on the x axis and therefore
%       out(1) is always 0 deg/rad.
%       
%   Revisions: None
%--------------------------------------------------------------------------------------------------%
    properties
        r;          % bond length inputs
        theta;      % the angles calculated
        rads;       % TRUE if output from run() and variable theta is in radians. Default: TRUE
    end
    methods
        function obj = RingOptimizer(r)     % r is a vector containing the lengths of the sides (must be >2 sides)
            if (length(r) < 3)
                throw(MException('RingOptimizer:NeedThreeSides','RingOptimizer: There must be >= 3 sides provided'));
            end
            obj.r = r;
            obj.rads = true;
        end
        function [in, out, lgm, exitflag] = run(obj, varargin)    % Run the non-linear system solver. Optional argument is string 'deg' for output in degrees
            % "in" is the internal angles, "out" is the external angles relative to the x axis,
            % "lgm" is the lagrangian multipliers
            % "exitflag", see documentation for fsolve(). 1 means a
            %       solution was found successfully
            if (~isempty(varargin) && strcmp(varargin{1}, 'deg'))
                obj.rads = false;
            end
                
            n = length(obj.r);
            options = optimset('Display','off');
            [res, ~, ~, exitflag] = fsolve(@obj.err, obj.guess, options);
            out = [0 res(1:n-1)];
            lgm = res(end-1:end);
            
            in = zeros(1,n);
            in(1) = pi - out(2);
            for i = 2:n-1
                in(i) = pi - out(i+1) + out(i);
            end
            in(end) = out(end) - pi;
            
            if (~obj.rads)
                out = rad2deg(out);
                in = rad2deg(in);
            end
            
            obj.theta = in;
        end
        function res = guess(obj)   % Produces the inital guess at the external angles (based on internal angles being those of a regular polygon)
            n = length(obj.r);
            res = zeros(1,n+1);
            for i=1:n-1
                res(i) = 2 *i * pi / n;
            end
            res(end-1) = 0.0;
            res(end) = 0.0;
        end
        function res = err(obj, par)    % Calculates the lagrangian function, returning the score values and the 2 LGM values
            n = length(obj.r);
            res = zeros(1,n+1);
            t = [0 par(1:n-1)];
            c = par(end-1:end);
            baseang = 2*pi / n;
            
            for i = 1:n-2
                res(i) = -2*(baseang - t(i+1) + t(i)) + 2*(baseang + t(i+1) - t(i+2)) - obj.r(i+1)*c(2)*cos(t(i+1))...
                    + obj.r(i+1)*c(1)*sin(t(i+1));
            end
            
            res(end-2) = -2*(baseang - t(end) + t(end-1)) + 2*(baseang - 2*pi + t(end)) - obj.r(end)*c(2)*cos(t(end))...
                + obj.r(end)*c(1)*sin(t(end));
            
            res(end-1) = obj.r(1);
            res(end) = 0;
            for i = 2:n
                res(end-1) = res(end-1) + obj.r(i)*cos(t(i));
                res(end) = res(end) + obj.r(i)*sin(t(i));
            end
        end
    end
end

