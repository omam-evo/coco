classdef Fitness < handle

    properties
        name
        %N
        f
        %y_hat
        f_stop
        f_eval
        f_best
        x_best
    end
    
    methods
        function obj = Fitness(name, params)
            obj.name = name;
            %obj.N = N;
            obj.f_eval = 0;
            obj.f_best = inf;
            if strcmpi(name, 'Sphere')
                obj.f = @(x) sphere(x);
            elseif strcmpi(name, 'Rastrigin')
                % [A,ALPHA] = [params(1),params(2)];
                obj.f = @(x) rastrigin(x, params(1), params(2));
            elseif strcmpi(name, 'Ellipsoid-H')
                obj.f = @(x) ellipsoidh(x); 
            elseif strcmpi(name, 'NoisySphere')
                obj.f = @(x) sphere(x) + randn(1,size(x,2)); 
            else
                error('Fitness not found')
            end
        end
        function [f_values] = get_f(obj, x)
            [N,LAMBDA] = size(x);
            f_values = obj.f(x);
            obj.f_eval = obj.f_eval + LAMBDA;
        end
        
    end
end

% vectorized fitness evaluation: size(x) = [N,LAMBDA]

function[f] = sphere(x)
  f = vecnorm(x, 2, 1).^2;
end

function [f] = rastrigin(x, a, alpha)
  f = sum(a - a * cos(alpha * x) + x.^2, 1);
end

function [f] = ellipsoidh(x)
  N = size(x,1);
  f = 10.^(6*((0:1:N-1))/(N-1)) * (x.*x);
end

