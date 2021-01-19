function [ Amu, Bmu, Cmu, Emu ] = evalTF( sys, mu )
%
%% PURPOSE:
%
% To evaluate a parameter-dependent transfer function of the form 
%
%   H(s;μ) := C(μ)*(sE(μ) − A(μ))^(−1)*B(μ),                            (1)
%
% at the parameter μ, where
%
%   A(μ) = a_1(μ)*A_1 + ... + a_{κ_A}(μ)*A_{κ_A}, 
%   B(μ) = b_1(μ)*B_1 + ... + b_{κ_B}(μ)*B_{κ_B},                       (2)
%   C(μ) = c_1(μ)*C_1 + ... + c_{κ_C}(μ)*C_{κ_C}, 
%   E(μ) = e_1(μ)*E_1 + ... + e_{κ_E}(μ)*E_{κ_E},
%
% and a_j, b_j, c_j, e_j are scalar-valued real-analytic functions and A_j,
% B_j, C_j, and E_j are fixed matrices. 
%
% ARGUMENTS
%
% Inputs:
%
% sys     : a struct containg the structural information of the function\
%           H(s;μ) specified as follows:
%    sys.A, sys.B, sys.C,  
%    sys.E                : arrays containing the matrices A_j, B_j, C_j,
%                           and E_j in (2), respectively. If κ_X is larger 
%                           than 1, i.e., the number of summands within the
%                           function X(μ) is larger than 1 and multiple 
%                           constant matrices X_1 to X_{κ_X} must be given,
%                           they must be stored in a cell array as 
%                           sys.X = { X_1, ..., X_{κ_X} }, where X is 
%                           either A, B, C, or E.  
%    sys.f.A, sys.f.B, 
%    sys.f.C, sys.f.E     : function handles specifying the functions 
%                           x_1, ..., x_{κ_X} stored in vector, i.e.,
%                           sys.f.X(k) = x_k, where x is either a, b, c, or
%                           e and X is either A, B, C, or E.
% mu      : the parameter value μ.
%
% Outputs:
%
% Amu, 
% Bmu,
% Cmu,
% Emu     : the realization matices A(μ), B(μ), C(μ), and E(μ),
%           respectively.
%
% REFERENCES:
%
% None.
%
% AUTHORS:
%
% Nicat Aliyev and Emre Mengi, Koc University, Istanbul, Turkey.
%
% Matthias Voigt, Technische Universitaet Berlin, Institut fuer Mathematik,
% Berlin, Germany.
%
% 07/11/2018.
%
% REVISIONS:
%
% -
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HINFOPT 1.0 Copyright (C) 2019 Nicat Aliyev, Emre Mengi, Matthias Voigt 
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% COMPUTATIONS.
%
% Evaluate the function handles at μ.
%
fA = sys.f.A( mu );
fB = sys.f.B( mu );
fC = sys.f.C( mu );
fE = sys.f.E( mu );
%
% Compute A(μ).
%
if iscell( sys.A )
    n = size( sys.A{ 1 }, 1 );
else
    n = size( sys.A, 1 );
end
%
Amu = sparse( n, n );
for k = 1:length( sys.A )
    Amu = Amu + fA( k )*sys.A{ k };
end
%
% Compute B(μ).
%
if iscell( sys.B )
    m = size( sys.B{ 1 }, 2 );
else
    m = size( sys.B, 2 );
end
%
Bmu = zeros( n, m );
for k = 1:length( sys.B )
    Bmu = Bmu + fB( k )*sys.B{ k };
end
%
% Compute C(μ).
%
if iscell( sys.C )
    p = size( sys.C{ 1 }, 1 );
else
    p = size( sys.C, 1 );
end
%
Cmu = zeros( p, n );
for k = 1:length( sys.C )
    Cmu = Cmu + fC( k )*sys.C{ k };
end
%
% Compute E(μ).
%
Emu = sparse( n, n );
for k = 1:length( sys.E )
    Emu = Emu + fE( k )*sys.E{ k };
end
%
return