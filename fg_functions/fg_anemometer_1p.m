function [ f, g ] = fg_anemometer_1p( mu, sys )
%
%% PURPOSE:
%
% To evaluate the L-infinity norm of the (reduced) anemometer_1p example 
% and its gradient at the parameter μ.
% The transfer function for this example is of the form
%
%   H(s;μ) := C*(sE − A(μ))^(−1)*B,                            
%
% where
%
%   A(μ) = A1 + μ*(A2 - A1),   μ in [0,1],
%
% ARGUMENTS
%
% Inputs:
%
% mu      : the parameter value μ.
% sys     : a struct containg the structual information of the function 
%           H(s;μ) specified as follows:
%    sys.A                : = {A1,A2}.
%    sys.B                : = B.
%    sys.C                : = C.
%    sys.E                : = E.             
%    sys.f.A              : = @(mu)[1-mu,mu].
%    sys.f.B              : = @(mu)1.
%    sys.f.C              : = @(mu)1.
%    sys.f.E              : = @(mu)1.
%    sys.method           : specifies the implementation of the 
%                           Boyd-Balakrishnan algorithm for computing the
%                           L-infinity norm as follows:
%                           = 0: use the FORTRAN routine AB13HD.f which is
%                                called by the gateway function generated
%                                by linorm_h.F;
%                           = 1: use the function 'norm' of the MATLAB
%                                Control System Toolbox (WARNING: This may
%                                not work, if the input descriptor matrix
%                                sys.E is singular at μ.)
%                           (default = 0).
%
% Outputs:
%
% f, g    : the L-infinity norm and its gradient at the point μ.
%
% REFERENCES:
%
% [1] N. Aliyev, P. Benner, E. Mengi, and M. Voigt.
%     A subspace framework for H-infinity norm minimization. Submitted for
%     publication, 2019.
%
% AUTHORS:
%
% Nicat Aliyev and Emre Mengi, Koc University, Istanbul, Turkey.
%
% Matthias Voigt, Technische Universitaet Berlin, Institut fuer Mathematik,
% Berlin, Germany.
%
% 07/12/2018.
%
% REVISIONS:
%
% Matthias Voigt, 03/2019.
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
%% SET OPTIONS.
%
% sys.method: specifies the implementation of the Boyd-Balakrishnan 
%             algorithm for computing the L-infinity norm
%
if isfield( sys, 'method' )
   method = sys.method;
else % default
   method = 0;
end
%
%% COMPUTATIONS.
%
A1 = sys.A{ 1 };
A2 = sys.A{ 2 };
%
[ Amu, Bmu, Cmu, Emu ] = evalTF( sys, mu );
%
% Evaluate the L-infinity norm.
%
if method == 0
    %
    % Use AB13HD.f/linorm_h.F.
    %
    [ f, z, infonorm ] = linorm_h( Amu, Emu, Bmu, Cmu, 0, 1, 0, 1, 0 );
    %
    if infonorm == 0
        if f( 2 ) == 0 || z( 2 ) == 0
            error( [ 'linorm_h: The reduced function is not in ', ...
                'L-infinity.' ] );
        else
            f = f( 1 )/f( 2 );
            z = z( 1 )/z( 2 );
        end
    else
        if infonorm == 7 || infonorm == 8
            error( [ 'linorm_h: The reduced function is not in ', ...
                'L-infinity.' ] );
        else
            %
            % AB13HD.f threw an error, the MATLAB function 'norm' is used
            % as a fallback option.
            %
            sysmu = dss( Amu, Bmu, Cmu, 0, Emu );
            [ f, z ] = norm( sysmu, Inf, 1e-6 );
        end
    end
else
    %
    % Use the MATLAB function 'norm'.
    %
    sysmu = dss( Amu, Bmu, Cmu, 0, Emu );
    [ f, z ] = norm( sysmu, Inf, 1e-6 );
end
%
% Evaluate the gradient using the analytic formula.
%
H = Cmu*( ( 1i*z*Emu - Amu )\Bmu );
H = full( H );
[ U, ~, V ] = svd( H );
u = U( :, 1 );
v = V( :, 1 );
%
dH = -Cmu/( 1i*z*Emu - Amu )*( A1 - A2 )*( ( 1i*z*Emu - Amu)\Bmu );
%
g = real( u'*dH*v );
%
return