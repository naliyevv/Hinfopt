function [ ci, cigrad ] = gransoConstr( mu, sys )
%
%% PURPOSE:
%
% To evaluate the inequality constraints needed for optimization using
% GRANSO in HINFOPT.
%
% The transfer function under consideration is defined parameters in a 
% closed box domain D. If the parameter μ is in D, then ci are elementwise
% nonpositive, otherwise they are not.
%
% ARGUMENTS
%
% Inputs:
%
% mu      : the parameter value μ.
% sys     : a struct containg the structural information of the as follows:       
% sys.bounds.lb, 
% sys.bounds.ub          : lower and upper bounds of the parameter box 
%                          domain D. Hereby, sys.bounds.lb(k) and
%                          sys.bounds.ub(k) are the bounds for the k-th
%                          parameter, respectively.
%
% Outputs:
%
% ci      : inequality constraints: if ci is elementwise nonpositive, then
%           the parameter value μ is in the domain D.
% cigrad  : gradient of the inequality constraints. 
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
% 07/16/2018.
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
%
%% COMPUTATIONS.
%
% Evaluate the constraint functions for the lower and upper bounds.
%
ci = [ sys.bounds.lb - mu; mu - sys.bounds.ub ];
%
% Evaluate the gradient of the constraint functions
%
lenlb = length( sys.bounds.lb );
cigrad = [ -eye( lenlb ), eye( lenlb ) ];
%
return