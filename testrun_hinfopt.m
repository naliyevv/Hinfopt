%
% Test script to run a set of example systems with HINFOPT v1.0.
%
% REFERENCES:
%
% [1] N. Aliyev, P. Benner, E. Mengi, and M. Voigt.
%     A subspace framework for H-infinity norm minimization. Submitted for
%     publication, 2019.
%
% AUTHOR:
%
% Matthias Voigt, Technische Universitaet Berlin, Institut fuer Mathematik,
% Berlin, Germany.
%
% 07/12/2018.
%
% REVISIONS:
%
% Matthias Voigt, 03/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HINFOPT 1.0 Copyright (C) 2018 Nicat Aliyev, Emre Mengi, Matthias Voigt 
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
delete diary.dia
diary( 'diary.dia' );
fprintf( '===========================================================\n' );
fprintf( 'TEST SERIES FOR hinfopt.m\n' );
fprintf( '===========================================================\n' );
fprintf( 'Date: %s.\n', date );
fprintf( 'Version: %s.\n', version );
fprintf( 'Computer: %s.\n', computer );
addpath( genpath( 'linorm_subsp_v1.3' ) );
addpath( genpath( 'fg_functions' ) );
addpath( genpath( 'GRANSO' ) );
%% EXAMPLE 1a
clear
load testexamples/anemometer_1p
options.prtlevel = 2;
options.initialParams = [ linspace( 0, 1, 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 1a: anemometer_1p using GRANSO\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 1b
clear
load testexamples/anemometer_1p
options.prtlevel = 2;
options.initialParams = [ linspace( 0, 1, 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 1;
options.gamma = -10000;
fprintf( '===========================================================\n' );
fprintf( 'Example 1b: anemometer_1p using eigopt\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 2a
clear
load testexamples/anemometer_3p
options.prtlevel = 2;
options.initialParams = [ linspace( 0, 1, 10 ); linspace( 1, 2, 10 ); ...
    linspace( 0.1, 2, 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 2a: anemometer_3p using GRANSO\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 2b
clear
load testexamples/anemometer_3p
options.prtlevel = 2;
options.initialParams = [ linspace( 0, 1, 10 ); linspace( 1, 2, 10 ); ...
    linspace( 0.1, 2, 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 1;
options.gamma = -10000;
fprintf( '===========================================================\n' );
fprintf( 'Example 2b: anemometer_3p using eigopt\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 3a
clear
load testexamples/SECM
options.prtlevel = 2;
options.initialParams = [ linspace( 1, exp( 2 ), 10 ); ...
    linspace( 1, exp( 2 ), 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 3a: SECM using GRANSO\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 3b
clear
load testexamples/SECM
options.prtlevel = 2;
options.initialParams = [ linspace( 1, exp( 2 ), 10 ); ...
    linspace( 1, exp( 2 ), 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 1;
options.gamma = -10000;
fprintf( '===========================================================\n' );
fprintf( 'Example 3b: SECM using eigopt\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 4a
clear
load testexamples/sin_membrane
options.prtlevel = 2;
options.initialParams = [ linspace( 2, 5, 10 ); ...
    linspace( 400, 750, 10 ); linspace( 3000, 3200, 10 ); ...
    linspace( 10, 12, 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 4a: sin_membrane using GRANSO\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 4b
clear
load testexamples/sin_membrane
options.prtlevel = 2;
options.initialParams = [ linspace( 2, 5, 10 ); ...
    linspace( 400, 750, 10 ); linspace( 3000, 3200, 10 ); ...
    linspace( 10, 12, 10 ); linspace( 0, 10, 10 ) ];
options.initialFreqs = linspace( 0, 10, 10 );
options.optmeth = 1;
options.gamma = -10000;
fprintf( '===========================================================\n' );
fprintf( 'Example 4b: sin_membrane using eigopt\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 5a
clear
load testexamples/synthetic
options.prtlevel = 2;
options.initialParams = [ linspace( 0.02, 1, 10 ); ...
    linspace( 0, 100, 10 ) ];
options.initialFreqs = linspace( 0, 100, 10 );
options.optmeth = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 5a: synthetic using GRANSO\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 5b
clear
load testexamples/synthetic
options.prtlevel = 2;
options.initialParams = [ linspace( 0.02, 1, 10 ); ...
    linspace( 0, 100, 10 ) ];
options.initialFreqs = linspace( 0, 100, 10 );
options.optmeth = 1;
options.gamma = -10000;
fprintf( '===========================================================\n' );
fprintf( 'Example 5b: synthetic using eigopt\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 6a
clear
load testexamples/T2DAL_BCI
options.prtlevel = 2;
options.initialParams = [ linspace( 1, 10000, 10 ); ...
    linspace( 1, 10000, 10 ); linspace( 1, 10000, 10 ); ...
    linspace( 0, 100, 10 ) ];
options.initialFreqs = linspace( 0, 100, 10 );
options.optmeth = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 6a: T2DAL_BCI using GRANSO\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%% EXAMPLE 6b
clear
load testexamples/T2DAL_BCI
options.prtlevel = 2;
options.initialParams = [ linspace( 1, 10000, 10 ); ...
    linspace( 1, 10000, 10 ); linspace( 1, 10000, 10 ); ...
    linspace( 0, 100, 10 ) ];
options.initialFreqs = linspace( 0, 100, 10 );
options.optmeth = 1;
options.gamma = -10000;
fprintf( '===========================================================\n' );
fprintf( 'Example 6b: T2DAL_BCI using eigopt\n' );
fprintf( '===========================================================\n' );
[ hinfnmin, mumin, info ] = hinfopt( sys, options );
%
diary off