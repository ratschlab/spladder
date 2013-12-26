function [engine, environment, basedir, mccdir, license_type] = determine_engine() ;
% [engine, environment, basedir, mccdir, license_type] = determine_engine() ;
% 
% returns either 'matlab' or 'octave' in variable engine
% and 'internal' or 'galaxy' in variable environment
% basedir is the directory of the matlab or octave instance
% mccdir is the directory of the matlab compiler (does not yet exist for octave)

global g_license_type ;

if isempty(g_license_type),
  g_license_type = 'academic' ;
end ;
license_type = g_license_type ;

if size(ver('Octave'), 1)
  engine='octave' ;
else
  engine='matlab' ;
end ;

environment='internal' ;
if isequal(whoami, 'galaxy'),
  environment='galaxy' ;
end ;

if isequal(engine, 'matlab') && isequal(environment,  'internal'),
  basedir = '/fml/ag-raetsch/share/software/matlab/' ;
elseif isequal(engine, 'matlab') && isequal(environment,  'galaxy'),
  basedir = '/home/galaxy/software/matlab-7.6/' ;
elseif isequal(engine, 'octave') && isequal(environment,  'internal'),
  basedir = '/fml/ag-raetsch/share/software/octave-3.0.3/' ;
elseif isequal(engine, 'octave') && isequal(environment,  'galaxy'),
  basedir = '/home/galaxy/software/octave' ;
end ;

if isequal(environment,  'internal'),
  mccdir = '/fml/ag-raetsch/share/software/matlab/' ;
elseif isequal(environment,  'galaxy'),
  mccdir = '/home/galaxy/software/matlab-7.6/' ;
end ;

return
