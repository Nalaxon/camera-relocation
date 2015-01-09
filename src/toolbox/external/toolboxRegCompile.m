% Compiles all mex routines that are part of toolbox.
%
% If you get warnings about "omp.h", your compiler may lack default OpenMp
% support (for info about OpenMP see http://en.wikipedia.org/wiki/OpenMP).
% To compile without OpenMP alter code below to 'useOmp(:)=0'; and re-run
% toolboxCompile. Note that this will disable parallelization and make some
% routines (in particular training certain classifier) much slower.
%
% USAGE
%  toolboxCompile
%
% INPUTS
%
% OUTPUTS
%
% EXAMPLE
%
% See also
%
% Piotr's Image&Video Toolbox      Version 3.25
% Copyright 2013 Piotr Dollar.  [pdollar-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see external/bsd.txt]

% compile options including openmp support for C++ files
opts = {'-output'};
if(exist('OCTAVE_VERSION','builtin')), opts={'-o'}; end
if( ispc ), optsOmp={'OPTIMFLAGS="$OPTIMFLAGS','/openmp"'}; else
  optsOmp={'CXXFLAGS="\$CXXFLAGS','-fopenmp"'};
  optsOmp=[optsOmp,'LDFLAGS="\$LDFLAGS','-fopenmp"'];
end
optsOmp=[optsOmp '-DUSEOMP'];

% list of files (missing /private/ part of directory)
fs={'classify/forestRegFindThr.cpp', 'classify/forestFindThr.cpp',...
  'classify/forestInds.cpp'};
n=length(fs); useOmp=zeros(1,n); useOmp([6 9])=1;

% compile every funciton in turn (special case for dijkstra)
disp('Compiling Piotr''s Toolbox.......................');
rd=fileparts(mfilename('fullpath')); rd=rd(1:end-9); tic;
errmsg=' -> COMPILE FAILURE: ''%s'' %s\n';
for i=1:n
  try %#ok<ALIGN>
    [d,f1,e]=fileparts(fs{i}); f=[rd '/' d '/private/' f1];
    if(useOmp(i)), optsi=[optsOmp opts]; else optsi=opts; end
    fprintf(' -> %s\n',[f e]); mex([f e],optsi{:},[f '.' mexext]);
  catch err, fprintf(errmsg,[f1 e],err.message); end
end
try %#ok<ALIGN>
  d=[rd '/matlab/private/']; fprintf(' -> %s\n',[d 'dijkstra1.cpp']);
  mex([d 'fibheap.cpp'], [d 'dijkstra1.cpp'], '-largeArrayDims', ...
    opts{:}, [d 'dijkstra1.' mexext]);
catch err, fprintf(errmsg,[f1 e],err.message); end
disp('..................................Done Compiling'); toc;