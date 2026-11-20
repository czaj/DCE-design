% RUN_TESTS  All tests of the DCE-design toolbox.
%   matlab -batch "cd('<toolbox>\tests'); run_tests"
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
dce_selftest;                 % generation, evaluation, Ngene translation, export
run(fullfile(fileparts(here), 'test_dce_block.m'));   % blocking
