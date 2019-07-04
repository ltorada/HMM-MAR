%
% Compile mex function
%
% Requirements: previous run of Makefile.m
%

mex -c -O ViterbiImpl.cpp
mex -O Viterbi.cpp ViterbiImpl.o InitData.o dCalc.o -R2018a






