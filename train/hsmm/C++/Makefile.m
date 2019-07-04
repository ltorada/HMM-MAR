%
% Compile mex function
%
% Other file dependencies: FBImpl.h, Mathe.h, error.h, matrix.h, cube.h, consts.h
%
%


mex -c -O InitData.cpp
mex -c -O dCalc.cpp
mex -c -O FBImpl.cpp InitData.o dCalc.o -R2018a
mex -O FB.cpp FBImpl.o InitData.o dCalc.o -R2018a

%mex -c -O lib.cpp
%mex -O FB.cpp lib.o -R2018a




