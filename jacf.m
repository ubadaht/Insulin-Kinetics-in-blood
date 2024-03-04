function [J] = jacf(t,X)

global kin k12 k21 k13 k31;


J = [-k13-k12 k21 k31 kin;
k12 -k21 0 0;
k13 0 -k31 0;
0 0 0 kin];


end