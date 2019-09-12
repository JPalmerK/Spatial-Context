function [spc,mu]=whiten_vec(sp)

mu=base3x(sp);
spc=sp-mu;
