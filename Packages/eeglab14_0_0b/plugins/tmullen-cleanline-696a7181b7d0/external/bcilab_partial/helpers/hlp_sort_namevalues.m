function res = hlp_sort_namevalues(varargin)
% Sort name-value pairs by name. 
% Sorted-NVPs = hlp_sort_namevalues(NVPs)
%
% Since the toolbox occasionally caches results of function calls,
% this is practical to uniformize arguments that lead to the same results, but look superficially different.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-03-28

names = varargin(1:2:end);
[names, I] = sort(names); I=([I*2-1; I*2]); I=I(:); %#ok<ASGLU>
res = varargin(I);
