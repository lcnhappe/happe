function res = hlp_homedir
% Get the user's home directory
res = char(java.lang.System.getProperty('user.home'));
