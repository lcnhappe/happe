function [inargs,outargs,h1] = hlp_inoutargs(filename,funcname)
% Retrieve the names of input/output arguments of the given m file.
% [In-Args,Out-Args,H1-Line] = hlp_inoutargs(Filename,Function-Name)
%
% In:
%   Filename        : full path to an m-file
%   Function-Name   : optional, name of a (possibly nested) function in the file, for which the args
%                     shall be retrieved
%
% Out:
%   In-Args  : cell array of input argument names for the function
%   Out-Args : cell array of output argument names for the function
%   H1-Line  : first line of help text
%
% Notes:
%   For heavy use, utl_fileversion is to be preferred (which caches the result of this function properly).
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-14

h1 = [];

if ~exist('funcname','var')
    funcname = []; end

if ~exist(filename,'file')
    error(['there is no file named ' filename]); end
f = fopen(filename,'r');
if f == -1
    error(['Cannot open the file named ' filename]); end
try
    while 1
        % read the file line by line
        l = fgetl(f);
        if ~ischar(l)
            error(['definition of the function ' funcname ' not found in ' filename]); end
        % remove commented part
        commpos = strfind(l,'%');
        if ~isempty(commpos)
            l = l(1:commpos(1)-1); end
        % see if it contains a function declaration
        funcpos = strfind(l,'function');
        if ~isempty(funcpos) && isempty(strtrim(l(1:funcpos-1)))            
            % we have the declaration line, find out where the input/output args are
            ostart = strfind(l,'[')+1;
            oend = strfind(l,']')-1;
            if isempty(ostart) || isempty(oend)
                ostart = funcpos+length('function');
                oend = strfind(l, '=')-1;
            end
            istart = strfind(l,'(')+1;
            iend = strfind(l,')')-1;
            % and read them in (note: both ranges could be empty)
            if ~(isempty(istart) || isempty(iend))
                inargs = hlp_split(l(istart:iend),', '); 
            else
                inargs = {};
            end
            if ~(isempty(ostart) || isempty(oend))
                outargs = hlp_split(l(ostart:oend),', '); 
            else
                outargs = {};
            end
            inargs = inargs(~cellfun(@isempty,inargs));
            outargs = outargs(~cellfun(@isempty,outargs));                
            if ~isempty(funcname)
                % and also make sure that we have the correct function name
                if ~isempty(oend)
                    nstart = strfind(l,'=')+1;
                else
                    nstart = ostart;
                end
                if ~isempty(istart)
                    nend = istart-2;
                else
                    nend = length(l);
                end                
                if strcmp(strtrim(l(nstart:nend)), funcname)
                    % match: we're done
                    break; 
                end
            else
                % we're done
                break;
            end
        end
    end
    try
        % check if we find a commented line... this would be the h1 line
        l = fgetl(f);
        comm = find(l=='%',1);
        if ~isempty(comm)
            h1 = strtrim(l(comm+1:end)); end
    catch
    end
    fclose(f);
catch e
    if f ~= -1
        fclose(f); end
    rethrow(e);
end
