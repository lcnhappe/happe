function hash = hlp_cryptohash(data,fromfile)
% Compute an MD5 hash of a file, string or generic data structure.
% Hash = hlp_cryptohash(Data,FromFile)
%
% In:
%   Data : data to be hashed; can be a filename, a string, or any other MATLAB data structure.
%
%   FromFile : if true, data is interpreted as a file name (default: false)
%
% Out:
%   Hash : MD5 hash (decimal, lowercase) of the Data
%
% Examples:
%   % calculate an md5 hash of a file
%   hlp_cryptohash('myscript.m',true);
%
%   % calculate an md5 hash of a data structure
%   hlp_cryptohash(lastmodel);
%
% See also:
%   hlp_fingerprint, hlp_serialize

%					      	     Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-10-10


if exist('fromfile','var') && fromfile
    % take data as a file name
    if ~ischar(data)
        error('To represent a file name, Data should be a string.'); end
    if ~exist(data,'file')
        error('The file %s does not exist.',data); end
    f = fopen(data,'r');
    try
        data = fread(f,Inf);
        fclose(f);
    catch e
        try 
            fclose(f);
        catch,end
        rethrow(e);
    end
else
    % take data literally
    if ~ischar(data)
        data = hlp_serialize(data); end
end

% use Java to hash the data (idea from Michael Kleder)
hasher = java.security.MessageDigest.getInstance('MD5');
hasher.update(uint8(data));
hash = dec2hex(typecast(hasher.digest,'uint8'),2)';
hash = lower(hash(:)');
