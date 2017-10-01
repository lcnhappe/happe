function X = arg_setdirect(X,value)
% Recursively set the arg_direct flag in a data structure to the given value.
% [Result] = arg_setdirect(Data,Value)
%

if nargin < 2
    error('please supply a value'); end
if isempty(X), return; end

if iscell(X)
    for k=1:numel(X)
        X{k} = arg_setdirect(X{k},value);
        if ischar(X{k}) && strcmpi(X{k},'arg_direct')
            X{k+1} = value;
        end
    end
elseif isstruct(X)
    if length(X) > 1
        for k=1:numel(X)
            X(k) = arg_setdirect(X(k),value); end
    else
        for fn=fieldnames(X)'
            X.(fn{1}) = arg_setdirect(X.(fn{1}),value); 
            if strcmp(fn{1},'arg_direct')
                X.(fn{1}) = value; end
        end
    end
end
