function params = arg_guidialog(func,varargin)
% Create an input dialog that displays input fields for a Function and Parameters.
% Parameters = arg_guidialog(Function, Options...)
%
% The Parameters that are passed to the function can be used to override some of its defaults. 
% The function must declare its arguments via arg_define. In addition, only a Subset of the function's specified arguments can be displayed.
%
% In:
%   Function : the function for which to display arguments
%
%   Options... : optional name-value pairs; possible names are:
%                 'Parameters' : cell array of parameters to the Function to override some of its defaults.
%
%                 'Subset' : Cell array of argument names to which the dialog shall be restricted; these arguments may contain . notation to index
%                            into arg_sub and the selected branch(es) of arg_subswitch/arg_subtoggle specifiers.
%                            Empty cells show up in the dialog as empty rows.
%
%                 'Title' : title of the dialog (by default: functionname())  
%
% Out:
%   Parameters : a struct that is a valid input to the Function.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-10-24

% parse arguments...
hlp_varargin2struct(varargin,{'params','Parameters'},{}, {'subset','Subset'},{}, {'dialogtitle','title','Title'}, [char(func) '()'], {'buttons','Buttons'},[]);

% obtain the argument specification for the function
rawspec = arg_report('rich', func, params); %#ok<*NODEF>
        

% extract a list of sub arguments...
[spec,subset] = obtain_items(rawspec,subset);

% create an inputgui() dialog...
geometry = repmat({[0.6 0.35]},1,length(spec)+length(buttons)/2);
geomvert = ones(1,length(spec)+length(buttons)/2);
% turn the spec into a UI list...
uilist = {};
for k = 1:length(spec)
     
    s = spec{k};
    if isempty(s)
        uilist(end+1:end+2) = {{} {}};
    else        
        if isempty(s.help)
            error(['Cannot display the argument ' subset{k} ' because it contains no description.']);
        else
            
            tag = subset{k};
            
            uilist{end+1} = {'Style','text', 'string',s.help{1}, 'fontweight','bold'};
            % depending on the type, we introduce different types of input widgets here...
            if iscell(s.range) && strcmp(s.type,'char')
                % string popup menu
                uilist{end+1} = {'Style','popupmenu', 'string',s.range,'value',find(strcmp(s.value,s.range)),'tag',tag};
            elseif strcmp(s.type,'logical')
                if length(s.range)>1
                    % multiselect
                    uilist{end+1} = {'Style','listbox', 'string',s.range, 'value',find(strcmp(s.value,s.range)),'tag',tag,'min',1,'max',100000}; 
                    geomvert(k) = min(3.5,length(s.range));
                else
                    % checkbox
                    uilist{end+1} = {'Style','checkbox', 'string','(set)', 'value',double(s.value),'tag',tag};
                end
            elseif strcmp(s.type,'char')
                % string edit
                uilist{end+1} = {'Style','edit', 'string', s.value,'tag',tag};
            else
                % expression edit
                if isinteger(s.value)
                    s.value = double(s.value); end
                uilist{end+1} = {'Style','edit', 'string', hlp_tostring(s.value),'tag',tag};
            end
            % append the tooltip string
            if length(s.help) > 1
                uilist{end} = [uilist{end} 'tooltipstring', regexprep(s.help{2},'\.\s+(?=[A-Z])','.\n')]; end
        end
    end

    if ~isempty(buttons) && k==buttons{1}
        % render a command button
        uilist(end+1:end+2) = {{} buttons{2}};
        buttons(1:2) = [];
    end
end

% invoke the GUI, obtaining a list of output values...
[outs,dummy,okpressed] = inputgui('geometry',geometry, 'uilist',uilist,'helpcom','disp(''coming soon...'')', 'title',dialogtitle,'geomvert',geomvert);

if ~isempty(okpressed)

    % remove blanks from the spec
    spec = spec(~cellfun('isempty',spec));
    subset = subset(~cellfun('isempty',subset));
    % turn the raw specification into a parameter struct (a non-direct one, since we will mess with it)
    params = arg_tovals(rawspec,false);

    % for each parameter produced by the GUI...
    for k = 1:length(outs)
        s = spec{k}; % current specifier
        v = outs{k}; % current raw value
        % do type conversion according to spec
        if iscell(s.range) && strcmp(s.type,'char')
            v = s.range{v};
        elseif strcmp(s.type,'expression')
            v = eval(v);
        elseif strcmp(s.type,'logical')
            if length(s.range)>1
                v = s.range(v);
            else
                v = logical(v);
            end
        elseif strcmp(s.type,'char')
            % noting to do
        else
            if ~isempty(v)
               v = eval(v); % convert back to numeric (or object, or cell) value
            end
        end
        % assign the converted value to params struct...
        params = assign(params,subset{k},v);
        
    end
    
    % now send the result through the function to check for errors and obtain a values structure...
    params = arg_report('rich',func,{params});
    params = arg_tovals(params,false);

else
    params = [];
end

% obtain a cell array of spec entries by name from the given specification
function [items,ids] = obtain_items(rawspec,requested,prefix)
if ~exist('prefix','var')
    prefix = ''; end
items = {}; 
ids = {};
% determine what subset of (possibly nested) items is requested
if isempty(requested)
    % look for a special argument/property arg_dialogsel, which defines the standard dialog representation for the given specification
    dialog_sel = find(cellfun(@(x)any(strcmp(x,'arg_dialogsel')),{rawspec.names}));
    if ~isempty(dialog_sel)
        requested = rawspec(dialog_sel).value; end
end
if isempty(requested)
    % empty means that all items are requested
    for k=1:length(rawspec)
        items{k} = rawspec(k);
        ids{k} = [prefix rawspec(k).names{1}];
    end
else
    % otherwise we need to obtain those items
    for k=1:length(requested)
        if ~isempty(requested{k})
            try
                items{k} = obtain(rawspec,requested{k});
            catch
                error(['The specified identifier (' prefix requested{k} ') could not be found in the function''s declared arguments.']);
            end
            ids{k} = [prefix requested{k}];
        end
    end    
end
% splice items that have children (recursively) into this list
for k = length(items):-1:1
    % special case: switch arguments are not spliced, but instead the argument that defines the option popupmenu will be retained
    if ~isempty(items{k}) && ~isempty(items{k}.children) && (~iscellstr(items{k}.range) || isempty(requested))
        [subitems, subids] = obtain_items(items{k}.children,{},[ids{k} '.']);
        if ~isempty(subitems)
            % and introduce blank rows around them
            items = [items(1:k-1) {{}} subitems {{}} items(k+1:end)];
            ids = [ids(1:k-1) {{}} subids {{}} ids(k+1:end)];
        end
    end
end

% remove items that cannot be displayed
retain = cellfun(@(x)isempty(x)||x.displayable,items);
items = items(retain);
ids = ids(retain);

% remove double blank rows
empties = cellfun('isempty',ids);
items(empties(1:end-1) & empties(2:end)) = [];
ids(empties(1:end-1) & empties(2:end)) = [];
 


% obtain a spec entry by name from the given specification
function item = obtain(rawspec,identifier)
% parse the . notation
dot = find(identifier=='.',1);
if ~isempty(dot)
    [head,rest] = deal(identifier(1:dot-1), identifier(dot+1:end));
else
    head = identifier;
    rest = [];
end
% search for the identifier at this level
names = {rawspec.names};
for k=1:length(names)
    if any(strcmp(names{k},head))
        % found a match!
        if isempty(rest)
            % return it
            item = rawspec(k);
        else
            % obtain the rest of the identifier
            item = obtain(rawspec(k).children,rest);
        end
        return;
    end
end
error(['The given identifier (' head ') was not found among the function''s declared arguments.']);



% assign a field with dot notation in a struct
function s = assign(s,id,v)
% parse the . notation
dot = find(id=='.',1);
if ~isempty(dot)
    [head,rest] = deal(id(1:dot-1), id(dot+1:end));
    if ~isfield(s,head)
        s.(head) = struct(); end
    s.(head) = assign(s.(head),rest,v);
else
    s.(id) = v;
end
