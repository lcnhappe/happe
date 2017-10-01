% turn off advanced expression facilities and load a dataset
exp_eval(Set(@disable_expressions,1))
eeg = io_loadset('data:/christian/Demo/imag.vhdr');


% run a filter on some data with some arguments using positional syntax (rest using defaults)
eeg = flt_fir(eeg, [], [5 10 25 30], 'bandpass');
% run a filter on some data with name-value pair syntax, leaving out a few parameters
eeg = flt_fir('signal',eeg, 'fspec',[5 10 25 30], 'fmode','bandpass', 'ftype','minimum-phase');
% run a filter on some data with name-value pair syntax, but using alternative (long) property names
eeg = flt_fir('signal',eeg, 'Frequencies',[5 10 25 30], 'Mode','bandpass', 'Type','minimum-phase');
% run a filter on some data with name-value pair syntax, and using a struct for good measure
eeg = flt_fir('signal',eeg, 'Frequencies',[5 10 25 30], struct('Mode','bandpass', 'Type','minimum-phase'), 'passripple',-20,'stopripple',-30);

% edit parameters for flt_fir, starting with some parameters
x=arg_guidialog(@flt_fir,'params',{'fspec',[5 10 25 30]})
% run the filter on some data using edited parameters, using the name-value-pair/struct syntax
% (any name-value pair may be replaced by a struct unambiguously)
eeg = flt_fir('signal',eeg,x);
% run the filter on some data using edited parameters x, but with overridden filter type
eeg = flt_fir('signal',eeg,x,'ftype','zero-phase');


% edit subset of flt_ica arguments for its default config
x=arg_guidialog(@flt_ica,'subset',{'variant.max_iter','variant.scheduler',[],'variant.pdftype','variant.num_models'})

% edit all ICA arguments for a particular config
x=arg_guidialog(@flt_ica,'params',{'variant','infomax', 'clean',{'strong', 'windows',{'window_len',200}}})

% edit all ICA arguments for the default config... 
x=arg_guidialog(@flt_ica)

% show an editable panel for ICA arguments (note that, e.g., the variant is switchable)
% the relative geometry has not been refined yet
h = arg_guipanel('Function',@flt_ica);
uiwait(gcf);
x = arg_tovals(h.GetPropertySpecification())


% show an editable panel for a simple BCI paradigm (note that many filters can be enabled and configured)
h = arg_guipanel('Function',@para_csp)

% show an editable panel for a particular variant of a BCI paradigm
h = arg_guipanel('Function',@para_speccsp,'Parameters',{'flt', {'iir',[7 30]}, 'ml',{'learner',{'svmperf','Loss','ROC'}}})

% show an editable panel for a custom paradigm (try replacing the learner by hkl (hierarchical kernel learning)
h = arg_guipanel('Function',@para_multimodel)


% obtain partially overridden defaults for some function (note the 
x = arg_report('vals',@ml_trainhkl,{'loss','squared', 'display',1})



% side feature: obtain the defaults set for a function 
% here one that has uncertain values for its hyper-parameters (here: cost & g)
x = arg_report('vals',@ml_trainsvmlight)
% form the cartesian product over parameter search ranges for these arguments, yielding a set of argument combinations to the function over which to search...
% (done automatically by the machine learning tools)
x = hlp_flattensearch(x)
x.parts{50}


% test
h = arg_guipanel('Function',@set_inject_events);
uiwait(gcf);
x = arg_tovals(h.GetPropertySpecification());x.segment

