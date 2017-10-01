function result = para_dataflow(varargin)
% Generic Signal Processing -> Feature Extraction -> Machine Learning BCI framework.
% Result = para_dataflow(FilterSetup, FeatureExtractionArguments, FeatureExtractionSetup, MachineLearningSetup, DialogSetup, ForwardedParameters...)
%
% Most BCI paradigms are implemented as a sequence of three major stages: Signal Processing, Feature Extraction and Machine Learning (see also bci_train).
% The Signal Processing stage operates on time series (plus meta-data), represented as EEGLAB datasets, and may contain multiple sub-stages, which together
% form a filter graph. The data passed from node to node in this graph is either continuous or epoched EEGLAB datasets, and may contain rich annotations, such as 
% channel locations, ICA decompositions, DIPFIT models, etc. The nodes are filter components which are shared among many paradigms, and most of them 
% are found in filters/flt_* and dataset_ops/set_*. For almost all paradigms, a default order of these stages can be defined, because several filters 
% can be arbitrarily ordered w.r.t. each other (linear operators), and most other filters make sense only when executed before or after certain other stages.
% The default signal processing pipeline is implemented in flt_pipeline; its parameters allow to selectively enable processing stages. Paradims derived from
% para_dataflow usually set their own defaults for flt_pipeline (i.e., enable and configure various stages by default), which can be further modified by the user.
%
% The simplest filter components are stateless (such as a surface laplacian filter) and operate on each sample individually, while other filters are 
% stateful (and have a certain "memory" of past data), such as FIR and IIR filters. Some of the stateful filters are time-variant (such as signal
% standardization) and some of those are adaptive (such as ICA). Some adaptive filters may be unsupervised, and others may depend on the target variable.
% The majority of filters is causal (i.e. does not need data from future samples to compute the transformed version of some sample) and can therefore be
% applied online, while some filters are non-causal (e.g., the signal envelope and zero-phase FIR filter), and can only be used for offline analyses 
% (e.g. for neuroscience). Finally, most filters operate on continuous data, while some filters operate on epoched/segmented data (such as time window 
% selection or fourier transform). All of these filter components are written against a unified framework.
%
% Following the Signal Processing stage, most paradigms implement the Feature Extraction stage (especially those paradigms which do not implement 
% adaptive statistical signal processing), in which signals are transformed into sets of feature vectors. At this stage, signal meta-data is 
% largely stripped off. The feature extraction performed by some paradigms is non-adaptive (such as windowed means or log-variance), while it is
% adaptive (and usually supervised) for others (e.g., CSP). Feature vectors can be viewed as (labeled) points in a high-dimensional space, the 
% feature space, which serves as the representation on which the last stage, the machine learning, operates.
%
% The machine learning stage defines a standardized computational framework: it is practically always adaptive, and thus involves a 'learn' case
% and a 'predict' case. In the learning case, labeled sets of feature vectors are received, processed & analyzed, their distribution w.r.t. the 
% target variables (labels) is estimated, and a predictive model (or prediction function) incorporating these relations is generated. In the prediction case, 
% the previously computed predictive model is applied to individual feature vectors to predict their label/target value (or a probability distribution 
% over possible label/target values).
%
% Likewise, a paradigm can be applied to data in a 'learn' mode, in which data is received, feature-extraction is possibly adapted, and a predictive model 
% is computed (which is an arbitrary data structure that incorporates the state of all adaptive stages), and a 'predict' mode, in which a previously
% computed predictive model is used to map data to a prediction by sending it through all of the paradigm's stages. Finally, a paradigm has a 'preprocess'
% mode, in which all the signal processing steps take place. Separating the preprocessing from the other two stages leaves more control to the framework
% (bci_train, onl_predict), for example to control the granularity (block size) of data that is fed through the processing stages, buffering, 
% caching of intermediate results, partitioning of datasets (for cross-validation, nested cross-validation and other resampling techniques) and 
% various optimizations (such as common subexpression elimination and lazy evaluation). This functionality is invisible to the paradigms.
% 
% The function para_dataflow represents a small sub-framework for the convenient implementation of paradigms that adhere to this overall three-stage system. 
% Paradigms may implement their functionality by calling into para_dataflow, setting some of its parameters in order to customize its standard system.
% Therefore, para_dataflow exposes a set of named parameters for each of the three stages. For signal processing, it exposes all the parameters of 
% flt_pipeline, the default signal processing pipeline, allowing paradims to enable various pipeline stages without having to care about their relative
% order or efficient execution. For feature extraction, it exposes the 'featureextract' and 'featureadapt' parameters, which are function handles which
% implement the feature extraction and feature adaption (if any) step of this processing phase; the 'featurevote' parameter specifies whether the 'featureadapt'
% stages requires a voting procedure in cases where more than two classes are present in the data. Very few constraints are imposed on the type of 
% inputs, outputs and internal processing of these functions, or on the type of data that is passed through it (EEGLAB datasets or STUDY sets, for example). 
% For the machine learning stage, the 'learner' parameter of ml_train is exposed, allowing to specify one of the ml_train*** / ml_predict*** functions 
% for learning and prediction, respectively. See ml_train for more explanations of the options.
%
% Paradigms making use of para_dataflow typically pass all user-specified parameters down to para_dataflow (so that the user has maximum control with no
% interference from the paradigm), and set up their own characteristic parameterization of para_dataflow as defaults for these user parameters.
%
% In:
%       Parameters... : parameters of the paradigm:
%                       * 'op' : one of the modes in which the paradigm can be applied to data:
%                                * 'preprocess', to pre-process the InputSet according to the paradigm (and parameters)
%                                * 'learn', to learn a predictive model from a pre-processed InputSet (and parameters)
%                                * 'predict', to make predictions given a pre-processed InputSet and a predictive model
%
%                       * 'data' : some data (usually an EEGLAB dataset)
%
%                       if op == 'preprocess':
%                       * all parameters of flt_pipeline can be supplied for preprocessing; (defaults: no defaults are imposed by para_dataflow itself)
%
%                       if op == 'learn'
%                       * 'featureadapt': the adaption stage of the feature extraction; function_handle, 
%                                          receives the preprocessed data and returns a model of feature-extraction parameters (default: returns [])
%                       * 'featureextract': the feature extraction stage; function_handle, receives the preprocessed data and 
%                                           the model from the featureadapt stage, and returns a NxF array of N feature vectors (for F features)
%                                           (default: vectorizes each data epoch)
%                       * 'featurevote': true if the 'featureadapt' function supports only two classes, so that voting is necessary when three
%                                        or more classes are in the data (default: false)
%                       * 'learner': parameter of ml_train; defines the machine-learning step that is applied to the features (default: 'lda')
%
%                       if op == 'predict'
%                       * 'model': the predictive model (as produced during the learning)
%
% Out:
%       Result : depends on the op; 
%                * if 'preprocess', this is the preprocessed dataset
%                * if 'learn', this is the learned model
%                * if 'predict', this is the predictions produced by the model, given the data
%
% Notes:
%   Pre-processing is usually a purely symbolic operation (i.e. symbolic data processing steps are added to the data expression); the framework
%   evaluates this expression (and potentially transforms it prior to that, for example for cross-validation) and passes the evaluated expression 
%   back to the paradigm for 'learning' and/or 'prediction' modes.
%
% Name:
%   Data-Flow Framework
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-29

if length(varargin) > 1 && iscell(varargin{1})
    % the paradigm is being invoked via a user function (which sets flt_defaults, etc.)
    [flt_defaults,fex_declaration,fex_defaults,ml_defaults,dialog_default] = deal(varargin{1:5});
    varargin = varargin(6:end);
else
    % the paradigm is being invoked directly
    [flt_defaults,fex_declaration,fex_defaults,ml_defaults,dialog_default] = deal({});
end

cfg = arg_define(varargin, ...
    ... % core arguments for the paradigm framework (passed by the framework)
    arg_norep('op',[],{'preprocess','learn','predict'},'Operation to execute on the data. Preprocess the raw data, learn a predictive model, or predict outputs given a model.'), ...
    arg_norep('data',[],[],'Data to be processed by the paradigm.'), ...
    arg_norep('model',[],[],'Model according to which to predict.'), ...
    ... % signal processing arguments (sourced from flt_pipeline)
    arg_sub({'flt','SignalProcessing'},flt_defaults,@flt_pipeline,'Signal processing stages. These can be enabled, disabled and configured for the given paradigm. The result of this stage flows into the feature extraction stage','cat','Signal Processing'), ...
    ... % arguments for the feature-extraction plugins (passed by the user paradigms)
    arg_sub({'fex','FeatureExtraction'},{},fex_declaration,'Parameters for the feature-extraction stage.','cat','Feature Extraction'), ...
    ... % feature-extraction plugin definitions (passed by the user paradigms)
    arg_sub({'plugs','PluginFunctions'},fex_defaults,{ ...
        arg({'adapt','FeatureAdaptor'},@default_feature_adapt,[],'The adaption function of the feature extraction. Function_handle, receives the preprocessed data, an options struct (with feature-extraction), and returns a model (which may just re-represent options).'),...
        arg({'extract','FeatureExtractor'},@default_feature_extract,[],'The feature extraction function. Function_handle, receives the preprocessed data and the model from the featureadapt stage, and returns a NxF array of N feature vectors (for F features).'), ...
        arg({'vote','FeatureAdaptorNeedsVoting'},false,[],'Feature-adaption function requires voting. Only relevant if the data contains three or more classes.') ...
        },'The feature-extraction functions','cat','Feature Extraction'), ...
    ... % machine learning arguments (sourced from ml_train)
    arg_sub({'ml','MachineLearning'},ml_defaults,@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.','cat','Machine Learning'), ...
    ... % configuration dialog layout
    arg({'arg_dialogsel','ConfigLayout'},dialog_default,[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));


% map all of cfg's fields into the function's workspace, for convenience
arg_toworkspace(cfg,true);

switch op
    case 'preprocess'
        % apply default signal processing
        result = flt_pipeline('signal',data,flt);
        
    case 'learn'
        classes = unique(set_gettarget(data));
        if ~(plugs.vote && numel(classes) > 2)
            % learn a model
            [result.featuremodel,result.predictivemodel] = learn_model(data,cfg);
        else
            % binary stage and more than two classes: learn 1-vs-1 models for voting
            result.classes = classes;
            for i=1:length(classes)
                for j=i+1:length(classes)
                    [result.voting{i,j}.featuremodel,result.voting{i,j}.predictivemodel] = learn_model(exp_eval(set_picktrials(data,'rank',{i,j})),cfg); end
            end
        end
        result.plugs = plugs;
        
    case 'predict'
        if ~isfield(model,'voting')
            % predict given the extracted features and the model
            result = ml_predict(model.plugs.extract(data,model.featuremodel), model.predictivemodel);
        else
            % 1-vs-1 voting is necessary, construct the aggregate result
            trialcount = exp_eval(set_partition(data,[]));
            result = {'disc' , zeros(trialcount,length(model.classes)), model.classes};
            % vote, adding up the probabilities from each vote
            for i=1:length(model.classes)
                for j=i+1:length(model.classes)
                    outcome = ml_predict(model.plugs.extract(data,model.voting{i,j}.featuremodel), model.voting{i,j}.predictivemodel);
                    result{2}(:,[i j]) = result{2}(:,[i j]) + outcome{2};
                end
            end
            % renormalize probabilities
            result{2} = result{2} ./ repmat(sum(result{2},2),1,size(result{2},2));
        end
end

function [featuremodel,predictivemodel] = learn_model(data,cfg)
% adapt the feature extractor
switch nargin(cfg.plugs.adapt)
    case 1
        featuremodel = cfg.plugs.adapt(data);
    case 2
        featuremodel = cfg.plugs.adapt(data,cfg.fex);
    otherwise
        featuremodel = cfg.plugs.adapt(data,cfg.fex,cfg);
end
% extract features & learn a predictive model
predictivemodel = ml_train('data',{cfg.plugs.extract(data,featuremodel),set_gettarget(data)},'learner',cfg.ml.learner);


function mdl = default_feature_adapt(data,args)
mdl = args;

function data = default_feature_extract(data,mdl)
data = squeeze(reshape(data.data,[],1,size(data.data,3)))';