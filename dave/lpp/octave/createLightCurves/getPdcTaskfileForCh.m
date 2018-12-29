function taskFile = getPdcTaskfileForCh(datapath, filename, chOrModOut, quarter, varargin)
%Return the appropriate taskfile for a given channel and quarter
%
%Inputs:
%dataPath:  	(string) Where is the data stored. e.g ~/proc
%filename:      (string) e.g pdc-inputs-0.mat 
%chOrModOut     Channel or modout. Integer is treated as a channel, mod-out
%               can be given as a decimal (e.g 3.1) or an array [3 1] 
%quarter:       (int) Which quarter to load data for
%
%Option inputs:
%ksop:                  (int) Use this ksop for the quarter. If not defined, use
%                       the official ksop number for that quarter
%isFullyQualifiedDataPath: If true, look for the data in dataPath/taskFileDir
%                       If false, look in
%                       dataPath/quarter/pipeline_results/lc/mpe_true.
%                       Default: false
%Returns:
%The structure stored in the .mat files called filename for the given channel

    taskDir = getPdcTaskdirForCh(datapath, chOrModOut, quarter, varargin{:});
    taskFilename = fullfile(taskDir, filename);
    
    tic;
    msg = sprintf('Loading %s', taskFilename);
    display(msg);
    retObj = load(taskFilename);
    toc;

    %Inelegant way of getting what we want from a load call
    namesCell = fieldnames(retObj);
    taskFile = retObj.(namesCell{1});
    
end
