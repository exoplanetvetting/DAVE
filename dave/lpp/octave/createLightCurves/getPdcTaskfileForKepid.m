function taskfile = getPdcTaskfileForKepid(datapath, filename, kepid, quarter, varargin)
%Load the appropriate taskfile for a given kepid
%$Id: getPdcTaskfileForKepid.m 50155 2013-01-17 00:53:14Z fmullall $
%$URL: svn+ssh://murzim.amn.nasa.gov/repo/so/trunk/Develop/jvc/taskfileIO/getPdcTaskfileForKepid.m $
%Fergal Mullally
%
%Inputs:
%dataPath:  	(string) Where is the data stored. e.g ~/proc
%filename:	(string) The taskfile of interest e.g pdc-outputs-0.mat
%kepid:         (int) which kepid to load the data for.
%quarter:       (int) Which quarter to load data for
%
%Option inputs:
%ksop:                  (int) Use this ksop for the quarter. If not defined, use
%                       the official ksop number for that quarter
%isFullyQualifiedDataPath: If true, look for the data in dataPath/taskFileDir
%                       If false, look in
%                       dataPath/quarter/pipeline_results/lc/mpe_true.
%                       Default: false
%Returns
%The structure stored in the approparte .mat file for the given kepid.
    ch = getChannelForKepid(kepid, quarter);
    taskfile = getPdcTaskfileForCh(datapath, filename, ch, quarter, varargin{:});
end
