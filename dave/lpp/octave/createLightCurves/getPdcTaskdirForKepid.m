function taskdir = getPdcTaskdirForKepid(datapath, kepid, quarter, varargin)
%Get the directory name that stores PDC output for a given kepid
%$Id: getPdcTaskdirForKepid.m 50263 2013-01-28 22:40:04Z fmullall $
%$URL: svn+ssh://murzim.amn.nasa.gov/repo/so/trunk/Develop/jvc/taskfileIO/getPdcTaskdirForKepid.m $
%Fergal Mullally
%
%Inputs:
%dataPath:  (string) Where is the data stored. e.g ~/proc
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
%A string for the path to the .mat files

    ch = getChannelForKepid(kepid, quarter);
    taskdir = getPdcTaskdirForCh(datapath, ch, quarter, varargin{:});
end
