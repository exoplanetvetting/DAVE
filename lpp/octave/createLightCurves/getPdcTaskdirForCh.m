function taskDir = getPdcTaskdirForCh(dataPath, chOrModOut, quarter, varargin)
%Find the directory name that stores a given channel and quarter
%$Id: getPdcTaskdirForCh.m 51755 2013-05-23 21:43:30Z fmullall $
%$URL: svn+ssh://murzim.amn.nasa.gov/repo/so/trunk/Develop/jvc/taskfileIO/getPdcTaskdirForCh.m $
%Fergal Mullally
%
%Inputs:
%dataPath:  (string) Where is the data stored. e.g ~/proc
%chOrModOut	Channel or modout. Integer is treated as a channel, mod-out
%		can be given as a decimal (e.g 3.1) or an array [3 1] 
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
%csci                   (str) Which csci to load. Defaults to pdc. Not
%                       well tested for other CSCIs.
%scMonth                (int, [1-3]) Return directory for the given short
%                       cadence month of *quarterly* processing
%monthlySc              (int, [1-3]) Return directory for the given short
%                       cadence month of the *monthly* processing
%monthly                (int, [1-3]) Use this to return directory for
%                       the given month of *long* cadence monthly processing
%Returns:
%A string giving the path to the .mat files


    %If modout is written 2.1, convert to array
    if chOrModOut ~= floor(chOrModOut)
        module = floor(chOrModOut);
        output = 10*(chOrModOut-module);
        chOrModOut = [module output];
    end

    %Process optional args
    optargin = size(varargin, 2);
    if mod(optargin, 2) ~=  0
        error('Odd number of optional arguments supplied');
    end
    
    %Set default values for opt args. KSOP is 0 by default.
    %If not set by input argument, we use a look up table, and only
    %if that value is needed (i.e not for monthly processings)
    isFullyQualifiedDataPath = false;
    cadenceType = 'LONG';
    ksop = 0;
    lcMonth = 0;
    scMonth = 0;
    monthlyProc = false;
    csci = 'pdc';

    for i = 1 : 2 : length(varargin)
        name = varargin{i};
        value = varargin{i+1};
        switch name
            case 'ksop'
                ksop = value;
            case 'isFullyQualifiedDataPath'
                isFullyQualifiedDataPath = value;
            case 'monthly'
                monthlyProc = true;
                lcMonth = value;
            case 'monthlySc';
                monthlyProc = true;
                scMonth = value;
            case 'scMonth'
                scMonth = value;
                cadenceType = 'SHORT';
            case 'csci'
                csci = value;
            otherwise
                error('Unrecognised option %s to getPdcTaskStruct', name);
        end
    end
 
    
    
    
    if ~isFullyQualifiedDataPath
        if monthlyProc
            quarterDir=sprintf('q%02i', quarter);
            ppath=fullfile(dataPath, quarterDir, 'pipeline_results', 'monthly');
            
            if lcMonth > 0
                p = sprintf('lc-m%i', lcMonth);
            else 
                p = sprintf('sc-m%i', scMonth);
            end
            dataPath = fullfile(ppath, p);
        else
            %Figure out the fully qualified data path
            %As I've said before, Matlab string handling is terrible
            quarterDir=sprintf('q%02i', quarter);
            ppath=fullfile(dataPath, quarterDir, 'pipeline_results', 'q*archive*ksop*');

            if ksop == 0
                ksop = getKsopForQuarter(quarter);
            end
            ppath = sprintf('%s%i*', ppath, ksop);
            archiveDir = dir(ppath);

            if isempty(archiveDir)
                error('Path not found: %s', ppath);
            end

            ppath=fullfile(dataPath, quarterDir, 'pipeline_results', archiveDir(1).name);

            if scMonth > 0
                scPath = sprintf('sc-m%1i', scMonth);
                dataPath = fullfile(ppath, scPath);
            else
                dataPath=fullfile(ppath,  '/lc/mpe_true/');
            end
        end
    end
    
    if isempty(dir(dataPath))
        error('Path not found: %s', dataPath)
    end
    
    csvFile = getTaskdirCsvFilename(dataPath, csci);
    
    %Call the SOC API
    cell = get_taskfiles_from_modout(csvFile, csci, chOrModOut, [], quarter, cadenceType);
    
    if isempty(cell)
        warning('taskfileIO:getPdcTaskdirForCh', 'No taskfile found');
        taskDir = 'NoTaskFileFound';
    else
        taskDir = cell{1};
    end
        
end
