function EXPT = context_expt(local)

    % creates EXPT structure for CCNL fMRI processing pipeline
    %
    % USAGE: EXPT = contextExpt()
    %
    % INPUTS:
    %   local (optional) - true if file path is on local computer, false if on NCF
    %
    % OUTPUTS:
    %   EXPT - experiment structure with fields
    %          .TR - repetition time
    %          .create_multi - function handle for creating multi structure
    %          .modeldir - where to put model results
    %          .subject(i).datadir - directory for subject data
    %          .subject(i).functional - .nii files for runs
    %          .subject(i).structural - .nii for structural scan
    %
    % Cody Kommers, July 2016
    % Momchil Tomov, Nov 2016
    
    % set default parameters
    %
    if nargin < 1
        [~, name] = system('hostname');
        if strfind(name, 'MacBook') % probably local           
            local = true;
        else
            local = false;
        end
    end
    
    % set main directory
    %
    if local
        exptdir = '/Users/memsql/Dropbox/research/context/'; % local group level
    else
        exptdir = '/ncf/gershman/Lab/ConLearn/'; % on CBS central server
    end
    
    % Load data from file with all subjects
    %
    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
    
    [allSubjects, subjdirs, nRuns] = context_getSubjectsDirsAndRuns();
    allSubjects
    metadata.allSubjects
    assert(isequal(allSubjects, metadata.allSubjects));
    
    for subj = 1:length(allSubjects)
        subjdir = [exptdir, 'subjects/', subjdirs{subj}, '/'];
        EXPT.subject(subj).datadir = [subjdir, 'preproc'];
        
        EXPT.subject(subj).structural = 'struct.nii';
        
        assert(nRuns{subj} == length(unique(data.runId(strcmp(data.participant, allSubjects{subj})))));
        for run = 1:nRuns{subj}
            EXPT.subject(subj).functional{run} = ['run',sprintf('%03d',run),'.nii'];
        end
        disp(EXPT.subject(subj));
    end
    
    % TR repetition time
    EXPT.TR = 2; %seconds
    % Function handle to create subject multi structure
    EXPT.create_multi = @context_create_multi;
    % Where you want model output data to live
    if local
        EXPT.modeldir = [exptdir, 'neural'];
    else
        EXPT.modeldir = [exptdir, 'glmOutput'];
    end
    
    % Where the data live, but not sure which data
    EXPT.datadir = [exptdir, 'testOutput'];


end