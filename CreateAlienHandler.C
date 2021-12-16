AliAnalysisGrid* CreateAlienHandler()
{
    // Check if user has a valid token, otherwise make one. This has limitations.
    // One can always follow the standard procedure of calling alien-token-init then
    //   source /tmp/gclient_env_$UID in the current shell.
    // if (!AliAnalysisGrid::CreateToken()) return NULL;
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode("test");
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    // plugin->SetROOTVersion("v5-34-02-1");
    //plugin->SetAliPhysicsVersion("vAN-20210708_ROOT6-1");//("vAN-20190301-1");
    plugin->SetAliPhysicsVersion("vAN-20211213_ROOT6-1");

    // plugin->SetAliPhysicsVersion("vAN-20210111");
     
    // Declare input data to be processed.
    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    
    //Location for debugging:
    // /alice/cern.ch/user/p/pwg_pp/JIRA/PWGPP-218/000246087/pass1
    
    // plugin->SetGridDataDir("/alice/data/2015/LHC15o"); 
    plugin->SetGridDataDir("/alice/sim/2020/LHC20i2b");


    //plugin->SetGridDataDir("/alice/data/2018/LHC18q");
    //plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
    
    //LHC17f2a_fast_fix
    //plugin->SetGridDataDir("/alice/sim/2017/LHC17f2a_fast_fix");
    
    //plugin->SetGridDataDir("/alice/data/2016/LHC1");
    //plugin->SetGridDataDir("/alice/data/2010/LHC10h");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Set data search pattern
    //plugin->SetDataPattern("/muon_calo_pass1/*ESDs.root");
    plugin->SetDataPattern("*ESDs.root");//("pass5_lowIR/*ESDs.root");
    //plugin->SetDataPattern("*ESDs.root");
    // Data pattern for reconstructed data
    //   plugin->SetDataPattern("*ESDs/pass4/*ESDs.root");
    plugin->SetRunPrefix("");   // real data: 000
    // plugin->SetUseSubmitPolicy();
    
   Long_t runnumbers[] = {288909};
    
    Int_t numberofruns = sizeof(runnumbers)/sizeof(Long_t);
    
    for(int atrun = 0;atrun<numberofruns;atrun++){
        plugin->AddRunNumber(runnumbers[atrun]);
    }
    
    //   plugin->AddRunNumber(104065);  // real data
    //   plugin->SetOutputSingleFolder("output");
    plugin->SetOutputToRunNo();
    // Method 2: Declare existing data files (raw collections, xml collections, root file)
    // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
    // XML collections added via this method can be combined with the first method if
    // the content is compatible (using or not tags)
    //   plugin->AddDataFile("tag.xml");
    //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
    // Define alien work directory where all files will be copied. Relative to alien $HOME
    plugin->SetGridWorkingDir("analysisV0_LHC20i2b");
    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    
    //uncomment for specific tasks
    //plugin->SetAnalysisSource("AliMultSelectionTask2.cxx");
    //plugin->SetAnalysisSource("AliAnalysisTaskSkeleton.cxx");
    
    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("libpythia6_4_21.so libOADB.so libPWGLFSTRANGENESS.so libRUN3.so");
    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    //   plugin->SetOutputFiles("lambda0recooutput.root");
    Bool_t lTrueFlag = kTRUE;
    plugin->SetDefaultOutputs(lTrueFlag);
    plugin->SetAdditionalLibs("AliAnalysisTaskQAMultistrangev2.h");
   // plugin->SetAnalysisSource("AliAnalysisTaskQAMultistrangev2.cxx ");
      
    // Optionally define the files to be archived.
    //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
    //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro("AliAnalysisTaskExtractCombined.C");
    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetSplitMaxInputFileNumber(150);
    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable("AliAnalysisTaskExtractCombined.sh");
    // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
    //   plugin->SetMaxInitFailed(5);
    // Optionally resubmit threshold.
    //   plugin->SetMasterResubmitThreshold(90);
    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(45000);
    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName("AliAnalysisTaskExtractCombined.jdl");
    // Optionally modify job price (default 1)
    plugin->SetPrice(1);
    // Optionally modify split mode (default 'se')
    plugin->SetSplitMode("se");
    //configure to merge files sequentially rather than all-at-once
    plugin->SetMergeViaJDL(kFALSE);
    
    plugin->SetNtestFiles(1);
    
    plugin->SetMaxMergeFiles(60);

    // plugin->SetMaxMergeStages(10);  // aimeric addition

    // plugin->WorkDirectorySize(10000); // default 5000MB

    return plugin;
}
