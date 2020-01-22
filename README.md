# P2O
Full system stock &amp; flow model of plastic use, environmental leakage, economics, jobs &amp; GHG emissions

The model was written in Matlab, version R2017b (9.3.0.713579), by Richard Bailey, and runs on Mac OS, Windows, Linux.

# Instructions for use
- Choose a numbered Archetype and Scenario (see below)
- Copy configuration files for chosen combination from ~/P2O/config_files_store in to ~/P2O/config_files 
- Open "run_P2O.m" in Matlab
- Set relevant choices within the setup sections
- Run script "run_P2O.m". Text output to the Command Window provides information as the model runs.
- Once finished, results are found in ~/P20/Output_files

The model can be run in 'single run' and 'batch' modes. 
For single runs (set 'batch_run=false()' within 'run_P2O.m'), the chosen archetype/scenario combination is specified explicitly within the code.
In batch mode (setting 'batch_run=true()' within 'run_P2O.m'), the archetype and scenario numbers are read from file 'arch_sen_to_run.csv'. Copies of the full directory structure (including all code and relevant configuration files), can be created by the user, along with appropriate 'arch_sen_to_run.csv' files, and set running without the need to change archetype/scenario specification within the code.
The included MS Excel sheet 'Input_Files_List.xlsx' provides details of the content of configuration files, so changes can be made to parameters/inputs if needed.

ARCHETYPES:               
1. HI_Urban      
2. HI_Rural 
3. UMI_Urban
4. UMI_Rural
5. LMI_Urban
6. LMI_Rural
7. LI_Urban
8. LI_Rural

SCENARIOS:
1. Baseline
2. Current commitments
3. Linear Scenario
4. Recycling scenario
5. R&S scenario
6. Scenario X (note 'Scenario X' is also referred to as the 'System Change Scenario, SCS' in associated publications)

