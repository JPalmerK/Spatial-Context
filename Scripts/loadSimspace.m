% Function for loading the worspace and required mat files on for PC and
% Linux locations

function whereAmI = loadSimspace()

whereAmI={};

try
    
    cd('/home/kpalmer/AnacondaProjects/Spatial-Context/Scripts');
    
    whereAmI(1) = {'/cache/kpalmer/quick_ssd/data/DCLDE_2013_10Channel/DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'};
    whereAmI(2) = {'DCLDE2013_RW_localizationsDCLDE_2013_10_Chan_all12_timed_July1919_localize_struct1_95.mat'};
    whereAmI(3) = {'DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_July1919_hyd1_95.mat'};
    whereAmI(4) = {'DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_July1919_calls1_95.mat'};
    whereAmI(5) = {'DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_July1919_array1_95.mat'};

    
        
    
catch ME
    %%
    cd('D:\Anaconda Projects\Spatial-Context\Scripts')
    
    whereAmI(1) = {'C:\Users\Kaitlin\Desktop\DCL2013_NEFSC_SBNMS_200903_metadata.xlsx'};
    whereAmI(2) = {'D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_array_struct1_95.mat'};
    whereAmI(3) = {'D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_localize_struct_95.mat'};
    whereAmI(4) = {'D:/data/SimStructures/DCLDE2013_RW_localizations_DCLDE_2013_10_Chan_all12_timed_Mar21_calls1_95.mat'}
    
    
end






end
