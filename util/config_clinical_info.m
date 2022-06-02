% get clinical info from csv file

if ispc
    FILE = string(ls(fullfile(PATH_DATA,PATIENTS{pat_i},'*.csv')));
    FILE_CLINIC = FILE(contains(FILE,"clinic"));
    FILE_MED = FILE(contains(FILE,"med"));
elseif ismac
    FILE = dir(fullfile(PATH_DATA,PATIENTS{pat_i},'*.csv'));
    FILE = {FILE.name};
    FILE_CLINIC = FILE{contains(FILE,"clinic")};
    FILE_MED = FILE{contains(FILE,"med")};

end

info_data.clinical_info = readtable(fullfile(PATH_DATA,PATIENTS{pat_i},FILE_CLINIC));
info_data.med_changes = readtable(fullfile(PATH_DATA,PATIENTS{pat_i},FILE_MED));
