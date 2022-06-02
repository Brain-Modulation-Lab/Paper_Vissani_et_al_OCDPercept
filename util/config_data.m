% get info on patients and sessions

PATIENTS = cellstr(ls(PATH_DATA));
PATIENTS = PATIENTS(~contains(PATIENTS,". "));

n_PATIENTS = numel(PATIENTS);

info_data = struct();

for pat_i = 1 : n_PATIENTS
    info_data(pat_i).PATIENT = PATIENTS{pat_i};
    if ispc
        FILES = cellstr(ls(fullfile(PATH_DATA,PATIENTS{pat_i},'*.json')));
        FILES = FILES(~contains(FILES,". "));
        info_data(pat_i).FILES = {FILES};
    elseif ismac
        FILES = dir(fullfile(PATH_DATA,PATIENTS{pat_i},'*.json'));
        FILES = {FILES.name};
        FILES = FILES(~contains(FILES,". "));
        info_data(pat_i).FILES = FILES;
    end
    info_data(pat_i).n_FILES = numel(FILES);
end

