function fileList = getFilenames(varargin)
% Recursively gets files within a directory
%   Either give just a dirName to get every file (excluding OS files)
%   Or give a dirName and a cell list of regexps that you want the files to match
  switch nargin
      case 1
          dirName = varargin{1};
      case 2
          dirName = varargin{1};
          includes = varargin{2};
  end

  excludes = {'.DS_Store'};
  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  
  % Include files
  if exist('includes')
      validFiles = cellfun(@(s) ~isempty(cell2mat(regexp(s,includes,'once'))) && isempty(cell2mat(regexp(s,excludes,'once'))),fileList);
  else
      validFiles = ~ismember(fileList,excludes);
  end
  
  % AND for exclude & include lists
      
  fileList = fileList(validFiles);
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..','Unprocessed','Unsorted'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    if exist('includes')
        fileList = [fileList; getFilenames(nextDir,includes)];  %# Recursively call getAllFiles
    else
        fileList = [fileList; getFilenames(nextDir)];
  end
end