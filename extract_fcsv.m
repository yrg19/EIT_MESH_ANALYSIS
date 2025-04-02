%% Function for extracting coordinates from fcsv file 

function coordinates = extract_fcsv(filename) 


% Open the file for reading
fileID = fopen(filename, 'r');

% Check if the file opened correctly
if fileID == -1
    error('Failed to open the file.');
end

% Read the file into a cell array line by line
rawData = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

% Extract lines
lines = rawData{1};

% Find the first non-header line (skipping lines starting with '#')
dataLines = lines(~startsWith(lines, '#'));

% Parse the data lines
numPoints = numel(dataLines);
coordinates = zeros(numPoints, 3); % Preallocate for x, y, z coordinates
labels = strings(numPoints, 1);    % Preallocate for labels

for i = 1:numPoints
    % Split the line by commas
    tokens = strsplit(dataLines{i}, ',');

    % Parse x, y, z coordinates
    coordinates(i, 1) = str2double(tokens{2}); % x
    coordinates(i, 2) = str2double(tokens{3}); % y
    coordinates(i, 3) = str2double(tokens{4}); % z

    % Parse label (usually the 12th column)
    if numel(tokens) >= 12
        labels(i) = tokens{12};
    end
end

end