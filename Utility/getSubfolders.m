function subfolders = getSubfolders(directory)
    % get the folder contents
    d = dir(directory);
    % remove all files (isdir property is 0)
    dfolders = d([d(:).isdir]); 
    % remove '.' and '..' 
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    % get only the names
    dfolders = {dfolders(:).name};
    for i = 1:length(dfolders)
        subfolders(i) = string(dfolders(i));
    end
end

