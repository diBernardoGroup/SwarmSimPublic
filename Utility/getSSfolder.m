function [SSdir] = getSSfolder()
    cur_dir = mfilename('fullpath');
    cur_dir = fileparts(cur_dir);
    while ~any(strcmp('Simulation', getSubfolders(cur_dir)))
    cur_dir = fileparts(cur_dir);
    end
    SSdir = cur_dir;
end

