function [out, filenames] = rotaball_importer(varargin)

n = 1;
for i = 2:2:length(varargin)
    refs(n).date = varargin{i-1};
    refs(n).expnumber = varargin{i};
    n = n+1;
end

all_dir = dir;

for e = 1:length(refs)
    ii = 1;
    for i = 1:length(all_dir)
        if ~isempty(strfind(all_dir(i).name, ['E' num2str(refs(e).expnumber) '_'])) && strcmp(all_dir(i).name(1), 'E')
            if strcmp(all_dir(i).name(end-3:end), '.txt')
                out{e}{ii} = rotaball_fileimporter([all_dir(i).folder '/' all_dir(i).name]);
                filenames{e}{ii,1} = all_dir(i).name;
                ii = ii + 1;
            end
        end
    end
end
        
