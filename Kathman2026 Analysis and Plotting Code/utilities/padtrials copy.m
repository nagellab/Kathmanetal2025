
function paddata = padtrials(data)

    fld = fieldnames(data);
    mxfrms = 0; 
    for t = 1:length(data)
        if max(data(t).frame) > mxfrms
            mxfrms = max(data(t).frame); 
        end %find length of longest trial to set dimensions of mean matrix
    end

    for f = 1:length(fld)
        if isa(data(1).(fld{f}), 'double') && ~strcmp(fld(f), 'frame')
            for t = 1:length(data)
                if ~isempty(data(t).dff_green)
                    uset = t; %just find a good trial (w dff data)
                end
            end
%             if f == 42
%                 keyboard; end
            if size(data(uset).(fld{f}), 1) > 1 && size(data(uset).(fld{f}), 2) == 1
                for t = 1:length(data)
                    if ~isempty(data(t).(fld{f})) && ~isempty(data(t).frame) && ...
                            length(data(t).(fld{f})) == length(data(t).frame)
                        paddata(t).(fld{f}) = nan(mxfrms, 1);
                        paddata(t).(fld{f})(data(t).frame) = data(t).(fld{f});
                    else
                        paddata(t).(fld{f}) = data(t).(fld{f});
                    end
                end
            elseif size(data(uset).(fld{f}), 1) > 1 && size(data(uset).(fld{f}), 2) == 2
                for t = 1:length(data)
                    if ~isempty(data(t).(fld{f})) && ~isempty(data(t).frame)
                        paddata(t).(fld{f}) = nan(mxfrms, 2);
                        paddata(t).(fld{f})(data(t).frame, :) = data(t).(fld{f});
                    else
                        paddata(t).(fld{f}) = data(t).(fld{f});
                    end
                end
            else
                for t = 1:length(data)
                    paddata(t).(fld{f}) = data(t).(fld{f});
                end
            end
        else
            for t = 1:length(data)
                paddata(t).(fld{f}) = data(t).(fld{f});
            end
        end
    end
% keyboard
% 
% clear all
% A = normpdf(1:8,3.5, 2);
% B = [1 2 6 4 7 5 2 1];
% C = [2 6 7 8 6 5 4 3];
% A = [A',B',C'];
% % Method 1, using mean
% x = 1 : size(A, 2); % Columns.
% y = 1 : size(A, 1); % Rows.
% [X, Y] = meshgrid(x, y)
% 
% 
% meanA = mean(A)
% % centerOfMassX = mean(A .* X) / meanA
% centerOfMassY = mean(A .* Y) ./ meanA
% 
% 
% % meanA = mean(A(:))
% % centerOfMassX = mean(A(:) .* X(:)) / meanA
% % centerOfMassY = mean(A(:) .* Y(:)) / meanA
% 
% 
% 
% % % Method 2: using regionprops
% % props = regionprops(true(size(A)), A, 'WeightedCentroid')