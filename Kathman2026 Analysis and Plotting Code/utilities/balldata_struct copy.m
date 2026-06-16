 function data = balldata_struct(inputdata, filenames)


%col1: timestamp
%col2: frame#
%col3: fictrac heading
%col4: fictrac speed
%col5: heading (?)
%col6: droll
%col7: dpitch
%col8: dyaw (raw)
%col9: dyaw (filtered (for voltage?))
%col10: motor voltage
%col11: motor sign
%col12: motor position
%col13: gain
%col14: light state/allo_shift
%col15: wind state
%col16: wind speed
%col17: wind start position
%col18: odor state
%col19: analog in value
%col20: plume t
%col21: plume x
%col22: plume dx
%col23: plume y
%col24: plume dy
%col25: plume calculated heading (unfiltered)
%col26: plume calculated heading (filtered)
%col27: plume calculated speed
%col28: plume odor value
%col29: plume compressed odor value
%col30: valve state
%col31: adaptive compression value (A)
%col32: 
%col33: 
%col34: 
%col35: 

% %Check filtered wind position
% figure(1); clf; hold on
% for t = 1:24
%     frame = inputdata{1}{t}(:,2);
%     fps = round(1/mean(diff(inputdata{1}{t}(:,1))));
%     ts = inputdata{1}{t}(:,1);
%     calc_ts = ts + diff([ts(1), frame(1)/fps]);
%     wp = inputdata{1}{t}(:,12)/2310*2*360;
%     plot(calc_ts, wp+t*1800)
% end


for t = 1:length(inputdata{1})
    trial.filename = filenames{1}{t};
    trial.fps = round(1/mean(diff(inputdata{1}{t}(:,1))));
    trial.wind_speed = floor(mode(inputdata{1}{t}(:,16))*1000/60/(2.5^2*pi)*100); %cm/s
    trial.winddir = mode(inputdata{1}{t}(:,17));
    
    
    trial.shift = round(max(abs(inputdata{1}{t}(:,14))));
    if ~isempty(trial.shift) && trial.shift > 0 && min(inputdata{1}{t}(:,14)) < 0
        trial.shift = trial.shift *-1;
    end

%     trial.light_on = max(inputdata{1}{t}(:,14));
    trial.wind_on = max(inputdata{1}{t}(:,15));
    trial.odor_on = max(inputdata{1}{t}(:,18));
    if isempty(strfind(filenames{1}{t}, 'CL'))
        trial.closed_loop = 0;
    else
        trial.closed_loop = 1;
    end
    trial.gain = mode(inputdata{1}{t}(:,13));
    trial.allo_pos = inputdata{1}{t}(:,14);
    trial.wind  = inputdata{1}{t}(:,15);
    trial.odor = inputdata{1}{t}(:,18);
    trial.ai = inputdata{1}{t}(:,19);

    
    trial.frame = inputdata{1}{t}(:,2);
    ts = inputdata{1}{t}(:,1);
    if isempty(ts)
        trial.calc_ts=NaN;
    else
        trial.calc_ts = ts + diff([ts(1), trial.frame(1)/trial.fps]); %offset ts so all trials aligned by frame#
    end
    
    head = unwrap(-inputdata{1}{t}(:,3))*180/pi; %fictrac heading
    if isempty(head)
        trial.calc_heading=NaN;
    else
        trial.calc_heading = wrapTo180(head-head(1)); %starts at zero
    end

    %%%%%% filter out data loss in motorposition data %%%%%%%%
    wp = inputdata{1}{t}(:,12)/2310*2*360;
    %wp = inputdata{1}{t}(:,12)/1000*360;
    trial.raw_windpos = wp;
    if isempty(wp)
        trial.calc_windpos=NaN;
    else
        %threshold filters derivatives for long durations of bad mposes (up to 50 samples)
        for ii = 1:60
            dwp = find(diff(wp) > 4); %deriv thresh
            if ~isempty(dwp)
                for i = 1:length(dwp)
                    if wp(dwp(i) + 1) > 32 %error pos thresh
                        wp(dwp(i)) = wp(dwp(i) + 1);
                    else
                        wp(dwp(i) + 1) = wp(dwp(i));
                    end
%                     plot(wp); pause(0.05)
                end
            end
        end
        trial.calc_windpos = wrapTo180(wp);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    if isempty(inputdata{1}{t}(:,4))
        trial.calc_speed=NaN;
    else
        trial.calc_speed = inputdata{1}{t}(:,4)*9.52/2*trial.fps;
    end
    if isempty(-inputdata{1}{t}(:,9))
        trial.calc_deltaz=NaN;
    else
        trial.calc_deltaz = -inputdata{1}{t}(:,9)*trial.fps;
    end
    if isempty(inputdata{1}{t}(:,7))
        trial.calc_deltapitch=NaN;
    else
        trial.calc_deltapitch = inputdata{1}{t}(:,7)*trial.fps;%*9.52/2*trial.fps; %forgot to add this initially,dont want to go back and reanalyze all old data
    end
    if isempty(inputdata{1}{t}(:,6))
        trial.calc_deltaroll=NaN;
    else
        trial.calc_deltaroll = inputdata{1}{t}(:,6)*trial.fps;
    end

    if ~all(isnan(inputdata{1}{t}(:,33)),1) %is there data in column 33 (which means its from old plume callback (i think just 220614?))

        trial.pybmt_heading = inputdata{1}{t}(:,26);
        trial.plume_t = inputdata{1}{t}(:,20);
        trial.plume_x = inputdata{1}{t}(:,21);
        trial.plume_dx = inputdata{1}{t}(:,22);
        trial.plume_y = inputdata{1}{t}(:,23);
        trial.plume_dy = inputdata{1}{t}(:,24);
        trial.plume_speed = inputdata{1}{t}(:,27);
        trial.plume_odor = inputdata{1}{t}(:,28);
        trial.plume_codor = inputdata{1}{t}(:,29);
        trial.plume_valve = inputdata{1}{t}(:,30);
        trial.plume_A = inputdata{1}{t}(:,31);
        trial.plume_beta = inputdata{1}{t}(:,32);
        trial.plume_spatgain = inputdata{1}{t}(:,33);
%         trial.plume_unfiltheading = inputdata{1}{t}(:,25);
    elseif ~all(isnan(inputdata{1}{t}(:,22)),1) 
        trial.pybmt_heading = inputdata{1}{t}(:,20);
        trial.reset_window = inputdata{1}{t}(:,21);
        trial.stim_col = mode(inputdata{1}{t}(:,22));
        trial.trig_state = inputdata{1}{t}(:,23);
        if ~all(isnan(inputdata{1}{t}(:,24)),1) 
            trial.stim_intensity = mode(inputdata{1}{t}(:,24));
            trial.stim_dwelltime = mode(inputdata{1}{t}(:,25));
        end

        trial.plume_t = inputdata{1}{t}(:,33);
        trial.plume_x = inputdata{1}{t}(:,34);
        trial.plume_dx = inputdata{1}{t}(:,35);
        trial.plume_y = inputdata{1}{t}(:,24);
        trial.plume_dy = inputdata{1}{t}(:,25);
        trial.plume_speed = inputdata{1}{t}(:,26);
        trial.plume_odor = inputdata{1}{t}(:,27);
        trial.plume_codor = inputdata{1}{t}(:,28);
        trial.plume_valve = inputdata{1}{t}(:,29);
        trial.plume_A = inputdata{1}{t}(:,30);
        trial.plume_beta = inputdata{1}{t}(:,31);
        trial.plume_spatgain = inputdata{1}{t}(:,32);
    else
        trial.pybmt_heading = inputdata{1}{t}(:,20);
        trial.plume_t = inputdata{1}{t}(:,21);
        trial.plume_x = inputdata{1}{t}(:,22);
        trial.plume_dx = inputdata{1}{t}(:,23);
        trial.plume_y = inputdata{1}{t}(:,24);
        trial.plume_dy = inputdata{1}{t}(:,25);
        trial.plume_speed = inputdata{1}{t}(:,26);
        trial.plume_odor = inputdata{1}{t}(:,27);
        trial.plume_codor = inputdata{1}{t}(:,28);
        trial.plume_valve = inputdata{1}{t}(:,29);
        trial.plume_A = inputdata{1}{t}(:,30);
        trial.plume_beta = inputdata{1}{t}(:,31);
        trial.plume_spatgain = inputdata{1}{t}(:,32);
    end

%if nicks code with reset to start trial or not
    if all(isnan(trial.plume_y))
        trial.allo_reset = round(median(inputdata{1}{t}(:,21)));
        trial.plume_t(:) = nan;
    else
        trial.allo_reset = nan;
    end
% if length(trial.allo_reset) ~= length(trial.plume_t)
%     keyboard
% end

%if will's code with offset
%     if all(isnan(trial.plume_y))
%         trial.allo_offset = inputdata{1}{t}(:,21);
%         trial.plume_x(:) = nan;
%     end

    %     out{t}.frame_num = data{1}{t}(:,2);
    %     out{t}.FTheading = data{1}{t}(:,3);
    %     out{t}.FTspeed = data{1}{t}(:,4);
    %     out{t}.heading = data{1}{t}(:,5);
    %     out{t}.delta_roll = data{1}{t}(:,6);
    %     out{t}.delta_pitch = data{1}{t}(:,7);
    %     out{t}.unfilt_delta_heading = data{1}{t}(:,8);
    %     out{t}.delta_heading = data{1}{t}(:,9);
    %     out{t}.motorvoltage = data{1}{t}(:,10);
    %     out{t}.motorsign = data{1}{t}(:,11);
    %     out{t}.motorposition = data{1}{t}(:,12);

    data(t)=trial;
    dsh = strfind(filenames{1}{t}, '-');
    file_ts(t,1) = str2num([filenames{1}{t}(dsh(end)-2:dsh(end)-1), filenames{1}{t}(dsh(end)+1:dsh(end)+6)]);

%     figure(t);clf; hold on; plot(inputdata{1}{t}(:,6)*trial.fps, 'b');plot(inputdata{1}{t}(:,6)*trial.fps+7, 'k');plot(inputdata{1}{t}(:,9)*trial.fps+14, 'r')
    
end
[x,i] = sort(file_ts); %sort by trial timestamp, not filename
data = data(i);

us = strfind(filenames{1}{1}, '_');
filename=[filenames{1}{1}(1:us(1)-1)];
save(filename,'data');
end
