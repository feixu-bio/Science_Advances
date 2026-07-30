%% calculate the distance to the boundary and record the distribution vs intensity
function[radial] = distance_boundary(T_alpha,spot,background_intensity)
    k=0;
    for i = 1:height(T_alpha)
        dotID = [];
        normalized_dis = [];
        max_intensity = [];
        mean_intensity = [];
        raw_dis = [];
        for j = 1:length(spot)
            %some alphashape is empty and some are only have few points
            %which need to be filter out.
            if ~isempty(T_alpha.alphaShape{i,1}.Points) && length(T_alpha.alphaShape{i,1}.Points)>100
                tf = inShape(T_alpha.alphaShape{i,1},spot(j).x,spot(j).y,spot(j).z);
                if tf == 1
                   I = nearestNeighbor(T_alpha.alphaShape{i,1},spot(j).x,spot(j).y,spot(j).z);
                   near_dot_x = T_alpha.alphaShape{i,1}.Points(I,1); %in um
                   % need to flip the y axis to match the orginal image and
                   % spot fitting.
                   near_dot_y = T_alpha.alphaShape{i,1}.Points(I,2); % in um
                   near_dot_z = T_alpha.alphaShape{i,1}.Points(I,3); % in um
                   % calculate the spot to boundary
                   spot_bound = ((spot(j).x - near_dot_x).^2+ ...
                                        (spot(j).y - near_dot_y).^2+ ... 
                                        (spot(j).z - near_dot_z).^2).^0.5;
                   r = T_alpha.radius(i);
                   norm_dis = spot_bound / r;
                   dotID = [dotID j];
                   raw_dis = [raw_dis spot_bound];
                   normalized_dis = [normalized_dis norm_dis];
                   max_intensity = [max_intensity spot(j).max_inty-background_intensity];
                   mean_intensity = [mean_intensity spot(j).mean_inty-background_intensity];
                   %normalize to the median intensity in the nucleus where the
                   %dots are.
                end
            end
        end
        if ~isempty(dotID) && length(dotID) <= 4
            % dot number in each nucleus should be less than 4. if not
            % filter out.
           k=k+1;
           radial(k).nucID = i;
           radial(k).dotID = dotID;
           radial(k).norm_dis = normalized_dis;
           radial(k).max_intensity = max_intensity;
           radial(k).mean_intensity = mean_intensity;
           radial(k).raw_dis = raw_dis;
        else
            continue
        end
    end
