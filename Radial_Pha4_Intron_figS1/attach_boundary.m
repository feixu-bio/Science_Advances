%% investigate the spots that are associated with the boundary.
function[attached] = attach_boundary(radial,LADs_thresh)
    indx = find([radial.raw_dis]<LADs_thresh);
    dotID = [radial.dotID];
    attached.spotID_attached = dotID(indx);
    attached.prop_attached = length(indx)/length([radial.raw_dis]); % how many proportion of the spots are attached
    intensity = [radial.max_intensity];
    attached.intensity_attached = intensity(indx); % extract the intensity of the attached spots.
