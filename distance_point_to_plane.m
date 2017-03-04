function distance=distance_point_to_plane(point,plane_formular)
% given a 3-D point's location, and a the coefficient of plane's formular
% computer the distance from the point to the plane
        distance= ([point;1]'*plane_formular)./(norm(plane_formular(1:3)));
end