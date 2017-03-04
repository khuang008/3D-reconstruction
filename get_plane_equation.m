function plane=get_plane_equation(points)
%given three points,compute the coefficient of plane's equation defined by these points

    plane_normal = cross(points(:,1)-points(:,2),points(:,3)-points(:,2));
    a=plane_normal(1);
    b=plane_normal(2);
    c=plane_normal(3);
    d=-points(:,1)'*plane_normal;
    plane=[a;b;c;d]./norm([a;b;c;d]);
    
    
    
end