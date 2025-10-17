function [angle] = angleBetweenVectors(v1,v2)
%ANGLEBETWEENVECTORS Summary of this function goes here
%   Detailed explanation goes here

    if size(v1,2)==2 && size(v1,1)~=2 && size(v2,2)==2 && size(v2,1)~=2
        v1=v1';
        v2=v2';
    end
    
    if size(v1,1)==2    % in 2D
        angle = atan2(dot(v1,[0 1; -1 0]*v2), dot(v1,v2));
    else                % in 3D
        angle = atan2(norm(cross(v1,v2)), dot(v1,v2));
    end
end

