function [p_traj] = plotTrajectory(xVec,thenDelete, color)
%
%
    p_traj = plot(xVec(:,:,1),xVec(:,:,2), 'color',color);
    
    if(thenDelete)
        drawnow
        delete(p_traj)
    end
end

