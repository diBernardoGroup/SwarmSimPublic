function [splot,lplot] = mysurfc(x_vec,y_vec,z,level,curveLengthTol)
%
%mysurfc surface polt with level set lines.
%   Used to plot 2D functions.
%   Short contour lines can be omitted from the plot.
%
%   [splot,lplot] = mysurfc(x_vec,y_vec,z,level,curveLengthTol)
%
%   Inputs:
%       x_vec are the x values (vector)
%       y_vec are the y values (vector)
%       z are the values of the surface to be plotted (matrix)
%       level is the value to draw the level line (scalar)
%       curveLengthTol fraction of the length of the longest level curve
%           to use as a threshold to draw more level curves. 0 draw all the level
%           curves, 1 draw only the longest level curve. (scalar [0,1])
%
%   Outputs:
%       splot is the plot of the surface
%       lplot is the plot of the level line
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    x_vec   (1,:)   double
    y_vec   (1,:)   double
    z               double 
    level           double = 0
    curveLengthTol  double = 0.6
end

levelset=contourc(x_vec, y_vec, z', [level level]);
nanindex=find(levelset(1,:)==level & floor(levelset(2,:)) == levelset(2,:));
levelset(1,nanindex)=nan;
[maxCurveLength,maxCurveindex]=max(levelset(2,nanindex));
cuveIndexes=find(levelset(2,nanindex)>=maxCurveLength*curveLengthTol);
curveLenghts=levelset(2,nanindex(cuveIndexes));

newlevelset=nan(2,0);
for i=1:length(curveLenghts)
    newlevelset=[newlevelset, levelset(:,nanindex(cuveIndexes(i)):nanindex(cuveIndexes(i))+curveLenghts(i))];
end

[x_mesh, y_mesh] = meshgrid(x_vec, y_vec);

splot=surf(x_mesh, y_mesh, z');
shading interp
view(2)
hold on
lplot=plot3(newlevelset(1,:),newlevelset(2,:),10*ones(length(newlevelset),1),'black-','LineWidth', 1);
caxis([0 1.5*level]);
colorbar;

if isempty(lplot); lplot=plot(nan,'black-','LineWidth', 1); end

end

