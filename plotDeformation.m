%-------------------------------------------------------------------------
%                      plotDeformation.m
%                      Michael Dodman
%                      January 2022
%
% This function plots the deformed and undeformed plate overlayed on each
% other with a 100X displacement magnification for clarity.
%-------------------------------------------------------------------------
% Variable Name      Definition
% ------------------------Input--------------------------------------------
% elem               Matrix of node indices of each element
% Nxy                [Y,X] matrix of node positions
% U                  Nodal Displacement Matrix
%-------------------------Output-------------------------------------------
% deformation_plot   Axes object for nodal deformation plot
%-------------------------Internal----------------------------------------
% row                Number of elements
% def                Nodal displacement values with 100X magnification
% pgon               Polyshape object for each element
% pgon_plot          Axes object for undeformed plot
% pgon_def           Polyshape object for each deformed element
% pgon_def_plot      Axes object for deformed plot
function [deformation_plot] = plotDeformation(elem,Nxy,U)

% Get total number of elements to iterate over
[row,~] = size(elem);
% Magnify deformation for plotting
def = Nxy+fliplr(U)*100;

deformation_plot = figure;

% Loop through and plot each deformed and undeformed element iteratively
for i = 1:row
    pgon(i) = polyshape(Nxy(elem(i,:),2),Nxy(elem(i,:),1));
    pgon_plot = plot(pgon(i),'linestyle','--','facealpha',0.1,'facecolor',[0.5 0.5 0.5],'linewidth',1,'alignvertexcenters','on');
    hold on
    scatter(Nxy(elem(i,:),2),Nxy(elem(i,:),1),200,'k','.');
    hold on
    pgon_def(i) = polyshape(def(elem(i,:),2),def(elem(i,:),1));
    pgon_def_plot = plot(pgon_def(i),'facealpha',0.05,'facecolor','r','edgecolor',[1 0 0],'linewidth',1);
    hold on
    scatter(def(elem(i,:),2),def(elem(i,:),1),200,'r','.');
    hold on
end

% Add legend and axes titles for plot
legend([pgon_plot pgon_def_plot],{'Undeformed','Deformed X100'})
xlabel('X (metres)')
ylabel('Y (metres)')
xlim([-5 55]*1e-3)
ylim([-5 40]*1e-3)
end