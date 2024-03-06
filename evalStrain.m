%-------------------------------------------------------------------------
%                      evalStrain.m
%                      Michael Dodman
%                      January 2022
%
% This function plots the deformed and undeformed plate overlayed on each
% other with a 100X displacement magnification for clarity.
%-------------------------------------------------------------------------
% Variable Name      Definition
% ------------------------Input--------------------------------------------
% elem               Matrix of node indices of each element
% Be                 [Y,X] matrix of node positions
% U                  Nodal Displacement Matrix
%-------------------------Output-------------------------------------------
% epsil              Element strain matrix
%-------------------------Internal----------------------------------------
% row                Number of elements
function [epsil] = evalStrain(elem,Be,U)

% Get number of elements in model
[row,~] = size(elem);
% Preallocate space for epsil matrix for efficiency
epsil = zeros(row,row);
% Calculate element strain matrices
for i = 1:row
    epsil(i,:) = Be(:,:,i)*reshape(U(elem(i,:),:).',[],1);
end
end