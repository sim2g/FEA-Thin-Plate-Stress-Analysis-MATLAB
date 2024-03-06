%-------------------------------------------------------------------------
%                      getKUmatrices.m
%                      Michael Dodman
%                      January 2022
%
%  This function calculates the global element stiffness matrices, nodal
%  displacements in X and Y, and the global strain matrix.
%-------------------------------------------------------------------------
% Variable Name      Definition
% ------------------------Input--------------------------------------------
% Nxy                [Y,X] matrix of node positions
% elem               Matrix of node indices of each element
% n_node             Number of nodes per element
% total_node         Total number of nodes in system
% load_matrix        [X,Y] matrix of load at each node
% load_nodes         Array containing indices of nodes with applied load
% H                  Length of applied uniform tension
%-------------------------Output-------------------------------------------
% K                  Global Element Stiffness Matrix
% U                  Nodal Displacement Matrix
% Be                 Global Strain Matrix
%--------------------------Internal----------------------------------------
% row                Number of rows in matrix elem
% col                Number of columns in matrix elem
% X                  Matrix of nodal X positions indexed with elem
% Y                  Matrix of nodal Y positions indexed with elem
% A                  Local strain matrix A
% B                  Local strain matrix B
% C                  Local strain matrix C
% J                  Jacobian Matrix
% BJ                 Strain Matrix with same dimensions as J
% Bx                 Reshaped Strain Matrix
% D                  Plane Stress Matrix
% L                  Jacobian Credit Points
% Hi                 Weight coefficients of each element
% XY                 Concatenation of X and Y into single matrix
% Ke                 Local Element Stiffness Matrix
% [p,q,r]            Loop counters for K matrix assembly
% [a,b,c,d]          K matrix assembly indices
function [K,U,Be] = getKUmatrices(Nxy,elem,n_node,total_node,load_matrix,load_nodes,H)

% Define global variables inside function
global E
global nu
global thickness
global DoF
global theta
global P

% Get elem matrix size
[row,col] = size(elem);

% Preallocate matrices which change size during loop for efficiency
X = zeros(row,col);
Y = zeros(row,col);
A = zeros(row,col);
B = zeros(row,col);
C = zeros(row,col);
Be = zeros(row,2*col,row);
Ke = zeros(col*2,col*2,row);
K = zeros(total_node*DoF,total_node*DoF);

J = zeros(DoF,n_node,row);
BJ = J;
Bx = zeros(row,DoF*n_node);

% Define plane stress matrix
D =[1 nu 0;
    nu 1 0;
    0 0 (1-nu)/2]*E/(1-nu^2);

% Define Jacobian credit points and weight coefficients
L = [0.5 0.5 0;
    0.5 0 0.5];
Hi = [1/6 1/6 1/6];

% Loop through K, Be and U calculation for each element
for i = 1: row
    
    % Get nodal X and Y coordinates for element
    for j = 1:col
        X(i,j) = Nxy(elem(i,j),2);
        Y(i,j) = Nxy(elem(i,j),1);
    end
    
    if n_node == 3
        
        % Calculate strain and stiffness matrices
        A(i,:) = [X(i,2)*Y(i,3)-X(i,3)*Y(i,2),X(i,3)*Y(i,1)-X(i,1)*Y(i,3),X(i,1)*Y(i,2)-X(i,2)*Y(i,1)];
        B(i,:) = [Y(i,2)-Y(i,3),Y(i,3)-Y(i,1),Y(i,1)-Y(i,2)];
        C(i,:) = [X(i,3)-X(i,2),X(i,1)-X(i,3),X(i,2)-X(i,1)];
        
        Be(:,:,i) = [B(i,1) 0 B(i,2) 0 B(i,3) 0;
            0 C(i,1) 0 C(i,2) 0 C(i,3);
            C(i,1) B(i,1) C(i,2) B(i,2) C(i,3) B(i,3)]/abs(sum(A(i,1:3)));
        
        
        Ke(:,:,i) = (Be(:,:,i)'*D*Be(:,:,i)*thickness)*(1/2*abs(sum(A(i,1:3))));
        
    elseif n_node == 6
        
        XY = [X(i,:);Y(i,:)]';
        
        % Calculate Jacobian, strain matrices and stiffness matrices
        % iteratively
        for k = 1:row
            J(:,:,k) = [4*L(1,k)-1,4*L(2,k),0,-4*L(2,k),4*L(1,k)+4*L(2,k)-3,4-4*L(2,k)-8*L(1,k);
                0,4*L(1,k),4*L(2,k)-1,4-8*L(2,k)-4*L(1,k),4*L(1,k)+4*L(2,k)-3,-4*L(1,k)];
            
            BJ(:,:,k) = inv(J(:,:,k)*XY)*J(:,:,k);
            
            Bx(:,:,k) = [reshape([BJ(1,:,k); zeros(1,n_node)],1,[]);
                reshape([zeros(1,n_node); BJ(2,:,k)],1,[]);
                reshape([BJ(2,:,k); BJ(1,:,k)],1,[])];
            
            Ke(:,:,i) = Ke(:,:,i) + thickness*Hi(i)*Bx(:,:,k)'*D*Bx(:,:,k)*det(J(:,:,k)*XY);
            
            Be = Be + Bx(:,:,i);
            
        end
        
    end
end

% Assemble global stiffness matrix from local stiffness matrices
for p = 1:size(elem(:,1))
    for q = 1:n_node
        for r = 1:n_node
            a = (q-1)*2+(1:2);
            b = (r-1)*2+(1:2);
            c = (elem(p,q)-1)*2+(1:2);
            d = (elem(p,r)-1)*2+(1:2);
            K(c,d) = K(c,d) + Ke(a,b,p);
        end
    end
end

% Boundary conditions for simple and hinged supports

K(1:3,:)=0;
K(:,1:3)=0;
K(1:3,1:3)=eye(3);

% Calculate load components in X and Y
if n_node == 3
    load_matrix(load_nodes,1) = H*thickness/2*P*sind(theta);
    load_matrix(load_nodes,2) = H*thickness/2*P*cosd(theta);
elseif n_node == 6
    load_matrix(load_nodes(1:2),1) = H*P*thickness/6*sind(theta);
    load_matrix(load_nodes(1:2),2) = H*P*thickness/6*cosd(theta);
    load_matrix(load_nodes(3),1) = H*P*thickness*2/3*sind(theta);
    load_matrix(load_nodes(3),2) = H*P*thickness*2/3*cosd(theta);
end

% Calculate nodal displacements and display as X and Y columns
U = inv(K)*reshape(load_matrix.',[total_node*DoF,1]);
U = reshape(U,[],total_node).';
end