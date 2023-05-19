%
%CrystalStability Generate random lattice configurations and perform local stability analysis.
%   Note: it does not run simulations.
%
%   See also: Launcher
%   
%   Authors:    Andrea Giusti
%   Date:       2023
%

close all
clear
clc

%% Parameters

% directory to save the results
outputDir='/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/stability of geometric lattices/simulations';

D=3;                        % number of dimensions [2 or 3]

Nagents = 25:100;           % range of the number of agents in the lattice
Ntimes=10;                  % How many simulations are launched

seed=0;                     % seed for random generator, if negative it is not set

LinkNumber=6*(D-1);         %number of links per agent in the lattice configuration (L)
                            %If D=2 then 6=triangular lattice, 4=square lattice, 3=hexagonal lattice
                            %If D=3 then 6=cubic lattice, 12=thetradic-octaedric lattice

Rmax= (sqrt(5-D)+1)/2;      % maximum lenght of a link (R_a). Must be in [1; Rnext]

% interaction function
syms x
a=0.5; c=24; f(x)= a/x^(c*2)-a/x^c %Lennard-Jones
%f(x)= -(x-1) %linear

epsilon = 10e-9;            % value approximated to zero

%% Preallocate variables
Nconfig=length(Nagents);
rigidity = nan(Nconfig,Ntimes);
equilibrium = nan(Nconfig,Ntimes);
stability = nan(Nconfig,Ntimes);
asim_stability = nan(Nconfig,Ntimes);
rigid_motion = nan(Nconfig,Ntimes);
hyperbolic = nan(Nconfig,Ntimes);
lambda= nan(Nconfig,Ntimes,7);
XVec= cell(Nconfig,Ntimes);

if seed>=0
    rng(seed,'twister'); % reproducible results
end
    
figure
tiledlayout(4,6,'TileSpacing','Compact','Padding','Compact');

for conf=1:Nconfig
    
    N=Nagents(conf);
    lattice_size=(floor(nthroot(N,D)+1))^D;
    
    disp(['Evaluating ',num2str(Ntimes),' configurations with N=',num2str(N),' and lattice_size=',num2str(lattice_size)])
    
    for rep =1:Ntimes
        
        % generic lattice
        X=perfectLactice(N, LinkNumber, D,true,true,lattice_size, 5000);
        XVec{conf,rep}=X;
        
        % build all links
        links=buildLinks(X, Rmax, false);
        
        %% Build incidence matrix
        n=size(X,1);
        m=size(links,1);
        B = buildIncidenceMatrixFromLinks(links,n);
        
        %% Define Matrices and Useful Quantities
        R=B'*X;                         %realative positions
        L=B*B';                         % Laplacian
        A=L-diag(diag(L));              % adjacency matrix
        % Le=B'*B;                        % edge laplacian
        % Ae=(Le~=0)-eye(m);              % adjacency matrix of the graph of the edges
        M = buildRigidityMatrix(X, B);  % rigidity matrix
        
        xvec = reshape(X',[],1);    % stack vector of positions
        rvec = reshape(R',[],1);    % stack vector of relative positions
        
        d=vecnorm(R,2,2);           %vector of distances
        Rhat=R./d;                  %unit vectors of relative positions
        
        df=diff(f,x);               % derivative of f
        
        %cosines=Rhat*Rhat'; % scalar vectors between unit vectors
        
        %% Analyse the configuration
        % infinitesimal rigidity
        rigidity(conf,rep) = rank(M)==D*N-D*(D+1)/2;
        
        equilibrium(conf,rep) = all(f(d)==0);
        
        %% Compute the Jacobian
        % Jacobian of the system of the distances
        %J=(Le*diag(df(d))).*cosines
        
        % Jacobian of the system of the positions
        %dG=nan(m,m,2*n);
        J=sparse(D*n,D*n);
        gamma = eval((df(d).*d - f(d)).*d.^-3);
        diag_gamma=sparse(diag(gamma));
        for k=1:(D*n)
            dG=sparse(diag_gamma * diag(M(:,k)));
            J(:,k)=kron(B * dG * B', eye(D)) * xvec;
        end
        
        %J=jordan(J)
        
        %% Analysis
        % analyse the Jacobian
        %[V,lambda] =eig(J);
        [V,lam]=eigs(J,7,'largestreal','Tolerance',epsilon/2,'MaxIterations',10^6);
        lam=diag(lam);
        
        % study equilibrium
        lam=vpa(subs(lam),3);
        lambda(conf,rep,:)=lam;
        hyperbolic(conf,rep)= all((real(lam) < -epsilon) | (real(lam) > epsilon));
        stability(conf,rep)= all(real(lam) < epsilon);
        asim_stability(conf,rep)= all(real(lam) < -epsilon);
        
        % get basis of the center manifold
        V0=V(:, -epsilon<real(lam) & real(lam)<epsilon);
        
        % check if the center manifold correspons to rigid motions
        if rigidity(conf,rep)
            rigid_motion(conf,rep) = norm(M*V0)<epsilon;
        end
        
        % check data are correct
        if ~rigidity(conf,rep) || ~equilibrium(conf,rep) || ~stability(conf,rep) || ~rigid_motion(conf,rep) || asim_stability(conf,rep)
            warning('Conditions not satisfied.')
        end
        
    end
    
    % print a sample config for various number of agents
    if ismember(conf, round(linspace(1,Nconfig,24)))
        nexttile
        plotSwarmInit(X,N,0,1.2)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        drawnow
    end
end

%% Display Results
disp(['rigidity=       ',num2str(all(rigidity,'all')), ' (expected 1)'])
disp(['equilibrium=    ',num2str(all(equilibrium,'all')), ' (expected 1)'])
disp(['stability=      ',num2str(all(stability,'all')), ' (expected 1)'])
disp(['asim_stability= ',num2str(any(asim_stability,'all')), ' (expected 0)'])
disp(['rigid_motion=   ',num2str(all(rigid_motion,'all')), ' (expected 1)'])

% create folder, save data and parameters
if outputDir
    counter=1;
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),'crystalStability_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),'crystalStability_',num2str(counter)]);
    mkdir(path)
    disp('Saving data in ' + string(path))
    save(fullfile(path, 'data'))
    
    fileID = fopen(fullfile(path, 'parameters.txt'),'wt');
    fprintf(fileID,'crystalStability\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Parameters: \n');
    fprintf(fileID,'Ntimes= %d\n',Ntimes);
    fprintf(fileID,'Nagents= %s\n',mat2str(reshape(Nagents,1,[])));
    fprintf(fileID,'D= %d\n',D);
    fprintf(fileID,'f= %s\n',char(f));
    fprintf(fileID,'seed= %d\n',seed);
    fprintf(fileID,'\nResults:\n');
    fprintf(fileID,['rigidity=       ',num2str(all(rigidity,'all')), ' (expected 1)\n']);
    fprintf(fileID,['equilibrium=    ',num2str(all(equilibrium,'all')), ' (expected 1)\n']);
    fprintf(fileID,['stability=      ',num2str(all(stability,'all')), ' (expected 1)\n']);
    fprintf(fileID,['asim_stability= ',num2str(any(asim_stability,'all')), ' (expected 0)\n']);
    fprintf(fileID,['rigid_motion=   ',num2str(all(rigid_motion,'all')), ' (expected 1)\n']);
    fclose(fileID);
end

%% Plots
% figure
% plot(graph(A),'XData',X(:,1),'YData',X(:,2), 'EdgeLabel', [1:m],'MarkerSize',10)
% axis equal
% title('lattice')

% figure
% plot(graph(Ae),'MarkerSize',10)
% axis equal
% title('edge graph')

% figure
% plot(real(lam),imag(lam),'o')
% hold on
% plot(0,0,'r*')
% plot([0,0],[-100,100],'r--')
% grid
% title('poles')

% figure
% fplot(f,[0,2])
% hold on
% scatter(1,0)
% ylim([-inf 1])
% grid
% title('f(x)')
% axis([0,2,-0.5, 1.1])


