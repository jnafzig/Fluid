%% This script solves fluid dynamics equations in one dimension.  Current
% setup is for the sod shock tube, which is a standard test for such
% solvers.  run plotscript.m afterwards to view results.

%% clear and start clock

USEPREVIOUS = false;
if USEPREVIOUS == false
    clearvars -except USEPREVIOUS;
end
tic;

%% Settings

% allow backflow
BACKFLOW = true;

% number of time steps to be recorded
NumRecTimeSteps = 100;
NumTimeStepsPerRec = 200;

% number of spatial steps to be recorded
NumRecSpaceSteps = 500;
NumSpaceStepsPerRec = 1;

% spatial limits
a = 0;
b = 1;

% physical constants
gamma = 1.4;

%% Grid Setup

% total number of grid points in time
Nt = NumRecTimeSteps*NumTimeStepsPerRec;

% total number of grid points in space
Nx = NumRecSpaceSteps*NumSpaceStepsPerRec;

% grid spacings
dx = (b-a)/Nx;
dt = dx/100;

% cell wall grid points
Xn = linspace(a,b,Nx+1)';

% cell centered grid points
X = (Xn(1:end-1) + Xn(2:end))/2;

% temporal grid points for calculation
T = linspace(0,dt*Nt,Nt);

%% Solar Geometry

% % area at each cell wall
% area = Xn.^2;
% 
% % volume of each shell
% volume = 1/2*(area(1:end-1)+area(2:end))*dx;
% 
% % gravitational potential (at the cell walls)
% phi = -1./Xn;

%% De Laval Nozzle Geometry

% intake = 1.0;
% throat = .7;
% outlet = 1.1;
% 
% env_area = linspace(intake,outlet,Nx+1);
% 
% x_nozzle = [linspace(0,pi/2,Nx/4),linspace(pi/2,pi,3*Nx/4+1)];
% 
% % area at each cell wall
% area = [env_area+(throat-env_area).*sin(x_nozzle).^2]';
% 
% % volume element of each cell
% volume = 1/2*(area(1:end-1)+area(2:end))*dx;
% 
% % gravitational potential (at the cell walls)
% phi = zeros(Nx+1,1);
% 

%% Sod Shock Tube Geometry

% area at each cell wall
area = ones(Nx+1,1);

% volume element of each cell
volume = 1/2*(area(1:end-1)+area(2:end))*dx;

% gravitational potential (at the cell walls)
phi = zeros(Nx+1,1);

%% Grid for recordings

% spatial grid points for recording
x = X(1:NumSpaceStepsPerRec:Nx);

% temporal grid points for recording
t = T(1:NumTimeStepsPerRec:Nt);

%% Initial and Boundary conditions
if USEPREVIOUS==false

% initialization
density = .35*ones(Nx,1); 
velocity = .5*ones(Nx,1);
pressure = .1*ones(Nx,1);

% sod shock tube: 
interface = .5;
density(X<=interface) = 1;
velocity(X<=interface) = 0;
pressure(X<=interface) = 1;

density(X>interface) = .125;
velocity(X>interface) = 0;
pressure(X>interface) = .1;

% %sod shock tube reversed:
% interface = .5;
% density(X<=interface) = .125;
% velocity(X<=interface) = 0;
% pressure(X<=interface) = .1;
% 
% density(X>interface) = 1;
% velocity(X>interface) = 0;
% pressure(X>interface) = 1;

% speed of sound 
cs = (gamma*pressure./density).^.5;
% maximum speed of sound
csmax = max(cs);

% momentum and energy densities
momentum = velocity.*density;
energy = pressure/(gamma-1) ...
    + 1/2*velocity.*momentum ...
    + density.*(phi(1:end-1)+phi(2:end))/2;

% left and right boundary conditions

% incoming flux
F1bcL = momentum(1);                                % mass flux
F2bcL = velocity(1).*momentum(1);                   % momentum flux
F3bcL = velocity(1).*(energy(1)+pressure(1));       % energy flux

% pressure and velocity at left cell
pressureL = pressure(1);
velocityL = velocity(1);

% incoming flux from right cell (only if we have infalling matter)
F1bcR = momentum(end);                                % mass flux
F2bcR = velocity(end).*momentum(end);                   % momentum flux
F3bcR = velocity(end).*(energy(end)+pressure(end));       % energy flux

% pressure and velocity at right boundary cell
pressureR = 0;
velocityR = velocity(end);

else
%% USEPREVIOUS=true --> so change incoming flux?

densityL = .5;
velocityL = 1.0;
pressureL = 1.0;

momentumL = velocityL.*densityL;
energyL = pressureL/(gamma-1) ...
    + 1/2*velocityL.*momentumL ...
    + densityL.*phi(1);

% incoming flux at left edge
F1bcL = momentumL;   
F2bcL = velocityL.*momentumL;                   % momentum flux
F3bcL = velocityL.*(energyL+pressureL);       % energy flux
    
end

%% Data storage allocation

% conserved quantities to be recorded
TotalMassRec = zeros(1,Nt);
TotalMomentumRec = zeros(1,Nt);
TotalEnergyRec = zeros(1,Nt);

% energies
TotalInternalRec = zeros(1,Nt);
TotalKineticRec = zeros(1,Nt);
TotalPotentialRec = zeros(1,Nt);

% conserved quantities in and out of the system
MassIN = zeros(1,Nt);
MassOUT = zeros(1,Nt);

MomentumIN = zeros(1,Nt);
MomentumOUT = zeros(1,Nt);
MomentumSource = zeros(1,Nt);

EnergyIN = zeros(1,Nt);
EnergyOUT = zeros(1,Nt);

% conserved variables to be recorded
densityRec = zeros(NumRecSpaceSteps,NumRecTimeSteps);
momentumRec = zeros(NumRecSpaceSteps,NumRecTimeSteps);
energyRec = zeros(NumRecSpaceSteps,NumRecTimeSteps);

% dependent quantities
pressureRec = zeros(NumRecSpaceSteps,NumRecTimeSteps);
velocityRec = zeros(NumRecSpaceSteps,NumRecTimeSteps);

% fluxes to be recorded
F1Rec = zeros(NumRecSpaceSteps+1,NumRecTimeSteps);
F2Rec = zeros(NumRecSpaceSteps+1,NumRecTimeSteps);
F3Rec = zeros(NumRecSpaceSteps+1,NumRecTimeSteps);

pressurebcRec = zeros(1,NumRecTimeSteps);
velocitybcRec = zeros(1,NumRecTimeSteps);
densitybcRec = zeros(1,NumRecTimeSteps);
momentumbcRec = zeros(1,NumRecTimeSteps);
energybcRec = zeros(1,NumRecTimeSteps);

%% Print output

fprintf('  Time Step   Percent Complete   \n');
fprintf('----------------------------------\n');

%% loop
for i = 1:Nt

    if (BACKFLOW==false)  % then use the standard scheme
%% Calculate fluxes     

    % flux of mass
    F1 = [F1bcL;velocity.*density]; 
    
    % flux of momentum
    F2 = [F2bcL;velocity.*momentum];
    
    % flux of energy
    F3 = [F3bcL;velocity.*(energy+pressure)];
    
    % note: total fluxes are on the cell walls so 
    % they end up size Nx+1
    
%% Calculate change in conserved variables

    % note: diff reduces vector dimension by 1 so result is size Nx
    % change in density
    D_density = - dt*diff(area.*F1)./volume;

    % change in momentum density: convective flux + force
    D_momentum = - dt*((diff(area.*F2))./volume ...
                        + diff([pressure;pressureR])/dx ...
                        + density.*diff(phi)/dx); 

    % change in energy density: convective flux + work done
    D_energy = - dt*(diff(area.*F3))./volume;
    
    else
% BACKFLOW = true, so break calculation of fluxes
% into right flowing and left flowing

%% Calculate fluxes for flow to the right     
% (depends on state of cell upwind)

    % flux of mass
    F1R = [F1bcL;velocity.*density]; 
    
    % flux of momentum
    F2R = [F2bcL;velocity.*momentum];
    
    % flux of energy
    F3R = [F3bcL;velocity.*(energy+pressure)];
    
    % note: total fluxes are on the cell walls so 
    % they end up size Nx+1
    
%% Calculate fluxes for flow to the left 
% (depends on state of cell upwind)

    % flux of mass
    F1L = [velocity.*density;F1bcR]; 
    
    % flux of momentum
    F2L = [velocity.*momentum;F2bcR];
    
    % flux of energy
    F3L = [velocity.*(energy+pressure);F3bcR];
    
    % note: total fluxes are on the cell walls so 
    % they end up size Nx+1
    
%% Combine left and right fluxes depending on velocity
    
    F1R([velocityL;velocity]<0) = 0; 
    F2R([velocityL;velocity]<0) = 0;
    F3R([velocityL;velocity]<0) = 0;  
    
    F1L([velocity;velocityR]>0) = 0; 
    F2L([velocity;velocityR]>0) = 0;
    F3L([velocity;velocityR]>0) = 0;  
    
%% For Reflective boundary conditions        
    F1R(end) = 0;
    F2R(end) = 0;
    F3R(end) = 0;
    
    F1L(1) = 0;
    F2L(1) = 0;
    F3L(1) = 0;
 
    F1 = F1R+F1L;
    F2 = F2R+F2L;
    F3 = F3R+F3L;
   
    % note: total fluxes are on the cell walls so 
    % they end up size Nx+1
    
%% Find the force on each cell. 
% (depends on pressure of cell downwind)

    % forces due to pressure gradients
    ForceR = diff([pressure;pressureR]); % force if cell is moving right 
    ForceL = diff([pressureL;pressure]); % force if cell is moving left
        
    ForceR(velocity<0) = 0; 
    ForceL(velocity>0) = 0;
 
    % total force is due to pressure gradient + gravity;
    Force = ForceR + ForceL + density.*diff(phi);
    
    % reflecting boundary conditions
    Force(end) = 0.5*(pressureR-pressure(end-1));
    Force(1) = 0.5*(pressure(2)-pressureL);
    
    % note: force acts on cell midpoint so size is Nx
 
%% Calculate change in conserved variables

    % note: diff reduces vector dimension by 1 so result is size Nx
    % change in density
    D_density = - dt*diff(area.*F1)./volume;

    % change in momentum density: convective flux + force
    D_momentum = - dt*((diff(area.*F2))./volume ...
                        + Force/dx); 

    % change in energy density: convective flux + work done
    D_energy = - dt*(diff(area.*F3))./volume;
          
    end

%% Record variables
    if mod(i-1,NumTimeStepsPerRec)==0
        % Print status output
        fprintf('%9.0f%15.1f\n',i,i/Nt*100);        

        index = ceil(i/NumTimeStepsPerRec);
        
        % conserved variables to be recorded
        densityRec(:,index) = density(1:NumSpaceStepsPerRec:Nx);
        momentumRec(:,index) = momentum(1:NumSpaceStepsPerRec:Nx);
        energyRec(:,index) = energy(1:NumSpaceStepsPerRec:Nx);

        % dependent quantities
        pressureRec(:,index) = pressure(1:NumSpaceStepsPerRec:Nx);
        velocityRec(:,index) = velocity(1:NumSpaceStepsPerRec:Nx);
        
        % fluxes
        F1Rec(:,index) = F1(1:NumSpaceStepsPerRec:Nx+1);
        F2Rec(:,index) = F2(1:NumSpaceStepsPerRec:Nx+1);
        F3Rec(:,index) = F3(1:NumSpaceStepsPerRec:Nx+1);

    end
    
    % conserved quantities to be recorded
    TotalMassRec(i) = sum(density.*volume);
    TotalMomentumRec(i) = sum(momentum.*volume);
    TotalEnergyRec(i) = sum(energy.*volume);
    
    % energy components
    TotalKineticRec(i) = 1/2*sum(velocity.*momentum.*volume);
    TotalInternalRec(i) = 1/(gamma-1)*sum(volume.*pressure);
    TotalPotentialRec(i) = sum(volume.*density...
                                .*(phi(1:end-1)+phi(2:end))/2);
    
    % inflow and outlfow
    MassIN(i) = dt*area(1)*F1(1);
    MassOUT(i) = dt*area(Nx+1)*F1(Nx+1);

    EnergyIN(i) = dt*area(1)*F3(1);
    EnergyOUT(i) = dt*area(Nx+1)*F3(Nx+1);
    
    MomentumIN(i) = dt*area(1)*F2(1);
    MomentumOUT(i) = dt*area(Nx+1)*F2(Nx+1);
    MomentumSource(i) = sum( - dt*( ...
                    + diff([pressure;pressureR])/dx ...
                    + density.*diff(phi)/dx).*volume );

%% Update conserved quantities

    density = density + D_density;
    momentum = momentum + D_momentum;
    energy = energy + D_energy;
    
%% Update boundary cell conditions

    % pressure and velocity at left cell
    pressureL = pressure(1);
    velocityL = velocity(1);

    % pressure and velocity at right boundary cell
    pressureR = pressure(end);
    velocityR = velocity(end);
   
%% Update other variables
    
    if (BACKFLOW == false)
    % then don't try to simulate infalling matter or it'll blow up.
    momentum(momentum<0)=0;
    end
    
    velocity = momentum./density;
    
    % equation of state:
    pressure = (gamma - 1)*(energy ...
        - 1/2*velocity.*momentum ...
        - density.*(phi(1:end-1)+phi(2:end))/2);
    
       
end

toc;
