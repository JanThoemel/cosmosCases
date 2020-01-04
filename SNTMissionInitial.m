function  [sstTemp,ns,altitude,panels,rho,Tatmos,v,radiusOfEarth,MeanMotion,mu,satelliteMass,panelSurface,...
          sstDesiredFunction,windOn,sunOn,deltaAngle,timetemp,totalTime,wakeAerodynamics,masterSatellite,...
          SSCoeff,inclination,SSParameters,meanAnomalyOffSet]=SNTMissionInitial(altitude,sunOn)
%% initial conditions for SNTMission Mission

%% simulation time constants
  totalTime         =120*90*60;   %% simulation period (approximate multiples of orbit periods),[s]
  compStep          =3;          %% computational step size in [s]
  lengthControlLoop =90;         %% in [s]
  timetemp  =0:compStep:lengthControlLoop; %% duration and interpolation timestep for each control loop
  %fprintf('\nnumber of float variables for each control loop: %d, size %f kbyte',size(timetemp,2)*(1+6+6+3)+6+1,(size(timetemp,2)*(1+6+6+3)+6+1)*4/1024); %% size of time vector, state vector, desired state vector and anglevector, C and MeanMotion of analytical solution
  wakeAerodynamics=0;             %% use of wake aerodynamics
  masterSatellite=0;              %% if 0 no master, 1 active master, 2 passive master,

  sstDesiredFunction=@SNTMissionDesired;
  windOn       =1;

  if not(altitude)
    sunOn       =1;
    altitude    =453e3;                 %% [m] 
  end
  J2On         =0; 

  if sunOn
    fprintf('\nonly dusk-dawn orbits for the time being! Arbitrary orbits require rotation of sunlight vector and accounting for eclipse');
  end
  deltaAngle   =30;                   %% roll,pitch,yaw angle resolution

  ns=2;
  satelliteMass=2;
  
  argumentOfPerigeeAtTe0=0;       %% not used yet
  trueAnomalyAtTe0=0;             %% not used yet
  ejectionVelocity=1;          %% [m/s]
  timeBetweenEjections=0.01;      %% [s]
  %startSecondPhase=8*2*pi/MeanMotion;  %% in s
  panelSurface=0.01;              %% [m^2]
  %panels=[1 1.5 4.5]; 
  panels=[0 0 3.5]; 

  %% other constants
  [rho,Tatmos,v,radiusOfEarth    ,mu,MeanMotion,SSOinclination,J2]=orbitalproperties(altitude);
  
  
  r0=radiusOfEarth+altitude;    %% in m

  inclination=SSOinclination;   %%
  %inclination=90               %% manual inclination if 90 or 0 deg SSCoeff is 1, no J2 effects are experienced
  
  SSCoeff=sqrt(1+3*J2*radiusOfEarth^2/8/r0^2*(1+3*cosd(2*inclination)))^J2On; %% J2On is 0, i.e. the effect shall not be taking into account then SSCoeff is zero
   
  %% initial conditions
  sstTemp=zeros(9,ns,size(timetemp,2));
  for i=1:ns
    sstTemp(1,i,1)=(i-1)*ejectionVelocity*timeBetweenEjections; %% x position
    sstTemp(4,i,1)=0;           %% u velocity
    sstTemp(7,i,1)=0;           %% alpha
    sstTemp(8,i,1)=0;           %% beta
    sstTemp(9,i,1)=0;           %% gamma
  end
  
  %% parameters for desired conditions  
  meanAnomalyOffSet=0;% %% for 0 satellite cross on the poles; for pi/2 satellite crosses at equator
  numberOfModes=10;
  SSParameters=zeros(8,ns,numberOfModes);
  AA=500; DD=0;

  %% generalized analytical solution, from Traub with variable renaming and accounting for different coordinate system, accelerations do not belong here
  SSParameters(1,2,1)=0;            %%
  SSParameters(2,2,1)=DD;        %%
  SSParameters(3,2,1)=0;            %%
  SSParameters(4,2,1)=0;            %%
  SSParameters(5,2,1)=0;            %%
  SSParameters(6,2,1)=0;            %%
  SSParameters(7,2,1)=0;            %%
  SSParameters(8,2,1)=0;            %%  

  SSParameters(1,3,1)=AA;         %% xmax (beta)
  SSParameters(2,3,1)=-DD;            %% xpermanent.only here, could assume any value
  SSParameters(3,3,1)=AA;         %% ymax
  SSParameters(4,3,1)=0;            %% ymaxdot0
  SSParameters(5,3,1)=0;            %% zmax (alpha)
  SSParameters(6,3,1)=0;            %% zpermanent, it is a permanent nadir-zenit movement should always be zero because it implies a permant ram-antiram movement. should this ever be non-zero, all has to be verified
  SSParameters(7,3,1)=0;            %% xz offset (not mentioned in Traub but in Ivanov)
  SSParameters(8,3,1)=0;            %% y offset  (not mentioned in Traub but in Ivanov)

end

