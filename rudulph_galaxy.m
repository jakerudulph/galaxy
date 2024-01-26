%A SIMPLE GALAXY

%SET UP INITIAL CONDITIONS AND CONSTANTS
clear 
n=input('Enter the number of stars in your galaxy. A wise number is between 50 and 500. ')
sm = 1.9891E30; % KILOGRAMS 
G = 6.672E-11; % Nm^2/kg^2 (GRAVITATIONAL CONSTANT)
sol = 299792458; % SPEED OF LIGHT, IN M/S 
secondsyear = 365*24*60*60;
ly = sol*secondsyear; %ONE LIGHTYEAR
parsec = 3.26*ly;
gr=15000; %GALAXY RADIUS IN parsecs
bhm=1.73E11; %BLACK HOLE MASS
bhm=bhm*sm; %SET IN UNITS OF SOLAR MASS
origin=0; %GALACTIC CENTER ORIGIN FOR BOTH X AND Y
tstep=2000000; %OUR STEP TIME IN MILLIONS OF YEARS
tinc = tstep*secondsyear; %THIS TIME IN MILLIONS OF YEARS*SECONDS/YR (SPEED OF LIGHT IN M/S)
by=5; %NUMBER OF BILLIONS OF YEARS TO RUN THE MODEL
ntime=by*10^9/tstep; %NUMBER OF TIME STEPS IN OUR BILLIONS OF YEARS

rO = gr*parsec; %THE OUTER BOUNDARY FOR GENERATING STARS IN PARSECS
rIn = (gr/6)*parsec; %IF OUR GALACTIC CENTER ENCOMPASSES STARS THAT ARE ORBITING FAST 
%NEAR THE CENTER, THEY MAY AS WELL BE CONSTITUENTS OF THE BLACK HOLE IN THE CENTER. THIS IS
%THE THRESHHOLD

%LET'S SET UP THE GALAXY CENTER. THE FIRST ROW WILL HOLD THE CENTER OF MASS
%VALUES: COLUMN 1 IS MASS, 2 IS X POSITION, 3 IS Y POSITION, 4 IS XV, 5 IS 
%YV, AND 6 IS DISTANCE FROM GALACTIC CENTER. WE'LL USE THIS MATRIX FOR
%STARS AS WELL.

masses = zeros(n,6);
masses(1,1) = bhm;
masses(1,2) = origin*parsec; %METERS
masses(1,3) = origin*parsec; %SAME
masses(1,4) = 0;
masses(1,5) = 0;
masses(1,6) = 0;


for i = 2:n, %FOR ALL MASSES WHICH ARE NOT THE BLACK HOLE (FIRST ROW)
    %CALCULATE JITTERED MASSES FOR STARS
    %DIDN'T HAVE THE STATISTICS TOOLBOX, OR WOULD HAVE DONE A NONCENTRAL
    %CHI^2 DISTRIBUTION OR POISSON
    mi = rand();
    m = (mi+.8)^3*sm; %NORMALIZE SO SOME STARS LARGER, SOME SMALLER THAN OUR SUN
    masses(i,1) = m;
    %CALCULATE DISTANCE FROM CENTER FOR STARS
    ri = rand();
    r = ri*(rO - rIn) + rIn;
    
    ang0 = 0; %CALCULATE ANGLE FOR GENERATION
    ang0 = rand();
    ang = ang0*2*pi;
    x = r*cos(ang);
    y = r*sin(ang);
    dx = cos(ang+pi/2);
    dy = sin(ang+pi/2);
    
    % COMPUTE GRAVITATIONAL ATTRACTION
    v = sqrt(G*bhm/r);
    masses(i,2) = x+origin*parsec; %METERS 
    masses(i,3) = y+origin*parsec; %SAME
    masses(i,4) = dx*v;
    masses(i,5) = dy*v;
    masses(i,6) = r;
    
    %MAKE A DUMMY SO WE CAN CHECK LATER
    masses1=masses;
    
end
subplot(2,2,1)
scatter(masses(:,2),masses(:,3),30,'b','filled') %CREATE SMALL BLUE STARS
hold on
scatter(masses(1,2),masses(1,3),'r')
axis([-2*gr*parsec 2*gr*parsec -2*gr*parsec 2*gr*parsec])
axis square
hold off
title('Initial Positions, Galactic Center in Red')
xlabel('X Plane, Meters')
ylabel('Y plane, Meters')

subplot(2,2,2)
hist(masses(2:n,1)/sm)
title('Histogram of Star Masses')
xlabel('Stellar Mass Equivalent')
ylabel('Frequency')

for j=1:tinc 
    %CALCULATING GRAVITATION ACCELERATION FOR EACH STAR
    %WE'LL START WITH EACH STAR TO ALL OTHER STARS
    
    for k=1:n
        
        %CALCULATE DISTANCES TO EACH OTHER STAR
        xd=masses(k,2); %EXTRACT X AND Y POSITIONS
        yd=masses(k,3);
        acc=[0 0]; %AN EMPTY MATRIX FOR XY ACCELERATIONS
        dr2 = xd^2+yd^2; 
        dstars=ones(n,2);
        massk=masses(1:n,1:3); %MASS AND XY COORDINATES FOR EACH STAR
        for q=2:n %USE THIS FOR LOOP TO CALCULATE AGGREGATE ACCELERATION FROM OTHER STARS
            massk(q,4)=massk(q,2)-xd;
            massk(q,5)=massk(q,3)-yd;
            massk(q,6)=(massk(q,4)/sqrt(dr2))*(G*(massk(q,1)))/dr2;
            massk(q,7)=(massk(q,5)/sqrt(dr2))*(G*(massk(q,1)))/dr2;
        end
        starG = [sum(massk(:,6)) sum(massk(:,7))];
        dcenter=[origin origin] - [xd yd]; %DISTANCE FROM GALACTIC CENTER
        
        acc = acc + starG + (dcenter/sqrt(dr2))*((G*(bhm+masses(k,1)))/dr2); 
        %GRAVITATIONAL ACCELERATION FOR X AND Y TO CENTER, MINDING AGGREGATE ACCELERATION FROM OTHER STARS
        masses(k,4)=masses(k,4)+acc(1)*tinc;
        masses(k,5)=masses(k,5)+acc(2)*tinc;
        masses(k,2) = masses(k,2) + masses(k,4)*tinc; 
        masses(k,3) = masses(k,3) + masses(k,5)*tinc;
        masses(k,6) = sqrt(dcenter(1)^2 + dcenter(2)^2);
    end

subplot(2,2,4) %PLOT MOTION
whitebg('k')
scatter(origin, origin,'r')
hold on
scatter(masses(:,2),masses(:,3),30,'w','filled') %CREATE SMALL WHITE STARS
title('Galactic Rotation')
xlabel('X Plane, Meters')
ylabel('Y plane, Meters')
axis([-2*gr*parsec 2*gr*parsec -2*gr*parsec 2*gr*parsec])
axis square
pause(0.05)
hold off

%PLOT GALAXY ROTATION CURVE
velocity=sqrt(masses(:,4).^2+(masses(:,5).^2));
%WHERE V = SQRT(ACCEL*DISTANCE FROM CENTER)
subplot(2,2,3)
scatter(masses(:,6)/ly,velocity)
title('Galaxy Rotation Curve')
xlabel('Distance from Center, Lightyears')  
ylabel('Velocity, m/s')
axis ([0.5E4 6E4 0 6E5])
hold off

end
