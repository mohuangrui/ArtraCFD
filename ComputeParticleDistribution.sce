//**************** Particle Distribution Configuration ***********

// the center and radius of the trace sphere
xCenter = 3;
yCenter = 3;
zCenter = 0;
R = 1;

// the radius of particles
// calculate radius by specifying the total number of particles
nParticles = 32;
r = R * sind(180 / nParticles);
// explicitly specify the radius instead, however non-zero gap
//r = 0.2;

//****************** Accumulation Calculation ******************
// this part is needed if accumulate new particles

nnParticles = 48;
// compute the new distribution radius
RR = (R + r) / (1 - sind(180 / nnParticles));
// compute the new particle radius
rr = RR * sind(180 / nnParticles);


nnnParticles = 96;
// compute the new distribution radius
RRR = (RR + rr) / (1 - sind(180 / nnnParticles));
// compute the new particle radius
rrr = RRR * sind(180 / nnnParticles);
// update variables
r = rrr;
R = RRR;

// ********************* Execute Part ***************************
// open the target file
FID=file('open',"Particle",'unknow');

// the angle gap of inclination angle between two contacted particles
phi = 2 * asind(r / R); // 3D
phi = 90; // 2D

// get total number of particles allowed on longitude circle
// note that inclination range is from 0 to 180 degree
nLongitude = int(180 / phi);

// adjust the phi according to n to ensure even spaced particles
phi = 180 / nLongitude;

// now generate the coordinates of particles for each lattitude circle
for n = 0:nLongitude
    if n == 0 | n == nLongitude then
        nLattitude = 1;
        xArray = xCenter + 0;
        yArray = yCenter + 0;
        zArray = zCenter + R * cosd(n * phi);
    else
        rLattitude = R * sind(n * phi);
        theta = 2 * asind(r / rLattitude);
        nLattitude = int(360 / theta);
        theta = 360 / nLattitude;
        // get azimuthal angle array
        thetaArray = (1:nLattitude)' * theta;
        // get the x, y, z coordinates array
        xArray = xCenter + rLattitude * cosd(thetaArray);
        yArray = yCenter + rLattitude * sind(thetaArray);
        zArray = zCenter + R * cosd(n * phi) * ones(nLattitude,1);
    end
    Particle = []; // initialize an empty matrix
    Particle(:,7) = zeros(nLattitude,1);
    Particle(:,6) = zeros(nLattitude,1);
    Particle(:,5) = zeros(nLattitude,1);
    Particle(:,4) = r * ones(nLattitude,1);
    Particle(:,3) = zArray;
    Particle(:,2) = yArray;
    Particle(:,1) = xArray;
    // write data out
    write(FID,Particle,'(50(f20.6))');
end

// close file
file('close',FID);

