//
//  main.m
//  Golf Ball Sim
//
//  Created by Devin Shelly on 12/28/13.
//  Copyright (c) 2013 Devin Shelly. All rights reserved.
//

#import <stdlib.h>

#pragma mark - Vector3D

typedef struct Vector3D
{
    double x;
    double y;
    double z;
} Vector3D;

Vector3D Vector3DMake(double x, double y, double z)
{
    Vector3D vector;
    vector.x = x;
    vector.y = y;
    vector.z = z;
    return vector;
}

Vector3D Vector3DMultiply(Vector3D vector, double scalar)
{
    return Vector3DMake(vector.x * scalar, vector.y * scalar, vector.z * scalar);
}

Vector3D Vector3DAdd(Vector3D vectorA, Vector3D vectorB)
{
    return Vector3DMake(vectorA.x + vectorB.x, vectorA.y + vectorB.y, vectorA.z + vectorB.z);
}

Vector3D Vector3DSubtract(Vector3D vectorA, Vector3D vectorB)
{
    return Vector3DAdd(vectorA, Vector3DMultiply(vectorB, -1.0));
}

Vector3D Vector3DCross(Vector3D vectorA, Vector3D vectorB)
{
    return Vector3DMake(vectorA.y * vectorB.z - vectorA.z * vectorB.y, vectorA.z * vectorB.x - vectorA.x * vectorB.z, vectorA.x * vectorB.y - vectorA.y * vectorB.x);
}

double Vector3DDot(Vector3D vectorA, Vector3D vectorB)
{
    return vectorA.x * vectorB.x + vectorA.y * vectorB.y + vectorA.z * vectorB.z;
}

double Vector3DLengthSq(Vector3D vector)
{
    return Vector3DDot(vector, vector);
}

double Vector3DLength(Vector3D vector)
{
    return sqrt(Vector3DLengthSq(vector));
}

Vector3D Vector3DNormalize(Vector3D vector)
{
    return Vector3DMultiply(vector, 1.0/Vector3DLength(vector));
}

Vector3D Vector3DRotateAroundAxisAngle(Vector3D vector, Vector3D axis, double angleInRadians)
{
    /* v is the vector, k is the axis to be rotated around */
    /* vrot = vcos0 + (kxv)sin0 + k(k*v)(1-cos0) */
    
    axis = Vector3DNormalize(axis);
    Vector3D vcos0 = Vector3DMultiply(vector, cos(angleInRadians));
    Vector3D kxvsin0 = Vector3DMultiply(Vector3DCross(axis, vector), sin(angleInRadians));
    Vector3D kkdotv1_cos0 = Vector3DMultiply(axis, Vector3DDot(axis, vector) * (1.0-cos(angleInRadians)));
    
    return Vector3DAdd(vcos0, Vector3DAdd(kxvsin0, kkdotv1_cos0));
}

double Vector3DAngleBetweenVectors(Vector3D vectorA, Vector3D vectorB)
{
    return acos(Vector3DDot(vectorA, vectorB) / Vector3DLength(vectorA) / Vector3DLength(vectorB));
}

BOOL Vector3DEqualToVector(Vector3D vectorA, Vector3D vectorB)
{
    return vectorA.x == vectorB.x && vectorA.y == vectorB.y && vectorA.z == vectorB.z;
}

Vector3D Vector3DProject(Vector3D vectorToProject, Vector3D vectorToProjectOnto)
{
    Vector3D normalizedVectorToProjectOnto = Vector3DNormalize(vectorToProjectOnto);
    double scalar = Vector3DDot(vectorToProject, normalizedVectorToProjectOnto);
    return Vector3DMultiply(normalizedVectorToProjectOnto, scalar);
}

#pragma mark - Statistical Functions

double drand48WithStandardNormalDistribution()
{
    double x1 = drand48();
    double x2 = drand48();
    return sqrt(-2*log(x1))*cos(2.0*M_PI*x2);
}

#pragma mark - Numerical Integration Functions

typedef struct State
{
    Vector3D x;
    Vector3D v;
} State;

State StateMake(Vector3D x, Vector3D v)
{
    State s;
    s.x = x;
    s.v = v;
    return s;
}

typedef struct Derivative
{
    Vector3D dx;
    Vector3D dv;
} Derivative;

Derivative DerivativeMake(Vector3D dx, Vector3D dv)
{
    Derivative d;
    d.dx = dx;
    d.dv = dv;
    return d;
}

Vector3D normalForce;
Vector3D gravitationalForce;
double rollingResistanceCoefficient;

Vector3D acceleration(State state)
{
    Vector3D velocityDirection = Vector3DNormalize(state.v);
    Vector3D rollingResistance = Vector3DMultiply(velocityDirection, rollingResistanceCoefficient * Vector3DLength(normalForce));
    Vector3D netGravitaionalForce = Vector3DAdd(normalForce, gravitationalForce);
    return Vector3DAdd(rollingResistance, netGravitaionalForce);
}

Derivative evaluate(State initial, double dt, Derivative derivative)
{
    Vector3D x = Vector3DAdd(initial.x, Vector3DMultiply(initial.v, dt));
    Vector3D v = Vector3DAdd(initial.v, Vector3DMultiply(derivative.dv, dt));
    State state = StateMake(x, v);
    
    Vector3D dv = acceleration(state);

    return DerivativeMake(state.v, dv);
}

State integrate(State state, double dt)
{
    Derivative zeroDerivative = DerivativeMake(Vector3DMake(0, 0, 0), Vector3DMake(0, 0, 0));
    Derivative a = evaluate(state, 0.0, zeroDerivative);
    Derivative b = evaluate(state, dt*0.5, a);
    Derivative c = evaluate(state, dt*0.5, b);
    Derivative d = evaluate(state, dt, c);
    
    Derivative weighteda = DerivativeMake(Vector3DMultiply(a.dx, 1.0/6.0), Vector3DMultiply(a.dv, 1.0/6.0));
    Derivative weightedb = DerivativeMake(Vector3DMultiply(b.dx, 2.0/6.0), Vector3DMultiply(b.dv, 2.0/6.0));
    Derivative weightedc = DerivativeMake(Vector3DMultiply(c.dx, 2.0/6.0), Vector3DMultiply(c.dv, 2.0/6.0));
    Derivative weightedd = DerivativeMake(Vector3DMultiply(d.dx, 1.0/6.0), Vector3DMultiply(d.dv, 1.0/6.0));
    
    Derivative derivative = DerivativeMake(Vector3DAdd(weighteda.dx, Vector3DAdd(weightedb.dx, Vector3DAdd(weightedc.dx, weightedd.dx))), Vector3DAdd(weighteda.dv, Vector3DAdd(weightedb.dv, Vector3DAdd(weightedc.dv, weightedd.dv))));
    
    return StateMake(Vector3DAdd(state.x, Vector3DMultiply(derivative.dx, dt)), Vector3DAdd(state.v, Vector3DMultiply(derivative.dv, dt)));
}
#pragma mark - Golf Ball Functions

double rollingResistanceCoefficientFromStimpmeter(double stimpmeterInFeet)
{
    /* Assume that air resistance is minimal at these speeds, just solve as though friction is the only force acting on the ball */
    /* Also assume that the coefficient of rolling resistance is a constant throughout the speeds a putted golf ball will travel */
    /* F = CrrN */
    
    
    double stimpmeterInMeters = stimpmeterInFeet * 0.3048;
    double initialVelocityInFtPerSec = 6.00;
    double initialVelocityInMPerSec = initialVelocityInFtPerSec * 0.3048;
    
    /* We have two equations and two unknowns: t and a */
    /* d = v0t + at^2/2 -> 0 = at^2/2 + v0t - d*/
    /* v = v0 + at */
    /* We know the distance at time tf is equal to the stimpmeter reading, and the velocity is zero */
    /* v(tf) = 0 -> tf = -v0/a */
    /* a(tf)^2/2 + v0tf = d */
    /* av0^2/2a^2 - v0^2/a = d */
    /* v0^2/(2a) - v0^2/a = d */
    /* 1/(2a) * (v0^2-2*v0^2) = d */
    /* 1/(2a) = d / (v0^2-2*v0^2) */
    /* a = (v0^2-2*v0^2)/d/2 */
    
    double acceleration = (initialVelocityInMPerSec * initialVelocityInMPerSec - 2 * initialVelocityInMPerSec * initialVelocityInMPerSec) / 2.0 / stimpmeterInMeters;
    double golfBallMassInKG = 0.0459262275;
    double golfBallWeight = golfBallMassInKG * 9.81;
    double rollingResistanceCoefficient = acceleration / golfBallWeight;
    
    return rollingResistanceCoefficient;
}

BOOL puttWillBeMade(State state, Vector3D holeLocation)
{
    
    
    return true;
}

Vector3D finalPositionAfterPutt(State initialState, Vector3D holeLocation)
{
    State state = initialState;
    
    double holeRadius = 0.1143; /* meters */
    double holeRadiusSq = holeRadius * holeRadius;
    while (Vector3DLengthSq(state.v) > 0.01)
    {
        state = integrate(state, 0.1);
        Vector3D ballToHole = Vector3DSubtract(state.x, holeLocation);
        double distanceToHoleSq = Vector3DLengthSq(ballToHole);
        
        if (distanceToHoleSq < holeRadiusSq)
        {
            if (puttWillBeMade(state, holeLocation))
            {
                state.x = holeLocation;
                state.v = Vector3DMake(0, 0, 0);
            }
        }
    }
    
    return state.x;
}

int main(int argc, const char * argv[])
{
    unsigned short seed = time(NULL);
    seed48(&seed);
    
    printf("Please input the stimpmeter value in feet: ");
    double stimpmeter;
    scanf("%lf", &stimpmeter);
    
    printf("The rolling resistance coefficient is %lf\n", rollingResistanceCoefficientFromStimpmeter(stimpmeter));
    
    printf("Please input the pitch of the green in degrees: ");
    double pitchInDegrees;
    scanf("%lf", &pitchInDegrees);
    double pitchInRadians = pitchInDegrees * M_PI/180.0;
    
    printf("Please input the roll of the green in degrees: ");
    double rollInDegrees;
    scanf("%lf", &rollInDegrees);
    double rollInRadians = rollInDegrees * M_PI/180.0;
    
    printf("Please input the initial speed of the ball in ft/s: ");
    double initialSpeedInFtPerS;
    scanf("%lf", &initialSpeedInFtPerS);
    double initialSpeedInMPerS = initialSpeedInFtPerS * 0.3048;
    
    printf("Please input the standard deviation of the ball's initial speed in ft/s: ");
    double initialSpeedStdDevInFtPerS;
    scanf("%lf", &initialSpeedStdDevInFtPerS);
    double initialSpeedStdDevInMPerS = initialSpeedStdDevInFtPerS * 0.3048;
    
    printf("Please input the standard deviation of target angle in degrees: ");
    double angleStdDevInDegrees;
    scanf("%lf", &angleStdDevInDegrees);
    double angleStdDevInRadians = angleStdDevInDegrees * M_PI/180.0;
    
    printf("Please input the number of trials to run: ");
    int numTrialsToRun;
    scanf("%d", &numTrialsToRun);
    
    Vector3D normalDirection = Vector3DMake(0, 1, 0);
    normalDirection = Vector3DRotateAroundAxisAngle(normalDirection, Vector3DMake(1, 0, 0), pitchInRadians);
    normalDirection = Vector3DRotateAroundAxisAngle(normalDirection, Vector3DMake(0, 0, 1), rollInRadians);
    
    Vector3D forwardDirection = Vector3DMake(0, 0, 1);
    forwardDirection = Vector3DRotateAroundAxisAngle(forwardDirection, Vector3DMake(1, 0, 0), pitchInRadians);
    forwardDirection = Vector3DRotateAroundAxisAngle(forwardDirection, Vector3DMake(0, 0, 1), rollInRadians);
    
    double golfBallMassInKG = 0.0459262275;
    
    gravitationalForce = Vector3DMake(0, -9.81 * golfBallMassInKG, 0);
    normalForce = Vector3DMultiply(Vector3DProject(gravitationalForce, normalDirection), -1);
    
    rollingResistanceCoefficient = rollingResistanceCoefficientFromStimpmeter(stimpmeter);
    
    Vector3D initialPosition = Vector3DMake(0, 0, 0);
    Vector3D initialVelocity = Vector3DMultiply(forwardDirection, initialSpeedInMPerS);
    
    Vector3D holeLocationInMeters = finalPositionAfterPutt(StateMake(initialPosition, initialVelocity), Vector3DMake(INFINITY, INFINITY, INFINITY));
    Vector3D holeLocationInFeet = Vector3DMultiply(holeLocationInMeters, 1.0/0.3048);
    
    printf("The hole is located at the following position (ft): %lf, %lf, %lf.\n", holeLocationInFeet.x, holeLocationInFeet.y, holeLocationInFeet.z);
    printf("The distance to the hole is %lf ft\n", Vector3DLength(holeLocationInFeet));
    int puttsMade = 0;
    double sumDistanceFromHoleInFeet = 0;
    for (int i = 0; i<numTrialsToRun; i++)
    {
        double angleVariation = drand48WithStandardNormalDistribution()*angleStdDevInRadians;
        double speedVariation = drand48WithStandardNormalDistribution()*initialSpeedStdDevInMPerS;
        double initialSpeed = initialSpeedInMPerS + speedVariation;
        Vector3D initialVelocity = Vector3DMultiply(Vector3DRotateAroundAxisAngle(forwardDirection, normalDirection, angleVariation), initialSpeed);
        Vector3D finalPositon = finalPositionAfterPutt(StateMake(initialPosition, initialVelocity), holeLocationInMeters);
        
        sumDistanceFromHoleInFeet += Vector3DLength(Vector3DSubtract(finalPositon, holeLocationInMeters))/0.3048;
        
        puttsMade += Vector3DEqualToVector(finalPositon, holeLocationInMeters);
    }
    
    printf("%d putts were made out of %d attempts.\n", puttsMade, numTrialsToRun);
    printf("%.1lf was the average distance from the hole in feet.\n", sumDistanceFromHoleInFeet/(double)numTrialsToRun);
    
    return 0;
}
