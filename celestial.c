//=====================================================================
//
//   ####  #####  ##      #####   ####  ######  ##    ###    ##
//  ##     ##     ##      ##     ##       ##    ##   ## ##   ##
//  ##     #####  ##      #####   ###     ##    ##  ##   ##  ##
//  ##     ##     ##      ##        ##    ##    ##  #######  ##
//   ####  #####  ######  #####  ####     ##    ##  ##   ##  ######
//
//=====================================================================

/*
    -Myvar:
    Coding conventions:
    1. Please add the thing your using in an include, with a comment
            ->
                //malloc, free
                #include <stdlib.h>

    2. Use Comment with one empty one above an one empty below
            ->
                //
                // Some Comment
                //

    3. Use Sections to make code navagation easyer
            ->
                //
                // [Helpers]
                //
                inline double rad_to_deg(double rad) { ... }

    4. use @ tags to mark things: (ensure tag index is upto date)
            ->
                // @Hack !!! this code is hacky
                // @Cleanup Code is messy the names need to be updated
                // @Speed !! malloc a big block insted of meany small malocs for cacheing

    5. prefix all procedures with a structs name if its related to that struct
            -> void orbital_elements_init(orbital_elements);
                insted of just init

 */

//
// Tag Index
//

/*
    * @Todo
    * @Speed
    * @Cleanup
    * @Refactoring
    * @Hack
    * @bug
 */

//
// Todo/Ideas
//

/*
    * Orbital Mechanics
    * Stellar classification
*/

//malloc, free
#include <stdlib.h>

//printf
#include <stdio.h>

//time_t
#include <time.h>

//M_PI, sqrt, abs, round
#include <math.h>

#include "celestial.h"

#ifdef TEST

//assert
#include <assert.h>

#endif

//
// [Implmentation]
//

//
// [Helpers]
//
static inline double rad_to_deg(double rad)
{
    // @Cleanup: vs code things M_PI does not exists but it does wtf
    return rad * (180.0 / M_PI);
}

static inline double deg_to_rad(double deg)
{
    // @Cleanup: vs code things M_PI does not exists but it does wtf
    return deg * (M_PI / 180.0);
}

// @Refactoring, @Cleanup : this can easly be confused with sqrT maby rename it
static inline double sqr(double n)
{
    return n * n;
}

const double GravConst_m_kg_s = 6.67408E-11;
const double GravConst_km_kg_s = 6.67408E-20;
//
// [Math]
//

//
// [Vector]
//

struct vector_t
{
    double x;
    double y;
    double z;
};

vector vector_create()
{
    return malloc(sizeof(struct vector_t));
}

void vector_free(vector vec)
{
    free(vec);
}

vector clone(vector vec)
{
    vector a = vector_create();
    a->x = vec->x;
    a->y = vec->y;
    a->z = vec->z;

    return a;
}

vector vector_init(double x, double y, double z)
{
    vector vec = vector_create();

    vec->x = x;
    vec->y = y;
    vec->z = z;

    return vec;
}

// maths
double vector_distance_to(vector a, vector b)
{
    double x = a->x - b->x;
    double y = a->y - b->y;
    double z = a->z - b->z;

    return sqrt(x * x + y * y + z * z);
}

void vector_abs(vector vec)
{
    vec->x = abs(vec->x);
    vec->y = abs(vec->y);
    vec->z = abs(vec->z);
}

double vector_length(vector vec)
{
    //
    return sqrt(sqr(vec->x) + sqr(vec->y) + sqr(vec->z));
}

double vector_normalize(vector vec)
{
    double l = vector_length(vec);

    vec->x /= l;
    vec->y /= l;
    vec->z /= l;

    return l;
}

double vector_dot(vector a, vector b)
{
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

vector vector_cross(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->y * b->z - a->z * b->y;
    vec->y = a->z * b->x - a->x * b->z;
    vec->z = a->x * b->y - a->y * b->x;

    return vec;
}

// operators
vector vector_add(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x + b->x;
    vec->y = a->y + b->y;
    vec->z = a->z + b->z;

    return vec;
}

vector vector_add_d(vector a, double b)
{
    vector vec = vector_create();

    vec->x = a->x + b;
    vec->y = a->y + b;
    vec->z = a->z + b;

    return vec;
}

vector vector_sub_d(vector a, double b)
{
    vector vec = vector_create();

    vec->x = a->x - b;
    vec->y = a->y - b;
    vec->z = a->z - b;

    return vec;
}

vector vector_div(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x / b->x;
    vec->y = a->y / b->y;
    vec->z = a->z / b->z;

    return vec;
}

vector vector_div_d(vector a, double b)
{
    vector vec = vector_create();
    vec->x = a->x / b;
    vec->y = a->y / b;
    vec->z = a->z / b;

    return vec;
}

vector vector_mul(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x * b->x;
    vec->y = a->y * b->y;
    vec->z = a->z * b->z;

    return vec;
}

vector vector_mul_d(vector a, double b)
{
    vector vec = vector_create();

    vec->x = a->x * b;
    vec->y = a->y * b;
    vec->z = a->z * b;

    return vec;
}

void vector_print(vector vec)
{
    printf("[X: %f, Y: %f, Z: %f]\n", vec->x, vec->y, vec->z);
}

//
// [Orbital_Mechanics]
//
struct orbital_elements_t
{
    double semi_major_axis_km;          //aka a
    double eccentricity;                //aka e
    double inclination_rad;
    double longitude_of_ascending_node_rad;
    double argument_of_periapsis_rad;
    double mean_anomaly_rad;            //aka M0 (at epoch)
    double epoch;                       //seconds past periapsis
    double gravitational_parameter_km;  //km^3 s^-2 (si units are m^3 s^-2
    double orbital_period_seconds;
    long mean_motion_rad;
    double apoapsis_km;
    double periapsis_km;
    double mass_kg;
    double parent_mass_kg;
};

struct state_vectors_t
{
    double parentMass;
    double myMass;
    vector velocity;
    vector position;
    double secondsPastEpochOfLastCalc;
};

// @Speed ! mallocing a big struct like this piecemeal is slow, it would be way faster to use an Array for Struct, with block mallocs
orbital_elements orbital_elements_create()
{
    return malloc(sizeof(struct orbital_elements_t));
}

void orbital_elements_free(orbital_elements elem) { free(elem); }

state_vectors state_vectors_create()
{
    return malloc(sizeof(struct state_vectors_t));
}

void state_vectors_free(state_vectors vecs)
{
    free(vecs);
}

void orbital_elements_init_from_major_planet(
    struct orbital_elements_t *elem,
    double parentMass,
    double myMass,
    double semiMajorAxis,
    double eccentricity,
    double inclination,
    double longitudeOfAscendingNode,
    double longitudeOfPeriapsis,
    double meanLongitude,
    double epoch)
{
}

void orbital_elements_calculate_extended_parameters(struct orbital_elements_t *elem)
{
    elem->gravitational_parameter_km = GravConst_km_kg_s * (elem->parent_mass_kg + elem->mass_kg); //km^3 s^-2
    elem->mean_motion_rad = sqrt(elem->gravitational_parameter_km / pow(elem->semi_major_axis_km, 3));

    elem->orbital_period_seconds = 2 * M_PI * sqrt(pow(elem->semi_major_axis_km, 3) / (elem->gravitational_parameter_km));

    elem->apoapsis_km = (1 + elem->eccentricity) * elem->semi_major_axis_km;
    elem->periapsis_km = (1 - elem->eccentricity) * elem->semi_major_axis_km;

}

void orbital_elements_init_from_vector(
    struct orbital_elements_t *elem,
    struct state_vectors_t state_vecs)
{
    double a;   //SemiMajorAxis
    double b;   //SemiMinorAxis
    double e;   //Eccentricity
    double i;   //inclination
    double loAN;//longditudeOfAccendingNode omega
    double aoP; //Argument of periapsis Omega
    double M0;  //MeanAnomaly at Epoch
    double E;   //EccentricAnomaly
    double sgp =  GravConst_m_kg_s * state_vecs.parentMass * state_vecs.parentMass; //standardGravParameter
    vector vel = state_vecs.velocity;
    vector pos = state_vecs.position;
    double p;   //


    //calculate angularVelocity and NodeVector
    vector angularVelocity = vector_cross(pos, vel);
    vector zed = vector_init(0, 0, 1);
    vector nodeVector = vector_cross(zed, angularVelocity);

    //calculate eccentricity
    vector eccentricityVector = vector_div_d(vector_cross(vel, angularVelocity),  sgp);
    eccentricityVector = eccentricityVector - vector_div_d(state_vecs.position , vector_length(state_vecs.position));
    e = vector_length(eccentricityVector);

    double specificOrbitalEnergy = pow(vector_length(vel),2) * vector_length(vel) * 0.5 - sgp / vector_length(pos);

    if (abs(e) > 1) //hypobola
    {
        a = -(-sgp / (2 * specificOrbitalEnergy)); //in this case the sma is negitive
        p = a * (1 - pow(e, 2));
    }
    else if (abs(e) < 1) //ellipse
    {
        a = -sgp / (2 * specificOrbitalEnergy);
        p = a * (1 - e * pow(e, 2));
    }
    else //parabola consider pushing this to just over parabola and using a hyperbola instead
    {
        p = pow(vector_length(angularVelocity), 2) / sgp;
        a = INFINITY;
    }
    b = a * sqrt(1 - pow(e, 2));
    double linierEccentricity = e * a;

    //inclination
    i = acos(angularVelocity->z / vector_length(angularVelocity));
    if(isnan(i))
        i = 0;
    //LonditudeOfAccendingNode
    double loANLen = nodeVector->x / vector_length(nodeVector);
    if(isnan(loANLen))
        loANLen = 0;
    else
    {
        loANLen = fmin(loANLen, -1);
        loANLen = fmax(loANLen, 1);
    }
    if(loANLen != 0)
        loAN = acos(loANLen);

    //calculate argumentOfPeriapsis
    if(loAN == 0)
    {
        aoP = atan2(eccentricityVector->y, eccentricityVector->x);
        if(vector_cross(pos, vel)->z < 0) //anticlockwise Orbit
            aoP = M_PI * 2 - aoP;
    }
    else
    {
        double aopLen = vector_dot(nodeVector, eccentricityVector);
        aopLen = aopLen / (vector_length(nodeVector) * e);
        aopLen = fmin(aopLen, 1);
        aopLen = fmax(aopLen, -1);
        aoP = acos(aopLen);
        if(eccentricityVector->z < 0) //anticlockwiseOrbit
            aoP = M_PI * 2 - aoP;
    }

    //calculate EccentricAnomaly
    double x = (pos->x * cos(-aoP)) - (pos->y * sin(-aoP));
    x = linierEccentricity + x;
    E = acos(x / a);
    //calculate MeanAnomaly at epoch
    M0 = e - e * sin(E);

    elem->semi_major_axis_km = a;
    elem->eccentricity = e;
    elem->apoapsis_km = (1 + e) * a;
    elem->periapsis_km  = (1 - e) * a;
    elem->longitude_of_ascending_node_rad = loAN;
    elem->argument_of_periapsis_rad = aoP;
    elem->inclination_rad = i;
    elem->mean_anomaly_rad = M0;
    elem->epoch = 0;
}

state_vectors orbital_elements_get_state_vectors(struct orbital_elements_t elem, struct state_vectors_t *state, double secondsSinceLastCalc)
{
    double trueAnomaly_rad;

    double secondsFromEpoch = state->secondsPastEpochOfLastCalc + secondsSinceLastCalc;

    //don't let secondsFromEpoch get too big
    while(secondsFromEpoch > elem.orbital_period_seconds)
    {
        secondsFromEpoch -= elem.orbital_period_seconds;
    }
    state->secondsPastEpochOfLastCalc = secondsFromEpoch;

    //calculate current mean anomaly:
    double currentMeanAnomaly = elem.mean_anomaly_rad;
    // Add nT
    currentMeanAnomaly += elem.mean_motion_rad * secondsFromEpoch;
    // Large nT can cause meanAnomaly to go past 2*Pi. Roll it down. It shouldn't, because timeSinceEpoch should be tapered above, but it has.
    currentMeanAnomaly = fmod(currentMeanAnomaly, (M_PI * 2));

    //calculate EccentricAnomaly:
    double eccentricAnomaly;
    short numIterations = 100;
    double e [numIterations];
    double epsilon = 1E-12; //amount of error we're happy with
    short i = 0;
    if(elem.eccentricity > 0.8)
    {
        e[i] = M_PI;
    }
    else
    {
        e[i] = currentMeanAnomaly;
    }

    do
    {
        e[i + 1] = e[i] - (e[i] - elem.eccentricity * sin(e[i]) - currentMeanAnomaly) / (1 - elem.eccentricity * cos(e[i]));
        i++;
    }
    while(abs(e[i] - e[i - 1]) > epsilon && i + 1 < numIterations);
    if(i + 1 >= numIterations)
    {
        eccentricAnomaly = 0; //failed to converge on answer;
    }
    else
    {
        eccentricAnomaly = e[i-1];
    }

    double x = cos(eccentricAnomaly) - elem.eccentricity;
    double y = sqrt(1 - elem.eccentricity * elem.eccentricity) * sin(eccentricAnomaly);

    trueAnomaly_rad = atan2(y,x);



    double radius_km = elem.semi_major_axis_km * (1 - elem.eccentricity * elem.eccentricity) / (1 + elem.eccentricity * cos(trueAnomaly_rad));
    double angleFromAN = trueAnomaly_rad + elem.argument_of_periapsis_rad;
    double lofAN = elem.longitude_of_ascending_node_rad;

    state->position->x = radius_km * (cos(lofAN) * cos(angleFromAN) - sin(lofAN) * sin(angleFromAN) * cos(elem.inclination_rad));
    state->position->y = radius_km * (sin(lofAN) * cos(angleFromAN) + cos(lofAN) * sin(angleFromAN) * cos(elem.inclination_rad));
    state->position->z = radius_km * (sin(elem.inclination_rad) * sin(angleFromAN));

    double positionLen = sqrt(state->position->x * state->position->x + state->position->y * state->position->y + state->position->z * state->position->z);
    double speed = sqrt(elem.gravitational_parameter_km * (2 / positionLen - 1 / elem.semi_major_axis_km));
    double linierEccentricity = elem.semi_major_axis_km * elem.eccentricity;


}
//
// Testing
//

#ifdef TEST

void vector_test()
{
    printf("Starting Vector Tests:\n");
    printf("--------------------\n");

    //init
    vector a = vector_init(10, 5, 7);
    vector b = vector_init(4, 3, 2);
    vector tmp = vector_init(-11, 5, 7);

    // vector_print(a);
    // vector_print(b);
    // vector_print(tmp);

    printf("vector_distance_to: ");
    assert(round(vector_distance_to(a, b)) == round(8.062258));
    printf("PASSED\n");

    printf("vector_abs: ");
    vector_abs(tmp);
    assert(tmp->x == 11);
    printf("PASSED\n");

    printf("vector_distance_to: ");
    assert(round(vector_length(a)) == round(13.190905));
    printf("PASSED\n");

    printf("vector_normalize: ");
    vector_normalize(tmp);
    assert(round(tmp->x) == round(0.787726));
    printf("PASSED\n");

    printf("vector_dot: ");
    assert(vector_dot(a, b) == 69);
    printf("PASSED\n");

    printf("vector_cross: ");
    vector_free(tmp);
    tmp = vector_cross(a, b);
    assert(tmp->y == 8);
    printf("PASSED\n");

    printf("vector_add: ");
    vector_free(tmp);
    tmp = vector_add(a, b);
    assert(tmp->x == a->x + b->x);
    printf("PASSED\n");

    printf("vector_sub: ");
    vector_free(tmp);
    tmp = vector_sub(a, b);
    assert(tmp->x == a->x - b->x);
    printf("PASSED\n");

    printf("vector_div: ");
    vector_free(tmp);
    tmp = vector_div(a, b);
    assert(tmp->x == a->x / b->x);
    printf("PASSED\n");

    printf("vector_mul: ");
    vector_free(tmp);
    tmp = vector_mul(a, b);
    assert(tmp->x == a->x * b->x);
    printf("PASSED\n");

    //free
    vector_free(a);
    vector_free(b);
    vector_free(tmp);
}

void orbital_elements_test()
{
    orbital_elements elm;
    elm = orbital_elements_create();

    orbital_elements_free(elm);
}

int main()
{
    vector_test();

    orbital_elements_test();

    return 0;
}

#endif
