#ifndef CELESTIAL_C_CELESTIAL_H
#define CELESTIAL_C_CELESTIAL_H

typedef int bool;
#define true 1
#define false 0

typedef struct orbital_elements_t *orbital_elements;

orbital_elements create_orbital_elements();

void free_orbital_elements(orbital_elements);

#endif //CELESTIAL_C_CELESTIAL_H

