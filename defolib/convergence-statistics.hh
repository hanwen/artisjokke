#ifndef CONVERGENCE_STATISTICS_HH
#define CONVERGENCE_STATISTICS_HH

#include "defo-proto.hh"

Deformation_state *read_reference_state (Mesh_connectivity*top) ;

extern long long vector_flop_count, element_flop_count;
extern int global_iteration_count;
extern int global_max_iteration_count;

void flop_stats (Deformation_state * def, Deformation_state *ref) ;

Real reference_state_energy_distance (Deformation_state * def , Deformation_state *ref);

long long flop_count ();
Real reference_state_distance (Deformation_state * def, Deformation_state * ref);

#endif
