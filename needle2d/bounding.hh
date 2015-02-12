#ifndef BOUNDING_HH
#define BOUNDING_HH
#include "proto.hh"

void nodes_bounding_sphere (Vector2 *center, Real * radius,
		       Link_array<Node> const &nodes, Deformation_state *def);

void nodes_bounding_box (Vector2 *min, Vector2 * max, Link_array<Node> const &nodes, Deformation_state*def);

#endif
