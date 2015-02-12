
#ifndef TEST_MESH_HH
#define TEST_MESH_HH
void test_scenario (char const *ch, Maubach_tree * , Deformation_state *, Needle_inserter * ins);
void test_mesh_connectivity (Mesh_connectivity*);
void test_mesh (Mesh_connectivity*);
Mesh_connectivity*generate_test_mesh ();


extern int mesh_size;

#endif
