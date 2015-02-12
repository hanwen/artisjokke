/*
  cube-generator.cc -- generate a tetrahedral mesh on a cubic grid.
 */

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

#include "mesh-connectivity.hh"
#include "defo-proto.hh"
#include "parray.hh"
#include "vector.hh"
#include "geometry.hh"
#include "deformation-state.hh"

Element *make_positive_tetrahedron (Mesh_connectivity *fem, Node ** ns, Deformation_state*);

struct Int_vector3
{
  int e_[3];
  int & operator [](int i) { return e_[i]; }
  Int_vector3 (int i, int j, int k)
    {
      e_[0] = i;
      e_[1] = j;
      e_[2] = k;
    }
  Int_vector3 ()
    {
      e_[0] = 0;
      e_[1] = 0;
      e_[2] = 0;
    }
};

/*
   Generate a tetrahedral model that consists of nodes in a uniform
   rectangular grid.  We generate one by one the cubelets: each
   cubelet contains six tetrahedrons. Triangles are created on the fly
   when needed.

   The cube is centered around the origin. 
 */
class Cube_generator
{
  Link_array<Node> nods_;
  Mesh_connectivity  *mesh_;
  Deformation_state* deformation_;
public:  
  Int_vector3 counts_;
  Real rand_fact_ ; 
  Real dims_[3];

  void generate_cubelet_volume (Int_vector3 idx);

  int index_3_to_1 (int i, int j, int k);
  int index_3_to_1 (Int_vector3 id);
  int try_generate_tetra (Int_vector3 p0,Int_vector3 p1, Int_vector3 p2,
			  Int_vector3 p3);

  bool valid_point (Int_vector3 p);
  void init ();  
  Cube_generator();
  void generate (  Mesh_connectivity *,Deformation_state*);
};

void
make_cube_mesh (Mesh_connectivity  * mt, Deformation_state*def,
		int cx, int cy, int cz,
		Real dx, Real dy, Real dz)
{

  Cube_generator cg;
  cg.counts_[0] = cx;
  cg.dims_[0] = dx;
  cg.counts_[1] = cy;
  cg.dims_[1] = dy;
  cg.counts_[2] = cz;
  cg.dims_[2] = dz;
  
  cg.generate (mt, def);
}

//Mesh_connectivity * generate_tetra (void);

int
Cube_generator::index_3_to_1 (int i, int j, int k)
{
  return counts_[1] *  counts_[2] * i + counts_[2] * j + k;
}

int
Cube_generator::index_3_to_1 (Int_vector3 id)
{
  return index_3_to_1 (id[0], id[1], id [2]);
}

bool
Cube_generator::valid_point (Int_vector3 p)
{
  for (int axis =0; axis < 3; axis++)
    if (p[axis] < 0 || p[axis] >= counts_[axis])
      return false;
  return true;    
}

Cube_generator::Cube_generator ()
{
  Int_vector3 c(2, 2, 2);
  Real d[3] = {1.0, 1.0, 1.0};
  counts_ = c;


  /*
    This random factor  is very confusing.
   */
  rand_fact_ = 0.00;
  for (int  i = 3; i--;)
    dims_[i] = d[i];
  init ();
}

void
Cube_generator::init ()
{
  mesh_ = 0;
  nods_.clear ();
}

int 
Cube_generator::try_generate_tetra (Int_vector3 p0, Int_vector3 p1, Int_vector3 p2,
				    Int_vector3 p3)
{
  Int_vector3 p[4];
  p[0] = p0;
  p[1] = p1;
  p[2] = p2;
  p[3] = p3;

  for (int i=0; i < 4; i++)
    if (!valid_point (p[i]))
      return 0;


  Node *nodes[4];
  for (int i=0; i < 4; i++)
    nodes[i] = nods_[index_3_to_1 (p[i])];

  make_positive_tetrahedron (mesh_, nodes, deformation_);

  return 1;
}

void
Cube_generator::generate_cubelet_volume (Int_vector3 idx)
{
  int x = idx [0];
  int y = idx [1];
  int z = idx [2];
  

  /*
    So sue me.  I'm too dense to figure out how to do this procedurally.
   */

  Int_vector3 a (x-1,y,  z-1);
  Int_vector3 b (x,  y,  z-1);
  Int_vector3 c (x,  y,  z);
  Int_vector3 d (x-1,y,  z);
  Int_vector3 e (x-1,y-1,z-1);
  Int_vector3 f (x,  y-1,z-1);
  Int_vector3 g (x,  y-1,z);
  Int_vector3 h (x-1,y-1,z);

  int tetras_generated=0;
  
  tetras_generated += try_generate_tetra (a, f, h, e);
  tetras_generated += try_generate_tetra (a, f, h, b);
  tetras_generated += try_generate_tetra (a, b, h, d);
  tetras_generated += try_generate_tetra (b, f, g, h);
  tetras_generated += try_generate_tetra (b, g, h, d);
  tetras_generated += try_generate_tetra (d, g, b, c);
}


/**
   Generate a cube of tetrahedrons.
 */
void
Cube_generator::generate (Mesh_connectivity  *m, Deformation_state *def)
{
  init ();
  mesh_ = m;
  for (int i=0; i <3; i++)
    assert (counts_[i] >=2);

  deformation_ =def;
  int nod_count = counts_[0] * counts_[1] * counts_[2];
  nods_.set_size (nod_count);
  for (int i=0; i < nods_.size (); i++)
    nods_[i] =0;

  fprintf (stdout, "generating cube of %d nodes (%d, %d, %d) ... ", nod_count,
	   counts_[0], counts_[1], counts_[2]

	   );
  fflush (stdout);
  
  Vector3 delta;
  for (int i=0; i < 3; i++)
    delta(i) = (counts_[i]-1) ? (dims_[i] / Real(counts_[i] -1.0)) : 0.0;
    
  /*
    generate the points
   */
  Int_vector3 idx;


  int c =0;

  for (idx[0]=0; idx[0] < counts_[0]; idx[0]++)
    for (idx[1]=0; idx[1] < counts_[1]; idx[1]++)
      for (idx[2]=0; idx[2] < counts_[2]; idx[2]++)
	{
	  Vector3 v;
	  for (int l = 0; l < 3; l++)
	    {
	      v(l) = idx[l] * delta (l) - dims_[l] * .5;
	      v(l) +=
		(dims_[l]/ counts_[l]) * rand_fact_ * random () / RAND_MAX * delta(l);
	    }

	  Node *n = mesh_->make_node ();
	  def->set_reference_location (n, v);	  
	  int k = index_3_to_1 (idx);
	  assert (!nods_[k]);
	  nods_[k] = n;

	  generate_cubelet_volume (idx);
	  if (!(++c %100))
	    {
	      printf ("[%d]", c);
	      fflush (stdout);
	    }
	}
  

  fprintf (stdout, "done\n");
}

/*
  Make a tetrahedron ensuring that volume is positive.
 */
Element *
make_positive_tetrahedron (Mesh_connectivity *fem, Node ** ns, Deformation_state*def)
{
  Simplex tetset(TETRAHEDRON, ns, 0);

  Real vol = simplex_volume3 (tetset, &reference_location3, def);

  if (fabs (vol) <  1e-18)
    {
      printf ( "Null volume tetrahedron!\n ");
      tetset.print();
      return 0;
    }
  
  if (vol < 0)
    tetset = tetset.get_mate ();

  Link_array<Element> elts;
  Array<Simplex> ss;
  ss.push (tetset);

  elts = fem->replace_elements (elts, ss);

  return elts[0];
}

