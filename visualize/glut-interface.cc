
#ifdef OPENGL

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <unistd.h>
#include <ctype.h>
#include <signal.h>

#include "debug.hh"
#include "opengl.hh"
#include "artisjokke-drawer.hh"
#include "mesh-connectivity.hh"
#include "glut-interface.hh"
#include "deformation-state.hh"
#include "setting.hh"

static int lastbutton;
static int lastx , lasty;
static int screen_height ;
static int screen_width ;
static bool shift_pressed;
static bool alt_pressed;
static bool ctrl_pressed;

Artisjokke_drawer * mesh_drawer_; 
static char ** argv_global;
static int argc_global;

Array<Vector3> try_wake;

void // GLUTCALLBACK
idle_func (void)
{
  glutPostRedisplay();
}

void
visible_func (int vis)
{
  if (vis == GLUT_VISIBLE)
    glutIdleFunc(idle_func);
  else
    glutIdleFunc(NULL);
}


static void
reshape(int w, int h)
{
  screen_width= w;
  screen_height = h;
  glViewport(0, 0, (GLint)w, (GLint)h);
}

static void
animate(void)
{
  mesh_drawer_->draw_frame();
  glutSwapBuffers();
  //  usleep (100000);  
}

static void
handle_key(unsigned char key, int , int )
{
  switch (key){
  case 'a':
    mesh_drawer_->axes_ = !mesh_drawer_->axes_;
    break;
  case '\e':
    exit (0);
    break ;

  case 'o':
    mesh_drawer_->original_ = !mesh_drawer_->original_;
    glutPostRedisplay();
    break ; 
    
  case 'b':
    mesh_drawer_->boundary_ = !mesh_drawer_->boundary_ ;
    glutPostRedisplay ();
    break;

  case '=':
    mesh_drawer_->print_stats ();
    break;     

  case '\f':
    mesh_drawer_->view_drag_. reset ();
    glutPostRedisplay();
    break;
    
  case 'e':
    mesh_drawer_->    edges_ = !mesh_drawer_->edges_;
    glutPostRedisplay();
    break;

  case 'd':
    mesh_drawer_->fog_ =     !mesh_drawer_->fog_ ;
    printf ("Depth cue %d\n", mesh_drawer_->fog_);
    
    mesh_drawer_->set_depth_cue ();
    glutPostRedisplay();
    break; 
    
  case 's':

    mesh_drawer_->solid_ = !mesh_drawer_->solid_;
    mesh_drawer_->lighting_ = mesh_drawer_->solid_;
    mesh_drawer_->set_lighting();

    printf  ("solid %d", mesh_drawer_->solid_);
    glutPostRedisplay();
    break;

  case 'l':
    mesh_drawer_->lighting_ = !mesh_drawer_->lighting_;
    printf ("Lighting set to  %d\n", mesh_drawer_->lighting_);
    break;
    
  case 'w':
    write_deformation_state (mesh_drawer_->deformation_, "out.model.state");
    break;

  case 'n':
    mesh_drawer_->nodes_ = !mesh_drawer_->nodes_;
    glutPostRedisplay();
    break;

  case 'c':
    glutPostRedisplay();
    break;

  case '?':
    print_settings ();
    break ;

  default:
    mesh_drawer_->process_key (key);
    break ; 
  }
}


static void
special_key (int key , int x, int y)
{
  int axis = -1;
  int dir = 0;

  switch (key)
    {
    case GLUT_KEY_LEFT:
      axis =0;
      dir = -1;
      break;
    case GLUT_KEY_UP:
      axis = 1;
      dir = 1;
      break;
    case GLUT_KEY_RIGHT:
      axis = 0;
      dir = 1;
      break;
    case GLUT_KEY_DOWN:
      axis = 1;
      dir = -1;
      break;
    default:
      break;
    }
}


static void
motion (int x,int y)
{
  Drag_state * current_drag  = &mesh_drawer_->view_drag_ ;
  if (shift_pressed)
    current_drag = &mesh_drawer_->needle_drag_;
  else if (ctrl_pressed )
    current_drag = &mesh_drawer_->light_drag_;
  else if (alt_pressed)
    current_drag = &mesh_drawer_->clip_drag_;
      
  if (lastbutton == 0)
    mesh_drawer_->rotate_drag (&current_drag->quat_,
                               screen_width, screen_height,
                               lastx, lasty,
                               x, y);
  else if (lastbutton == 1)
    mesh_drawer_->scale_drag (&current_drag->scale_,
                               screen_width, screen_height,
                               lastx, lasty,
                               x, y);
  else if (lastbutton == -1)
    mesh_drawer_->translate_drag (&current_drag->translate_,
                               screen_width, screen_height,
                               lastx, lasty,
                               x, y);
  
  lastx = x;
  lasty = y;
  glutPostRedisplay();    
}

static void
mouse_event (int button, int state, int x, int y)
{
  lastx = x;
  lasty = y;
  
  int b = 0;
  //  int s = 0;
  switch (button) {
  case GLUT_MIDDLE_BUTTON:
    b = 0;break;
  case GLUT_LEFT_BUTTON:
    b = -1; break;
  case GLUT_RIGHT_BUTTON:
    b = 1; break;
  }
  lastbutton = b;

  shift_pressed = (glutGetModifiers () & GLUT_ACTIVE_SHIFT);
  alt_pressed = (glutGetModifiers () & GLUT_ACTIVE_ALT);
  ctrl_pressed =   (glutGetModifiers () & GLUT_ACTIVE_CTRL);
}


void
glut_setup (Artisjokke_drawer *md, int argc, char * argv[])
{
  printf ("\n\nSetting up GLUT OpenGL display.\n");

  argv_global = argv;
  argc_global = argc;
  
  glutInit(&argc, argv);
  int type = GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE;
  glutInitDisplayMode(type);
  
  int wp = 0;
  int ws = 512;

  glutInitWindowPosition (wp, wp);
  glutInitWindowSize (ws, ws);
  
  if (glutCreateWindow("Artisjokke") == GL_FALSE)
    {
      exit(1);
    }

  mesh_drawer_ = md;

  glutReshapeFunc(&reshape);
  glutKeyboardFunc(&handle_key);
  glutDisplayFunc(&animate);
  glutMouseFunc(&mouse_event);
  glutMotionFunc(&motion);
  glutSpecialFunc (&special_key);
  glutVisibilityFunc (&visible_func);
  //  glutTimerFunc(100, timer, 0);
  glutIdleFunc(idle_func);

  md->setup ();
  glutMainLoop();
}
#endif 
