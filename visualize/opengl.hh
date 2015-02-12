
#ifdef OPENGL
#ifndef OPENGL_HH
#define OPENGL_HH

/*
Windows include MUST be before GL includes, otherwise
the compiler generates the wrong calling convention (sigh)
*/

#if defined(_WIN32)
#include <windows.h>
#endif


#include "GL/glut.h"

#endif
#endif
