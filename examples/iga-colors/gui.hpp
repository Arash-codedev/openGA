// This library is free and distributed under
// Mozilla Public License Version 2.0.

#include <GL/glut.h>
#include <GL/freeglut.h>
#include <SDL2/SDL.h>
#include <thread>

bool to_reopen_window=false;
double gui_subject_R,gui_subject_G,gui_subject_B;

void timer(int value )
{
	(void) value;

	if(to_reopen_window)
	{
		glutLeaveMainLoop();
	}

	if(!glutGetWindow())
		return ;
	
	glutPostRedisplay();
	glutTimerFunc(30, timer, 1);

}

void display()
{
	if(!glutGetWindow())
		return ;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	float r=float(gui_subject_R/255.0);
	float g=float(gui_subject_G/255.0);
	float b=float(gui_subject_B/255.0);
	glClearColor(r,g,b,1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, 100, 0, 100, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glColor3ub( 255, 255, 255 );
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_COLOR_ARRAY );
	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_COLOR_ARRAY );

	glFlush();
	glutSwapBuffers();
}

void reshape(int w, int h)
{
	if(!glutGetWindow())
		return ;
	glViewport(0, 0, w, h);
}

void init_gui()
{
	int argc=1;
	char *argv[1] = {(char*)"something"};
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);

	glutInitWindowSize(150,150);

	glutCreateWindow("IGA");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutTimerFunc(30, timer, 1);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
		GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	
	glutMainLoopEvent();
}

void refresh_gui()
{
	if(to_reopen_window)
	{
		if(!glutGetWindow())
		{
			glutMainLoop();
			init_gui();
		}
		to_reopen_window=false;
	}
	if(glutGetWindow())
		glutPostRedisplay();
	glutMainLoopEvent();
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	// boost::this_thread::sleep(boost::posix_time::milliseconds(100));
}
