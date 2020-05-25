#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "point.h"        // for point_t
#include "tab_point.h"    // for tabular_point_free, tabular_point_add, tabu...
#include "convex_hull.h"
//#define FPRINT 1

double tab_x[6] = {0,0.3,0.5,0.8,0.9,0.95};
double tab_y[6] = {0,3,6,12,18,23};

int cmp_double (double a, double b)
{
    double diff = a - b;
    double precision = 0.0000000001;

    if ( diff < precision && diff > -1*precision)
    	return 0;
    else if (diff < precision)
    	return -1;
    else
    	return 1;
}

int check_points( tabular_point_t* t)
{
    if (t->index != 6)
	return 0;
    for (int i = 0; i < 6; i++)
	{
	    if (cmp_double (t->tab[i]->x, tab_x[i]) != 0 ||
		cmp_double (t->tab[i]->y, tab_y[i]) != 0)
		return 0;
	}
    return 1;
}

int main(int argc, char**argv)
{
    srand(time(NULL));
#ifdef FPRINT
    FILE *graph_file;
    char init[] = "/tmp/point1";
    graph_file = fopen (init, "w");
#endif
    tabular_point_t* list = tabular_point_create ();

    //to test the sensibility of our selection!!
    int CST = 1;
    double PREC = 0.1;
    if (argc > 2)
	{
	    CST = atoi(argv[1]);
	    PREC = strtod (argv[2], NULL);
	}
    //add several bad points around our good point, to test if our
    //convex hull doesn't make a mistake!
    for (int i = 0; i < 6; i++)
	{
	    tab_y[i]*= CST;
	    double x = tab_x[i];
	    double y = tab_y[i];
	    y = tab_y[i]*(1+PREC);
	    tabular_point_add (list, i, x, y);
#ifdef FPRINT
	    fprintf(graph_file,  "%f \t %f\n",x,y);
#endif	    
	    x = tab_x[i]*(1-PREC);	    
	    tabular_point_add (list, i, x, y);
#ifdef FPRINT
	    fprintf(graph_file,  "%f \t %f\n",x,y);
#endif
	    tabular_point_add (list, i, tab_x[i], tab_y[i]);
#ifdef FPRINT
	    fprintf(graph_file,  "%f \t %f\n",tab_x[i],tab_y[i]);
#endif	    

	}
#ifdef FPRINT
    fclose(graph_file);
#endif
    
#ifdef FPRINT
    char init2[] = "/tmp/point2";
    graph_file = fopen (init2, "w");
#endif

    tabular_point_t* res = convex_hull (list);

#ifdef FPRINT    
    int len = tabular_point_get_index (res);
    point_t *el;
    for (int i=0; i < len; i++)
	{
	    el = tabular_point_get_point(res, i);
	    fprintf(graph_file,  "%f \t %f\n", el->x, el->y);
	}
    fclose(graph_file);
#endif
    //check our result!  test if the convex is exactly our tab_x and tab_y!
    //if PREC < 0 then you insert good points, so the
    //convex_hull isn't tab_x:tab_y!
    if ((check_points (res) && PREC < 0) && (!check_points (res) && PREC > 0))
	{
	    fprintf (stderr, "error with our convex hull.\n If you"
		     " want to see what's happen: use the define FPRINT!\n");
	    tabular_point_free (res);
	    tabular_point_free (list);
	    return EXIT_FAILURE;
	}

    //free
    tabular_point_free (res);
    tabular_point_free (list);


    //constant point
    return EXIT_SUCCESS;
}
