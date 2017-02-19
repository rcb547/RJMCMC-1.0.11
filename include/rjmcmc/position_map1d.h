#ifndef position_map1d_h
#define position_map1d_h

typedef struct _position_map1d position_map1d_t;

position_map1d_t *
position_map1d_create(int max_partitions, double minx, double maxx);

void
position_map1d_destroy(position_map1d_t *p);

void
position_map1d_clone(const position_map1d_t *src,
		     position_map1d_t *dst);

int
position_map1d_npartitions(const position_map1d_t *p);

int
position_map1d_insert(position_map1d_t *p,
		      double x,
		      int iy);

int 
position_map1d_delete(position_map1d_t *p,
		      double x,
		      int iy);

int
position_map1d_move(position_map1d_t *p,
		    double x,
		    double new_x);

int 
position_map1d_small_move(position_map1d_t *p,
			  double x,
			  double new_x);

int 
position_map1d_nearest(position_map1d_t *p,
		       double x,
		       double *nx);

double 
position_map1d_position_of_index(position_map1d_t *p,
				 int iy);

int 
position_map1d_next_index(position_map1d_t *p,
			  double x);

double
position_map1d_next_position(position_map1d_t *p,
			     double x);

int
position_map1d_predecessor_of_point(position_map1d_t *p,
				    double x);

int
position_map1d_predecessor_of_index(position_map1d_t *p,
				    int iy);

int 
position_map1d_successor_of_point(position_map1d_t *p,
				  double x);

int 
position_map1d_successor_of_index(position_map1d_t *p,
				  int iy);

int
position_map1d_traverse_intervals(position_map1d_t *p,
				  int (*interval_cb)(void *user_arg,
						     double xmin,
						     double xmax,
						     int iy,
						     int riy),
				  void *user_arg);

int 
position_map1d_fill_list(position_map1d_t *p,
			 double *positions,
			 int *npartitions);

void
position_map1d_dump(position_map1d_t *p);

int
position_map1d_valid(position_map1d_t *p);

#endif /* position_map1d_h */
