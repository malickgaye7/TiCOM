struct ADVECTION {
  int ADVECTION;
 
  float gamma;
  float timestep;
  float fine_tune_dt;
  float dt_reduced;
  float fixed_timestep;
  float max_dimensionless_time;

  int markers,markers_per_ele;
	int marker_type[2];
	int marker_type_prev[2];
	int *element[2];
	int *element_prev[2];
  float marker_maxdist,marker_mindist;
  int markerIX,markerIZ;

 
  int min_timesteps;  
  int max_timesteps;
  int max_total_timesteps;
  int timesteps;
  int total_timesteps;
  int temp_iterations;
  int max_substeps;
  int sub_iterations;
  int last_sub_iterations; 

  float vel_substep_aggression;
  float temp_updatedness;
  float visc_updatedness;

  float lid_defining_velocity;
  float sub_layer_sample_level;

 
 } advection;


