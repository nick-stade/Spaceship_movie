#include "drawing3d.c"

int main(){
  int start = 0, end = 200;
  int idf;
  int w,h;
  w = 800; h = 800;
  
  char frame[100];
  renspeed = 1;

  idf = create_new_xwd_map (w,h) ;// returns -1 on error, 1 if ok
      if (idf == -1) { printf("failure\n") ;  exit(0) ; }

  int earth = init_xwd_map_from_file ("earthJ.xwd") ;// returns -1 on error, 1 if ok

  int cab = init_xwd_map_from_file ("cab2J.xwd") ;// returns -1 on error, 1 if ok

	int dstar = init_xwd_map_from_file ("deathstarJ.xwd") ;// returns -1 on error, 1 if ok

  int cab1 = init_xwd_map_from_file ("cab1J.xwd") ;// returns -1 on error, 1 if ok

  int ship = init_xwd_map_from_file ("space2J.xwd") ;// returns -1 on error, 1 if ok

  int stars = init_xwd_map_from_file ("starsJ.xwd") ;// returns -1 on error, 1 if ok


  double zbuff[800][800];


	//for light:
	double xx[3],yy[3],zz[3];
  double light[3];
	double rgb[3];
	double v[4][4],vi[4][4],vt[4][4],vti[4][4];
	int Tn;
	int Ttypelist[100];double Tvlist[100];
	double us[2][2];
	double t;int fnum = start;
	double eye[3], coi[3], up[3] ;
	char key;
  double theta, step,shipstep;
  
	while (fnum <= end) {
		for(int i = 0; i < 800; i++){for(int j = 0 ; j< 800; j++){zbuff[i][j] = 1000000;}}
	    t = 0.01*fnum ;

      shipstep = .015-.013*t;
      if(t > 1.1){shipstep  = 10;}
      eye[0] = -6*t + 10 - pow(t,2.5) - pow(t-.5,10);
      eye[1] = -6*t + 10 - pow(t,2.5) - pow(t-.5,10);
      eye[2] = t*40 - 40;

      coi[0] = t*10;
      coi[1] = pow(t,2.5) + pow(t-.7,15);
      coi[2] = 200 + t*10;

	    up[0]  = eye[0]; 
	    up[1]  = eye[1] + 1 ;
	    up[2]  = eye[2];
	    M3d_view(v,vi,eye,coi,up);
      theta = t*90;

      clear_xwd_map(idf,0,0,0);
//coi
      // Tn = 0 ; 
      // Ttypelist[Tn] = TX ; Tvlist[Tn] =  t*10; ; Tn++ ;
      // Ttypelist[Tn] = TY ; Tvlist[Tn] =  pow(t,2.5) + pow(t-.7,15);  ; Tn++ ;
      // Ttypelist[Tn] = TZ ; Tvlist[Tn] =  200+t*10 ;  Tn++ ;
      // M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
      // M3d_mat_mult(v,v,vt);
      // g[0] = spherex; g[1] = spherey; g[2] = spherez;
      // us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = M_PI;
      // step = .003;
      // rgb[0] = .5; rgb[1] = 0; rgb[2] = .5;
      // plot(zbuff,g,v,us,rgb,idf,step);
      // M3d_mat_mult(v,v,vti);

      //earth
      Tn = 0 ; 
      Ttypelist[Tn] = SX ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = SY ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = RX ; Tvlist[Tn] =  60.0 ; Tn++ ;
      Ttypelist[Tn] = RZ ; Tvlist[Tn] =  30.0 ; Tn++ ;
      Ttypelist[Tn] = RY ; Tvlist[Tn] =  theta/4 ; Tn++ ;
      Ttypelist[Tn] = TZ ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = TY ; Tvlist[Tn] =  800.0 ; Tn++ ;
      Ttypelist[Tn] = TX ; Tvlist[Tn] =  800.0 ; Tn++ ;
      M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
      M3d_mat_mult(v,v,vt);
      g[0] = spherex; g[1] = spherey; g[2] = spherez;
      us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = M_PI;
      step = .00025;
      plot_map_nospec(zbuff,g,v,us,earth,idf,step);
      M3d_mat_mult(v,v,vti);


      //deathstar
      Tn = 0 ; 
      Ttypelist[Tn] = SX ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = SY ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = RX ; Tvlist[Tn] =  90.0 ; Tn++ ;
      Ttypelist[Tn] = RY ; Tvlist[Tn] =  -theta/2 ; Tn++ ;
      Ttypelist[Tn] = TZ ; Tvlist[Tn] =  10000.0 ; Tn++ ;
      Ttypelist[Tn] = TY ; Tvlist[Tn] =  1000.0 ; Tn++ ;
      Ttypelist[Tn] = TX ; Tvlist[Tn] = -1500.0 ; Tn++ ;
      Ttypelist[Tn] = RY ; Tvlist[Tn] =  -theta/20 + 16 ; Tn++ ;
      Ttypelist[Tn] = RX ; Tvlist[Tn] =  theta/15 - 16 ; Tn++ ;
      M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
      M3d_mat_mult(v,v,vt);
      g[0] = spherex; g[1] = spherey; g[2] = spherez;
      us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = M_PI;
      step = .005;
      plot_map_nospec(zbuff,g,v,us,dstar,idf,step);
      M3d_mat_mult(v,v,vti);


      //stars
      Tn = 0 ; 
      Ttypelist[Tn] = SX ; Tvlist[Tn] =  10000.0 ; Tn++ ;
      Ttypelist[Tn] = SY ; Tvlist[Tn] =  10000.0 ; Tn++ ;
      Ttypelist[Tn] = SZ ; Tvlist[Tn] =  10000.0 ; Tn++ ;
      Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -theta/16 ; Tn++ ;
      Ttypelist[Tn] = RY ; Tvlist[Tn] =  theta/16 ; Tn++ ;
      Ttypelist[Tn] = RX ; Tvlist[Tn] =  theta/24; Tn++ ;
      M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
      M3d_mat_mult(v,v,vt);
      g[0] = spherex; g[1] = spherey; g[2] = spherez;
      us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = M_PI;
      step = .001;
      plot_map_nospec(zbuff,g,v,us,stars,idf,step);
      M3d_mat_mult(v,v,vti);

      if(fnum < 110){
        //cone
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  2/2.8 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  2/2.8 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = conex; g[1] = coney; g[2] = retv;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = -2.8; us[1][1] = 2.8;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep);  
        M3d_mat_mult(v,v,vti);
        //main cylinder
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1.0 ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = cosu; g[1] = sinu; g[2] = retv;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = -2.8; us[1][1] = 2.8;
        rgb[0] = .7; rgb[1] = .7; rgb[2] = .7;
        step = .01;
        plot(zbuff,g,v,us,rgb,idf,shipstep);
        M3d_mat_mult(v,v,vti);


        //positive half
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = washerx; g[1] = washery; g[2] = washerz;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = 8;
        step = .005;
        plot_map(zbuff,g,v,us,cab,idf,shipstep*.1);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  4.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] = -9.0 ; Tn++ ;
        Ttypelist[Tn] = RY ; Tvlist[Tn] =  90.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = collumnx; g[1] = collumny; g[2] = collumnz;
        us[0][0] = 0; us[0][1] = 8; us[1][0] = 0; us[1][1] = 4.5;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep/4);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  4.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] = -9.0 ; Tn++ ;
        Ttypelist[Tn] = RX ; Tvlist[Tn] =  90.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = collumnx; g[1] = collumny; g[2] = collumnz;
        us[0][0] = 0; us[0][1] = 8; us[1][0] = 0; us[1][1] = 4.5;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep/4);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = cosu; g[1] = sinu; g[2] = retv;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = -.7; us[1][1] = .7;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  .5 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  4.2 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = spherex; g[1] = spherey; g[2] = spherez;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = M_PI;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep);
        M3d_mat_mult(v,v,vti);



        //negative half
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = washerx; g[1] = washery; g[2] = washerz;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = 8;
        step = .005;
        plot_map(zbuff,g,v,us,cab1,idf,shipstep*.1);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  4.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] = -9.0 ; Tn++ ;
        Ttypelist[Tn] = RY ; Tvlist[Tn] =  90.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = collumnx; g[1] = collumny; g[2] = collumnz;
        us[0][0] = 0; us[0][1] = 8; us[1][0] = 0; us[1][1] = 4.5;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep/4);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.7 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  4.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] = -9.0 ; Tn++ ;
        Ttypelist[Tn] = RX ; Tvlist[Tn] =  90.0 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = collumnx; g[1] = collumny; g[2] = collumnz;
        us[0][0] = 0; us[0][1] = 8; us[1][0] = 0; us[1][1] = 4.5;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep/4);      
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -3.5 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = cosu; g[1] = sinu; g[2] = retv;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = -.7; us[1][1] = .7;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep);
        M3d_mat_mult(v,v,vti);
        Tn = 0 ; 
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  2 ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  .5 ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -4.2 ; Tn++ ;
        Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -theta  ; Tn++ ;
        M3d_make_movement_sequence_matrix(vt,vti,Tn,Ttypelist,Tvlist);
        M3d_mat_mult(v,v,vt);
        g[0] = spherex; g[1] = spherey; g[2] = spherez;
        us[0][0] = 0; us[0][1] = 2*M_PI; us[1][0] = 0; us[1][1] = M_PI;
        step = .01;
        plot_map(zbuff,g,v,us,ship,idf,shipstep);
        M3d_mat_mult(v,v,vti);
      }
	  	
	  	
      sprintf(frame,"im%04d.xwd",fnum);
      xwd_map_to_named_xwd_file(idf, frame) ;
      printf("frame %d rendered\n",fnum);
      fnum++;
	  }
}