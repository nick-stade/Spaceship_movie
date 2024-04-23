#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
char type[100];
double params[100][2][2];
int    num_objects ;



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////





void Draw_ellipsoid (int onum)
{
  int n,i,j ;
  double p,t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    for(j = 0; j < n; j++){
      t = i*2*M_PI/n ;
      p = j*M_PI/n;
      xyz[0] = cos(t)*sin(p) ;
      xyz[1] = sin(t)*sin(p) ;
      xyz[2] = cos(p);
      M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
      x = xyz[0] ;
      y = xyz[1] ;
      G_point(x,y) ;
    }
  }

}




void Draw_the_scene()
{
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    Draw_ellipsoid(onum) ;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



double intersect_hyperbola(double dx, double dy, double S[], double uv[],double params[]){
  double a,b,c,x,y;
  a = pow(dx,2) - pow(dy,2);
  b = 2*S[0]*dx - 2*S[1]*dy;
  c = pow(S[0],2) - pow(S[1],2) - 1;
  double t1,t2;
  if(pow(b,2) - 4*a*c < 0){return -1;}
  t1 = (-b + sqrt(pow(b,2) - 4*a*c))/(2*a);
  t2 = (-b - sqrt(pow(b,2) - 4*a*c))/(2*a);
  if(t1 > t2 && t2 > 0){
    x = S[0] + t2*dx;
    y = S[1] + t2*dy;
    uv[0] = atan2(y,x);
    return t2;
  }else if(t1 > 0){
    x = S[0] + t1*dx;
    y = S[1] + t1*dy;
    uv[0] = atan2(y,x);
    return t1;
  }else{return -1;}
}



//return t
double intersect_sphere(double dx, double dy,double dz,double S[],double uv[]){
  double a,b,c,x,y,z;
  a = pow(dx,2) + pow(dy,2) + pow(dz,2);
  b = 2*S[0]*dx + 2*S[1]*dy + 2*S[2]*dz;
  c = pow(S[0],2) + pow(S[1],2) + pow(S[2],2) - 1;
  double t1,t2;
  if(pow(b,2) - 4*a*c < 0){return -1;}
  t1 = (-b + sqrt(pow(b,2) - 4*a*c))/(2*a);
  t2 = (-b - sqrt(pow(b,2) - 4*a*c))/(2*a);
  if(t1 > t2 && t2 > 0){
    x = S[0] + t2*dx;
    y = S[1] + t2*dy;
    z = S[2] + t2*dz;
    uv[0] = atan2(y,x);
    uv[1] = atan2(sqrt(pow(x,2) + pow(y,2)),z);
    return t2;
  }else if(t1 > 0){
    x = S[0] + t1*dx;
    y = S[1] + t1*dy;
    z = S[2] + t1*dz;
    uv[0] = atan2(y,x);
    uv[1] = atan2(sqrt(pow(x,2) + pow(y,2)),z);
    return t1;
  }else{return -1;}
}

void normal_usphere(double uv[], int onum,double result[],double obpt[]){
  double dzdx,dzdy,dzdp[3];
  obpt[0] = cos(uv[0])*sin(uv[1]); obpt[1] = sin(uv[0])*sin(uv[1]); obpt[2] = cos(uv[1]);
  dzdx = 2*obpt[0]; dzdy = 2*obpt[1]; dzdx = 2*obpt[2];
  dzdp[0] = obinv[onum][0][0]*dzdx + obinv[onum][1][0]*dzdy;
  dzdp[1] = obinv[onum][0][1]*dzdx + obinv[onum][1][1]*dzdy;
  dzdp[2] = 0;
  M3d_normalize(dzdp,dzdp);
  result[0] = dzdp[0];
  result[1] = dzdp[1];
  result[2] = dzdp[2];
}

double intersect(double dx, double dy,double dz, double S[], double uv[],int onum){
  if(type[onum] == 'S'){
    return intersect_sphere(dx,dy,dz,S,uv);
  }
}

void rayint (double Rsource[3], double Rtip[3]){
  //pRres is the resulting parameters that define a point on an object. in 2d this is one parameter. in 3d this is u and v
  double dx,dy,dz,pRes[num_objects][2],point[3],t[num_objects],tfinal;
  int closest;
  double nsource[3], ntip[3];
  tfinal = 10000;closest = -1;
  for(int i = 0; i < num_objects; i++){
    M3d_mat_mult_pt(nsource,obinv[i],Rsource);
    M3d_mat_mult_pt(ntip,obinv[i],Rtip);
    dx = ntip[0] - nsource[0];
    dy = ntip[1] - nsource[1];
    dz = ntip[2] - nsource[2];
    t[i] = intersect(dx,dy,dz,nsource,pRes[i],i);
    if(t[i] < tfinal && t[i] != -1){tfinal = t[i];closest = i;}
  }
  double uv[2];
  double normalvector[3];
  if(closest != -1){
    uv[0] = pRes[closest][0]; uv[1] = pRes[closest][1];
    // normal_usphere(pRes[closest],closest,normalvector,point);
    // normalvector[0] *= 30;
    // normalvector[1] *= 30;

    point[0] = cos(uv[0])*sin(uv[1]); point[1] = sin(uv[0])*sin(uv[1]); point[2] = cos(uv[1]);
    M3d_mat_mult_pt(point,obmat[closest],point);
    G_rgb(color[closest][0],color[closest][1],color[closest][2]) ;
    //G_line(Rtip[0],Rtip[1],point[0],point[1]);
    //G_line(point[0],point[1],point[0]+dzdp[0],point[1]+dzdp[1]);
    //G_fill_circle(Rtip[1],Rtip[2],2);
    G_point(Rtip[0],Rtip[1]);
    printf("%lf, %lf\n",Rtip[0],Rtip[1]);
  }
}



int test01()
{
  double vm[4][4], vi[4][4];
  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  double Rsource[3];
  double Rtip[3];
  double argb[3] ;

    //////////////////////////////////////////////////////////////////////
    M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
    //////////////////////////////////////////////////////////////////////

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.0 ;
    color[num_objects][1] = 0.8 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  70   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   25   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  200   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  1000   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;
    type[num_objects] = 'S';

    num_objects++ ; // don't forget to do this

    /*//////////////////////////////////////////////////////////////
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  180   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   40   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  400   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  550   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.3 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   75   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   35   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  150   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  360   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  500   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  130   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   30   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -15   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  700   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this        
    /////////////////////////////////////////////////////////////*/

    

    G_rgb(0,0,0) ;
    G_clear() ;

    
    
    Rsource[0] =  0 ;  Rsource[1] =  0 ;  Rsource[2] = -400/tan(30) ;   
    
    double ytip,xtip;
    for (ytip = -400 ; ytip <= 400 ; ytip++) {
      for (xtip = -400; xtip <= 400 ; xtip++) {
        Rtip[0]    = xtip ;  Rtip[1]  = ytip ;  Rtip[2]   =  0  ;    
        rayint(Rsource,Rtip);
      }
    }

    G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
    while (G_wait_key() != 'q') ;
    G_save_image_to_file("2d_Simple_Raytracer.xwd") ;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(800,800);
  test01() ;
}