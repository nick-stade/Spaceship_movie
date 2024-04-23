#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"


// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;
//to speed up rendering. 1 is normal, 10 is faster
int renspeed = 1;


int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}






void light_model (double irgb[3],
                  double *xx, double *yy, double *zz, 
                  double argb[3])
// irgb == inherent color of object (input to this function)
// xx[],yy[],zz[] are points in the polygon
// argb == actual color of object (output of this function)
{
  double Eye[3] ;
  Eye[0] = 0 ; Eye[1] = 0 ; Eye[2] = 0 ; 

  double P[3]  ;
  P[0] = xx[0] ;  P[1] = yy[0] ;  P[2] = zz[0] ;

  double a[3] ;
  a[0] = xx[1] - xx[0] ;  a[1] = yy[1] - yy[0] ;  a[2] = zz[1] - zz[0] ;

  double b[3] ;
  b[0] = xx[2] - xx[0] ;  b[1] = yy[2] - yy[0] ;  b[2] = zz[2] - zz[0] ;
 
  double N[3] ;
  M3d_x_product (N, a,b) ;

  Light_Model (irgb, Eye, P, N, argb) ;
}



double (*g[3]) (double u,double v) ;

void plot(double zbuff[][800], double (*f[3])(double u, double v),double transform[][4],
          double params[][2],double rgb[],int idf,double step){
    double P[3];
    double hither = .1;
    double H = tan(M_PI/3);
    double xx[3],yy[3],zz[3];
    double argb[3];
    double b = .005;
    int count = 0;
    int xt,yt;
    int e;
    MAX_DIFFUSE = 0.5;
    for(double u = params[0][0]; u <= params[0][1]; u+=step*renspeed){
        for(double v = params[1][0];v <= params[1][1];v+=step*renspeed){
          P[0] = f[0](u,v); P[1] = f[1](u,v); P[2] = f[2](u,v);
          M3d_mat_mult_pt(P,transform,P);
          if(P[2] == 0){P[2] = .0001;}
          xt = (int)(400 + 400*P[0]*(H/P[2])); 
          yt = (int)(400 + 400*P[1]*(H/P[2]));
          if(xt >= 0 && xt < 800 && yt >= 0 && yt < 800){
            if(P[2] < zbuff[yt][xt] && P[2] > hither){
                xx[0] = f[0](u,v);   yy[0] = f[1](u,v);   zz[0] = f[2](u,v);
                xx[1] = f[0](u+b,v); yy[1] = f[1](u+b,v); zz[1] = f[2](u+b,v);
                xx[2] = f[0](u,v+b); yy[2] = f[1](u,v+b); zz[2] = f[2](u,v+b);
                M3d_mat_mult_points(xx,yy,zz,transform,xx,yy,zz,3);
                light_model(rgb,xx,yy,zz,argb);
                e = set_xwd_map_color(idf, xt,yt, argb[0],argb[1],argb[2]) ;
                zbuff[yt][xt] = zz[0];
            }
          }
        }
    }
}

void plot_map(double zbuff[][800], double (*f[3])(double u, double v),double transform[][4],
              double params[][2],int id,int idf, double step){
    double P[3];
    double hither = .1;
    double H = tan(M_PI/3);
    double xx[3],yy[3],zz[3];
    double argb[3],rgb[3];
    double b = .005;
    int count = 0;
    int xt,yt;

    int ud,vd;
    int e,d[2] ;
    int width,height ;
    e = get_xwd_map_dimensions(id, d) ;
    if (e == -1) { printf("failure\n") ;  exit(0) ; }
    width = d[0] ; height = d[1] ;
    for(double u = params[0][0]; u <= params[0][1]; u+=step*renspeed){
        for(double v = params[1][0];v <= params[1][1];v+=step*renspeed){
          P[0] = f[0](u,v); P[1] = f[1](u,v); P[2] = f[2](u,v);
          M3d_mat_mult_pt(P,transform,P);
          if(P[2] == 0){P[2] = .0001;}
          xt = (int)(400 + 400*P[0]*(H/P[2]));
          yt = (int)(400 + 400*P[1]*(H/P[2]));
          if(xt >= 0 && xt < 800 && yt >= 0 && yt < 800){
            if(P[2] < zbuff[yt][xt] && P[2] > hither){
                xx[0] = f[0](u,v);   yy[0] = f[1](u,v);   zz[0] = f[2](u,v);
                xx[1] = f[0](u+b,v); yy[1] = f[1](u+b,v); zz[1] = f[2](u+b,v);
                xx[2] = f[0](u,v+b); yy[2] = f[1](u,v+b); zz[2] = f[2](u,v+b);
                M3d_mat_mult_points(xx,yy,zz,transform,xx,yy,zz,3);\

                ud = (int)((u-params[0][0])*width/(params[0][1]-params[0][0])); 
                vd = (int)((v-params[1][0])*height/(params[1][1]-params[1][0]));
                e = get_xwd_map_color(id, ud,vd,rgb) ; // returns -1 on error, 1 if ok
                if (e == -1) { continue ; }

                light_model(rgb,xx,yy,zz,argb);
                e = set_xwd_map_color(idf, xt,yt, argb[0],argb[1],argb[2]) ;
                zbuff[yt][xt] = zz[0];
            }
          }
        }
    }
}

void plot_map_nospec(double zbuff[][800], double (*f[3])(double u, double v),double transform[][4],
              double params[][2],int id,int idf, double step){
    double P[3];
    double hither = .1;
    double H = tan(M_PI/3);
    double xx[3],yy[3],zz[3];
    double argb[3],rgb[3];
    double b = .005;
    int count = 0;
    int xt,yt;
    MAX_DIFFUSE = .8;
    int ud,vd;
    int e,d[2] ;
    int width,height ;
    e = get_xwd_map_dimensions(id, d) ;
    if (e == -1) { printf("failure\n") ;  exit(0) ; }
    width = d[0] ; height = d[1] ;
    for(double u = params[0][0]; u <= params[0][1]; u+=step*renspeed){
        for(double v = params[1][0];v <= params[1][1];v+=step*renspeed){
          P[0] = f[0](u,v); P[1] = f[1](u,v); P[2] = f[2](u,v);
          M3d_mat_mult_pt(P,transform,P);
          if(P[2] == 0){P[2] = .0001;}
          xt = (int)(400 + 400*P[0]*(H/P[2]));
          yt = (int)(400 + 400*P[1]*(H/P[2]));
          if(xt >= 0 && xt < 800 && yt >= 0 && yt < 800){
            if(P[2] < zbuff[yt][xt] && P[2] > hither){
                xx[0] = f[0](u,v);   yy[0] = f[1](u,v);   zz[0] = f[2](u,v);
                xx[1] = f[0](u+b,v); yy[1] = f[1](u+b,v); zz[1] = f[2](u+b,v);
                xx[2] = f[0](u,v+b); yy[2] = f[1](u,v+b); zz[2] = f[2](u,v+b);
                M3d_mat_mult_points(xx,yy,zz,transform,xx,yy,zz,3);\

                ud = (int)((u-params[0][0])*width/(params[0][1]-params[0][0])); 
                vd = (int)((v-params[1][0])*height/(params[1][1]-params[1][0]));
                e = get_xwd_map_color(id, ud,vd,rgb) ; // returns -1 on error, 1 if ok
                if (e == -1) { continue ; }

                light_model(rgb,xx,yy,zz,argb);
                e = set_xwd_map_color(idf, xt,yt, argb[0],argb[1],argb[2]) ;
                zbuff[yt][xt] = zz[0];
            }
          }
        }
    }
    MAX_DIFFUSE = .5;
}

double spherex(double u,double v){return sin(v)*cos(u);}
double spherey(double u,double v){return sin(v)*sin(u);}
double spherez(double u,double v){return cos(v);}
double torx(double u, double v){return cos(u)*(3+cos(v));}
double tory(double u, double v){return sin(u)*(3+cos(v));}
double torz(double u, double v){return sin(v);}

double squarex(double u){
  if(u >= 7*M_PI/4 || u <= M_PI/4){
    return 1;
  }else if(u >= M_PI/4 && u <= 3*M_PI/4){
    return cos(2*(u-(M_PI/4)));
  }else if(u >= 5*M_PI/4 && u <= 7*M_PI/4){
    return cos(2*(u-(5*M_PI/4)));
  }else return -1;
}
double squarey(double u){
    if(u >= 7*M_PI/4 || u <= M_PI/4){
    return cos(2*(u+M_PI/4));
  }else if(u >= M_PI/4 && u <= 3*M_PI/4){
    return 1;
  }else if(u >= 3*M_PI/4 && u <= 5*M_PI/4){
    return cos(2*(u-(3*M_PI/4)));
  }else return -1;
}

double absqx(double u){
  if(u > 8){u = fmod(u,8);}
  if(u >= 0 && u < 2){
    return 1 - u;
  }else if(u >= 2 && u < 4){
    return -1;
  }else if(u >= 4 && u < 6){
    return (u - 4) - 1;
  }else if(u >= 6 && u <= 8){
    return 1;
  }
}

double absqy(double u){
  if(u > 8){u = fmod(u,8);}
  if(u >= 0 && u < 2){
    return 1;
  }else if(u >= 2 && u < 4){
    return 1 - (u - 2);
  }else if(u >= 4 && u < 6){
    return -1;
  }else if(u >= 6 && u <= 8){
    return u - 7;
  }
}


double washerx(double u, double v){return cos(u)*(10+absqx(v));}
double washery(double u, double v){return sin(u)*(10+absqx(v));}
double washerz(double u, double v){return absqy(v);}

double collumnx(double u, double v){return absqx(u);}
double collumny(double u, double v){return absqy(u);}
double collumnz(double u, double v){return v;}

double retu(double u, double v){return u;}
double retv(double u, double v){return v;}
double ret0(double u, double v){return 0;}
double cosu(double u, double v){return cos(u);}
double sinu(double u, double v){return sin(u);}
double cosv(double u, double v){return cos(v);}
double sinv(double u, double v){return sin(v);}

double circx(double u, double v){return v*cos(u);}
double circy(double u, double v){return v*sin(u);}

double conex(double u, double v){return cos(u)*v;}
double coney(double u, double v){return sin(u)*v;}


