#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

double (*g[2]) (double u) ;

void plot(double (*fx)(double u),double (*fy)(double),double transform[][4],double params[]){
	double x[10000],y[10000],z[10000];
	int count = 0;
	for(double i = params[0];i <= params[1];i+=.01){
		x[count] = fx(i);y[count] = fy(i);z[count]=0;
		count++;
	}
	G_rgb(1,1,1);
	M3d_mat_mult_points(x,y,z,transform,x,y,z,count);
	for(int i = 0; i < count; i++){
		G_point(x[i],y[i]);
	}
}

void plotnormal(double (*fx)(double u),double (*fy)(double),double transform[][4],double params[]){
    double x[10000],y[10000],z[10000];
    double dx,dy,len;
    int count = 0;
    for(double i = params[0];i <= params[1];i+=.01){
        x[count] = fx(i);y[count] = fy(i);z[count]=0;
        count++;
    }
    G_rgb(1,1,1);
    M3d_mat_mult_points(x,y,z,transform,x,y,z,count);
    for(int i = 0; i < count; i++){
        G_point(x[i],y[i]);
        if(i%10 == 0 && i != 0){
            dx = x[i] - x[i-1];
            dy = y[i] - y[i-1];
            len = sqrt(pow(dx,2)+pow(dy,2));
            G_line(x[i],y[i],x[i]+((20*dy)/len),y[i]-((20*dx)/len));
        }
    }
}

double sgn(double u){
    if(u < 0){return -1;}
    else if(u > 0){return 1;}
    return 0;
}
double sum4(double u){return pow(1-u*u*u*u, 0.25);}
double retur(double u){return u;}
double sgx(double u){return sgn(cos(u))*pow((cos(u)),2);}
double sgy(double u){return sgn(sin(u))*pow((sin(u)),2);}
double asgx(double u){return sgn(cos(u))*pow(cos(u),4);}
double asgy(double u){return sgn(sin(u))*pow(sin(u),4);}
double lemx(double u){return pow(cos(u),3);}
double archx(double u){return u+sin(u);}
double sq(double u){return u*u;}

double absqx(double u){
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


int main(){
	G_init_graphics(800,800);
	G_rgb(0,0,0);
	G_clear();
	double v[4][4],vi[4][4];
	int Tn;
	int Ttypelist[100];double Tvlist[100];
	double us[2];
	G_rgb(1,1,1);

	//circle
	Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = cos; g[1] = sin;
  us[0] = .25*M_PI;us[1] = 1.5*M_PI;
  plotnormal(g[0],g[1],v,us);

  //sum4
  Tn = 0;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  170.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = retur;g[1] = sum4;
  us[0] = -1;us[1] = 1;
  plotnormal(g[0],g[1],v,us);

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  670.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = absqx;g[1] = absqy;
  us[0] = 0; us[1] = 8;
  plotnormal(g[0],g[1],v,us);

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  460.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = sgx;g[1] = sgy;
  us[0] = 0; us[1] = 2*M_PI;
  plotnormal(g[0],g[1],v,us);

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   80.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   45.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  130.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  650.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = asgx;g[1] = asgy;
  us[0] = 0;us[1] = 2*M_PI;
  plotnormal(g[0],g[1],v,us);

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = cosh;g[1] = sinh;
  us[0] = -1;us[1] = 1.5;
  plotnormal(g[0],g[1],v,us);

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  125.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  125.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  620.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  210.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = lemx;g[1] = sin;
  us[0] = 0;us[1] = 2*M_PI;
  plotnormal(g[0],g[1],v,us);

Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  140.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  200.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = retur;g[1] = sq;
  us[0] = -1; us[1] = 2;
    plotnormal(g[0],g[1],v,us);

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  15.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  15.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  130.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  50.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix(v,vi,Tn,Ttypelist,Tvlist);
  g[0] = archx; g[1] = cos;
  us[0] = M_PI; us[1] = 7*M_PI;
  plotnormal(g[0],g[1],v,us);

  G_wait_key();
}



