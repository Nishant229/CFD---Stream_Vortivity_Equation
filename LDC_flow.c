#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

int Velocity_profile(int M,int N,double U,double Re,double dx,double dy,double u[N+1][M+1],double v[N+1][M+1],double psi[N+1][M+1],double omega[N+1][M+1]);

int main() {
	//Reading Problem data
	double L=1.0, H=1.0;

	//Taking Inputs from user
	int M,N;
	/*printf("Please enter the Number of Grid Points on X-axis: ");
	scanf("%d",&M);
	printf("Please enter the Number of Grid Points on X-axis: ");
	scanf("%d",&N);*/
	M=129, N=129;
	double dx=L/(M-1), dy=H/(N-1);
	double U1=1, U2=0,U3=0,U4=0,U=U1;
	double V1=0, V2=0,V3=0,V4=0;
	double Re[]={100,400};
	int n=sizeof(Re)/sizeof(Re[0]);
	
	double u[N+1][M+1];
	double v[N+1][M+1];
	double psi[N+1][M+1];
	double omega[N+1][M+1];
	
	// Initialization
	for (int j=2;j<=N-1;j++) {
		for (int i=2;i<=M-1;i++) {
			u[j][i] = 0;
			v[j][i] = 0;
			psi[j][i] = 0;
			omega[j][i] = 0;
		}
	}
	
	//// Boundary Conditions
	//Left Boundary
	for (int j=1;j<=N;j++) {
		u[j][1] = U4;
		v[j][1] = V4;
		psi[j][1] = 0;
		omega[j][1] = -2*(psi[j][2] - psi[j][1])/(dx*dx);
	}
	//Right Boundary
	for (int j=1;j<=N;j++) {
		u[j][M] = U2;
		v[j][M] = V2;
		psi[j][M] = 0;
		omega[j][M] = -2*(psi[j][M-1] - psi[j][M])/(dx*dx);
	}
	//Top Boundary
	for (int i=1;i<=M;i++) {
		u[N][i] = U1;
		v[N][i] = V1;
		psi[N][i] = 0;
		omega[N][i] =  -2*(psi[N-1][i] - psi[N][i] + U*dy)/(dy*dy);
	}
	//Bottom Boundary
	for (int i=1;i<=M;i++) {
		u[1][i] = U3;
		v[1][i] = V3;
		psi[1][i] = 0;
		omega[1][i] =  -2*(psi[2][i] - psi[1][i])/(dy*dy);
	}
	
	clock_t t;
	double time_taken[n];
	for(int i=0;i<n;i++) {
		t = clock();
		//Stream Vorticity Method
		int iterations = Velocity_profile(M,N,U,Re[i],dx,dy,u,v,psi,omega);
		t = clock() - t;
		printf("Re = %lf\n",Re[i]);
		time_taken[i] = ((double)t)/CLOCKS_PER_SEC;
		printf("The program took %lf seconds to execute\n",time_taken[i]);
	}
	for(int i=0;i<n;i++) {
		printf("The program took %lf seconds to execute for Re = %lf\n",time_taken[i],Re[i]);
	}
	
	return 0;
}

int Velocity_profile(int M,int N,double U,double Re,double dx,double dy,double u[N+1][M+1],double v[N+1][M+1],double psi[N+1][M+1],double omega[N+1][M+1]) {
	double psi_old[N+1][M+1], omega_old[N+1][M+1];
	// Initialization
	for (int j=1;j<=N;j++) {
		for (int i=1;i<=M;i++) {
			psi_old[j][i] = psi[j][i];
			omega_old[j][i] = omega[j][i];
		}
	}
	
	FILE *fp1;
	char fname[100] = "Stream_Vorticity_error_RE_";
	char lname[100] = ".dat";
	char str[100];
	sprintf(str,"%.2f",Re);
	strcat(fname,str);
	strcat(fname,lname);
	fp1=fopen(fname,"w");
	fprintf(fp1,"NOI \t psi_error \t omega_error\n\n");
	
	double psi_error  = 1;
	double omega_error  = 1;
	double e = pow(10,-6);
	double beta=dx/dy;
	int NOI=0;
	while (psi_error>e || omega_error>e) {
		omega_error = 0;
		psi_error = 0;
		
		//Calculating omega
		for (int j=2;j<=N-1;j++) {
			for(int i=2;i<=M-1;i++) {
				omega[j][i] = (0.5/(1+pow(beta,2))) * ( (1-(psi[j+1][i]-psi[j-1][i])*(beta*Re/4))*(omega[j][i+1]) + (1+(psi[j+1][i]-psi[j-1][i])*(beta*Re/4))*(omega[j][i-1]) + (1+(psi[j][i+1]-psi[j][i-1])*(Re/(4*beta)))*(pow(beta,2)*omega[j+1][i]) + (1-(psi[j][i+1]-psi[j][i-1])*(Re/(4*beta)))*(pow(beta,2)*omega[j-1][i]) );
				omega_error = omega_error + pow((omega[j][i] - omega_old[j][i]),2 );
			}
		}
		
		//calculating psi
		for(int j=2;j<=N-1;j++) {
			for (int i=2;i<=M-1;i++) {
				psi[j][i] =  (0.5/(1+pow(beta,2))) * ( omega[j][i]*(dx*dx) + pow(beta,2)*(psi[j+1][i] + psi[j-1][i]) + (psi[j][i+1] + psi[j][i-1]) );
				psi_error = psi_error + pow((psi[j][i] - psi_old[j][i]),2 );
			}
		}
		
		//Updating Boundary Conditions
		//Left Boundary
		for (int j=1;j<=N;j++) {
			omega[j][1] = -2*(psi[j][2] - psi[j][1])/(dx*dx);
		}
		//Right Boundary
		for (int j=1;j<=N;j++) {
			omega[j][M] = -2*(psi[j][M-1] - psi[j][M])/(dx*dx);
		}
		//Top Boundary
		for (int i=1;i<=M;i++) {
			omega[N][i] =  -2*(psi[N-1][i] - psi[N][i] + U*dy)/(dy*dy);
		}
		//Bottom Boundary
		for (int i=1;i<=M;i++) {
			omega[1][i] =  -2*(psi[2][i] - psi[1][i])/(dy*dy);
		}
		
		//updating the old loop
		for (int j=1;j<=N;j++) {
			for (int i=1;i<=M;i++) {
				psi_old[j][i] = psi[j][i];
				omega_old[j][i] = omega[j][i];
			}
		}
		
		NOI++;
		
		//Error Calculation
		psi_error=pow(psi_error/((M-2)*(N-2)),0.5);
		omega_error=pow(omega_error/((M-2)*(N-2)),0.5);
		printf("%d -\t %.8f \t %.8f\n",NOI,psi_error,omega_error);
		fprintf(fp1,"%d \t %.10lf \t %.10lf\n",NOI, psi_error,omega_error);
		
		
	}
	
	for (int j=2;j<=N-1;j++) {
		for (int i=2;i<=M-1;i++) {
			u[j][i] = (psi[j+1][i] - psi[j-1][i] ) / (2*dy);
			v[j][i] = -(psi[j][i+1] - psi[j][i-1] ) / (2*dx);
		}
	}
	
	////	Correcting the corner points
	u[1][N]=0.5*(u[1][N-1]+u[2][N]); // Topmost left corner
	u[N][M]=0.5*(u[N][M-1]+u[N-1][M]); // topmost right corner
	
	//Printing the final Iterations and Error
	printf("Stream Vorticity Equation: \n");
	printf("No. Of Iterations - %d \n psi_Error %lf\n omega_Error: %lf\n",NOI,psi_error,omega_error);
	
	
	
	//Vertical Centerline: u_velocity
	FILE *fp2;
	char f1name[100] = "u_centerline_Velocity_RE_";
	char l1name[100] = ".dat";
	strcat(f1name,str);
	strcat(f1name,l1name);
	fp2=fopen(f1name,"w");
	for (int j=1;j<=N;j++) {
		fprintf(fp2,"%lf\t%lf\n",u[j][51],(j-1)*dy);
	}
	fclose(fp2);
	
	//Vertical Centerline: u_velocity
	FILE *fp3;
	char f2name[100] = "v_centerline_Velocity_RE_";
	char l2name[100] = ".dat";
	strcat(f2name,str);
	strcat(f2name,l2name);
	fp3=fopen(f2name,"w");
	for (int i=1;i<=M;i++) {
		fprintf(fp3,"%lf\t%lf\n",(i-1)*dx,v[51][i]);
	}
	fclose(fp3);
	
	
	//Results
	double x,y;
	FILE *fp4;
	char f3name[100] = "Result_RE_";
	char l3name[100] = ".dat";
	strcat(f3name,str);
	strcat(f3name,l3name);
	fp4=fopen(f3name,"w");
	fprintf(fp4,"ZONE I=%d, J=%d\n",M,N);
	for(int j=1;j<=N;j++) {
		y=(j-1)*(dy);
		for(int i=1;i<=M;i++) {
			x=(i-1)*dx;
			fprintf(fp4,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[j][i],v[j][i],psi[j][i],omega[j][i]);
		}
	}
	fclose(fp4);
	
	return NOI;
	
	
}
