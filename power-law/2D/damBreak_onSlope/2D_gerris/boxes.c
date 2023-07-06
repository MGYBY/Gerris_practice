#include<stdio.h>
 /*Creation of boxes a*b*c
COMMAND:
gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt
*/
void processorDist(int xInd, int yInd, int totalInd);

int main()
{
   int a, b, c, box,i ,j, connection;
   double x, y, lx, nd, nv, ly;
   a = 10;
   b = 3;
   x = 0.0;
   y = 0.0;
   nd = 0.380937;
   lx = nd*75.0*5.0*a/b;
   ly = nd*75.0*5.0;

//    scanf("%d%d%d",&a,&b,&c);


connection=0;

for ( j = 1; j <= b; j=j+1 ) {
	for ( i = 1; i <= a; i=i+1 ) {
		x = i*lx/a;
		y = j*ly/b;
		// processor 1 sweep -- bottom
		if (j==1) {
			// left BC
			if (i==1) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t left = Boundary { \n \
				\t     BcDirichlet P depth_bc(t) \n \
				\t     BcDirichlet U velocity_bc(t) \n \
				\t     BcDirichlet V 0.0 \n \
				\t } \n \
				\t bottom = Boundary \n \
				} \n");
			}
			// right BC
			else if (i==a) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t right = Boundary { \n \
				\t     BcNeumann P 0.0 \n \
				\t     BcNeumann U 0.0 \n \
				\t     BcNeumann V 0.0 \n \
				\t } \n \
				\t bottom = Boundary \n \
				} \n");
			}
			// middle boxes
			else {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t bottom = Boundary \n \
				} \n");
			}
		}

		//processor 2 sweep -- middle
		if (j==2) {
			// left BC
			if (i==1) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t left = Boundary { \n \
				\t     BcDirichlet P depth_bc(t) \n \
				\t     BcDirichlet U velocity_bc(t) \n \
				\t     BcDirichlet V 0.0 \n \
				\t } \n \
				} \n");
			}
			// right BC
			else if (i==a) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t right = Boundary { \n \
				\t     BcNeumann P 0.0 \n \
				\t     BcNeumann U 0.0 \n \
				\t     BcNeumann V 0.0 \n \
				\t } \n \
				} \n");
			}
			// middle boxes
			else {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				} \n");
			}
		}

		// processor 3 -- top
		if (j==3) {
			// left BC
			if (i==1) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t left = Boundary { \n \
				\t     BcDirichlet P depth_bc(t) \n \
				\t     BcDirichlet U velocity_bc(t) \n \
				\t     BcDirichlet V 0.0 \n \
				\t } \n \
				\t top = Boundary \n \
				} \n");
			}
			// right BC
			else if (i==a) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t right = Boundary { \n \
				\t     BcNeumann P 0.0 \n \
				\t     BcNeumann U 0.0 \n \
				\t     BcNeumann V 0.0 \n \
				\t } \n \
				\t top = Boundary \n \
				} \n");
			}
			// middle boxes
			else {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t top = Boundary \n \
				} \n");
			}
		}


	}
}
// }
// horizontal sweep
for( i = 1; i <= (a-1); i = i + 1 ) {
	for( j = 1; j <= b; j = j + 1 ) {
	box = ((j-1)*a+i);

	printf("%d %d right\n",box,box+1);
	connection=connection+1;

	}
}

// vertical sweep
for( i = 1; i <= a; i = i + 1 ) {
	for( j = 1; j <= (b-1); j = j + 1 ) {
	box = ((j-1)*a+i);

	printf("%d %d top\n",box,box+a);
	connection=connection+1;

	}
}

//  printf("\n\n\n\nThe Number of connection is: %d",connection);
   return 0;
}

void processorDist(int xInd, int yInd, int totalInd)
{
	// for now, the processor partition is only x-dependent
	double lc1, lc2;
	lc1 = 9.50/24.0*totalInd;
	lc2 = 19.0/24.0*totalInd;
	if (xInd<=lc1)
	{
		printf("\t pid = 0 \n");
	}
	else if (xInd<=lc2 && xInd>lc1)
	{
		printf("\t pid = 1 \n");
	}
	else
	{
		printf("\t pid = 2 \n");
	}
}
