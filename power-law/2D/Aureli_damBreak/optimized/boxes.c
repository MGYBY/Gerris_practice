#include<stdio.h>
 /*Creation of boxes a*b*c
COMMAND:
gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt
*/

void processorDist(int xInd, int yInd, int totalInd);

int main()
{
   int a, b, c, box,i ,j, connection, processorNum;
   a = 13;
   b = 6;
   processorNum = 3;
//    scanf("%d%d%d",&a,&b,&c);


connection=0;

for ( j = 1; j <= b; j=j+1 ) {
	for ( i = 1; i <= a; i=i+1 ) {
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
				\t left = Boundary \n \
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
				\t right = Boundary \n \
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
		if (j==2 || j==3 || j==4 || j==5) {
			// left BC
			if (i==1) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t left = Boundary \n \
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
				\t right = Boundary \n \
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
		if (j==6) {
			// left BC
			if (i==1) {
				printf("# No. %d box \n", ((j-1)*a+i));
				printf("GfsBox {\n");
				processorDist(i, j, a);
// 				printf(" \
// 				\t pid = 0 \n");
				printf("\
				\t left = Boundary \n \
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
				\t right = Boundary \n \
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
	lc1 = 9.0/24.0*totalInd;
	lc2 = 18.0/24.0*totalInd;
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
