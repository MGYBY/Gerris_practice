#include<stdio.h>
 /*Creation of boxes a*b*c
COMMAND:
gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt
*/
int main()
{
   int a, b, c, box,i ,j, connection, processorNum;
   a = 5;
   b = 3;
   processorNum = 3; // a dummy variable here
//    scanf("%d%d%d",&a,&b,&c);


connection=0;

for ( i = 0; i < a*b; i=i+1 ) {
	if(i<b)
	{
		printf("# No. %d box \n", (i+1));
		printf("GfsBox {\n");
		if (i%b == 0) {
		printf(" \
\t pid = 0 \n");
		}
		else if (i%b == 1) {
		printf(" \
\t pid = 1 \n");
		}
		else{
			printf(" \
\t pid = 2 \n");
		}
		printf("\
\t left = Boundary { \n \
\t     BcDirichlet P depth_bc(t) \n \
\t     BcDirichlet U velocity_bc(t) \n \
\t     BcDirichlet V 0.0 \n \
\t } \n \
\t top = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t } \n \
\t bottom = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t} \n \
} \n");
	}
	else if (i>(a*b-b-1)) {
		printf("# No. %d box \n", (i+1));
		printf("GfsBox {\n");
		if (i%b == 0) {
		printf(" \
\t pid = 0 \n");
		}
		else if (i%b == 1) {
		printf(" \
\t pid = 1 \n");
		}
		else{
			printf(" \
\t pid = 2 \n");
		}
		printf("\
\t right = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t } \n \
\t top = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t } \n \
\t bottom = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t} \n \
} \n");
	}
	else
	{
		printf("# No. %d box \n", (i+1));
		printf("GfsBox {\n");
		if (i%b == 0) {
		printf(" \
\t pid = 0 \n");
		}
		else if (i%b == 1) {
		printf(" \
\t pid = 1 \n");
		}
		else{
			printf(" \
\t pid = 2 \n");
		}
		printf("\
\t top = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t } \n \
\t bottom = Boundary { \n \
\t     BcNeumann P 0.0 \n \
\t     BcNeumann U 0.0 \n \
\t     BcNeumann V 0.0 \n \
\t} \n \
} \n");
	}
}
// }
// horizontal sweep
for( i = 1; i <= (a-1); i = i + 1 ) {
	for( j = 1; j <= b; j = j + 1 ) {
	box = b*(i-1)+j;

	printf("%d %d right\n",box,box+b);
	connection=connection+1;

	}
}

// vertical sweep
for( i = 1; i <= a; i = i + 1 ) {
	for( j = 1; j <= (b-1); j = j + 1 ) {
	box = b*(i-1)+j;

	printf("%d %d top\n",box,box+1);
	connection=connection+1;

	}
}

//  printf("\n\n\n\nThe Number of connection is: %d",connection);
   return 0;
}
