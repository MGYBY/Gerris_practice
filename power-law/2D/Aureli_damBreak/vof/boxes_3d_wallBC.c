#include<stdio.h>
 /*Creation of boxes a*b*c
COMMAND:
gcc boxes.c -o boxes.out && ./boxes.out - > boxes.txt
*/
int main()
{
int a, b, c, box,i ,j, k, connection, processorNum;
a = 104;
b = 9;
c = 24;
processorNum = 3;
//    scanf("%d%d%d",&a,&b,&c);


connection=0;

for ( j = 1; j <= b; j=j+1 ) {
	for (k = 1; k <= c; k=k+1) {
		for ( i = 1; i <= a; i=i+1 ) {
			// bottom plane
			if (j==1) {
				// bottom, back line
				if (k==1) {
					// bottom, back, left box
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t back = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					// bottom, back, right box
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t back = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t back = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
				}
				// bottom plane, middle line series
				else if (k>1 && k<=6) {
					// bottom, back, left box
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
		// 				printf(" \
		// 				\t pid = 0 \n");
						printf("\
						\t left = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						} \n");
					}
				}
				// bottom plane, side walls, exclude front boxes
				else if (k>6 && k<c) {
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
		// 				printf(" \
		// 				\t pid = 0 \n");
						printf("\
						\t left = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else if (i>1 && i<32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						} \n");
					}
					else if (i==32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t right = Boundary \n \
						} \n");
					}
					else if (i==33) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t right = Boundary \n \
						\t left = Boundary \n \
						} \n");
					}
					else if (i==34) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t left = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						} \n");
					}
				}
				// bottom plane, front boxes
				else if (k==c) {
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
		// 				printf(" \
		// 				\t pid = 0 \n");
						printf("\
						\t left = Boundary \n \
						\t front = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t front = Boundary \n \
						\t bottom = Boundary \n \
						} \n");
					}
					else if (i>1 && i<32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t right = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==33) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t right = Boundary \n \
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==34) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t bottom = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
				}
			}

			// top plane
			else if (j==b) {
				// top, back line
				if (k==1) {
					// top, back, left box
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t back = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					// top, back, right box
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t back = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t back = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
				}
				// top plane, middle line series
				else if (k>1 && k<=6) {
					// top, back, left box
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
		// 				printf(" \
		// 				\t pid = 0 \n");
						printf("\
						\t left = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						} \n");
					}
				}
				// top plane, side walls, exclude front boxes
				else if (k>6 && k<c) {
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else if (i>1 && i<32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						} \n");
					}
					else if (i==32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t right = Boundary \n \
						} \n");
					}
					else if (i==33) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t right = Boundary \n \
						\t left = Boundary \n \
						} \n");
					}
					else if (i==34) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t left = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						} \n");
					}
				}
				// top plane, front boxes
				else if (k==c) {
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t front = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t front = Boundary \n \
						\t top = Boundary \n \
						} \n");
					}
					else if (i>1 && i<32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t right = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==33) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t right = Boundary \n \
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==34) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t top = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
				}
			}

			// middle plane
			else {
				// middle, back line
				if (k==1) {
					// middle, back, left box
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t back = Boundary \n \
						} \n");
					}
					// middle, back, right box
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t back = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t back = Boundary \n \
						} \n");
					}
				}
				// middle plane, middle line series
				else if (k>1 && k<=6) {
					// middle, back, left box
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
		// 				printf(" \
		// 				\t pid = 0 \n");
						printf("\
						\t left = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						} \n");
					}
				}
				// middle plane, side walls, exclude front boxes
				else if (k>6 && k<c) {
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						} \n");
					}
					else if (i>1 && i<32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						} \n");
					}
					else if (i==32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						} \n");
					}
					else if (i==33) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t left = Boundary \n \
						} \n");
					}
					else if (i==34) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						} \n");
					}
				}
				// middle plane, front boxes
				else if (k==c) {
					if (i==1) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==a) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i>1 && i<32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t front = Boundary \n \
						} \n");
					}
					else if (i==32) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==33) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t right = Boundary \n \
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else if (i==34) {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t left = Boundary \n \
						\t front = Boundary \n \
						} \n");
					}
					else {
						printf("# No. %d box \n", (((k-1)*a+i)+(j-1)*a*c));
						printf("GfsBox {\n");
						printf("\
						\t front = Boundary \n \
						} \n");
					}
				}
			}
		}
	}
}

// connectivity
// longitudinal sweep
for( j = 1; j <= b; j = j + 1 ) {
	for (k = 1; k <= c; k=k+1 ) {
		for( i = 1; i <= (a-1); i = i + 1 ) {
			box = ((k-1)*a+i)+(j-1)*a*c;
			printf("%d %d right\n", box, box+1);
			connection=connection+1;
		}
	}
}

// transverse sweep
for( j = 1; j <= b; j = j + 1 ) {
	for( i = 1; i <= a; i = i + 1 ) {
		for (k = 1; k <= (c-1); k = k + 1 ) {
			box = ((k-1)*a+i)+(j-1)*a*c;
			printf("%d %d front\n", box, box+a);
			connection=connection+1;
		}
	}
}

// vertical sweep
for( k = 1; k <= c; k = k + 1 ) {
	for( i = 1; i <= a; i = i + 1 ) {
		for (j = 1; j <= (b-1); j = j + 1 ) {
			box = ((k-1)*a+i)+(j-1)*a*c;
			printf("%d %d top\n", box, box+(a*c));
			connection=connection+1;
		}
	}
}

// connectivity summary, to be deleted
printf("\n\n\n\nThe Number of connection is: %d",connection);
return 0;
}
