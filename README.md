# Gerris practice

Some Gerris (Popinet 2003) codes.

Built-from-source code change:
Line 416 of "init.c":
change to:
```
MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
```
