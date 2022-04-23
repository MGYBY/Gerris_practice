# Gerris practice

Some Gerris (Popinet 2003) codes.

Some tips and tricks:
1. Built-from-source code change:
Line 416 of "init.c":
change to:
```
MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
```
2. Gerris-style syntax-highlighting:
C-syntax is good enough. Follow this thread:

*https://unix.stackexchange.com/questions/137989/how-to-change-default-syntax-highlighting-for-header-file-in-kate*
