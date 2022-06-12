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
C-syntax is good enough. Follow these threads:

*https://unix.stackexchange.com/questions/137989/how-to-change-default-syntax-highlighting-for-header-file-in-kate*

*https://gitlab.com/-/snippets/4499*
3. Possible error in some distros:
```
gerris: file `column.gfs' is not a valid simulation file
column.gfs:212:12: cannot load module: /usr/local/lib/gerris/libgfsview2D.so: cannot open shared object file: No such file or directory
```
Solution: make a symbolic link to the correct path:
```
sudo ln /usr/lib/x86_64-linux-gnu/gerris/libgfsview2D-0.0.1.so libgfsview2D.so

sudo ln /usr/lib/x86_64-linux-gnu/gerris/libgfsview3D-0.0.1.so libgfsview3D.so
```
