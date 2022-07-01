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
**Gfsview cannot be installed in Ubuntu 22.04 from source?**

4. Gfsview compilation error for Ubuntu 22.04:
```
/usr/bin/ld: gfsview2D-main.o: undefined reference to symbol 'XSync'
/usr/bin/ld: /lib/x86_64-linux-gnu/libX11.so.6: error adding symbols: DSO missing from command line
collect2: error: ld returned 1 exit status
make[4]: *** [Makefile:395: gfsview2D] Error 1
```
Solution: include necessary libraries
```
./configure LIBS=-lX11
```

5. Some good reads:

*https://resources.oreilly.com/examples/9781565922259/*
