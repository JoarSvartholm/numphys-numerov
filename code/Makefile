CC = gcc
LOADLIBES =  -L/usr/local/lib -lgsl -lgslcblas -lm
CFLAGS = -g -I/usr/local/include
OBJS =
DEPS = 

PIC = -fPIC

#gaussian: gaussian.c
#	gcc $(CFLAGS) gaussian.c -o gaussian.o $(LOADLIBES)


harm: harm.o $(OBJS)
	$(CC) $(CFLAGS)  $(PIC) -Wall $(INCLUDE)  $^ -o $@ $(LOADLIBES)


%:%.o $(OBJS)#	 $(DEPS)
	$(CC) $(CFLAGS)  $(PIC) -Wall $(INCLUDE)  $^   -o $@ $(LOADLIBES)

%.o:%.c $(DEPS)
		$(CC) $(CFLAGS) -c  $(PIC) -Wall $(INCLUDE)  $< -o $@ #$(LOADLIBES)

#$(LIBODE): $(OBJS)
#	ar -cr $@  $?

clean:
	rm -rf *.so *.dylib *.o *~ *.a
