LIBNAME=tripack
#LD=ld
#FFLAGS=-g -fpic
#FFLAGS=-g -O -fpic
#SHLIBLDFLAGS=-shared

OBJS=	tripack.o inhull.o voronoi.o


$(LIBNAME):	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS) #-lf2c

clean:
	@rm -f *.o *.so tritest fort.2

realclean:
	@rm -f Makefile *.o *.so

test: tritest 

tritest: tritest.f tripack.o
	$(F77)  $(FFLAGS) tritest.f tripack.o -o tritest


# compilation with f77 -g -fpic fails on Linux with g77 0.5.19.1!
# (assembler complains about unknown i386 instructions, so it seems to be
#  an g77 internal problem)
# But it works when the PICFLAG is omitted:
#%.o: %.f
#	$(F77) -g -c $< -o $@
# compilation with f77 -O2 -fpic works, but now interpp() (it calls idbvip)
# crashes R!

