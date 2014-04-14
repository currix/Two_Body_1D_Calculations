BINPATH = ../bin/
#
#     HO
2body_1D_HO_wsaxon_SRC = nrtype.f90 constants.f90 pot_param.f90	\
var_definition.f90 module_1body_ho.f90 module_2body.f90 module_2body_ho.f90 \
wsaxon_pot.f90 ddeltak.f90 x_ho_matel.f90 x2_ho_matel.f90 funsub.f90 \
e_gstate_ho.f90 compute_kmin.f90 res_int_2body_HO.f90
#
2body_1D_HO_wsaxon_OBJ = $(2body_1D_HO_wsaxon_SRC:.f90=.o)
#
#     THO
2body_1D_THO_wsaxon_SRC = nrtype.f90 constants.f90 pot_param.f90	\
var_definition.f90 module_1body_tho.f90 module_2body.f90 module_2body_tho.f90 \
wsaxon_pot.f90 ddeltak.f90 \
e_gstate_tho.f90 compute_amin_tho.f90 res_int_2body_THO.f90
#
2body_1D_THO_wsaxon_OBJ = $(2body_1D_THO_wsaxon_SRC:.f90=.o)
#
#     ISQW
2body_1D_ISQW_wsaxon_SRC = nrtype.f90 constants.f90 pot_param.f90	\
var_definition.f90 module_1body_isqw.f90 module_2body.f90 module_2body_isqw.f90 \
wsaxon_pot.f90 ddeltak.f90 x_isqw_matel.f90 x2_isqw_matel.f90 \
 res_int_2body_ISQW.f90
#
2body_1D_ISQW_wsaxon_OBJ = $(2body_1D_ISQW_wsaxon_SRC:.f90=.o)
#
#
FC = gfortran
FOPT	= -c 
FOPT   += -Wall
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
FOPT   += $(MODLIB)
#
##LAPACK95_PATH = /usr/local/lib/lapack95
##LAPACK95 = $(LAPACK95_PATH)/lapack95.a
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
# Cluster
LAPACK77 = -llapack
BLAS     = -lblas
# NAG
NAGLIB = -L/usr/local/lib/NAG -lnag
#
LIBS     = $(NAGLIB) $(LAPACK95) $(LAPACK77) $(BLAS)

all: HO THO ISQW
all_static: HO_static THO_static ISQW_static
#
HO: wsaxon_HO_2body_1D
HO_static: wsaxon_HO_2body_1D_static
#
THO: wsaxon_THO_2body_1D
THO_static: wsaxon_THO_2body_1D_static
#
ISQW: wsaxon_ISQW_2body_1D
ISQW_static: wsaxon_ISQW_2body_1D_static
#
.PHONY : all
#
# Standard compilation rules
#
%.o: %.f90
	$(FC) $(FOPT) -c $<
#
%.o: %.F90
	$(FC) $(FOPT) -c $<
#
%.o: %.f
	$(FC) $(FOPT) $<
#
#####################################################################
wsaxon_HO_2body_1D: $(2body_1D_HO_wsaxon_OBJ)
	$(FC) -o $(BINPATH)/wsaxon_2body_1D_HO  $(2body_1D_HO_wsaxon_OBJ) $(LIBS) 
wsaxon_HO_2body_1D_static: $(2body_1D_HO_wsaxon_OBJ)
	$(FC) -static -o $(BINPATH)/wsaxon_2body_1D_HO  $(2body_1D_HO_wsaxon_OBJ) $(LIBS) 
#####################################################################
wsaxon_THO_2body_1D: $(2body_1D_THO_wsaxon_OBJ)
	$(FC) -o $(BINPATH)/wsaxon_2body_1D_THO  $(2body_1D_THO_wsaxon_OBJ) $(LIBS) 
wsaxon_THO_2body_1D_static: $(2body_1D_THO_wsaxon_OBJ)
	$(FC) -static -o $(BINPATH)/wsaxon_2body_1D_THO  $(2body_1D_THO_wsaxon_OBJ) $(LIBS) 
####################################################################
wsaxon_ISQW_2body_1D: $(2body_1D_ISQW_wsaxon_OBJ)
	$(FC) -o $(BINPATH)/wsaxon_2body_1D_ISQW  $(2body_1D_ISQW_wsaxon_OBJ) $(LIBS) 
wsaxon_ISQW_2body_1D_static: $(2body_1D_ISQW_wsaxon_OBJ)
	$(FC) -static -o $(BINPATH)/wsaxon_2body_1D_ISQW  $(2body_1D_ISQW_wsaxon_OBJ) $(LIBS) 
#####################################################################
clean:
	@rm -f $(2body_1D_HO_wsaxon_OBJ) $(2body_1D_THO_wsaxon_OBJ) $(2body_1D_ISQW_wsaxon_OBJ)
#####################################################################



