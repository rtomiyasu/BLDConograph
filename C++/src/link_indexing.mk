# object files linked for building exe files

# ROOT
OBJS_INDEXING += \
object/ControlFile.o \
object/ControlParam.o \
object/main_indexing.o \
object/p_out.o \

# bravais_type
OBJS_INDEXING += \
object/bravais_type/BravaisType.o \
object/bravais_type/BravaisLattice.o \


# centring_type
OBJS_INDEXING += \
object/centring_type/CentringType.o \
object/centring_type/bravais_lat.o \


# lattice_symmetry
OBJS_INDEXING += \
object/lattice_symmetry/lattice_symmetry.o \
object/lattice_symmetry/LatticeFigureOfMerit.o \
object/lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.o \
object/lattice_symmetry/ReducedLatticeToCheckBravais.o \


# point_group
OBJS_INDEXING += \
object/point_group/coset_representative_data.o \
object/point_group/point_gp_data.o \


# symmetric_operation
OBJS_INDEXING += \
object/symmetric_operation/MillerIndex.o \
object/symmetric_operation/S1.o \
object/symmetric_operation/SymmetricOperation.o \
object/symmetric_operation/translation_vector.o \
object/symmetric_operation/StringS1.o \


# utility_data_structure
OBJS_INDEXING += \
object/utility_data_structure/VCData.o \


# utility_func
OBJS_INDEXING += \
object/utility_func/gcd.o \
object/utility_func/lattice_constant.o \
object/utility_func/zstring.o \


# utility_rw_param
OBJS_INDEXING += \
object/utility_rw_param/I_ReadData.o \
object/utility_rw_param/RWParam_void.o \


#zerror_type
OBJS_INDEXING += \
object/zerror_type/error_mes.o \


# zlog
OBJS_INDEXING += \
object/zlog/rlog.o \
object/zlog/zlog.o \

