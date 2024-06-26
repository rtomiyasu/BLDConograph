# object files linked for building exe files

# ROOT
OBJS_INDEXING += \
object/ControlFile.o \
object/ControlParam.o \
object/main_indexing.o \

# utility_lattice_reduction
OBJS_INDEXING += \
object/utility_lattice_reduction/error_stable_bravais_3D.o \
object/utility_lattice_reduction/mlist.o \

# utility_func
OBJS_INDEXING += \
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

