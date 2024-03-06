
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/lattice_symmetry/lattice_symmetry.cc \
$(ROOT)/lattice_symmetry/LatticeFigureOfMerit.cc \
$(ROOT)/lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.cc \
$(ROOT)/lattice_symmetry/ReducedLatticeToCheckBravais.cc \

OBJS += \
object/lattice_symmetry/lattice_symmetry.o \
object/lattice_symmetry/LatticeFigureOfMerit.o \
object/lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.o \
object/lattice_symmetry/ReducedLatticeToCheckBravais.o \

DEPS += \
${addprefix object/lattice_symmetry/, \
lattice_symmetry.d \
LatticeFigureOfMerit.d \
LatticeFigureOfMeritToCheckSymmetry.d \
ReducedLatticeToCheckBravais.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/lattice_symmetry/%.o: $(ROOT)/lattice_symmetry/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '


