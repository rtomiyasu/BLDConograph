
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \

OBJS += \

DEPS += \
${addprefix object/utility_lattice_reduction/, \
}


# Each subdirectory must supply rules for building sources it contributes
object/utility_lattice_reduction/%.o: $(ROOT)/utility_lattice_reduction/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '


