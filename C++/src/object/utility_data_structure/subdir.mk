
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/utility_data_structure/VCData.cc \

OBJS += \
object/utility_data_structure/VCData.o \

DEPS += \
${addprefix object/utility_data_structure/, \
VCData.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/utility_data_structure/%.o: $(ROOT)/utility_data_structure/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '


