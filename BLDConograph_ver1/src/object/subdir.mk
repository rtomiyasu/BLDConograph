
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/ControlFile.cc \
$(ROOT)/ControlParam.cc \
$(ROOT)/p_out.cc \
$(ROOT)/main_indexing.cc \

OBJS += \
object/ControlFile.o \
object/ControlParam.o \
object/p_out.o \
object/main_indexing.o \

DEPS += \
${addprefix object/, \
ControlFile.d \
ControlParam.d \
p_out.d \
main_indexing.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/%.o: $(ROOT)/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '
