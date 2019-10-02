# IMPOTRANT NOTE: Unless set in the environment,
# values of the following variables need to be supplied via the command line
# (1) $(BOOST_DIR)  # default is empty, path to boost library (only the headers will be needed)
# (2) $(PROCESSOR)  # may take values 0 or 1: to set use of the pseudo-random number generator
 					# from boost or (in case of cuda availablility) curand respectively; default is 1
# You may also need to supply explicitely executables for $(AR), $(AS) and $(CC)
# to override or supplement the environment variables
# e.g., the call from the directory where the Makefile is, may look like:
# make PROCESSOR=0 -B all BOOST_DIR=/path/to/boost_includes AS=/usr/bin/as AR=/usr/bin/ar CC=/usr/local/gcc-9.1.0/bin/gcc-9.1.0
# Alternatively, hardcode them directly here.
#
# The executable is ../app/mitosim_X, where "X" is either "CPU" or "CUDA", depending on the $(PROCESSOR)
#
# To check the actual value of a variable $(VAR):
# make print-VAR

print-%  : ; @echo $* = $($*)

UTILS_DIR := ../external/

MITOSIM_DIR := ../
MITOSIM_APP_DIR := $(MITOSIM_DIR)app/

SRCS_MAIN := main.cpp
SRCS_MAIN_FULL := $(patsubst %,$(MITOSIM_APP_DIR)%,$(SRCS_MAIN))
SRCS_FULL += $(SRCS_MAIN_FULL)

OBJS_DIR := ./obj/
OBJS += $(SRCS_MAIN:.cpp=.o)
OBJS_FULL = $(patsubst %,$(OBJS_DIR)%,$(OBJS))

ASMS_DIR := asm/
ASMS += $(SRCS_MAIN:.cpp=.s)
ASMS_FULL = $(patsubst %,$(ASMS_DIR)%,$(ASMS))

MITOSIM_INC_DIR := $(MITOSIM_DIR)include/
REACTS_INC_DIR := $(MITOSIM_DIR)include/reactions/
INCDIR = -I$(BOOST_DIR) -I$(UTILS_DIR) -I$(MITOSIM_INC_DIR) -I$(REACTS_INC_DIR)

TARGET_DIR := ../app/
TARGET_CORE := mitosim

PROCESSOR_CPU := 0
PROCESSOR_CUDA := 1

ifeq ($(PROCESSOR), $(PROCESSOR_CPU))
	CPP_VARS += -DPROCESSOR=0
	TARGET := $(TARGET_CORE)_CPU
else
	CPP_VARS += -DPROCESSOR=1
	TARGET := $(TARGET_CORE)_CUDA
endif
TARGET := $(TARGET_DIR)$(TARGET)
all : $(TARGET)

CC := gcc
CFLAGS = -c -Wall -O3 -msse4 -march=native $(INCDIR) $(CPP_VARS) -std=c++17
LFLAGS = -Wall -L$(UTILS_DIR)lib -lpthread -lstdc++ -lutils

$(TARGET) : $(OBJS_FULL)
	$(CC) $(LFLAGS) $^ -o $@

$(OBJS_DIR)main.o : $(ASMS_DIR)main.s
	$(AS) $< -o $@

$(ASMS_DIR)main.s : $(MITOSIM_APP_DIR)main.cpp \
	$(REACTS_INC_DIR)ntw_fission.h \
	$(REACTS_INC_DIR)ntw_fusion11.h \
	$(REACTS_INC_DIR)ntw_fusion12.h \
	$(REACTS_INC_DIR)ntw_fusion1u.h \
	$(REACTS_INC_DIR)reaction.h \
	$(REACTS_INC_DIR)fission.h \
	$(REACTS_INC_DIR)fusion.h \
	$(REACTS_INC_DIR)fusion11.h \
	$(REACTS_INC_DIR)fusion12.h \
	$(REACTS_INC_DIR)fusion1u.h \
	$(MITOSIM_INC_DIR)fusion_candidates.h \
	$(MITOSIM_INC_DIR)edge.h \
	$(MITOSIM_INC_DIR)segment.h \
	$(MITOSIM_INC_DIR)structure.h \
	$(MITOSIM_INC_DIR)core_transformer.h \
	$(MITOSIM_INC_DIR)ability_for_fission.h \
	$(MITOSIM_INC_DIR)ability_for_fusion.h \
	$(MITOSIM_INC_DIR)network.h \
	$(MITOSIM_INC_DIR)edge.h \
	$(MITOSIM_INC_DIR)segment.h \
	$(MITOSIM_INC_DIR)simulation.h \
	$(MITOSIM_INC_DIR)network.h \
	$(MITOSIM_INC_DIR)config.h
	$(CC) $(CFLAGS) -S $< -o $@

#clean:

.PHONY: clean

clean:
	rm -rf $(OBJS_DIR)*.o $(ASMS_DIR)*.s $(TARGET)
