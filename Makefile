CXX = g++
CXXFLAGS += -O2 -Wall -Werror -Wextra
ROOTFLAGS := `root-config --cflags --libs`

BUILDDIR = ./build

SRCS = $(wildcard *.C)
EXES = $(patsubst %.C,%.exe,$(SRCS))
DEPS = $(patsubst %.C,$(BUILDDIR)/%.d,$(SRCS))

.PHONY: all clean

all: $(EXES)

%.exe: %.C
	@mkdir -p $(BUILDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -MMD -MF $(BUILDDIR)/$(@D)/$(*F).d $< -o $@

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -rf $(BUILDDIR)/*

-include $(DEPS)
