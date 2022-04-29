#
# Makefile for the program tirade
#
# Turn on debugging
TARGET=tirade.x

SILENCE_FLAG +=
CC=gcc
CFLAGS=-g $(SILENCE_FLAG)
TEX=latex
BIB=bibtex

SRCDIR=src
OBJDIR=obj
BINDIR=bin
DOCDIR=doc/LaTeX

LIBS=-lm

SRC=$(SRCDIR)/tirade.c $(SRCDIR)/tidal_module.c $(SRCDIR)/conduction_module.c $(SRCDIR)/utility.c
#OBJ=$(OBJDIR)/tirade.o $(OBJDIR)/tidal_module.o $(OBJDIR)/conduction_module.o $(OBJDIR)/utility.o
OBJ=$(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
DOC=tidal_doc

all: $(TARGET)

$(TARGET): $(OBJ)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $(BINDIR)/$(TARGET) $(OBJ) $(LIBS)
	@echo "Target built!"

$(OBJ): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# Option to compile the LaTeX documentation.
# Not really a need for this on the user end since the pdf is included in the package
# Will need access to the tex file, the bibliography file, the bib. style file and 
# the figures. 
#docs:
#	cd $(DOCDIR); pwd; \
#	$(TEX) $(DOC).tex; \
#	$(BIB) $(DOC); \
#	$(TEX) $(DOC).tex; \
#	$(TEX) $(DOC).tex; \
#	dvipdf $(DOC).dvi ../pdf/$(DOC).pdf

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(TARGET)
	@echo "All clean!"
