FiguresGP = $(wildcard *.gp)
FiguresPDF = $(patsubst %.gp, %.pdf, $(FiguresGP))

all: clean $(FiguresPDF)

clean:
	$(RM) -rv $(FiguresPDF)

%.pdf: %.gp
	gnuplot $^
	pdflatex $(patsubst %.pdf, %.tex, $@) > dump.txt