
MANPAGES = $(basename $(wildcard  *.1.md))

all: $(MANPAGES)

%.1: %.1.md
	pandoc -s --from=markdown --to=man $< > $@

clean:
	$(RM) $(MANPAGES)

.PHONY: clean
