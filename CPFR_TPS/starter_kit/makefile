default: compile link


compile:
#	rm *.o
	$(FLUPRO)/flutil/fff  source_mod.f


link:
	$(FLUPRO)/flutil/ldpm3qmd source_mod.o -o fluka_EF_9MeV.exe

clean:
	rm *.o
