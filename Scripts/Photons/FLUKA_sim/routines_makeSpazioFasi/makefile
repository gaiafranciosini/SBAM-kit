default: compile link


compile:
#	rm *.o
	$(FLUPRO)/flutil/fff  UpdateCurrentParticle.f
	$(FLUPRO)/flutil/fff  mgdraw.f
	$(FLUPRO)/flutil/fff  mgdraw_lib.f
	$(FLUPRO)/flutil/fff  usrini.f
	$(FLUPRO)/flutil/fff  usrein.f
	$(FLUPRO)/flutil/fff  usreou.f
	$(FLUPRO)/flutil/fff  usrout.f
	$(FLUPRO)/flutil/fff  source.f	

link:
	$(FLUPRO)/flutil/ldpm3qmd usrini.o usrout.o usrein.o usreou.o mgdraw_lib.o mgdraw.o source.o UpdateCurrentParticle.o -o SAFEST_10MeV.exe
