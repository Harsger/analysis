.PHONY : all
all:fitNclust investigateCRF track postprocessor sampleEvents converter studyCRF qaqc moduleRawTracker

fitNclust: fitter.C analysis.h dicvecvec.o 
	g++ -std=c++11 -Ofast fitter.C dicvecvec.o -o fitNclust `pkg-config opencv --cflags` `pkg-config opencv --libs` `root-config --glibs --cflags`

investigateCRF: investigator.C analysis.h dicvecvec.o 
	g++ -std=c++11 -Ofast investigator.C dicvecvec.o -o investigateCRF `pkg-config opencv --cflags` `pkg-config opencv --libs` `root-config --glibs --cflags`

track: tracker.C analysis.h dicvecvec.o 
	g++ -std=c++11 -Ofast tracker.C dicvecvec.o -o track `pkg-config opencv --cflags` `pkg-config opencv --libs` `root-config --glibs --cflags`

postprocessor: postprocessor.C analysis.h dicvecvec.o
	g++ -std=c++11 -Ofast postprocessor.C dicvecvec.o -o postprocessor `pkg-config opencv --cflags` `pkg-config opencv --libs` `root-config --glibs --cflags`

sampleEvents: sampleEvents.C analysis.h dicvecvec.o
	g++ -std=c++11 -Ofast sampleEvents.C dicvecvec.o -o sampleEvents `pkg-config opencv --cflags` `pkg-config opencv --libs` `root-config --glibs --cflags`

converter: converter.C dicvecvec.o
	g++ -std=c++11 -Ofast converter.C dicvecvec.o -o converter `pkg-config opencv --cflags` `pkg-config opencv --libs` `root-config --glibs --cflags`

dicvecvec.o: headvecvec.h analysis.h
	rootcint -f dicvecvec.cxx -c headvecvec.h
	g++ -std=c++11 `root-config --glibs --cflags` -c dicvecvec.cxx

studyCRF: studyCRF.C
	g++ -std=c++11 -Ofast studyCRF.C -o studyCRF `root-config --glibs --cflags`

qaqc: qaqc.C
	g++ -std=c++11 -Ofast qaqc.C -o qaqc `root-config --glibs --cflags`

moduleRawTracker: moduleRawTracker.C dicvecvec.o
	g++ -std=c++11 -Ofast moduleRawTracker.C dicvecvec.o -o moduleRawTracker `root-config --glibs --cflags`

clean:
	rm fitNclust investigateCRF track postprocessor sampleEvents converter dicvecvec.cxx dicvecvec.o studyCRF qaqc moduleRawTracker
