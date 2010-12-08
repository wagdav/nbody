NBody: NBody.chpl
	chpl -o nbody NBody.chpl

clean:
	rm -f nbody
