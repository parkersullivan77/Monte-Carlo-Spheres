raytracer: main.o raytracer.o cam.o light.o sphere.o ray.o triangle.o
	g++ -w -O3 main.o raytracer.o cam.o light.o sphere.o ray.o triangle.o -o raytracer
main.o: main.cc
	g++ -w -O3 -c main.cc 

raytracer.o: raytracer.cc raytracer.h
	g++ -w -O3 -c raytracer.cc raytracer.h

cam.o: cam.cc cam.h
	g++ -w -O3 -c cam.cc cam.h

light.o: light.cc light.h
	g++ -w -O3 -c light.cc light.h
	

sphere.o:sphere.cc sphere.h
	g++ -w -O3 -c sphere.cc sphere.h

ray.o:ray.cc ray.h
	g++ -w -O3 -c ray.cc ray.h

triangle.o:triangle.cc triangle.h
	g++ -w -O3 -c triangle.cc triangle.h
	
clean:
	rm *.o raytracer


